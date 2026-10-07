#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

#include "nuX_M1_closure.hxx"
#include "nuX_M1_sources.hxx"
#include "nuX_baryon_mass.hxx"
#include "nuX_utils.hxx"

namespace nuX_M1 {

using namespace std;
using namespace Loop;
using namespace nuX_Utils;

// Fallback if header didn't set it
#ifndef MAX_GROUPSPECIES
#define MAX_GROUPSPECIES 3
#endif

// Map closure string to a device-safe closure selector.
static inline closure_t pick_closure(const char *name) {
  if (CCTK_Equals(name, "Eddington"))
    return CLOSURE_EDDINGTON;
  if (CCTK_Equals(name, "Kershaw"))
    return CLOSURE_KERSHAW;
  if (CCTK_Equals(name, "Minerbo"))
    return CLOSURE_MINERBO;
  if (CCTK_Equals(name, "thin"))
    return CLOSURE_THIN;
  char msg[BUFSIZ];
  snprintf(msg, BUFSIZ, "Unknown closure \"%s\"", name);
  CCTK_ERROR(msg);
  return CLOSURE_EDDINGTON; // not reached
}

void launch_source_increment_kernel(
    GridDescBaseDevice const &grid, GF3D2layout const layout_cc,
    tensor::slicing_geometry_const const geom,
    tensor::fluid_velocity_field_const const fidu, cGH const *const cctkGH,
    closure_t const closure_fun, CCTK_REAL const dt, CCTK_INT const ngroups,
    CCTK_INT const nspecies, CCTK_REAL const *const nuX_m1_mask,
    CCTK_REAL const *const fidu_w_lorentz, CCTK_REAL const *const fidu_velx,
    CCTK_REAL const *const fidu_vely, CCTK_REAL const *const fidu_velz,
    CCTK_REAL *const rN, CCTK_REAL *const rE, CCTK_REAL *const rFx,
    CCTK_REAL *const rFy, CCTK_REAL *const rFz, CCTK_REAL *const chi,
    CCTK_REAL const *const eta_0, CCTK_REAL const *const eta_1,
    CCTK_REAL const *const abs_0, CCTK_REAL const *const abs_1,
    CCTK_REAL const *const scat_1, CCTK_REAL const *const nueave,
    CCTK_REAL *const source_update_delta_N,
    CCTK_REAL *const source_update_delta_E,
    CCTK_REAL *const source_update_delta_Fx,
    CCTK_REAL *const source_update_delta_Fy,
    CCTK_REAL *const source_update_delta_Fz, CCTK_REAL const rad_N_floor,
    CCTK_REAL const rad_E_floor, CCTK_REAL const rad_eps,
    CCTK_REAL const closure_epsilon, CCTK_INT const closure_maxiter,
    CCTK_INT const use_fallback, CCTK_INT const source_force_nonlinear,
    CCTK_REAL const source_thick_limit, CCTK_REAL const source_scat_limit,
    CCTK_REAL const source_therm_limit, CCTK_INT const source_maxiter,
    CCTK_REAL const source_epsabs, CCTK_REAL const source_epsrel) {
  grid.loop_int_device<1, 1, 1>(grid.nghostzones, [=] CCTK_DEVICE(
                                                      const PointDesc &p) {
    const int ijk = layout_cc.linear(p.i, p.j, p.k);
    int const groupspec = ngroups * nspecies;

    if (nuX_m1_mask[ijk]) {
      for (int ig = 0; ig < groupspec; ++ig) {
        int const i4D = layout_cc.linear(p.i, p.j, p.k, ig);
        source_update_delta_N[i4D] = 0.0;
        source_update_delta_E[i4D] = 0.0;
        source_update_delta_Fx[i4D] = 0.0;
        source_update_delta_Fy[i4D] = 0.0;
        source_update_delta_Fz[i4D] = 0.0;
      }
      return;
    }

    tensor::generic<CCTK_REAL, 4, 1> beta_u;
    geom.get_shift_vec(p, &beta_u);
    const CCTK_REAL alp_ijk = geom.get_lapse(p);
    const CCTK_REAL betax_ijk = beta_u(1);
    const CCTK_REAL betay_ijk = beta_u(2);
    const CCTK_REAL betaz_ijk = beta_u(3);
    const CCTK_REAL W_ijk = fidu_w_lorentz[ijk];

    tensor::metric<4> g_dd;
    tensor::inv_metric<4> g_uu;
    tensor::generic<CCTK_REAL, 4, 1> n_u, n_d;
    tensor::generic<CCTK_REAL, 4, 2> gamma_ud;
    geom.get_metric(p, &g_dd);
    geom.get_inv_metric(p, &g_uu);
    geom.get_normal(p, &n_u);
    geom.get_normal_form(p, &n_d);
    geom.get_space_proj(p, &gamma_ud);
    const CCTK_REAL volform_ijk = sqrt(
        nuX_Utils::metric::spatial_det(g_dd(1, 1), g_dd(1, 2), g_dd(1, 3),
                                       g_dd(2, 2), g_dd(2, 3), g_dd(3, 3)));

    tensor::generic<CCTK_REAL, 4, 1> u_u, u_d;
    tensor::generic<CCTK_REAL, 4, 2> proj_ud;
    fidu.get(p, &u_u);
    tensor::contract(g_dd, u_u, &u_d);
    calc_proj(u_d, u_u, &proj_ud);

    tensor::generic<CCTK_REAL, 4, 1> v_u, v_d;
    pack_v_u(fidu_velx[ijk], fidu_vely[ijk], fidu_velz[ijk], &v_u);
    tensor::contract(g_dd, v_u, &v_d);

    for (int ig = 0; ig < groupspec; ++ig) {
      int const i4D = layout_cc.linear(p.i, p.j, p.k, ig);
      assert(isfinite(rN[i4D]));
      assert(isfinite(rE[i4D]));
      assert(isfinite(rFx[i4D]));
      assert(isfinite(rFy[i4D]));
      assert(isfinite(rFz[i4D]));

      CCTK_REAL Estar = rE[i4D];
      tensor::generic<CCTK_REAL, 4, 1> Fstar_d;
      pack_F_d(betax_ijk, betay_ijk, betaz_ijk, rFx[i4D], rFy[i4D], rFz[i4D],
               &Fstar_d);
      apply_floor(g_uu, &Estar, &Fstar_d, rad_E_floor, rad_eps);
      assert(isfinite(Estar));
      assert(isfinite(Fstar_d(1)));
      assert(isfinite(Fstar_d(2)));
      assert(isfinite(Fstar_d(3)));

      if (source_rates_are_zero(eta_0[i4D], eta_1[i4D], abs_0[i4D], abs_1[i4D],
                                scat_1[i4D])) {
        tensor::symmetric2<CCTK_REAL, 4, 2> P_dd;
        calc_closure(cctkGH, p.i, p.j, p.k, ig, closure_fun, g_dd, g_uu, n_d,
                     W_ijk, u_u, v_d, proj_ud, Estar, Fstar_d, &chi[i4D], &P_dd,
                     closure_epsilon, closure_maxiter, use_fallback != 0);
        source_update_delta_N[i4D] = 0.0;
        source_update_delta_E[i4D] = 0.0;
        source_update_delta_Fx[i4D] = 0.0;
        source_update_delta_Fy[i4D] = 0.0;
        source_update_delta_Fz[i4D] = 0.0;
        continue;
      }

      CCTK_REAL Nstar = std::max(rN[i4D], rad_N_floor);
      CCTK_REAL Enew = Estar;
      tensor::generic<CCTK_REAL, 4, 1> Fnew_d;

      tensor::symmetric2<CCTK_REAL, 4, 2> P_dd;
      calc_closure(cctkGH, p.i, p.j, p.k, ig, closure_fun, g_dd, g_uu, n_d,
                   W_ijk, u_u, v_d, proj_ud, Estar, Fstar_d, &chi[i4D], &P_dd,
                   closure_epsilon, closure_maxiter, use_fallback != 0);

      tensor::symmetric2<CCTK_REAL, 4, 2> rT_dd;
      assemble_rT(n_d, Estar, Fstar_d, P_dd, &rT_dd);

      CCTK_REAL const Jstar = calc_J_from_rT(rT_dd, u_u);
      tensor::generic<CCTK_REAL, 4, 1> Hstar_d;
      calc_H_from_rT(rT_dd, u_u, proj_ud, &Hstar_d);

      const CCTK_REAL dtau = alp_ijk * dt / W_ijk;
      CCTK_REAL Jnew =
          (Jstar + dtau * eta_1[i4D] * volform_ijk) / (1 + dtau * abs_1[i4D]);

      CCTK_REAL const khat = abs_1[i4D] + scat_1[i4D];
      tensor::generic<CCTK_REAL, 4, 1> Hnew_d;
      for (int a = 1; a < 4; ++a)
        Hnew_d(a) = Hstar_d(a) / (1 + dtau * khat);
      Hnew_d(0) = 0.0;
      for (int a = 1; a < 4; ++a)
        Hnew_d(0) -= Hnew_d(a) * (u_u(a) / u_u(0));

      CCTK_REAL const H2 = tensor::dot(g_uu, Hnew_d, Hnew_d);

      chi[i4D] = CCTK_REAL(1.0 / 3.0);

      CCTK_REAL const dthick = 3.0 * (1.0 - chi[i4D]) / 2.0;
      CCTK_REAL const dthin = 1.0 - dthick;

      for (int a = 0; a < 4; ++a)
        for (int b = a; b < 4; ++b) {
          rT_dd(a, b) =
              Jnew * u_d(a) * u_d(b) + Hnew_d(a) * u_d(b) + Hnew_d(b) * u_d(a) +
              dthin * Jnew * (H2 > 0 ? Hnew_d(a) * Hnew_d(b) / H2 : 0.0) +
              dthick * Jnew * (g_dd(a, b) + u_d(a) * u_d(b)) / 3.0;
        }

      Enew = calc_J_from_rT(rT_dd, n_u);
      calc_H_from_rT(rT_dd, n_u, gamma_ud, &Fnew_d);
      apply_floor(g_uu, &Enew, &Fnew_d, rad_E_floor, rad_eps);

      SourceUpdateContext source_ctx(
          cctkGH, p.i, p.j, p.k, ig, closure_epsilon, closure_maxiter,
          use_fallback != 0, dt, alp_ijk, g_dd, g_uu, n_d, n_u, gamma_ud, u_d,
          u_u, v_d, v_u, proj_ud, W_ijk, Estar, Fstar_d, Estar, Fstar_d,
          volform_ijk * eta_1[i4D], abs_1[i4D], scat_1[i4D]);
      const int source_result = source_update(
          source_ctx, closure_fun, &chi[i4D], &Enew, &Fnew_d,
          source_force_nonlinear, source_thick_limit, source_scat_limit,
          source_maxiter, source_epsabs, source_epsrel);
      if (source_result == NUX_M1_SOURCE_FAIL)
        source_solver_abort();

      assert(isfinite(Enew));
      assert(isfinite(Fnew_d(1)));
      assert(isfinite(Fnew_d(2)));
      assert(isfinite(Fnew_d(3)));
      assert(isfinite(chi[i4D]));
      if (!source_state_is_realizable(g_uu, Enew, Fnew_d))
        source_solver_abort();

      // The nonlinear solve returns the physical collision state. Record
      // that increment before any numerical floor/realizability repair.
      source_update_delta_E[i4D] = Enew - Estar;
      source_update_delta_Fx[i4D] = Fnew_d(1) - Fstar_d(1);
      source_update_delta_Fy[i4D] = Fnew_d(2) - Fstar_d(2);
      source_update_delta_Fz[i4D] = Fnew_d(3) - Fstar_d(3);

      // Use the accepted physical source state consistently for the
      // energy closure and number source. The optional transport-margin
      // repair is applied only when committing the evolved moments.
      apply_closure(g_dd, g_uu, n_d, W_ijk, u_u, v_d, proj_ud, Enew, Fnew_d,
                    chi[i4D], &P_dd);

      tensor::symmetric2<CCTK_REAL, 4, 2> T_dd;
      assemble_rT(n_d, Enew, Fnew_d, P_dd, &T_dd);
      Jnew = calc_J_from_rT(T_dd, u_u);
      assert(isfinite(Jnew));

      CCTK_REAL Gamma =
          compute_Gamma(W_ijk, v_u, Jnew, Enew, Fnew_d, rad_E_floor, rad_eps);
      assert(isfinite(Gamma));

      if (!source_uses_thermalized_number_limit(dtau, abs_0[i4D],
                                                source_therm_limit)) {
        source_update_delta_N[i4D] =
            (Nstar + dt * alp_ijk * volform_ijk * eta_0[i4D]) /
                (1 + dt * alp_ijk * abs_0[i4D] / Gamma) -
            Nstar;
      } else {
        source_update_delta_N[i4D] =
            (nueave[i4D] > 0 ? (Gamma * Jnew) / nueave[i4D] - Nstar : 0.0);
      }
      assert(isfinite(source_update_delta_E[i4D]));
      assert(isfinite(source_update_delta_Fx[i4D]));
      assert(isfinite(source_update_delta_Fy[i4D]));
      assert(isfinite(source_update_delta_Fz[i4D]));
      assert(isfinite(source_update_delta_N[i4D]));
    }
  });
}

void launch_apply_source_kernel(
    GridDescBaseDevice const &grid, GF3D2layout const layout_cc,
    tensor::slicing_geometry_const const geom, CCTK_INT const ngroups,
    CCTK_INT const nspecies, CCTK_REAL const *const nuX_m1_mask,
    CCTK_REAL *const netabs, CCTK_REAL *const netheat, CCTK_REAL *const rN,
    CCTK_REAL *const rE, CCTK_REAL *const rFx, CCTK_REAL *const rFy,
    CCTK_REAL *const rFz, CCTK_REAL *const momx, CCTK_REAL *const momy,
    CCTK_REAL *const momz, CCTK_REAL *const tau, CCTK_REAL *const DYe,
    CCTK_REAL const *const source_limiter_theta,
    CCTK_REAL *const source_update_delta_N,
    CCTK_REAL *const source_update_delta_E,
    CCTK_REAL *const source_update_delta_Fx,
    CCTK_REAL *const source_update_delta_Fy,
    CCTK_REAL *const source_update_delta_Fz, CCTK_REAL const dt,
    CCTK_REAL const mb, CCTK_INT const backreact,
    bool const implicit_backreaction, CCTK_REAL const rad_N_floor,
    CCTK_REAL const rad_E_floor, CCTK_REAL const rad_eps) {
  assert(dt > 0.0);
  const CCTK_REAL inv_dt = 1.0 / dt;
  grid.loop_int_device<1, 1, 1>(grid.nghostzones, [=] CCTK_DEVICE(
                                                      const PointDesc &p) {
    const int ijk = layout_cc.linear(p.i, p.j, p.k);
    netabs[ijk] = 0;
    netheat[ijk] = 0;
    if (nuX_m1_mask[ijk])
      return;

    tensor::generic<CCTK_REAL, 4, 1> beta_u;
    geom.get_shift_vec(p, &beta_u);
    const CCTK_REAL betax_ijk = beta_u(1);
    const CCTK_REAL betay_ijk = beta_u(2);
    const CCTK_REAL betaz_ijk = beta_u(3);
    tensor::inv_metric<4> g_uu;
    geom.get_inv_metric(p, &g_uu);

    int const groupspec = ngroups * nspecies;
    const CCTK_REAL theta = source_limiter_theta[ijk];
    for (int ig = 0; ig < groupspec; ++ig) {
      int const i4D = layout_cc.linear(p.i, p.j, p.k, ig);

      assert(isfinite(theta));
      assert(isfinite(rE[i4D]));
      assert(isfinite(rFx[i4D]));
      assert(isfinite(rFy[i4D]));
      assert(isfinite(rFz[i4D]));
      assert(isfinite(source_update_delta_E[i4D]));
      assert(isfinite(source_update_delta_Fx[i4D]));
      assert(isfinite(source_update_delta_Fy[i4D]));
      assert(isfinite(source_update_delta_Fz[i4D]));
      assert(isfinite(source_update_delta_N[i4D]));

      CCTK_REAL E = rE[i4D];
      tensor::generic<CCTK_REAL, 4, 1> F_d;
      pack_F_d(betax_ijk, betay_ijk, betaz_ijk, rFx[i4D], rFy[i4D], rFz[i4D],
               &F_d);
      apply_floor(g_uu, &E, &F_d, rad_E_floor, rad_eps);

      const CCTK_REAL collision_delta_E = theta * source_update_delta_E[i4D];
      const CCTK_REAL collision_delta_Fx = theta * source_update_delta_Fx[i4D];
      const CCTK_REAL collision_delta_Fy = theta * source_update_delta_Fy[i4D];
      const CCTK_REAL collision_delta_Fz = theta * source_update_delta_Fz[i4D];
      const CCTK_REAL collision_delta_N = theta * source_update_delta_N[i4D];

      E += collision_delta_E;
      const CCTK_REAL Fx_new = F_d(1) + collision_delta_Fx;
      const CCTK_REAL Fy_new = F_d(2) + collision_delta_Fy;
      const CCTK_REAL Fz_new = F_d(3) + collision_delta_Fz;

      assert(isfinite(E));
      assert(isfinite(Fx_new));
      assert(isfinite(Fy_new));
      assert(isfinite(Fz_new));

      pack_F_d(betax_ijk, betay_ijk, betaz_ijk, Fx_new, Fy_new, Fz_new, &F_d);
      apply_floor(g_uu, &E, &F_d, rad_E_floor, rad_eps);

      CCTK_REAL N = max(rN[i4D], rad_N_floor) + collision_delta_N;
      N = max(N, rad_N_floor);

      const CCTK_REAL collision_delta_DYe =
          -mb * ((ig == 0 ? collision_delta_N : 0.0) -
                 (ig == 1 ? collision_delta_N : 0.0));

      if (backreact) {
        netabs[ijk] += collision_delta_DYe;
        netheat[ijk] -= collision_delta_E;
      }

      // ODESolvers recovers the implicit RHS from the accepted diagonal-stage
      // increment. Couple only the physical collision increment here; the
      // radiation floor and realizability repairs below are numerical and
      // must not exchange energy, momentum, or lepton number with matter.
      if (implicit_backreaction) {
        momx[ijk] -= collision_delta_Fx;
        momy[ijk] -= collision_delta_Fy;
        momz[ijk] -= collision_delta_Fz;
        tau[ijk] -= collision_delta_E;
        DYe[ijk] += collision_delta_DYe;

        assert(isfinite(momx[ijk]));
        assert(isfinite(momy[ijk]));
        assert(isfinite(momz[ijk]));
        assert(isfinite(tau[ijk]));
        assert(isfinite(DYe[ijk]));
      }

      // Reuse the stage workspace after committing the update. The legacy
      // semi-implicit backreaction consumes this accepted collision rate.
      source_update_delta_N[i4D] = collision_delta_N * inv_dt;
      source_update_delta_E[i4D] = collision_delta_E * inv_dt;
      source_update_delta_Fx[i4D] = collision_delta_Fx * inv_dt;
      source_update_delta_Fy[i4D] = collision_delta_Fy * inv_dt;
      source_update_delta_Fz[i4D] = collision_delta_Fz * inv_dt;

      rE[i4D] = E;
      unpack_F_d(F_d, &rFx[i4D], &rFy[i4D], &rFz[i4D]);
      rN[i4D] = N;

      assert(isfinite(rN[i4D]));
      assert(isfinite(rE[i4D]));
      assert(isfinite(rFx[i4D]));
      assert(isfinite(rFy[i4D]));
      assert(isfinite(rFz[i4D]));
    }
  });
}

extern "C" void nuX_M1_CalcUpdate(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_M1_CalcUpdate;
  DECLARE_CCTK_PARAMETERS;

  if (verbose) {
    CCTK_INFO("nuX_M1_CalcUpdate");
  }

  closure_t const closure_fun = pick_closure(closure);

  const bool semi_implicit = CCTK_Equals(method, "semi-implicit");
  const bool implicit_backreaction =
      backreact && !semi_implicit &&
      CCTK_Equals(imex_backreaction, "implicit");
  const int semi_implicit_stage_value =
      semi_implicit ? int(*semi_implicit_stage) : 1;
  if (semi_implicit && semi_implicit_stage_value != 1 &&
      semi_implicit_stage_value != 2)
    CCTK_VERROR("Invalid nuX semi-implicit stage %d",
                semi_implicit_stage_value);
  CCTK_REAL const dt =
      CCTK_DELTA_TIME / CCTK_REAL(semi_implicit_stage_value);

  if (verbose) {
    CCTK_VINFO("Applying implicit source step at time %e with dt %e",
               cctkGH->cctk_time, dt);
  }

  const GridDescBaseDevice grid(cctkGH);
  const GF3D2layout layout_cc(cctkGH, {1, 1, 1});
  const GF3D2layout layout_vc(cctkGH, {0, 0, 0});

  tensor::slicing_geometry_const geom(layout_vc, layout_cc, alp, betax, betay,
                                      betaz, gxx, gxy, gxz, gyy, gyz, gzz, kxx,
                                      kxy, kxz, kyy, kyz, kzz);
  tensor::fluid_velocity_field_const fidu(layout_vc, layout_cc, alp, betax,
                                          betay, betaz, fidu_w_lorentz,
                                          fidu_velx, fidu_vely, fidu_velz);

  // particle_mass is in MeV
  CCTK_REAL const mb = nuX_Utils::AverageBaryonMass(particle_mass);
  const bool limit_backreaction = backreact;

  if (verbose) {
    CCTK_INFO("nuX_M1_CalcUpdate 1");
  }

  // Step 1: compute source increments into grid functions. This keeps the
  // per-thread state for limiter/application kernels small. The analytic
  // thick-limit state constructed below also seeds the nonlinear solve when
  // source_force_nonlinear is enabled.
  launch_source_increment_kernel(
      grid, layout_cc, geom, fidu, cctkGH, closure_fun, dt, ngroups, nspecies,
      nuX_m1_mask, fidu_w_lorentz, fidu_velx, fidu_vely, fidu_velz, rN, rE,
      rFx, rFy, rFz, chi, eta_0, eta_1, abs_0, abs_1, scat_1, nueave,
      source_update_delta_N, source_update_delta_E, source_update_delta_Fx,
      source_update_delta_Fy, source_update_delta_Fz, rad_N_floor, rad_E_floor,
      rad_eps, closure_epsilon, closure_maxiter, use_fallback,
      source_force_nonlinear, source_thick_limit, source_scat_limit,
      source_therm_limit, source_maxiter, source_epsabs, source_epsrel);

  // Step 2: compute the per-cell source limiter from the scratch increments.
  grid.loop_int_device<1, 1, 1>(
      grid.nghostzones, [=] CCTK_DEVICE(const PointDesc &p) {
        const int ijk = layout_cc.linear(p.i, p.j, p.k);
        int const groupspec = ngroups * nspecies;
        if (nuX_m1_mask[ijk]) {
          source_limiter_theta[ijk] = 0.0;
          return;
        }

        tensor::generic<CCTK_REAL, 4, 1> beta_u;
        geom.get_shift_vec(p, &beta_u);
        const CCTK_REAL betax_ijk = beta_u(1);
        const CCTK_REAL betay_ijk = beta_u(2);
        const CCTK_REAL betaz_ijk = beta_u(3);
        tensor::inv_metric<4> g_uu;
        geom.get_inv_metric(p, &g_uu);

        CCTK_REAL theta = 1.0;
        if (source_limiter >= 0) {
          CCTK_REAL delta_tau = 0.0;
          CCTK_REAL delta_DYe = 0.0;
          for (int ig = 0; ig < groupspec; ++ig) {
            int const i4D = layout_cc.linear(p.i, p.j, p.k, ig);

            CCTK_REAL Estar = rE[i4D];
            tensor::generic<CCTK_REAL, 4, 1> Fstar_d;
            pack_F_d(betax_ijk, betay_ijk, betaz_ijk, rFx[i4D], rFy[i4D],
                     rFz[i4D], &Fstar_d);
            apply_floor(g_uu, &Estar, &Fstar_d, rad_E_floor, rad_eps);
            if (source_update_delta_E[i4D] < 0) {
              theta = min(theta, -source_limiter * max(Estar, 0.0) /
                                     source_update_delta_E[i4D]);
            }

            CCTK_REAL Nstar = max(rN[i4D], rad_N_floor);
            if (source_update_delta_N[i4D] < 0) {
              theta = min(theta, -source_limiter * max(Nstar, 0.0) /
                                     source_update_delta_N[i4D]);
            }

            if (limit_backreaction) {
              delta_tau -= source_update_delta_E[i4D];
              delta_DYe +=
                  -mb * ((ig == 0 ? source_update_delta_N[i4D] : 0.0) -
                         (ig == 1 ? source_update_delta_N[i4D] : 0.0));
            }
          }

          if (limit_backreaction) {
            if (delta_tau < 0.0)
              theta = min(theta, -source_limiter * max(tau[ijk], 0.0) /
                                     delta_tau);

            if (dens[ijk] > 0.0) {
              const CCTK_REAL delta_Ye = delta_DYe / dens[ijk];
              if (delta_Ye > 0.0)
                theta = min(theta,
                            source_limiter *
                                max(source_Ye_max - Ye[ijk], 0.0) / delta_Ye);
              else if (delta_Ye < 0.0)
                theta = min(theta,
                            source_limiter *
                                min(source_Ye_min - Ye[ijk], 0.0) / delta_Ye);
            }
          }
          theta = max(CCTK_REAL(0), theta);
        }
        source_limiter_theta[ijk] = theta;
      });

  // Step 3: apply the limited update.
  launch_apply_source_kernel(
      grid, layout_cc, geom, ngroups, nspecies, nuX_m1_mask, netabs, netheat,
      rN, rE, rFx, rFy, rFz, momx, momy, momz, tau, DYe,
      source_limiter_theta, source_update_delta_N, source_update_delta_E,
      source_update_delta_Fx, source_update_delta_Fy, source_update_delta_Fz,
      dt, mb, backreact, implicit_backreaction, rad_N_floor, rad_E_floor,
      rad_eps);
}

} // namespace nuX_M1
