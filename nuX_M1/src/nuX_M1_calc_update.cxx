#include <cassert>
#include <algorithm>
#include <cmath>
#include <cstdio>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "nuX_M1_macro.hxx"
#include "nuX_M1_closure.hxx"
#include "nuX_utils.hxx"
#include "nuX_M1_sources.hxx"
#include "avg_baryon_mass.hpp"

namespace nuX_M1 {

using namespace std;
using namespace Loop;
using namespace nuX_Utils;

// Fallback if header didn't set it
#ifndef MAX_GROUPSPECIES
#define MAX_GROUPSPECIES 3
#endif

// Map closure string → function pointer
static inline closure_t pick_closure_fun(const char *name) {
  if (CCTK_Equals(name, "Eddington"))
    return eddington;
  if (CCTK_Equals(name, "Kershaw"))
    return kershaw;
  if (CCTK_Equals(name, "Minerbo"))
    return minerbo;
  if (CCTK_Equals(name, "thin"))
    return thin;
  char msg[BUFSIZ];
  snprintf(msg, BUFSIZ, "Unknown closure \"%s\"", name);
  CCTK_ERROR(msg);
  return eddington; // not reached
}

extern "C" void nuX_M1_CalcUpdate(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_M1_CalcUpdate;
  DECLARE_CCTK_PARAMETERS;

  if (verbose) {
    CCTK_INFO("nuX_M1_CalcUpdate");
  }

  closure_t const closure_fun = pick_closure_fun(closure);

  CCTK_REAL const dt =
      CCTK_DELTA_TIME / static_cast<CCTK_REAL>(*TimeIntegratorStage);
  const bool apply_backreact = backreact && (1 == *TimeIntegratorStage);

  if (verbose) {
    CCTK_VINFO("Integrated to time, dt, TimeIntegratorStage: %e, %e, %e",
               cctkGH->cctk_time, dt,
               static_cast<CCTK_REAL>(*TimeIntegratorStage));
  }

  const GridDescBaseDevice grid(cctkGH);
  const GF3D2layout layout_cc(cctkGH, {1, 1, 1});
  const GF3D2layout layout_vc(cctkGH, {0, 0, 0});

  // Geometry & velocity field accessors
  tensor::slicing_geometry_const geom(layout_vc, layout_cc, alp, betax, betay,
                                      betaz, gxx, gxy, gxz, gyy, gyz, gzz, kxx,
                                      kxy, kxz, kyy, kyz, kzz);
  tensor::fluid_velocity_field_const fidu(layout_vc, layout_cc, alp, betax,
                                          betay, betaz, fidu_w_lorentz,
                                          fidu_velx, fidu_vely, fidu_velz);

  // particle_mass is in MeV
  CCTK_REAL const mb = AverageBaryonMass(particle_mass);

  if (verbose) {
    CCTK_INFO("nuX_M1_CalcUpdate 1");
  }

  // Steps
  // 1. F^m   = F^k + dt/2 [ A[F^k] + S[F^m]   ]
  // 2. F^k+1 = F^k + dt   [ A[F^m] + S[F^k+1] ]
  // At each step we solve an implicit problem in the form
  //    F = F^* + cdt S[F]
  // Where F^* = F^k + cdt A
  grid.loop_int_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        const int ijk = layout_cc.linear(p.i, p.j, p.k);
        netabs[ijk] = 0;
        netheat[ijk] = 0;
        if (nuX_m1_mask[ijk]) {
          return;
        }
        tensor::generic<CCTK_REAL, 4, 1> beta_u;
        geom.get_shift_vec(p, &beta_u);
        const CCTK_REAL alp_ijk = geom.get_lapse(p);
        const CCTK_REAL betax_ijk = beta_u(1);
        const CCTK_REAL betay_ijk = beta_u(2);
        const CCTK_REAL betaz_ijk = beta_u(3);
        const CCTK_REAL W_ijk = fidu_w_lorentz[ijk];
        const CCTK_REAL dens_ijk = dens[ijk];
        // Metric, normal, projectors
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

        // Fluid 4-velocity, projector
        tensor::generic<CCTK_REAL, 4, 1> u_u, u_d;
        tensor::generic<CCTK_REAL, 4, 2> proj_ud;
        fidu.get(p, &u_u);
        tensor::contract(g_dd, u_u, &u_d);
        calc_proj(u_d, u_u, &proj_ud);

        // Fiducial 3-velocity
        tensor::generic<CCTK_REAL, 4, 1> v_u, v_d;
        pack_v_u(fidu_velx[ijk], fidu_vely[ijk], fidu_velz[ijk], &v_u);
        tensor::contract(g_dd, v_u, &v_d);

        // Per-(group×species) accumulators
        int const groupspec = ngroups * nspecies;

        // Source RHS are stored here
        // TODO: to use ngroups * nspecies instead of MAX_GROUPSPECIES
        CCTK_REAL DrE[MAX_GROUPSPECIES];
        CCTK_REAL DrFx[MAX_GROUPSPECIES];
        CCTK_REAL DrFy[MAX_GROUPSPECIES];
        CCTK_REAL DrFz[MAX_GROUPSPECIES];
        CCTK_REAL DrN[MAX_GROUPSPECIES];
        // --------------------------
        // Step 1 — compute sources
        // --------------------------
        for (int ig = 0; ig < groupspec; ++ig) {
          int const i4D = layout_cc.linear(p.i, p.j, p.k, ig);
          assert(isfinite(rN[i4D]));
          assert(isfinite(rE[i4D]));
          assert(isfinite(rFx[i4D]));
          assert(isfinite(rFy[i4D]));
          assert(isfinite(rFz[i4D]));
          assert(isfinite(rN_p[i4D]));
          assert(isfinite(rE_p[i4D]));
          assert(isfinite(rFx_p[i4D]));
          assert(isfinite(rFy_p[i4D]));
          assert(isfinite(rFz_p[i4D]));
#if (NUX_M1_SRC_METHOD == NUX_M1_SRC_EXPL)

          // Explicit sources (THC EXPL)

          // Radiation fields
          CCTK_REAL E = rE[i4D];
          tensor::generic<CCTK_REAL, 4, 1> F_d;
          pack_F_d(betax_ijk, betay_ijk, betaz_ijk, rFx[i4D], rFy[i4D],
                   rFz[i4D], &F_d);
          tensor::generic<CCTK_REAL, 4, 1> F_u, S_d, tS_d;
          tensor::contract(g_uu, F_d, &F_u);

          // Compute radiation quantities in the fluid frame
          CCTK_REAL J = rJ[i4D];
          CCTK_REAL Gamma =
              compute_Gamma(W_ijk, v_u, J, E, F_d, rad_E_floor, rad_eps);
          tensor::generic<CCTK_REAL, 4, 1> H_d;
          pack_H_d(rHt[i4D], rHx[i4D], rHy[i4D], rHz[i4D], &H_d);

          // Compute radiation sources
          calc_rad_sources(eta_1[i4D] * volform_ijk, abs_1[i4D], scat_1[i4D],
                           u_d, J, H_d, &S_d);
          DrE[ig] = dt * calc_rE_source(alp_ijk, n_u, S_d);
          calc_rF_source(alp_ijk, gamma_ud, S_d, &tS_d);
          DrFx[ig] = dt * tS_d(1);
          DrFy[ig] = dt * tS_d(2);
          DrFz[ig] = dt * tS_d(3);
          DrN[ig] = dt * alp_ijk *
                    (volform_ijk * eta_0[i4D] - abs_0[i4D] * rN[i4D] / Gamma);
          assert(isfinite(DrE[ig]));
          assert(isfinite(DrFx[ig]));
          assert(isfinite(DrFy[ig]));
          assert(isfinite(DrFz[ig]));
          assert(isfinite(DrN[ig]));

#else // NUX_M1_SRC_METHOD != NUX_M1_SRC_EXPL

          // Predictor (advect radiation)
          CCTK_REAL Estar = rE_p[i4D] + dt * rE_rhs[i4D];
          tensor::generic<CCTK_REAL, 4, 1> Fstar_d;
          pack_F_d(betax_ijk, betay_ijk, betaz_ijk,
                   rFx_p[i4D] + dt * rFx_rhs[i4D],
                   rFy_p[i4D] + dt * rFy_rhs[i4D],
                   rFz_p[i4D] + dt * rFz_rhs[i4D], &Fstar_d);
          apply_floor(g_uu, &Estar, &Fstar_d, rad_E_floor, rad_eps);
          CCTK_REAL Nstar = std::max(rN_p[i4D] + dt * rN_rhs[i4D], rad_N_floor);
          CCTK_REAL Enew = Estar;
          tensor::generic<CCTK_REAL, 4, 1> Fnew_d;

          tensor::symmetric2<CCTK_REAL, 4, 2> P_dd;
          calc_closure(cctkGH, p.i, p.j, p.k, ig, closure_fun, g_dd, g_uu, n_d,
                       W_ijk, u_u, v_d, proj_ud, Estar, Fstar_d, &chi[i4D],
                       &P_dd, closure_epsilon, closure_maxiter);

          tensor::symmetric2<CCTK_REAL, 4, 2> rT_dd;
          assemble_rT(n_d, Estar, Fstar_d, P_dd, &rT_dd);

          CCTK_REAL const Jstar = calc_J_from_rT(rT_dd, u_u);
          tensor::generic<CCTK_REAL, 4, 1> Hstar_d;
          calc_H_from_rT(rT_dd, u_u, proj_ud, &Hstar_d);

          CCTK_REAL const dtau = dt / W_ijk;
          CCTK_REAL Jnew = (Jstar + dtau * eta_1[i4D] * volform_ijk) /
                           (1 + dtau * abs_1[i4D]);

          CCTK_REAL const khat = (abs_1[i4D] + scat_1[i4D]);
          tensor::generic<CCTK_REAL, 4, 1> Hnew_d;
          for (int a = 1; a < 4; ++a)
            Hnew_d(a) = Hstar_d(a) / (1 + dtau * khat);
          Hnew_d(0) = 0.0;
          for (int a = 1; a < 4; ++a)
            Hnew_d(0) -= Hnew_d(a) * (u_u(a) / u_u(0));

          CCTK_REAL const H2 = tensor::dot(g_uu, Hnew_d, Hnew_d);

#if (NUX_M1_SRC_METHOD == NUX_M1_SRC_BOOST)
          CCTK_REAL const xi = (Jnew > rad_E_floor ? sqrt(H2) / Jnew : 0.0);
          chi[i4D] = closure_fun(xi);
#else
          chi[i4D] = (CCTK_REAL)(1.0 / 3.0);
#endif

          CCTK_REAL const dthick = 3.0 * (1.0 - chi[i4D]) / 2.0;
          CCTK_REAL const dthin = 1.0 - dthick;

          for (int a = 0; a < 4; ++a)
            for (int b = a; b < 4; ++b) {
              rT_dd(a, b) =
                  Jnew * u_d(a) * u_d(b) + Hnew_d(a) * u_d(b) +
                  Hnew_d(b) * u_d(a) +
                  dthin * Jnew * (H2 > 0 ? Hnew_d(a) * Hnew_d(b) / H2 : 0.0) +
                  dthick * Jnew * (g_dd(a, b) + u_d(a) * u_d(b)) / 3.0;
            }

          Enew = calc_J_from_rT(rT_dd, n_u);
          calc_H_from_rT(rT_dd, n_u, gamma_ud, &Fnew_d);
          apply_floor(g_uu, &Enew, &Fnew_d, rad_E_floor, rad_eps);

#if (NUX_M1_SRC_METHOD == NUX_M1_SRC_IMPL)
          (void)source_update(
              cctkGH, p.i, p.j, p.k, ig, closure_fun, closure_epsilon,
              closure_maxiter, dt, alp_ijk, g_dd, g_uu, n_d, n_u, gamma_ud, u_d,
              u_u, v_d, v_u, proj_ud, W_ijk, Estar, Fstar_d, Estar, Fstar_d,
              volform_ijk * eta_1[i4D], abs_1[i4D], scat_1[i4D], &chi[i4D],
              &Enew, &Fnew_d, source_thick_limit, source_scat_limit,
              source_maxiter, source_epsabs, source_epsrel);

          apply_floor(g_uu, &Enew, &Fnew_d, rad_E_floor, rad_eps);

          apply_closure(g_dd, g_uu, n_d, W_ijk, u_u, v_d, proj_ud, Enew, Fnew_d,
                        chi[i4D], &P_dd);

          tensor::symmetric2<CCTK_REAL, 4, 2> T_dd;
          assemble_rT(n_d, Enew, Fnew_d, P_dd, &T_dd);
          Jnew = calc_J_from_rT(T_dd, u_u);
#endif

          DrE[ig] = Enew - Estar;
          DrFx[ig] = Fnew_d(1) - Fstar_d(1);
          DrFy[ig] = Fnew_d(2) - Fstar_d(2);
          DrFz[ig] = Fnew_d(3) - Fstar_d(3);

          CCTK_REAL Gamma = compute_Gamma(W_ijk, v_u, Jnew, Enew, Fnew_d,
                                          rad_E_floor, rad_eps);

          if (source_therm_limit < 0 || dt * abs_0[i4D] < source_therm_limit) {
            DrN[ig] = (Nstar + dt * alp_ijk * volform_ijk * eta_0[i4D]) /
                          (1 + dt * alp_ijk * abs_0[i4D] / Gamma) -
                      Nstar;
          } else {
            DrN[ig] =
                (nueave[i4D] > 0 ? (Gamma * Jnew) / nueave[i4D] - Nstar : 0.0);
          }

#endif // NUX_M1_SRC_METHOD
        } // ig

        // --------------------------
        // Step 2 — source limiter
        // --------------------------
        CCTK_REAL theta = 1.0;
        if (source_limiter >= 0) {
          CCTK_REAL DTau_sum = 0.0;
          CCTK_REAL DDxp_sum = 0.0;
          for (int ig = 0; ig < groupspec; ++ig) {
            int const i4D = layout_cc.linear(p.i, p.j, p.k, ig);

            CCTK_REAL Estar = rE_p[i4D] + dt * rE_rhs[i4D];
            if (DrE[ig] < 0) {
              theta = min(theta, -source_limiter * max(Estar, 0.0) / DrE[ig]);
            }
            DTau_sum -= DrE[ig];

            CCTK_REAL Nstar = rN_p[i4D] + dt * rN_rhs[i4D];
            if (DrN[ig] < 0) {
              theta = min(theta, -source_limiter * max(Nstar, 0.0) / DrN[ig]);
            }
            const CCTK_REAL DDxp_ig =
                -mb * ((ig == 0 ? DrN[ig] : 0.0) - (ig == 1 ? DrN[ig] : 0.0));
            DDxp_sum += DDxp_ig;
          }
          CCTK_REAL const DYe = DDxp_sum / dens_ijk;
          if (DTau_sum < 0) {
            theta = min(theta, -source_limiter * max(tau[ijk], 0.0) / DTau_sum);
          }
          if (DYe > 0) {
            theta = min(theta, source_limiter *
                                   max(source_Ye_max - Ye[ijk], 0.0) / DYe);
          } else if (DYe < 0) {
            theta = min(theta, source_limiter *
                                   min(source_Ye_min - Ye[ijk], 0.0) / DYe);
          }
          theta = max((CCTK_REAL)0.0, theta);
        }

        // --------------------------
        // Step 3 — apply updates
        // --------------------------
        for (int ig = 0; ig < groupspec; ++ig) {
          int const i4D = layout_cc.linear(p.i, p.j, p.k, ig);

          CCTK_REAL E = rE_p[i4D] + dt * rE_rhs[i4D] + theta * DrE[ig];
          const CCTK_REAL Fx_new =
              rFx_p[i4D] + dt * rFx_rhs[i4D] + theta * DrFx[ig];
          const CCTK_REAL Fy_new =
              rFy_p[i4D] + dt * rFy_rhs[i4D] + theta * DrFy[ig];
          const CCTK_REAL Fz_new =
              rFz_p[i4D] + dt * rFz_rhs[i4D] + theta * DrFz[ig];
          tensor::generic<CCTK_REAL, 4, 1> F_d;
          pack_F_d(betax_ijk, betay_ijk, betaz_ijk, Fx_new, Fy_new, Fz_new,
                   &F_d);
          apply_floor(g_uu, &E, &F_d, rad_E_floor, rad_eps);

          CCTK_REAL N = rN_p[i4D] + dt * rN_rhs[i4D] + theta * DrN[ig];
          N = max(N, rad_N_floor);
          const CCTK_REAL DDxp_ig =
              -mb * ((ig == 0 ? DrN[ig] : 0.0) - (ig == 1 ? DrN[ig] : 0.0));

          if (apply_backreact) {
            assert(ngroups == 1);
            assert(nspecies == 3);
            momx[ijk] -= theta * DrFx[ig];
            momy[ijk] -= theta * DrFy[ig];
            momz[ijk] -= theta * DrFz[ig];
            tau[ijk] -= theta * DrE[ig];
            DYe[ijk] += theta * DDxp_ig;
            netabs[ijk] += theta * DDxp_ig;
            netheat[ijk] -= theta * DrE[ig];
          }

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

} // namespace nuX_M1
