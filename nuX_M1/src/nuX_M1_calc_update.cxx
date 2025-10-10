#include <cassert>
#include <algorithm>
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

  // Stage dt
  CCTK_REAL const dt =
      CCTK_DELTA_TIME / static_cast<CCTK_REAL>(*TimeIntegratorStage);

  if (verbose) {
    CCTK_VINFO("Integrated to time, dt, TimeIntegratorStage: %e, %e, %e", cctkGH->cctk_time, dt, static_cast<CCTK_REAL>(*TimeIntegratorStage));
  }
   
  --(*TimeIntegratorStage);
  // CCTK_REAL const dt = CCTK_DELTA_TIME;

  // Geometry & velocity field accessors
  tensor::slicing_geometry_const geom(alp, betax, betay, betaz, gxx, gxy, gxz,
                                      gyy, gyz, gzz, kxx, kxy, kxz, kyy, kyz,
                                      kzz, volform);
  tensor::fluid_velocity_field_const fidu(alp, betax, betay, betaz,
                                          fidu_w_lorentz, fidu_velx, fidu_vely,
                                          fidu_velz);

  // particle_mass is in MeV
  CCTK_REAL const mb = AverageBaryonMass(particle_mass);

  const GridDescBaseDevice grid(cctkGH);
  const GF3D2layout layout2(cctkGH, {1, 1, 1});

  if (verbose) { CCTK_INFO("nuX_M1_CalcUpdate 1");}

  // Steps
  // 1. F^m   = F^k + dt/2 [ A[F^k] + S[F^m]   ]
  // 2. F^k+1 = F^k + dt   [ A[F^m] + S[F^k+1] ]
  // At each step we solve an implicit problem in the form
  //    F = F^* + cdt S[F]
  // Where F^* = F^k + cdt A
  grid.loop_int_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        const int ijk = layout2.linear(p.i, p.j, p.k);
        netabs[ijk] = 0;
        netheat[ijk] = 0;
        if (nuX_m1_mask[ijk]) {
          return;
        }
        if (verbose) { CCTK_INFO("nuX_M1_CalcUpdate 2");}
        // Metric, normal, projectors
        tensor::metric<4> g_dd;
        tensor::inv_metric<4> g_uu;
        tensor::generic<CCTK_REAL, 4, 1> n_u, n_d;
        tensor::generic<CCTK_REAL, 4, 2> gamma_ud;
        geom.get_metric(ijk, &g_dd);
        geom.get_inv_metric(ijk, &g_uu);
        geom.get_normal(ijk, &n_u);
        geom.get_normal_form(ijk, &n_d);
        geom.get_space_proj(ijk, &gamma_ud);

        // Fluid 4-velocity, projector
        tensor::generic<CCTK_REAL, 4, 1> u_u, u_d;
        tensor::generic<CCTK_REAL, 4, 2> proj_ud;
        fidu.get(ijk, &u_u);
        tensor::contract(g_dd, u_u, &u_d);
        calc_proj(u_d, u_u, &proj_ud);

        // Fiducial 3-velocity
        tensor::generic<CCTK_REAL, 4, 1> v_u, v_d;
        pack_v_u(fidu_velx[ijk], fidu_vely[ijk], fidu_velz[ijk], &v_u);
        tensor::contract(g_dd, v_u, &v_d);

        // Work vectors
        tensor::generic<CCTK_REAL, 4, 1> F_d, Fstar_d, Fnew_d;

        // Per-(group×species) accumulators
        int const groupspec = ngroups * nspecies;

        // Source RHS are stored here
        // TODO: to use ngroups * nspecies instead of MAX_GROUPSPECIES
        CCTK_REAL DrE[MAX_GROUPSPECIES];
        CCTK_REAL DrFx[MAX_GROUPSPECIES];
        CCTK_REAL DrFy[MAX_GROUPSPECIES];
        CCTK_REAL DrFz[MAX_GROUPSPECIES];
        CCTK_REAL DrN[MAX_GROUPSPECIES];
        CCTK_REAL DDxp[MAX_GROUPSPECIES];
        if (verbose) { CCTK_INFO("nuX_M1_CalcUpdate 3");}
        // --------------------------
        // Step 1 — compute sources
        // --------------------------
        for (int ig = 0; ig < groupspec; ++ig) {
          int const i4D = layout2.linear(p.i, p.j, p.k, ig);
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
          pack_F_d(betax[ijk], betay[ijk], betaz[ijk], rFx[i4D], rFy[i4D],
                   rFz[i4D], &F_d);
          tensor::generic<CCTK_REAL, 4, 1> F_u, S_d, tS_d;
          tensor::contract(g_uu, F_d, &F_u);

          // Compute radiation quantities in the fluid frame
          CCTK_REAL J = rJ[i4D];
          CCTK_REAL const Gamma = compute_Gamma(fidu_w_lorentz[ijk], v_u, J, E,
                                                F_d, rad_E_floor, rad_eps);
	  if (!isfinite(Gamma) || Gamma < 1.0) Gamma = 1.0; // JayMOD
	  tensor::generic<CCTK_REAL, 4, 1> H_d;
          pack_H_d(rHt[i4D], rHx[i4D], rHy[i4D], rHz[i4D], &H_d);

          // Compute radiation sources
          calc_rad_sources(eta_1[i4D] * volform[ijk], abs_1[i4D], scat_1[i4D],
                           u_d, J, H_d, &S_d);
          DrE[ig] = dt * calc_rE_source(alp[ijk], n_u, S_d);
          calc_rF_source(alp[ijk], gamma_ud, S_d, &tS_d);
          DrFx[ig] = dt * tS_d(1);
          DrFy[ig] = dt * tS_d(2);
          DrFz[ig] = dt * tS_d(3);
          DrN[ig] = dt * alp[ijk] *
                    (volform[ijk] * eta_0[i4D] - abs_0[i4D] * rN[i4D] / Gamma);
          if (verbose) { CCTK_INFO("nuX_M1_CalcUpdate 4");}
	  assert(isfinite(DrE[ig]));
	  assert(isfinite(DrFx[ig]));
	  assert(isfinite(DrFy[ig]));
	  assert(isfinite(DrFz[ig]));
	  assert(isfinite(DrN[ig]));

#else // NUX_M1_SRC_METHOD != NUX_M1_SRC_EXPL

          // Here we boost to the fluid frame, compute fluid matter
          // interaction, and boost back. These values are used as
          // initial guess for the implicit solve.

          // Predictor (advect radiation)
          CCTK_REAL Estar = rE_p[i4D] + dt * rE_rhs[i4D];
          pack_F_d(betax[ijk], betay[ijk], betaz[ijk],
                   rFx_p[i4D] + dt * rFx_rhs[i4D],
                   rFy_p[i4D] + dt * rFy_rhs[i4D],
                   rFz_p[i4D] + dt * rFz_rhs[i4D], &Fstar_d);
          apply_floor(g_uu, &Estar, &Fstar_d, rad_E_floor, rad_eps);
          CCTK_REAL Nstar = std::max(rN_p[i4D] + dt * rN_rhs[i4D], rad_N_floor);
          CCTK_REAL Enew = Estar;

          // Closure (GPU-safe; gsl pointer ignored)
          // Compute quantities in the fluid frame
          tensor::symmetric2<CCTK_REAL, 4, 2> P_dd;
          calc_closure(cctkGH, p.i, p.j, p.k, ig, closure_fun, g_dd, g_uu, n_d,
                       fidu_w_lorentz[ijk], u_u, v_d, proj_ud, Estar, Fstar_d,
                       &chi[i4D], &P_dd, closure_epsilon, closure_maxiter);

          // Build T^{μν} in normal frame
          tensor::symmetric2<CCTK_REAL, 4, 2> rT_dd;
          assemble_rT(n_d, Estar, Fstar_d, P_dd, &rT_dd);

          CCTK_REAL const Jstar = calc_J_from_rT(rT_dd, u_u);
          tensor::generic<CCTK_REAL, 4, 1> Hstar_d;
          calc_H_from_rT(rT_dd, u_u, proj_ud, &Hstar_d);

          // Matter interaction estimate (fluid frame)
          CCTK_REAL const dtau = dt / fidu_w_lorentz[ijk];
          CCTK_REAL Jnew = (Jstar + dtau * eta_1[i4D] * volform[ijk]) /
                           (1 + dtau * abs_1[i4D]);

          // Only three components of H^a are independent H^0 is found by
          // requiring H^a u_a = 0
          CCTK_REAL const khat = (abs_1[i4D] + scat_1[i4D]);
          tensor::generic<CCTK_REAL, 4, 1> Hnew_d;
          for (int a = 1; a < 4; ++a)
            Hnew_d(a) = Hstar_d(a) / (1 + dtau * khat);
          Hnew_d(0) = 0.0;
          for (int a = 1; a < 4; ++a)
            Hnew_d(0) -= Hnew_d(a) * (u_u(a) / u_u(0));

          // Update T^{μν} pieces
          CCTK_REAL const H2 = tensor::dot(g_uu, Hnew_d, Hnew_d);

// TODO: Boost library is not GPU compatible, so the first condition is never
// true. Hence we have the second "(NUX_M1_SRC_METHOD == NUX_M1_SRC_IMPL)"
// condition
#if (NUX_M1_SRC_METHOD == NUX_M1_SRC_BOOST) ||                                 \
    (NUX_M1_SRC_METHOD == NUX_M1_SRC_IMPL)
          // BOOST (and IMPlicit fallback here): compute xi and use chosen
          // closure
          CCTK_REAL const xi = (Jnew > rad_E_floor ? sqrt(H2) / Jnew : 0.0);
          chi[i4D] = closure_fun(xi);
      // chi[i4D] = closure_fun ? closure_fun(xi) : (CCTK_REAL)(1.0 / 3.0);
#else
          // Thick-limit default
          chi[i4D] = (CCTK_REAL)(1.0 / 3.0);
#endif

          CCTK_REAL const dthick = 3.0 * (1.0 - chi[i4D]) / 2.0;
          CCTK_REAL const dthin = 1.0 - dthick;
          if (verbose) { CCTK_INFO("nuX_M1_CalcUpdate 5");}
          for (int a = 0; a < 4; ++a)
            for (int b = a; b < 4; ++b) {
              rT_dd(a, b) =
                  Jnew * u_d(a) * u_d(b) + Hnew_d(a) * u_d(b) +
                  Hnew_d(b) * u_d(a) +
                  dthin * Jnew * (H2 > 0 ? Hnew_d(a) * Hnew_d(b) / H2 : 0.0) +
                  dthick * Jnew * (g_dd(a, b) + u_d(a) * u_d(b)) / 3.0;
            }

          // Boost back to lab frame
          Enew = calc_J_from_rT(rT_dd, n_u);
          calc_H_from_rT(rT_dd, n_u, gamma_ud, &Fnew_d);
          apply_floor(g_uu, &Enew, &Fnew_d, rad_E_floor, rad_eps);

#if (NUX_M1_SRC_METHOD == NUX_M1_SRC_IMPL)
          // Compute interaction with matter
          (void)source_update(
              cctkGH, p.i, p.j, p.k, ig, closure_fun, closure_epsilon,
              closure_maxiter, dt, alp[ijk], g_dd, g_uu, n_d, n_u, gamma_ud,
              u_d, u_u, v_d, v_u, proj_ud, fidu_w_lorentz[ijk], Estar, Fstar_d,
              Estar, Fstar_d, volform[ijk] * eta_1[i4D], abs_1[i4D],
              scat_1[i4D], &chi[i4D], &Enew, &Fnew_d, source_thick_limit,
              source_scat_limit, source_maxiter);
          apply_floor(g_uu, &Enew, &Fnew_d, rad_E_floor, rad_eps);

          // Update closure
          apply_closure(g_dd, g_uu, n_d, fidu_w_lorentz[ijk], u_u, v_d, proj_ud,
                        Enew, Fnew_d, chi[i4D], &P_dd);
          if (verbose) { CCTK_INFO("nuX_M1_CalcUpdate 6");}
          // Compute new radiation energy density in the fluid frame
          tensor::symmetric2<CCTK_REAL, 4, 2> T_dd;
          assemble_rT(n_d, Enew, Fnew_d, P_dd, &T_dd);
          Jnew = calc_J_from_rT(T_dd, u_u);
#endif

          // Changes (predictor → updated)
          DrE[ig] = Enew - Estar;
          DrFx[ig] = Fnew_d(1) - Fstar_d(1);
          DrFy[ig] = Fnew_d(2) - Fstar_d(2);
          DrFz[ig] = Fnew_d(3) - Fstar_d(3);

          // Updated Gamma (lab ↔ fluid relation)
          (void)compute_Gamma(fidu_w_lorentz[ijk], v_u, Jnew, Enew, Fnew_d,
                              rad_E_floor, rad_eps);

          // Number density update
          // N^k+1 = N^* + dt ( eta - abs N^k+1 )
          if (source_therm_limit < 0 || dt * abs_0[i4D] < source_therm_limit) {
            DrN[ig] =
                (Nstar + dt * alp[ijk] * volform[ijk] * eta_0[i4D]) /
                    (1 + dt * alp[ijk] * abs_0[i4D] /
                             (fidu_w_lorentz[ijk] > 0 ? fidu_w_lorentz[ijk]
                                                      : 1.0)) -
                Nstar;
            // The neutrino number density is updated assuming the neutrino
            // average energies are those of the equilibrium
          } else {
            DrN[ig] = (nueave[i4D] > 0
                           ? (fidu_w_lorentz[ijk] * Jnew) / nueave[i4D] - Nstar
                           : 0.0);
          }
#endif // NUX_M1_SRC_METHOD
           
          // Fluid lepton sources (ν_e – \barν_e)
          DDxp[ig] = -mb * (DrN[ig] * (ig == 0) - DrN[ig] * (ig == 1));
        } // ig
        if (verbose) { CCTK_INFO("nuX_M1_CalcUpdate 7");}
        // --------------------------
        // Step 2 — source limiter
        // --------------------------
        CCTK_REAL theta = 1.0;
        if (source_limiter >= 0) {
          CCTK_REAL DTau_sum = 0.0;
          CCTK_REAL DDxp_sum = 0.0;
          for (int ig = 0; ig < groupspec; ++ig) {
            int const i4D = layout2.linear(p.i, p.j, p.k, ig);

            CCTK_REAL Estar = rE_p[i4D] + dt * rE_rhs[i4D];
            if (DrE[ig] < 0) {
              theta = min(theta, -source_limiter * max(Estar, 0.0) / DrE[ig]);
            }
            DTau_sum -= DrE[ig];

            CCTK_REAL Nstar = rN_p[i4D] + dt * rN_rhs[i4D];
            if (DrN[ig] < 0) {
              theta = min(theta, -source_limiter * max(Nstar, 0.0) / DrN[ig]);
            }
            DDxp_sum += DDxp[ig];
          }
          CCTK_REAL const DYe = DDxp_sum / dens[ijk];
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
        if (verbose) { CCTK_INFO("nuX_M1_CalcUpdate 8");}
        // --------------------------
        // Step 3 — apply updates
        // --------------------------
        for (int ig = 0; ig < groupspec; ++ig) {
          int const i4D = layout2.linear(p.i, p.j, p.k, ig);

	  assert(isfinite(DrE[ig]));
          assert(isfinite(DrFx[ig]));
          assert(isfinite(DrFy[ig]));
          assert(isfinite(DrFz[ig]));
          assert(isfinite(DrN[ig]));
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
          // Update radiation quantities
          CCTK_REAL E = rE_p[i4D] + dt * rE_rhs[i4D] + theta * DrE[ig];
          F_d(1) = rFx_p[i4D] + dt * rFx_rhs[i4D] + theta * DrFx[ig];
          F_d(2) = rFy_p[i4D] + dt * rFy_rhs[i4D] + theta * DrFy[ig];
          F_d(3) = rFz_p[i4D] + dt * rFz_rhs[i4D] + theta * DrFz[ig];
          apply_floor(g_uu, &E, &F_d, rad_E_floor, rad_eps);

          CCTK_REAL N = rN_p[i4D] + dt * rN_rhs[i4D] + theta * DrN[ig];
          N = max(N, rad_N_floor);
	  assert(isfinite(momx[ijk]));
          // Compute back reaction on the fluid
          // NOTE: fluid backreaction is only needed at the last substep
          if (backreact && 0 == *TimeIntegratorStage) {
            assert(ngroups == 1);
            assert(nspecies == 3);
	    assert(isfinite(momx[ijk]));
	    printf("momx is finite: assert before modification is true!");
            momx[ijk] -= theta * DrFx[ig];
	    assert(isfinite(momx[ijk]));
            momy[ijk] -= theta * DrFy[ig];
            momz[ijk] -= theta * DrFz[ig];
            tau[ijk] -= theta * DrE[ig];
            DYe[ijk] += theta * DDxp[ig];
            // densxn[ijk] -= theta * DDxp[ig];
            netabs[ijk] += theta * DDxp[ig];
            netheat[ijk] -= theta * DrE[ig];
          }
          if (verbose) { CCTK_INFO("nuX_M1_CalcUpdate 9");}
          // Save updated results into grid functions
          rE[i4D] = E;
          unpack_F_d(F_d, &rFx[i4D], &rFy[i4D], &rFz[i4D]);
          rN[i4D] = N;
	  if (!isfinite(rN[i4D]) || !isfinite(rE[i4D]) ||
            !isfinite(rFx[i4D]) || !isfinite(rFy[i4D]) || !isfinite(rFz[i4D])) {
            printf("Non-finite before sync at (i,j,k)=(%d,%d,%d), ig=%d\n", p.i,p.j,p.k,ig);
          }

	  assert(isfinite(rN[i4D]));
	  assert(isfinite(rE[i4D]));
	  assert(isfinite(rFx[i4D]));
	  assert(isfinite(rFy[i4D]));
	  assert(isfinite(rFz[i4D]));
	  if (verbose) { CCTK_INFO("nuX_M1_CalcUpdate 10");}
        }
      });
}

} // namespace nuX_M1
