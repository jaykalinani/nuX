#include <loop_device.hxx>

#include <algorithm>
#include <cassert>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

#include "nuX_M1_diagnostics.hxx"
#include "nuX_M1_opacity_utils.hxx"
#include "nuX_M1_weak_equil.hxx"
#include "nuX_rate_units.hxx"
#include "nuX_utils.hxx"
#include "nuX_weakrates.hxx"
#include "setup_eos.hxx"

namespace nuX_M1 {
using namespace std;
using namespace Loop;
using namespace nuX_Utils;
using namespace EOSX;

#ifndef MAX_GROUPSPECIES
#define MAX_GROUPSPECIES 3
#endif

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline bool
weak_rate_is_valid(CCTK_REAL const value, CCTK_REAL const max_abs) {
  return isfinite(value) && value >= CCTK_REAL(0) &&
         (max_abs < CCTK_REAL(0) || value <= max_abs);
}

template <bool store, typename Layout, typename EOS>
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
record_equilibrium_diagnostics(const EquilibriumDiagnostics<store> &diagnostics,
                               const Layout &layout, const int i, const int j,
                               const int k, const int ng,
                               const int equilibrium_status,
                               const CCTK_REAL rho, const CCTK_REAL nb,
                               const CCTK_REAL particle_mass,
                               const CCTK_REAL eps_total, const CCTK_REAL ylep,
                               const CCTK_REAL temperature, const CCTK_REAL ye,
                               EOS *const eos) {
  if constexpr (store) {
    CCTK_REAL residual[2] = {CCTK_REAL(1.0e300), CCTK_REAL(1.0e300)};
    const CCTK_REAL state[2] = {temperature, ye};
    if (func_eq_weak(rho, nb, particle_mass, eps_total, ylep, state, residual,
                     eos) != 0) {
      residual[0] = CCTK_REAL(1.0e300);
      residual[1] = CCTK_REAL(1.0e300);
    }
    diagnostics.record_all(layout, i, j, k, ng, equilibrium_status, residual[0],
                           residual[1]);
  }
}

template <bool store_source_diagnostics>
void CalcOpacityWeakRates(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_M1_CalcOpacityWeakRates;
  DECLARE_CCTK_PARAMETERS;

  // Match the legacy dev integrator: evaluate opacities on the half-step
  // predictor pass and hold them fixed for the full-step corrector.
  if (CCTK_Equals(method, "semi-implicit") && *semi_implicit_stage == 1)
    return;

  if (verbose) {
    CCTK_INFO("nuX_M1_CalcOpacityWeakRates");
  }

  const GridDescBaseDevice grid(cctkGH);
  const GF3D2layout layout_cc(cctkGH, {1, 1, 1});
  const GF3D2layout layout_vc(cctkGH, {0, 0, 0});
  const GF3D2<const CCTK_REAL> gf_gxx(layout_vc, gxx);
  const GF3D2<const CCTK_REAL> gf_gxy(layout_vc, gxy);
  const GF3D2<const CCTK_REAL> gf_gxz(layout_vc, gxz);
  const GF3D2<const CCTK_REAL> gf_gyy(layout_vc, gyy);
  const GF3D2<const CCTK_REAL> gf_gyz(layout_vc, gyz);
  const GF3D2<const CCTK_REAL> gf_gzz(layout_vc, gzz);
  const GF3D2<const CCTK_REAL> gf_alp(layout_vc, alp);

  // Opacity trapping is a macro-step decision. ODESolvers temporarily changes
  // CCTK_DELTA_TIME for diagonal implicit source solves, so use the saved step
  // dt.
  const CCTK_REAL step_delta_time = ODESolvers_GetStepDeltaTime();
  CCTK_REAL const dt =
      step_delta_time > 0.0 ? step_delta_time : CCTK_DELTA_TIME;

  auto eos_3p = global_eos_3p_tab3d;
  auto weakrates = nuX_WeakRates::global_weakrates;
  if (!eos_3p)
    CCTK_ERROR("nuX_M1_CalcOpacityWeakRates requires a tabulated EOS");
  if (!weakrates)
    CCTK_ERROR("nuX_M1_CalcOpacityWeakRates requires nuX_WeakRates");
  const EquilibriumDiagnostics<store_source_diagnostics> diagnostics(
      equilibrium_status, equilibrium_lepton_residual,
      equilibrium_energy_residual);

  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        const int ijk = layout_cc.linear(p.i, p.j, p.k);

        if (nuX_m1_mask[ijk]) {
          for (int ig = 0; ig < nspecies * ngroups; ++ig) {
            int const i4D = layout_cc.linear(p.i, p.j, p.k, ig);
            abs_0[i4D] = 0.0;
            abs_1[i4D] = 0.0;
            eta_0[i4D] = 0.0;
            eta_1[i4D] = 0.0;
            scat_1[i4D] = 0.0;
            nueave[i4D] = 0.0;
            diagnostics.clear(i4D);
          }
          return;
        }

        assert(nspecies == 3);
        assert(ngroups == 1);
        const int ng = nspecies * ngroups;
        for (int ig = 0; ig < ng; ++ig)
          diagnostics.clear(layout_cc.linear(p.i, p.j, p.k, ig));

        const CCTK_REAL rho_raw = rho[ijk];
        const CCTK_REAL temp_raw = temperature[ijk];
        const CCTK_REAL ye_raw = Ye[ijk];
        if (!isfinite(rho_raw) || rho_raw <= CCTK_REAL(0) ||
            !isfinite(temp_raw) || temp_raw <= CCTK_REAL(0) ||
            !isfinite(ye_raw) || ye_raw < CCTK_REAL(0) ||
            ye_raw > CCTK_REAL(1) || !isfinite(particle_mass) ||
            particle_mass <= CCTK_REAL(0)) {
          for (int ig = 0; ig < ng; ++ig) {
            const int i4D = layout_cc.linear(p.i, p.j, p.k, ig);
            abs_0[i4D] = 0.0;
            abs_1[i4D] = 0.0;
            eta_0[i4D] = 0.0;
            eta_1[i4D] = 0.0;
            scat_1[i4D] = 0.0;
            nueave[i4D] = 0.0;
            diagnostics.invalid(i4D);
          }
          return;
        }
        const CCTK_REAL rhoL =
            fmin(fmax(rho_raw, eos_3p->rgrho.min), eos_3p->rgrho.max);
        const CCTK_REAL tempL =
            fmin(fmax(temp_raw, eos_3p->rgtemp.min), eos_3p->rgtemp.max);
        const CCTK_REAL yeL =
            fmin(fmax(ye_raw, eos_3p->rgye.min), eos_3p->rgye.max);
        const CCTK_REAL nb_nm3 = rhoL * rate_units::code_density_to_g_nm3 /
                                 (particle_mass * rate_units::mev_mass_to_g);
        const CCTK_REAL nb_transport =
            nb_nm3 / rate_units::transport_number_density_to_nm3;
        const CCTK_REAL nb_fm3 =
            nb_nm3 / rate_units::physical_number_density_fm3_to_nm3;

        const CCTK_REAL gxx_cc = tensor::interp_v2c(gf_gxx, p);
        const CCTK_REAL gxy_cc = tensor::interp_v2c(gf_gxy, p);
        const CCTK_REAL gxz_cc = tensor::interp_v2c(gf_gxz, p);
        const CCTK_REAL gyy_cc = tensor::interp_v2c(gf_gyy, p);
        const CCTK_REAL gyz_cc = tensor::interp_v2c(gf_gyz, p);
        const CCTK_REAL gzz_cc = tensor::interp_v2c(gf_gzz, p);
        const CCTK_REAL alphaL = tensor::interp_v2c(gf_alp, p);
        const CCTK_REAL spatial_det = nuX_Utils::metric::spatial_det(
            gxx_cc, gxy_cc, gxz_cc, gyy_cc, gyz_cc, gzz_cc);
        const CCTK_REAL wL = fidu_w_lorentz[ijk];
        if (!isfinite(spatial_det) || spatial_det <= CCTK_REAL(0) ||
            !isfinite(alphaL) || alphaL <= CCTK_REAL(0) || !isfinite(wL) ||
            wL < CCTK_REAL(1) || !isfinite(dt) || dt < CCTK_REAL(0)) {
          for (int ig = 0; ig < ng; ++ig) {
            const int i4D = layout_cc.linear(p.i, p.j, p.k, ig);
            abs_0[i4D] = 0.0;
            abs_1[i4D] = 0.0;
            eta_0[i4D] = 0.0;
            eta_1[i4D] = 0.0;
            scat_1[i4D] = 0.0;
            nueave[i4D] = 0.0;
            diagnostics.invalid(i4D);
          }
          return;
        }
        const CCTK_REAL volformL = sqrt(spatial_det);
        const CCTK_REAL proper_dt = alphaL * dt / wL;

        CCTK_REAL nudens_0[4], nudens_1[4];
        for (int ig = 0; ig < ng; ++ig) {
          const int i4D = layout_cc.linear(p.i, p.j, p.k, ig);
          const CCTK_REAL in_fac = (ig == 2 && ng == 3) ? 0.25 : 1.0;
          nudens_0[ig] = in_fac * rnnu[i4D] / volformL;
          nudens_1[ig] = in_fac * rJ[i4D] / volformL;
          if (ig == 2 && ng == 3) {
            nudens_0[3] = nudens_0[ig];
            nudens_1[3] = nudens_1[ig];
          }
        }

        using tabulated_eos = EOSX::eos_3p_tabulated3d;
        const CCTK_REAL lr = log(rhoL);
        const CCTK_REAL lt = log(tempL);
        const auto eos_state = eos_3p->interptable->interpolate<
            tabulated_eos::EV::MU_E, tabulated_eos::EV::MU_P,
            tabulated_eos::EV::MU_N, tabulated_eos::EV::XA,
            tabulated_eos::EV::XH, tabulated_eos::EV::XN, tabulated_eos::EV::XP,
            tabulated_eos::EV::ABAR, tabulated_eos::EV::ZBAR>(lr, lt, yeL);
        const nuX_WeakRates::EOSState weak_eos = {
            rhoL * rate_units::code_density_to_g_nm3 /
                rate_units::per_cm3_to_per_nm3,
            tempL,
            yeL,
            particle_mass,
            eos_state[0],
            eos_state[1],
            eos_state[2],
            eos_state[3],
            eos_state[4],
            eos_state[5],
            eos_state[6],
            eos_state[7],
            eos_state[8]};
        const auto coeffs = weakrates->compute_rates(weak_eos);

        CCTK_REAL kappa_1_loc[MAX_GROUPSPECIES];
        CCTK_REAL abs_0_loc[MAX_GROUPSPECIES], abs_1_loc[MAX_GROUPSPECIES];
        CCTK_REAL scat_1_loc[MAX_GROUPSPECIES];
        CCTK_REAL eta_0_loc[MAX_GROUPSPECIES], eta_1_loc[MAX_GROUPSPECIES];

        for (int ig = 0; ig < ng; ++ig) {
          const CCTK_REAL kappa_0_a =
              coeffs.kappa_0_a[ig] * rate_units::per_cm_to_per_nm;
          const CCTK_REAL kappa_a =
              coeffs.kappa_a[ig] * rate_units::per_cm_to_per_nm;
          const CCTK_REAL kappa_s =
              coeffs.kappa_s[ig] * rate_units::per_cm_to_per_nm;
          const CCTK_REAL emissivity_0 =
              coeffs.eta_0[ig] * rate_units::per_cm3_to_per_nm3;
          const CCTK_REAL emissivity =
              coeffs.eta[ig] * rate_units::per_cm3_to_per_nm3;

          abs_0_loc[ig] = kappa_0_a * rate_units::code_length_to_nm;
          abs_1_loc[ig] = kappa_a * rate_units::code_length_to_nm;
          scat_1_loc[ig] = kappa_s * rate_units::code_length_to_nm;
          kappa_1_loc[ig] = abs_1_loc[ig] + scat_1_loc[ig];
          eta_0_loc[ig] = emissivity_0 /
                          rate_units::transport_number_density_to_nm3 *
                          rate_units::code_time_to_s;
          eta_1_loc[ig] = emissivity /
                          rate_units::code_energy_density_to_mev_nm3 *
                          rate_units::code_time_to_s;

          if (!weak_rate_is_valid(abs_0_loc[ig], opacity_rate_max) ||
              !weak_rate_is_valid(abs_1_loc[ig], opacity_rate_max) ||
              !weak_rate_is_valid(scat_1_loc[ig], opacity_rate_max) ||
              !weak_rate_is_valid(eta_0_loc[ig], opacity_rate_max) ||
              !weak_rate_is_valid(eta_1_loc[ig], opacity_rate_max)) {
            abs_0_loc[ig] = 0.0;
            abs_1_loc[ig] = 0.0;
            scat_1_loc[ig] = 0.0;
            kappa_1_loc[ig] = 0.0;
            eta_0_loc[ig] = 0.0;
            eta_1_loc[ig] = 0.0;
          }
        }

        const CCTK_REAL tau = min(sqrt(abs_1_loc[0] * kappa_1_loc[0]),
                                  sqrt(abs_1_loc[1] * kappa_1_loc[1])) *
                              proper_dt;

        // Compute the optically thin state first so it can be used as a
        // deterministic fallback for an invalid trapped equilibrium.
        CCTK_REAL nudens_0_thin[MAX_GROUPSPECIES];
        CCTK_REAL nudens_1_thin[MAX_GROUPSPECIES];
        const auto thin_densities = weakrates->equilibrium_densities(weak_eos);
        for (int ig = 0; ig < ng; ++ig) {
          nudens_0_thin[ig] = thin_densities.number[ig] *
                              rate_units::per_cm3_to_per_nm3 /
                              rate_units::transport_number_density_to_nm3;
          nudens_1_thin[ig] = thin_densities.energy[ig] *
                              rate_units::per_cm3_to_per_nm3 /
                              rate_units::code_energy_density_to_mev_nm3;
          fallback_equilibrium_moments(nudens_0_thin[ig], nudens_1_thin[ig],
                                       CCTK_REAL(0), CCTK_REAL(0));
        }

        CCTK_REAL nudens_0_trap[MAX_GROUPSPECIES] = {CCTK_REAL(0)};
        CCTK_REAL nudens_1_trap[MAX_GROUPSPECIES] = {CCTK_REAL(0)};
        if (opacity_tau_trap >= 0 && tau > opacity_tau_trap) {
          const CCTK_REAL epsL = eos_3p->eps_from_rho_temp_ye(rhoL, tempL, yeL);
          CCTK_REAL eps_total = epsL;
          for (int ig = 0; ig < ng; ++ig) {
            eps_total +=
                rJ[layout_cc.linear(p.i, p.j, p.k, ig)] / (volformL * rhoL);
          }

          const CCTK_REAL ylep_e =
              yeL + (nudens_0[0] - nudens_0[1]) / nb_transport;
          CCTK_REAL temp_trap = tempL;
          CCTK_REAL ye_trap = yeL;
          int ierr = BetaEquilibriumTrapped(rhoL, nb_fm3, particle_mass,
                                            eps_total, ylep_e, temp_trap,
                                            ye_trap, tempL, yeL, eos_3p);
          int equilibrium_status_local = ierr == 0 ? 1 : 0;
          if (ierr) {
            ierr =
                BetaEquilibriumTrapped(rhoL, nb_fm3, particle_mass, epsL, yeL,
                                       temp_trap, ye_trap, tempL, yeL, eos_3p);
            equilibrium_status_local = ierr == 0 ? 2 : -1;
          }

          record_equilibrium_diagnostics(diagnostics, layout_cc, p.i, p.j, p.k,
                                         ng, equilibrium_status_local, rhoL,
                                         nb_fm3, particle_mass, eps_total,
                                         ylep_e, temp_trap, ye_trap, eos_3p);

          if (ierr == 0) {
            CCTK_REAL mu_p_trap, mu_n_trap, mu_e_trap;
            eos_3p->mu_pne_from_rho_temp_ye(rhoL, temp_trap, ye_trap,
                                            mu_p_trap, mu_n_trap, mu_e_trap);
            const nuX_WeakRates::EOSState weak_trap_eos = {
                weak_eos.rho, temp_trap, ye_trap, particle_mass, mu_e_trap,
                mu_p_trap,    mu_n_trap, 0.0,       0.0,           0.0,
                0.0,          0.0,       0.0};
            const auto trap_densities =
                weakrates->equilibrium_densities(weak_trap_eos);
            for (int ig = 0; ig < ng; ++ig) {
              nudens_0_trap[ig] = trap_densities.number[ig] *
                                  rate_units::per_cm3_to_per_nm3 /
                                  rate_units::transport_number_density_to_nm3;
              nudens_1_trap[ig] = trap_densities.energy[ig] *
                                  rate_units::per_cm3_to_per_nm3 /
                                  rate_units::code_energy_density_to_mev_nm3;
              fallback_equilibrium_moments(
                  nudens_0_trap[ig], nudens_1_trap[ig], nudens_0_thin[ig],
                  nudens_1_thin[ig]);
            }
          } else {
            // A failed nonlinear solve does not define a trapped state.  Do
            // not turn its last iterate into opacities merely because the
            // resulting moments happen to be finite.
            for (int ig = 0; ig < ng; ++ig) {
              nudens_0_trap[ig] = nudens_0_thin[ig];
              nudens_1_trap[ig] = nudens_1_thin[ig];
            }
          }
        }

        for (int ig = 0; ig < ng; ++ig) {
          const int i4D = layout_cc.linear(p.i, p.j, p.k, ig);
          CCTK_REAL nudens_0, nudens_1;
          if (opacity_tau_trap < 0 || tau <= opacity_tau_trap) {
            nudens_0 = nudens_0_thin[ig];
            nudens_1 = nudens_1_thin[ig];
          } else if (tau > opacity_tau_trap + opacity_tau_delta) {
            nudens_0 = nudens_0_trap[ig];
            nudens_1 = nudens_1_trap[ig];
          } else {
            const CCTK_REAL lam = (tau - opacity_tau_trap) / opacity_tau_delta;
            nudens_0 = lam * nudens_0_trap[ig] + (1 - lam) * nudens_0_thin[ig];
            nudens_1 = lam * nudens_1_trap[ig] + (1 - lam) * nudens_1_thin[ig];
          }

          nueave[i4D] = nudens_0 > 0.0 ? nudens_1 / nudens_0 : 0.0;

          const CCTK_REAL corr_fac = opacity_mean_energy_correction(
              rnnu[i4D], rJ[i4D], nudens_0, nudens_1,
              opacity_corr_fac_max);

          scat_1[i4D] = corr_fac * scat_1_loc[ig];
          if (ig == 2) {
            eta_0[i4D] = corr_fac * eta_0_loc[ig];
            eta_1[i4D] = corr_fac * eta_1_loc[ig];
            abs_0[i4D] = nudens_0 > rad_N_floor ? eta_0[i4D] / nudens_0 : 0.0;
            abs_1[i4D] = nudens_1 > rad_E_floor ? eta_1[i4D] / nudens_1 : 0.0;
          } else {
            abs_0[i4D] = corr_fac * abs_0_loc[ig];
            abs_1[i4D] = corr_fac * abs_1_loc[ig];
            eta_0[i4D] = abs_0[i4D] * nudens_0;
            eta_1[i4D] = abs_1[i4D] * nudens_1;
          }

          if (!weak_rate_is_valid(abs_0[i4D], opacity_rate_max) ||
              !weak_rate_is_valid(abs_1[i4D], opacity_rate_max) ||
              !weak_rate_is_valid(scat_1[i4D], opacity_rate_max) ||
              !weak_rate_is_valid(eta_0[i4D], opacity_rate_max) ||
              !weak_rate_is_valid(eta_1[i4D], opacity_rate_max) ||
              !isfinite(nueave[i4D]) || nueave[i4D] < CCTK_REAL(0)) {
            abs_0[i4D] = 0.0;
            abs_1[i4D] = 0.0;
            eta_0[i4D] = 0.0;
            eta_1[i4D] = 0.0;
            scat_1[i4D] = 0.0;
            nueave[i4D] = 0.0;
          }
        }
      });
}

extern "C" void nuX_M1_CalcOpacityWeakRates(CCTK_ARGUMENTS) {
  DECLARE_CCTK_PARAMETERS;
  if (source_diagnostics)
    CalcOpacityWeakRates<true>(cctkGH);
  else
    CalcOpacityWeakRates<false>(cctkGH);
}

} // namespace nuX_M1
