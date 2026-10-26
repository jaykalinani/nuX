#include <loop_device.hxx>

#include <algorithm>
#include <cassert>
#include <limits>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

#include "m1_opacities.hpp"
#include "nuX_M1_closure.hxx"
#include "nuX_M1_nurates_cell.hxx"
#include "nuX_M1_nurates.hxx"
#include "nuX_M1_nurates_equilibrium.hxx"
#include "nuX_M1_opacity_utils.hxx"
#include "nuX_M1_weak_equil.hxx"
#include "nuX_rate_units.hxx"
#include "nuX_utils.hxx"
#include "setup_eos.hxx"

namespace nuX_M1 {
using namespace std;
using namespace Loop;
using namespace nuX_NuRates;
using namespace nuX_Utils;
using namespace EOSX;

#ifndef MAX_GROUPSPECIES
#define MAX_GROUPSPECIES 3
#endif

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline bool
rate_is_valid(CCTK_REAL const value, CCTK_REAL const max_abs) {
  return isfinite(value) && value >= CCTK_REAL(0) &&
         (max_abs < CCTK_REAL(0) || value <= max_abs);
}

// Keep these failures active in release GPU builds too. Ordinary assert (and
// some AMReX device assertions) can disappear when NDEBUG is defined.
CCTK_HOST CCTK_DEVICE inline void
require_nurates(bool const valid, char const *message) {
  if (valid)
    return;
#if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
  printf("%s\n", message);
  device_abort();
#else
  CCTK_ERROR(message);
#endif
}

void CalcOpacityNuRates(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_M1_CalcOpacityNuRates;
  DECLARE_CCTK_PARAMETERS;

  // NuRates depends on the radiation population, not only on matter.
  // Refresh at the predictor state before the corrector transport evaluation.

  if (verbose) {
    CCTK_INFO("nuX_M1_CalcOpacityNuRates");
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

  // NuRates Setup
  // Init structs for nurates calls
  MyQuadrature my_quad = {.type = kGauleg,
                          .alpha = -42.,
                          .dim = 1,
                          .nx = 10,
                          .ny = 1,
                          .nz = 1,
                          .x1 = 0.,
                          .x2 = 1.,
                          .y1 = -42.,
                          .y2 = -42.,
                          .z1 = -42.,
                          .z2 = -42.,
                          .points = {0},
                          .w = {0}};
  GaussLegendre(&my_quad);

  // Use the configured interaction flags and corrections.
  OpacityFlags opacity_flags = global_opac_flags;
  OpacityParams opacity_pars = global_opac_params;

  // Setup EOS
  auto eos_3p = global_eos_3p_tab3d;
  if (!eos_3p)
    CCTK_ERROR("nuX_M1_CalcOpacityNuRates requires a tabulated EOS");

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
          }
          return;
        }

        assert(nspecies == 3);
        assert(ngroups == 1);
        const int ng = nspecies * ngroups;
        /*---------------- vvv NuRates boilerplate vvv -------------*/
        // Init GreyOpacs struct
        GreyOpacityParams my_grey_opacity_params = {};

        // Neutrino reactions
        my_grey_opacity_params.opacity_flags = opacity_flags;

        // Opacity parameters
        my_grey_opacity_params.opacity_pars = opacity_pars;

        // Convert Thermodynamic Data to nurates
        const CCTK_REAL rho_raw = rho[ijk];
        const CCTK_REAL temp_raw = temperature[ijk];
        const CCTK_REAL ye_raw = Ye[ijk];
        CCTK_REAL rhoL;
        require_nurates(
            prepare_nurates_matter(eos_3p, rho_raw, temp_raw, ye_raw,
                                    particle_mass, rhoL,
                                    my_grey_opacity_params.eos_pars),
            "nuX_M1 NuRates: invalid matter/EOS state; rates cannot be evaluated");
        const CCTK_REAL tempL = my_grey_opacity_params.eos_pars.temp;

        // Convert M1 Data to nurates
        CCTK_REAL const gxx_cc = tensor::interp_v2c(gf_gxx, p);
        CCTK_REAL const gxy_cc = tensor::interp_v2c(gf_gxy, p);
        CCTK_REAL const gxz_cc = tensor::interp_v2c(gf_gxz, p);
        CCTK_REAL const gyy_cc = tensor::interp_v2c(gf_gyy, p);
        CCTK_REAL const gyz_cc = tensor::interp_v2c(gf_gyz, p);
        CCTK_REAL const gzz_cc = tensor::interp_v2c(gf_gzz, p);
        CCTK_REAL const alphaL = tensor::interp_v2c(gf_alp, p);
        const CCTK_REAL spatial_det = nuX_Utils::metric::spatial_det(
            gxx_cc, gxy_cc, gxz_cc, gyy_cc, gyz_cc, gzz_cc);
        CCTK_REAL const wL = fidu_w_lorentz[ijk];
        require_nurates(
            isfinite(spatial_det) && spatial_det > CCTK_REAL(0) &&
                isfinite(alphaL) && alphaL > CCTK_REAL(0) &&
                isfinite(wL) && wL >= CCTK_REAL(1) && isfinite(dt) &&
                dt >= CCTK_REAL(0),
            "nuX_M1 NuRates: invalid metric, observer, or timestep");
        const CCTK_REAL volformL = sqrt(spatial_det);
        CCTK_REAL const proper_dt = alphaL * dt / wL;
        CCTK_REAL nudens_0[total_num_species];
        CCTK_REAL nudens_1[total_num_species];
        for (int ig = 0; ig < ngroups * nspecies; ++ig) {
          const int i4D = layout_cc.linear(p.i, p.j, p.k, ig);
          const CCTK_REAL in_fac = (ig == 2 && ng == 3) ? 0.25 : 1.0;

          nudens_0[ig] = in_fac * rnnu[i4D] / volformL;
          nudens_1[ig] = in_fac * rJ[i4D] / volformL;
          my_grey_opacity_params.m1_pars.n[ig] =
              nudens_0[ig] * nuX_ndens_conv; // transport unit to nm^-3
          my_grey_opacity_params.m1_pars.J[ig] =
              nudens_1[ig] * nuX_edens_conv; // CU to MeV nm^-3
          const CCTK_REAL n_nr = my_grey_opacity_params.m1_pars.n[ig];
          const CCTK_REAL J_nr = my_grey_opacity_params.m1_pars.J[ig];
          constexpr CCTK_REAL chi_slack =
              64 * numeric_limits<CCTK_REAL>::epsilon();
          require_nurates(
              equilibrium_moments_are_valid(n_nr, J_nr) &&
                  isfinite(chi[i4D]) &&
                  chi[i4D] >= CCTK_REAL(1) / 3 - chi_slack &&
                  chi[i4D] <= CCTK_REAL(1) + chi_slack,
              "nuX_M1 NuRates: invalid radiation moments for spectral reconstruction");
          my_grey_opacity_params.m1_pars.chi[ig] =
              min(max(chi[i4D], CCTK_REAL(1) / 3), CCTK_REAL(1));

          // Fill data for anti-heavy neutrinos if only 3 species are evolved.
          if (ig == 2 && ng == 3) {
            nudens_0[3] = in_fac * rnnu[i4D] / volformL;
            nudens_1[3] = in_fac * rJ[i4D] / volformL;
            my_grey_opacity_params.m1_pars.n[3] =
                my_grey_opacity_params.m1_pars.n[ig];
            my_grey_opacity_params.m1_pars.J[3] =
                my_grey_opacity_params.m1_pars.J[ig];
            my_grey_opacity_params.m1_pars.chi[3] =
                my_grey_opacity_params.m1_pars.chi[ig];
          }
        }

        // Keep the original moments for the equilibrium energy/lepton
        // budget; the temporary floor spectrum is only a rate reference.
        M1Quantities budget_moments = my_grey_opacity_params.m1_pars;
        for (int sp = 0; sp < 3; ++sp) {
          budget_moments.n[sp] = nudens_0[sp] * nuX_ndens_conv;
          budget_moments.J[sp] = nudens_1[sp] * nuX_edens_conv;
        }
        bool floor_only[3];
        nurates_rate_reference(my_grey_opacity_params.m1_pars, volformL,
                                tempL, rad_N_floor, rad_E_floor, floor_only);

        // Distribution parameters
        my_grey_opacity_params.distr_pars = CalculateDistrParamsFromM1(
            &my_grey_opacity_params.m1_pars, &my_grey_opacity_params.eos_pars);
        for (int ig = 0; ig < total_num_species; ++ig) {
          const auto &spectrum = my_grey_opacity_params.distr_pars;
          require_nurates(
              isfinite(spectrum.temp_t[ig]) && spectrum.temp_t[ig] > 0 &&
                  isfinite(spectrum.eta_t[ig]) &&
                  isfinite(spectrum.temp_f[ig]) && spectrum.temp_f[ig] > 0 &&
                  isfinite(spectrum.beta_f[ig]) && spectrum.beta_f[ig] >= 0,
              "nuX_M1 NuRates: nonfinite reconstructed radiation spectrum");
          const auto &radiation = my_grey_opacity_params.m1_pars;
          if (spectrum.w_t[ig] > 0 && radiation.n[ig] > 0) {
            if (!ThickDistributionMatchesM1(
                    radiation.n[ig], radiation.J[ig], spectrum.temp_t[ig],
                    spectrum.eta_t[ig])) {
              // Failure-only context: no diagnostic grid storage or extra
              // work on successful rate evaluations.
              printf("nuX_M1 NuRates: spectrum mismatch at (%e,%e,%e), "
                     "species %d: n=%e J=%e Tnu=%e eta=%e chi=%e "
                     "Tmatter=%e rho=%e step_dt=%e log_y=%e\n",
                     p.x, p.y, p.z, ig, radiation.n[ig], radiation.J[ig],
                     spectrum.temp_t[ig], spectrum.eta_t[ig],
                     radiation.chi[ig], tempL, rhoL, dt,
                     4 * log(radiation.n[ig]) - 3 * log(radiation.J[ig]) -
                         log(kBS_FourPi_hc3));
              require_nurates(
                  false,
                  "nuX_M1 NuRates: reconstructed FD spectrum does not match moments");
            }
          }
        }

        // Set up quadrature on GPU
        MyQuadrature gpu_quad;
        gpu_quad.nx = my_quad.nx;
        for (int idx = 0; idx < gpu_quad.nx; idx++) {
          gpu_quad.w[idx] = my_quad.w[idx];
          gpu_quad.points[idx] = my_quad.points[idx];
        }

        const M1OpacitiesNonThermalSeparated coeffs =
            ComputeM1OpacitiesNonThermalSeparated(
                &gpu_quad, &gpu_quad, &my_grey_opacity_params);

        // Convert emissivities, opacities from nurates
        nurates_raw_opacities_t raw[MAX_GROUPSPECIES];

        for (int ig = 0; ig < ngroups * nspecies; ++ig) {
          const CCTK_REAL out_fac = (ig == 2 && ng == 3) ? 4.0 : 1.0;

          raw[ig] = convert_nurates_opacities(coeffs, ig, out_fac);
          const auto &rates = raw[ig];
          if (!rate_is_valid(rates.abs_0, opacity_rate_max) ||
              !rate_is_valid(rates.abs_1_thermal, opacity_rate_max) ||
              !rate_is_valid(rates.abs_1_nonthermal, opacity_rate_max) ||
              !rate_is_valid(rates.scat_1, opacity_rate_max) ||
              !rate_is_valid(rates.eta_0, opacity_rate_max) ||
              !rate_is_valid(rates.eta_1_thermal, opacity_rate_max) ||
              !rate_is_valid(rates.eta_1_nonthermal, opacity_rate_max)) {
            printf("nuX_M1 NuRates: invalid rates at (%d,%d,%d), species %d: "
                   "absN=%e absTh=%e absNT=%e scat=%e etaN=%e etaTh=%e "
                   "etaNT=%e (opacity_rate_max=%e)\n",
                   p.i, p.j, p.k, ig, rates.abs_0, rates.abs_1_thermal,
                   rates.abs_1_nonthermal, rates.scat_1, rates.eta_0,
                   rates.eta_1_thermal, rates.eta_1_nonthermal, opacity_rate_max);
            require_nurates(
                false, "nuX_M1 NuRates: refusing to silently disable interactions");
          }
        }
        /*---------------- ^^^ NuRates boilerplate ^^^ -------------*/

        CCTK_REAL equilibrium_number[3], equilibrium_energy[3];
        nurates_equilibrium_targets(
            eos_3p, rhoL, my_grey_opacity_params.eos_pars, budget_moments, raw,
            particle_mass, proper_dt, opacity_tau_trap, opacity_tau_delta,
            equilibrium_number, equilibrium_energy);

        // Correct cross-sections for incoming neutrino energy
        for (int ig = 0; ig < ngroups * nspecies; ++ig) {
          int const i4D = layout_cc.linear(p.i, p.j, p.k, ig);

          // Set the neutrino black body function
          const CCTK_REAL nudens_0 = equilibrium_number[ig];
          const CCTK_REAL nudens_1 = equilibrium_energy[ig];

          // Set the neutrino energies
          nueave[i4D] = nudens_0 > 0.0 ? nudens_1 / nudens_0 : 0.0;

          // Correct absorption opacities for non-LTE effects
          // (kappa ~ E_nu^2)
          const auto &radiation = my_grey_opacity_params.m1_pars;
          const auto rates = correct_nurates_opacities(
              raw[ig], ig == 2, floor_only[ig], rnnu[i4D], rJ[i4D],
              radiation.n[ig] > THRESHOLD_N, radiation.J[ig] > THRESHOLD_J,
              nudens_0, nudens_1,
              opacity_corr_fac_max, rad_N_floor, rad_E_floor);
          abs_0[i4D] = rates.abs_0;
          abs_1[i4D] = rates.abs_1;
          scat_1[i4D] = rates.scat_1;
          eta_0[i4D] = rates.eta_0;
          eta_1[i4D] = rates.eta_1;

          if (!rate_is_valid(abs_0[i4D], opacity_rate_max) ||
              !rate_is_valid(abs_1[i4D], opacity_rate_max) ||
              !rate_is_valid(scat_1[i4D], opacity_rate_max) ||
              !rate_is_valid(eta_0[i4D], opacity_rate_max) ||
              !rate_is_valid(eta_1[i4D], opacity_rate_max) ||
              !isfinite(nueave[i4D]) || nueave[i4D] < CCTK_REAL(0)) {
            printf("nuX_M1 NuRates: invalid corrected rates at (%d,%d,%d), "
                   "species %d: absN=%e absE=%e scat=%e etaN=%e etaE=%e "
                   "mean=%e (opacity_rate_max=%e)\n",
                   p.i, p.j, p.k, ig, abs_0[i4D], abs_1[i4D], scat_1[i4D],
                   eta_0[i4D], eta_1[i4D], nueave[i4D], opacity_rate_max);
            require_nurates(
                false, "nuX_M1 NuRates: invalid equilibrium-corrected coefficients");
          }
        }
      });
}

extern "C" void nuX_M1_CalcOpacityNuRates(CCTK_ARGUMENTS) {
  CalcOpacityNuRates(cctkGH);
}

} // namespace nuX_M1
