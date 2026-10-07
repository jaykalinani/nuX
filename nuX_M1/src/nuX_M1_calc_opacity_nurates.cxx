#include <loop_device.hxx>

#include <algorithm>
#include <cassert>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

#include "m1_opacities.hpp"
#include "nuX_M1_opacity_utils.hxx"
#include "nuX_M1_weak_equil.hxx"
#include "nuX_rate_units.hxx"
#include "nuX_utils.hxx"
#include "setup_eos.hxx"

namespace nuX_M1 {
// using namespace thc;
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

void CalcOpacityNuRates(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_M1_CalcOpacityNuRates;
  DECLARE_CCTK_PARAMETERS;

  // Match the legacy dev integrator: evaluate opacities on the half-step
  // predictor pass and hold them fixed for the full-step corrector.
  if (CCTK_Equals(method, "semi-implicit") && *semi_implicit_stage == 1)
    return;

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

  // Opacity flags
  OpacityFlags opacity_flags = global_opac_flags;

  // Opacity parameters (corrections all switched off)
  OpacityParams opacity_pars = global_opac_params;

  // Setup EOS
  auto eos_3p = global_eos_3p_tab3d;
  if (!eos_3p)
    CCTK_ERROR("nuX_M1_CalcOpacityNuRates requires a tabulated EOS");
  // Setup Printer
  // thc::Printer::start(
  //         "[INFO|THC|THC_M1_CalcOpacity]: ",
  //         "[WARN|THC|THC_M1_CalcOpacity]: ",
  //         "[ERR|THC|THC_M1_CalcOpacity]: ",
  //         m1_max_num_msg, m1_max_num_msg);

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
          }
          return;
        }
        const CCTK_REAL rhoL =
            fmin(fmax(rho_raw, eos_3p->rgrho.min), eos_3p->rgrho.max);
        const CCTK_REAL tempL =
            fmin(fmax(temp_raw, eos_3p->rgtemp.min), eos_3p->rgtemp.max);
        const CCTK_REAL yeL =
            fmin(fmax(ye_raw, eos_3p->rgye.min), eos_3p->rgye.max);
        CCTK_REAL nb_nr =
            rhoL * nuX_dens_conv / (particle_mass * kBS_MeVtog); // CU to nm^-3
        const CCTK_REAL nb_transport = nb_nr / nuX_ndens_conv;
        const CCTK_REAL nb_fm3 =
            nb_nr / rate_units::physical_number_density_fm3_to_nm3;
        my_grey_opacity_params.eos_pars.nb = nb_nr;
        my_grey_opacity_params.eos_pars.temp = tempL;
        my_grey_opacity_params.eos_pars.ye = yeL;

        // NuRates uses yn/yp as the abundances of free neutrons/protons in
        // beta reactions, nucleon scattering, and the HR98/BRT06
        // bremsstrahlung kernels.  Ye and 1-Ye are total charge fractions and
        // include nucleons bound in nuclei, so obtain the free fractions from
        // the tabulated composition instead.
        using tabulated_eos = EOSX::eos_3p_tabulated3d;
        const auto eos_state = eos_3p->interptable->interpolate<
            tabulated_eos::EV::MU_P, tabulated_eos::EV::MU_N,
            tabulated_eos::EV::MU_E, tabulated_eos::EV::XN,
            tabulated_eos::EV::XP>(log(rhoL), log(tempL), yeL);
        const CCTK_REAL mu_pL = eos_state[0];
        const CCTK_REAL mu_nL = eos_state[1];
        const CCTK_REAL mu_eL = eos_state[2];
        const CCTK_REAL xnL = eos_state[3];
        const CCTK_REAL xpL = eos_state[4];
        if (!isfinite(mu_pL) || !isfinite(mu_nL) || !isfinite(mu_eL)) {
          for (int ig = 0; ig < ng; ++ig) {
            const int i4D = layout_cc.linear(p.i, p.j, p.k, ig);
            abs_0[i4D] = 0.0;
            abs_1[i4D] = 0.0;
            eta_0[i4D] = 0.0;
            eta_1[i4D] = 0.0;
            scat_1[i4D] = 0.0;
            nueave[i4D] = 0.0;
          }
          return;
        }
        my_grey_opacity_params.eos_pars.yn =
            isfinite(xnL) ? min(max(xnL, CCTK_REAL(0)), CCTK_REAL(1))
                          : CCTK_REAL(1) - yeL;
        my_grey_opacity_params.eos_pars.yp =
            isfinite(xpL) ? min(max(xpL, CCTK_REAL(0)), CCTK_REAL(1)) : yeL;

        my_grey_opacity_params.eos_pars.mu_p = mu_pL;
        my_grey_opacity_params.eos_pars.mu_n = mu_nL;
        my_grey_opacity_params.eos_pars.mu_e = mu_eL;

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
          }
          return;
        }
        const CCTK_REAL volformL = sqrt(spatial_det);
        CCTK_REAL const proper_dt = alphaL * dt / wL;
        CCTK_REAL nudens_0[4],
            nudens_1[4]; // force this to be 4 b/c nurates expects 4
        for (int ig = 0; ig < ngroups * nspecies; ++ig) {
          const int i4D = layout_cc.linear(p.i, p.j, p.k, ig);
          const CCTK_REAL in_fac = (ig == 2 && ng == 3) ? 0.25 : 1.0;

          nudens_0[ig] = in_fac * rnnu[i4D] / volformL;
          nudens_1[ig] = in_fac * rJ[i4D] / volformL;
          my_grey_opacity_params.m1_pars.n[ig] =
              nudens_0[ig] * nuX_ndens_conv; // transport unit to nm^-3
          my_grey_opacity_params.m1_pars.J[ig] =
              nudens_1[ig] * nuX_edens_conv; // CU to MeV nm^-3
          my_grey_opacity_params.m1_pars.chi[ig] = chi[i4D];

          // Fill data for anti-heavy neutrinos if only 3 species are evolved.
          if (ig == 2 && ng == 3) {
            nudens_0[3] = in_fac * rnnu[i4D] / volformL;
            nudens_1[3] = in_fac * rJ[i4D] / volformL;
            my_grey_opacity_params.m1_pars.n[3] =
                nudens_0[3] * nuX_ndens_conv; // transport unit to nm^-3
            my_grey_opacity_params.m1_pars.J[3] =
                nudens_1[3] * nuX_edens_conv; // CU to MeV nm^-3
            my_grey_opacity_params.m1_pars.chi[3] = chi[i4D];
          }
        }

        // Distribution parameters
        my_grey_opacity_params.distr_pars = CalculateDistrParamsFromM1(
            &my_grey_opacity_params.m1_pars, &my_grey_opacity_params.eos_pars);

        // Set up quadrature on GPU
        MyQuadrature gpu_quad;
        gpu_quad.nx = my_quad.nx;
        for (int idx = 0; idx < gpu_quad.nx; idx++) {
          gpu_quad.w[idx] = my_quad.w[idx];
          gpu_quad.points[idx] = my_quad.points[idx];
        }

        M1Opacities coeffs =
            ComputeM1Opacities(&gpu_quad, &gpu_quad, &my_grey_opacity_params);

        // Convert emissivities, opacities from nurates
        CCTK_REAL kappa_1_loc[MAX_GROUPSPECIES];
        CCTK_REAL abs_0_loc[MAX_GROUPSPECIES], abs_1_loc[MAX_GROUPSPECIES];
        CCTK_REAL scat_1_loc[MAX_GROUPSPECIES];
        CCTK_REAL eta_0_loc[MAX_GROUPSPECIES], eta_1_loc[MAX_GROUPSPECIES];

        for (int ig = 0; ig < ngroups * nspecies; ++ig) {
          const int i4D = layout_cc.linear(p.i, p.j, p.k, ig);
          const CCTK_REAL out_fac = (ig == 2 && ng == 3) ? 4.0 : 1.0;

          abs_0_loc[ig] = coeffs.kappa_0_a[ig] * nuX_length_conv;
          abs_1_loc[ig] = coeffs.kappa_a[ig] * nuX_length_conv;
          scat_1_loc[ig] = coeffs.kappa_s[ig] * nuX_length_conv;
          kappa_1_loc[ig] = abs_1_loc[ig] + scat_1_loc[ig];

          eta_0_loc[ig] =
              coeffs.eta_0[ig] / nuX_ndens_conv * nuX_time_conv * out_fac;
          eta_1_loc[ig] =
              coeffs.eta[ig] / nuX_edens_conv * nuX_time_conv * out_fac;

          if (!rate_is_valid(abs_0_loc[ig], opacity_rate_max) ||
              !rate_is_valid(abs_1_loc[ig], opacity_rate_max) ||
              !rate_is_valid(scat_1_loc[ig], opacity_rate_max) ||
              !rate_is_valid(eta_0_loc[ig], opacity_rate_max) ||
              !rate_is_valid(eta_1_loc[ig], opacity_rate_max)) {
            abs_0_loc[ig] = 0.0;
            abs_1_loc[ig] = 0.0;
            scat_1_loc[ig] = 0.0;
            kappa_1_loc[ig] = 0.0;
            eta_0_loc[ig] = 0.0;
            eta_1_loc[ig] = 0.0;
          }
        }
        /*---------------- ^^^ NuRates boilerplate ^^^ -------------*/

        // An effective optical depth used to decide whether to compute
        // the black body function for neutrinos assuming neutrino trapping
        // or at a fixed temperature and Ye
        CCTK_REAL const tau = min(sqrt(abs_1_loc[0] * kappa_1_loc[0]),
                                  sqrt(abs_1_loc[1] * kappa_1_loc[1])) *
                              proper_dt;

        // Compute the optically thin equilibrium first. It is also the
        // deterministic fallback if the trapped equilibrium or its moments
        // cannot be evaluated.
        CCTK_REAL nudens_0_thin[MAX_GROUPSPECIES];
        CCTK_REAL nudens_1_thin[MAX_GROUPSPECIES];
        NeutrinoDens(mu_nL, mu_pL, mu_eL, tempL, nudens_0_thin[0],
                     nudens_0_thin[1], nudens_0_thin[2], nudens_1_thin[0],
                     nudens_1_thin[1], nudens_1_thin[2]);
        for (int ig = 0; ig < ng; ++ig) {
          fallback_equilibrium_moments(nudens_0_thin[ig], nudens_1_thin[ig],
                                       CCTK_REAL(0), CCTK_REAL(0));
        }

        // Compute the neutrino black-body functions for trapped neutrinos.
        CCTK_REAL nudens_0_trap[MAX_GROUPSPECIES] = {CCTK_REAL(0)};
        CCTK_REAL nudens_1_trap[MAX_GROUPSPECIES] = {CCTK_REAL(0)};
        if (opacity_tau_trap >= 0 && tau > opacity_tau_trap) {

          const CCTK_REAL epsL = eos_3p->eps_from_rho_temp_ye(rhoL, tempL, yeL);
          CCTK_REAL eps_total = epsL;

          for (int ig = 0; ig < ngroups * nspecies; ++ig) {
            // rJ is densitized energy density.  Convert it to a physical
            // energy density and then to the same specific-energy units as
            // the EOS before forming the conserved total.
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
          // ierr = WeakEquilibrium(
          //         rho[ijk], temperature[ijk], Y_e[ijk],
          //         nudens_0[0], nudens_0[1], nudens_0[2],
          //         nudens_1[0], nudens_1[1], nudens_1[2],
          //         &temperature_trap, &Y_e_trap,
          //         &nudens_0_trap[0], &nudens_0_trap[1], &nudens_0_trap[2],
          //         &nudens_1_trap[0], &nudens_1_trap[1], &nudens_1_trap[2]);
          if (ierr) {
            // Try to recompute the weak equilibrium using neglecting
            // current neutrino data
            ierr =
                BetaEquilibriumTrapped(rhoL, nb_fm3, particle_mass, epsL, yeL,
                                       temp_trap, ye_trap, tempL, yeL, eos_3p);
            // ierr = WeakEquilibrium(
            //         rho[ijk], temperature[ijk], Y_e[ijk],
            //         0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            //         &temperature_trap, &Y_e_trap,
            //         &nudens_0_trap[0], &nudens_0_trap[1], &nudens_0_trap[2],
            //         &nudens_1_trap[0], &nudens_1_trap[1], &nudens_1_trap[2]);
            if (ierr) {
              // Keep the initialized fallback state (tempL, yeL) in this lane.
            }
          }

          if (ierr == 0) {
            CCTK_REAL mu_p_trap, mu_n_trap, mu_e_trap;
            eos_3p->mu_pne_from_rho_temp_ye(rhoL, temp_trap, ye_trap,
                                            mu_p_trap, mu_n_trap, mu_e_trap);

            NeutrinoDens(mu_n_trap, mu_p_trap, mu_e_trap, temp_trap,
                         nudens_0_trap[0], nudens_0_trap[1], nudens_0_trap[2],
                         nudens_1_trap[0], nudens_1_trap[1],
                         nudens_1_trap[2]);

            for (int ig = 0; ig < ng; ++ig) {
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

        // ierr = NeutrinoDensity(
        //         rho[ijk], temperature[ijk], Y_e[ijk],
        //         &nudens_0_thin[0], &nudens_0_thin[1], &nudens_0_thin[2],
        //         &nudens_1_thin[0], &nudens_1_thin[1], &nudens_1_thin[2]);
        // assert(!ierr);

        // Correct cross-sections for incoming neutrino energy
        for (int ig = 0; ig < ngroups * nspecies; ++ig) {
          int const i4D = layout_cc.linear(p.i, p.j, p.k, ig);

          // Set the neutrino black body function
          CCTK_REAL nudens_0, nudens_1;
          if (opacity_tau_trap < 0 || tau <= opacity_tau_trap) {
            nudens_0 = nudens_0_thin[ig];
            nudens_1 = nudens_1_thin[ig];
          } else if (tau > opacity_tau_trap + opacity_tau_delta) {
            nudens_0 = nudens_0_trap[ig];
            nudens_1 = nudens_1_trap[ig];
          } else {
            CCTK_REAL const lam = (tau - opacity_tau_trap) / opacity_tau_delta;
            nudens_0 = lam * nudens_0_trap[ig] + (1 - lam) * nudens_0_thin[ig];
            nudens_1 = lam * nudens_1_trap[ig] + (1 - lam) * nudens_1_thin[ig];
          }

          // Set the neutrino energies
          nueave[i4D] = nudens_0 > 0.0 ? nudens_1 / nudens_0 : 0.0;

          // Correct absorption opacities for non-LTE effects
          // (kappa ~ E_nu^2)
          const CCTK_REAL corr_fac = opacity_mean_energy_correction(
              rnnu[i4D], rJ[i4D], nudens_0, nudens_1,
              opacity_corr_fac_max);

          // Extract scattering opacity
          // scat_1[i4D] = corr_fac*(kappa_1_loc[ig] - abs_1_loc[ig]);
          scat_1[i4D] = corr_fac * scat_1_loc[ig];

          // Enforce Kirchhoff's laws.
          // . For the heavy lepton neutrinos this is implemented by
          //   changing the opacities.
          // . For the electron type neutrinos this is implemented by
          //   changing the emissivities.
          // It would be better to have emissivities and absorptivities
          // that satisfy Kirchhoff's law.
          if (ig == 2) {
            eta_0[i4D] = corr_fac * eta_0_loc[ig];
            eta_1[i4D] = corr_fac * eta_1_loc[ig];
            abs_0[i4D] = (nudens_0 > rad_N_floor ? eta_0[i4D] / nudens_0 : 0);
            abs_1[i4D] = (nudens_1 > rad_E_floor ? eta_1[i4D] / nudens_1 : 0);
          } else {
            abs_0[i4D] = corr_fac * abs_0_loc[ig];
            abs_1[i4D] = corr_fac * abs_1_loc[ig];
            eta_0[i4D] = abs_0[i4D] * nudens_0;
            eta_1[i4D] = abs_1[i4D] * nudens_1;
          }

          if (!rate_is_valid(abs_0[i4D], opacity_rate_max) ||
              !rate_is_valid(abs_1[i4D], opacity_rate_max) ||
              !rate_is_valid(scat_1[i4D], opacity_rate_max) ||
              !rate_is_valid(eta_0[i4D], opacity_rate_max) ||
              !rate_is_valid(eta_1[i4D], opacity_rate_max) ||
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

extern "C" void nuX_M1_CalcOpacityNuRates(CCTK_ARGUMENTS) {
  CalcOpacityNuRates(cctkGH);
}

} // namespace nuX_M1
