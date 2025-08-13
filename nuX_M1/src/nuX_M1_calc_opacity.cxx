//  Templated Hydrodynamics Code: an hydro code built on top of HRSCCore
//  Copyright (C) 2020, David Radice <david.radice@psu.edu>
//
//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include <loop_device.hxx>

#include <algorithm>
#include <cassert>
#include <sstream>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "m1_opacities.hpp"
#include "nuX_M1_weak_equil.hxx"
#include "setup_eos.hxx"

// #include "thc_printer.hh"
// #include "thc_M1_macro.h"
// 
// #include "utils.hh"

namespace nuX_M1 {
// using namespace thc;
using namespace std;
using namespace Loop;
using namespace nuX_Rates;
using namespace EOSX;

extern "C" void nuX_M1_CalcOpacity(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS_nuX_M1_CalcOpacity;
    DECLARE_CCTK_PARAMETERS;

    // Opacities are constant throught the timestep
    if (*TimeIntegratorStage != 2) {
        return;
    }

    if (verbose) {
        CCTK_INFO("nuX_M1_CalcOpacity");
    }

    const GridDescBaseDevice grid(cctkGH);
    const GF3D2layout layout2(cctkGH, {1, 1, 1});

    CCTK_REAL const dt = CCTK_DELTA_TIME;

    // NuRates Setup
    // Init structs for nurates calls
    MyQuadrature my_quad = {.type   = kGauleg,
      .alpha  = -42.,
      .dim    = 1,
      .nx     = 10,
      .ny     = 1,
      .nz     = 1,
      .x1     = 0.,
      .x2     = 1.,
      .y1     = -42.,
      .y2     = -42.,
      .z1     = -42.,
      .z2     = -42.,
      .points = {0},
      .w      = {0}};
    GaussLegendre(&my_quad);
    
    // Opacity flags
    OpacityFlags opacity_flags = global_opac_flags;
 
    // Opacity parameters (corrections all switched off)
    OpacityParams opacity_pars = global_opac_params;
 
    // Setup EOS
    auto eos_3p = global_eos_3p_tab3d;

    // Setup Printer
    // thc::Printer::start(
    //         "[INFO|THC|THC_M1_CalcOpacity]: ",
    //         "[WARN|THC|THC_M1_CalcOpacity]: ",
    //         "[ERR|THC|THC_M1_CalcOpacity]: ",
    //         m1_max_num_msg, m1_max_num_msg);

  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
            const int ijk = layout2.linear(p.i, p.j, p.k);

            if (nuX_m1_mask[ijk]) {
                for (int ig = 0; ig < nspecies*ngroups; ++ig) {
                    int const i4D = layout2.linear(p.i, p.j, p.k, ig);
                    abs_0[i4D] = 0.0;
                    abs_1[i4D] = 0.0;
                    eta_0[i4D] = 0.0;
                    eta_1[i4D] = 0.0;
                    scat_1[i4D] = 0.0;
                }
                return;
            }

            assert(nspecies == 3);
            assert(ngroups == 1);
            const int ng = nspecies * ngroups;

            /*---------------- vvv NuRates boilerplate vvv -------------*/
            // Init GreyOpacs struct
            GreyOpacityParams my_grey_opacity_params;

            // Neutrino reactions
            my_grey_opacity_params.opacity_flags = opacity_flags;

            // Opacity parameters
            my_grey_opacity_params.opacity_pars = opacity_pars;

            // Convert Thermodynamic Data to nurates
            CCTK_REAL rhoL = rho[ijk];
            CCTK_REAL tempL = temperature[ijk];
            CCTK_REAL yeL = Ye[ijk];
            CCTK_REAL nb_nr = rhoL * nuX_dens_conv / (particle_mass * kBS_MeVtog); // CU to nm^-3
            CCTK_REAL nbL = nb_nr / nuX_ndens_conv;
            my_grey_opacity_params.eos_pars.nb   = nb_nr;
            my_grey_opacity_params.eos_pars.temp = tempL;
            my_grey_opacity_params.eos_pars.yp   = yeL;
            my_grey_opacity_params.eos_pars.yn   = 1.0 - yeL;

            CCTK_REAL mu_pL, mu_nL, mu_eL;
            eos_3p->mu_pne_from_valid_rho_temp_ye(rhoL, tempL, yeL, mu_pL, mu_nL, mu_eL);
            my_grey_opacity_params.eos_pars.mu_p = mu_pL;
            my_grey_opacity_params.eos_pars.mu_n = mu_nL;
            my_grey_opacity_params.eos_pars.mu_e = mu_eL;

            // Convert M1 Data to nurates
            CCTK_REAL volformL = volform[ijk];
            CCTK_REAL nudens_0[4], nudens_1[4]; // force this to be 4 b/c nurates expects 4
            for (int ig = 0; ig < ngroups * nspecies; ++ig) {
              const int i4D = layout2.linear(p.i, p.j, p.k, ig);
              const CCTK_REAL dup_fac = 1.0 / ((1.0 + (ig > 1)) * (1.0 + (ng == 3)));

              nudens_0[ig] = dup_fac * rnnu[i4D] / volformL;
              nudens_1[ig] = dup_fac * rJ[i4D] / volformL;
              my_grey_opacity_params.m1_pars.n[ig]   = nudens_0[ig] * nuX_ndens_conv; // fm^-3 to nm^-3
              my_grey_opacity_params.m1_pars.J[ig]   = nudens_1[ig] * nuX_edens_conv; // CU to MeV nm^-3
              my_grey_opacity_params.m1_pars.chi[ig] = chi[i4D];

              // Fill data for anti-heavy neutrinos if only 3 species are evolved.
              if (ig == 2 && ng == 3) {
                nudens_0[3] = dup_fac * rnnu[i4D] / volformL;
                nudens_1[3] = dup_fac * rJ[i4D] / volformL;
                my_grey_opacity_params.m1_pars.n[3]   = nudens_0[3] * nuX_ndens_conv; // fm^-3 to nm^-3
                my_grey_opacity_params.m1_pars.J[3]   = nudens_1[3] * nuX_edens_conv; // CU to MeV nm^-3
                my_grey_opacity_params.m1_pars.chi[3] = chi[i4D];
              }
            }

            // Distribution parameters
            my_grey_opacity_params.distr_pars = CalculateDistrParamsFromM1(
              &my_grey_opacity_params.m1_pars,
              &my_grey_opacity_params.eos_pars);

            // Set up quadrature on GPU
            MyQuadrature gpu_quad;
            gpu_quad.nx = my_quad.nx;
            for (int idx = 0; idx < gpu_quad.nx; idx++)
            {
              gpu_quad.w[idx]      = my_quad.w[idx];
              gpu_quad.points[idx] = my_quad.points[idx];
            }

            M1Opacities coeffs = ComputeM1Opacities(&gpu_quad, &gpu_quad,
                &my_grey_opacity_params);

            // Convert emissivities, opacities from nurates
            CCTK_REAL kappa_0_loc[ng], kappa_1_loc[ng];
            CCTK_REAL abs_0_loc[ng], abs_1_loc[ng];
            CCTK_REAL scat_0_loc[ng], scat_1_loc[ng];
            CCTK_REAL eta_0_loc[ng], eta_1_loc[ng];

            for (int ig = 0; ig < ngroups * nspecies; ++ig) {
              const int i4D = layout2.linear(p.i, p.j, p.k, ig);
              const CCTK_REAL out_fac = (1.0 + (ig > 1)) * (1.0 + (ng == 3));

              abs_0_loc[ig] = coeffs.kappa_0_a[ig] * nuX_length_conv;
              abs_1_loc[ig] = coeffs.kappa_a[ig] * nuX_length_conv;
              scat_0_loc[ig] = 0.0;
              scat_1_loc[ig] = coeffs.kappa_s[ig] * nuX_length_conv;
              kappa_0_loc[ig] = abs_0_loc[ig] + scat_0_loc[ig];
              kappa_1_loc[ig] = abs_1_loc[ig] + scat_1_loc[ig];

              eta_0_loc[ig] = coeffs.eta_0[ig] / nuX_ndens_conv * nuX_time_conv * out_fac;
              eta_1_loc[ig] = coeffs.eta[ig] / nuX_edens_conv * nuX_time_conv * out_fac;
            }
            /*---------------- ^^^ NuRates boilerplate ^^^ -------------*/
      
            // An effective optical depth used to decide whether to compute
            // the black body function for neutrinos assuming neutrino trapping
            // or at a fixed temperature and Ye
            CCTK_REAL const tau = min(
                    sqrt(abs_1_loc[0]*kappa_1_loc[0]),
                    sqrt(abs_1_loc[1]*kappa_1_loc[1]))*dt;

            // Compute the neutrino black body functions assuming trapped neutrinos
            CCTK_REAL nudens_0_trap[ng], nudens_1_trap[ng];
            if (opacity_tau_trap >= 0 && tau > opacity_tau_trap) {

                CCTK_REAL epsL = eos_3p->eps_from_valid_rho_temp_ye(rhoL, tempL, yeL);
                CCTK_REAL etot = epsL;

                for (int ig = 0; ig < ngroups * nspecies; ++ig) {
                  etot += rJ[layout2.linear(p.i, p.j, p.k, ig)];
                }

                // CCTK_REAL const nudens_0[3] = {
                //     rnnu[layout2.linear(i, j, k, 0)]/volform[ijk],
                //     rnnu[layout2.linear(i, j, k, 1)]/volform[ijk],
                //     rnnu[layout2.linear(i, j, k, 2)]/volform[ijk],
                // };
                // CCTK_REAL const nudens_1[3] = {
                //     rJ[layout2.linear(i, j, k, 0)]/volform[ijk],
                //     rJ[layout2.linear(i, j, k, 1)]/volform[ijk],
                //     rJ[layout2.linear(i, j, k, 2)]/volform[ijk],
                // };

                // TODO: Change BetaEq call to accept more lepton fractions if 4 species are evolved
                CCTK_REAL ylep_e = yeL - (nudens_0[0] - nudens_0[1]) / nbL;
                CCTK_REAL temp_trap, ye_trap;
                int ierr = BetaEquilibriumTrapped(rhoL, nbL, etot, ylep_e, temp_trap, ye_trap, tempL, yeL, eos_3p);
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
                    ierr = BetaEquilibriumTrapped(rhoL, nbL, epsL, yeL, temp_trap, ye_trap, tempL, yeL, eos_3p);
                    // ierr = WeakEquilibrium(
                    //         rho[ijk], temperature[ijk], Y_e[ijk],
                    //         0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                    //         &temperature_trap, &Y_e_trap,
                    //         &nudens_0_trap[0], &nudens_0_trap[1], &nudens_0_trap[2],
                    //         &nudens_1_trap[0], &nudens_1_trap[1], &nudens_1_trap[2]);
                    if (ierr) {
                        printf("Could not find the weak equilibrium!");
                        // printf("Reflevel = " << ilogb(cctkGH->cctk_levfac[0]) << endl;
                        printf("Iteration = %i\n", cctk_iteration);
                        printf("(i, j, k) = (%i, %i, %i)\n", p.i, p.j, p.k);
                        printf("(x, y, z) = (%e, %e, %e)\n", p.x, p.y, p.z);
                        printf("rho = %e\n", rhoL);
                        printf("temperature = %e\n", tempL);
                        printf("Y_e = %e\n", yeL);
                        printf("nudens_0 = %e, %e, %e, %e\n", nudens_0[0], nudens_0[1], nudens_0[2], nudens_0[3]);
                        printf("nudens_1 = %e, %e, %e, %e\n", nudens_1[0], nudens_1[1], nudens_1[2], nudens_1[3]);

                        // ostringstream ss;
                        // ss << "Could not find the weak equilibrium!" << endl;
                        // ss << "Reflevel = " << ilogb(cctkGH->cctk_levfac[0]) << endl;
                        // ss << "Iteration = " << cctk_iteration << endl;
                        // ss << "(i, j, k) = (" << i << ", " << j << ", " << k << ")\n";
                        // ss << "(x, y, z) = (" << x[ijk] << ", " << y[ijk] << ", "
                        //                       << z[ijk] << ")\n";
                        // ss << "rho = " << rho[ijk] << endl;
                        // ss << "temperature = " << temperature[ijk] << endl;
                        // ss << "Y_e = " << Y_e[ijk] << endl;
                        // ss << "alp = " << alp[ijk] << endl;
                        // ss << "nudens_0 = " << nudens_0[0] << " " << nudens_0[1]
                        //                     << " " << nudens_0[2] << endl;
                        // ss << "nudens_1 = " << nudens_1[0] << " " << nudens_1[1]
                        //                     << " " << nudens_1[2] << endl;
                        // Printer::print_warn(ss.str());
                    }
                }

                CCTK_REAL mu_p_trap, mu_n_trap, mu_e_trap;
                eos_3p->mu_pne_from_valid_rho_temp_ye(rhoL, temp_trap, ye_trap, mu_p_trap, mu_n_trap, mu_e_trap);

                NeutrinoDens(mu_n_trap, mu_p_trap, mu_e_trap, temp_trap, nudens_0_trap[0], nudens_0_trap[1],
                  nudens_0_trap[2], nudens_1_trap[0], nudens_1_trap[1], nudens_1_trap[2]);

                if (ng == 4) {
                  nudens_0_trap[2] *= 0.5;
                  nudens_1_trap[2] *= 0.5;
                  nudens_0_trap[3] = nudens_0_trap[2];
                  nudens_1_trap[3] = nudens_1_trap[2];
                }

                assert(isfinite(nudens_0_trap[0]));
                assert(isfinite(nudens_0_trap[1]));
                assert(isfinite(nudens_0_trap[2]));
                assert(isfinite(nudens_1_trap[0]));
                assert(isfinite(nudens_1_trap[1]));
                assert(isfinite(nudens_1_trap[2]));
            }

            // Compute the neutrino black body function assuming fixed temperature and Y_e
            CCTK_REAL nudens_0_thin[3], nudens_1_thin[3];
            NeutrinoDens(mu_nL, mu_pL, mu_eL, tempL, nudens_0_thin[0], nudens_0_thin[1],
              nudens_0_thin[2], nudens_1_thin[0], nudens_1_thin[1], nudens_1_thin[2]);

            // NeutrinoDens assumes 3 species transport. Split heavy density if 4 species are used
            if (ng == 4) {
              nudens_0_thin[2] *= 0.5;
              nudens_1_thin[2] *= 0.5;
              nudens_0_thin[3] = nudens_0_thin[2];
              nudens_1_thin[3] = nudens_1_thin[2];
            }

            // ierr = NeutrinoDensity(
            //         rho[ijk], temperature[ijk], Y_e[ijk],
            //         &nudens_0_thin[0], &nudens_0_thin[1], &nudens_0_thin[2],
            //         &nudens_1_thin[0], &nudens_1_thin[1], &nudens_1_thin[2]);
            // assert(!ierr);

            // Correct cross-sections for incoming neutrino energy
            for (int ig = 0; ig < ngroups*nspecies; ++ig) {
                int const i4D = layout2.linear(p.i, p.j, p.k, ig);

                // Set the neutrino black body function
                CCTK_REAL nudens_0, nudens_1;
                if (opacity_tau_trap < 0 || tau <= opacity_tau_trap) {
                    nudens_0 = nudens_0_thin[ig];
                    nudens_1 = nudens_1_thin[ig];
                }
                else if (tau > opacity_tau_trap + opacity_tau_delta) {
                    nudens_0 = nudens_0_trap[ig];
                    nudens_1 = nudens_1_trap[ig];
                }
                else {
                    CCTK_REAL const lam = (tau - opacity_tau_trap)/opacity_tau_delta;
                    nudens_0 = lam*nudens_0_trap[ig] + (1 - lam)*nudens_0_thin[ig];
                    nudens_1 = lam*nudens_1_trap[ig] + (1 - lam)*nudens_1_thin[ig];
                }

                // Set the neutrino energies
                nueave[i4D] = nudens_1/nudens_0;

                // Correct absorption opacities for non-LTE effects
                // (kappa ~ E_nu^2)
                CCTK_REAL corr_fac = 1.0;
                corr_fac = (rJ[i4D]/rnnu[i4D])*(nudens_0/nudens_1);
                if (!isfinite(corr_fac)) {
                    corr_fac = 1.0;
                }
                corr_fac *= corr_fac;
                corr_fac = max(1.0/opacity_corr_fac_max, min(corr_fac, opacity_corr_fac_max));

                // Extract scattering opacity
                // scat_1[i4D] = corr_fac*(kappa_1_loc[ig] - abs_1_loc[ig]);
                scat_1[i4D] = corr_fac*scat_1_loc[ig];

                // Enforce Kirchhoff's laws.
                // . For the heavy lepton neutrinos this is implemented by
                //   changing the opacities.
                // . For the electron type neutrinos this is implemented by
                //   changing the emissivities.
                // It would be better to have emissivities and absorptivities
                // that satisfy Kirchhoff's law.
                if (ig == 2) {
                    eta_0[i4D] = corr_fac*eta_0_loc[ig];
                    eta_1[i4D] = corr_fac*eta_1_loc[ig];
                    abs_0[i4D] = (nudens_0 > rad_N_floor ? eta_0[i4D]/nudens_0 : 0);
                    abs_1[i4D] = (nudens_1 > rad_E_floor ? eta_1[i4D]/nudens_1 : 0);
                }
                else {
                    abs_0[i4D] = corr_fac*abs_0_loc[ig];
                    abs_1[i4D] = corr_fac*abs_1_loc[ig];
                    eta_0[i4D] = abs_0[i4D]*nudens_0;
                    eta_1[i4D] = abs_1[i4D]*nudens_1;
                }
            }
        }); // UTILS_ENDLOOP3(thc_m1_calc_opacity);
    // Done with printing
    // thc::Printer::stop();
}

} // namespace 
