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

#include "nuX_fakerates.hxx"
#include "setup_eos.hxx"

namespace nuX_M1 {
using namespace std;
using namespace Loop;
using namespace nuX_FakeRates;
using namespace EOSX;

#ifndef MAX_GROUPSPECIES
#define MAX_GROUPSPECIES 3
#endif

extern "C" void nuX_M1_CalcFakeOpacity(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_M1_CalcFakeOpacity;
  DECLARE_CCTK_PARAMETERS;

  // Opacities are constant throught the timestep
  if (*TimeIntegratorStage != 2) {
    return;
  }

  if (verbose) {
    CCTK_INFO("nuX_M1_CalcFakeOpacity");
  }

  const GridDescBaseDevice grid(cctkGH);
  const GF3D2layout layout2(cctkGH, {1, 1, 1});

  CCTK_REAL const dt = CCTK_DELTA_TIME;

  FakeRatesDef *myfakerates = global_fakerates;

  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        const int ijk = layout2.linear(p.i, p.j, p.k);

        if (nuX_m1_mask[ijk]) {
          for (int ig = 0; ig < nspecies * ngroups; ++ig) {
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
        // TODO: currently use MAX_GROUPSPECIES instead of ng for array
        // initialization

        CCTK_REAL rhoL = rho[ijk];
        CCTK_REAL tempL = temperature[ijk];
        CCTK_REAL yeL = Ye[ijk];
        CCTK_REAL nb_nr =
            rhoL * nuX_dens_conv / (particle_mass * kBS_MeVtog); // CU to nm^-3
        CCTK_REAL nbL = nb_nr / nuX_ndens_conv;
        M1Opacities coeffs = myfakerates->ComputeFakeOpacities(rhoL);

        // Convert M1 Data to nurates
        CCTK_REAL volformL = volform[ijk];
        CCTK_REAL nudens_0[4],
            nudens_1[4]; // force this to be 4 b/c nurates expects 4
        for (int ig = 0; ig < ngroups * nspecies; ++ig) {
          const int i4D = layout2.linear(p.i, p.j, p.k, ig);
          const CCTK_REAL dup_fac =
              1.0 / ((1.0 + (ig > 1)) * (1.0 + (ng == 3)));

          nudens_0[ig] = dup_fac * rnnu[i4D] / volformL;
          nudens_1[ig] = dup_fac * rJ[i4D] / volformL;

          // Fill data for anti-heavy neutrinos if only 3 species are evolved.
          if (ig == 2 && ng == 3) {
            nudens_0[3] = dup_fac * rnnu[i4D] / volformL;
            nudens_1[3] = dup_fac * rJ[i4D] / volformL;
          }
        }

        // Convert emissivities, opacities from nurates
        CCTK_REAL kappa_0_loc[MAX_GROUPSPECIES], kappa_1_loc[MAX_GROUPSPECIES];
        CCTK_REAL abs_0_loc[MAX_GROUPSPECIES], abs_1_loc[MAX_GROUPSPECIES];
        CCTK_REAL scat_0_loc[MAX_GROUPSPECIES], scat_1_loc[MAX_GROUPSPECIES];
        CCTK_REAL eta_0_loc[MAX_GROUPSPECIES], eta_1_loc[MAX_GROUPSPECIES];

        for (int ig = 0; ig < ngroups * nspecies; ++ig) {
          const int i4D = layout2.linear(p.i, p.j, p.k, ig);
          const CCTK_REAL out_fac = (1.0 + (ig > 1)) * (1.0 + (ng == 3));

          abs_0_loc[ig] = coeffs.kappa_0_a[ig];
          abs_1_loc[ig] = coeffs.kappa_a[ig];
          scat_0_loc[ig] = 0.0;
          scat_1_loc[ig] = coeffs.kappa_s[ig];
          kappa_0_loc[ig] = abs_0_loc[ig];
          kappa_1_loc[ig] = abs_1_loc[ig];

          eta_0_loc[ig] = coeffs.eta_0[ig];
          eta_1_loc[ig] = coeffs.eta[ig];
        }

        // An effective optical depth used to decide whether to compute
        // the black body function for neutrinos assuming neutrino trapping
        // or at a fixed temperature and Ye
        CCTK_REAL const tau = min(sqrt(abs_1_loc[0] * kappa_1_loc[0]),
                                  sqrt(abs_1_loc[1] * kappa_1_loc[1])) *
                              dt;

        // Compute the neutrino black body functions assuming trapped neutrinos
        CCTK_REAL nudens_0_trap[MAX_GROUPSPECIES],
            nudens_1_trap[MAX_GROUPSPECIES];
        if (opacity_tau_trap >= 0 && tau > opacity_tau_trap) {

          CCTK_REAL epsL = eps[ijk];
          CCTK_REAL etot = epsL;

          for (int ig = 0; ig < ngroups * nspecies; ++ig) {
            etot += rJ[layout2.linear(p.i, p.j, p.k, ig)];
          }

          // TODO: Change BetaEq call to accept more lepton fractions if 4
          // species are evolved
          CCTK_REAL ylep_e = yeL - (nudens_0[0] - nudens_0[1]) / nbL;
          CCTK_REAL temp_trap = tempL;
          CCTK_REAL ye_trap = ylep_e;

          CCTK_REAL mu_p_trap, mu_n_trap, mu_e_trap;
          // eos_3p->mu_pne_from_valid_rho_temp_ye(rhoL, temp_trap, ye_trap,
          // mu_p_trap, mu_n_trap, mu_e_trap);

          myfakerates->FakeNeutrinoDens(
              rhoL, nudens_0_trap[0], nudens_0_trap[1], nudens_0_trap[2],
              nudens_1_trap[0], nudens_1_trap[1], nudens_1_trap[2]);

          // NOTE: the block below will never be run because ng is always
          // assumed to be 3
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

        // Compute the neutrino black body function assuming fixed temperature
        // and Y_e
        CCTK_REAL nudens_0_thin[3], nudens_1_thin[3];
        myfakerates->FakeNeutrinoDens(rhoL, nudens_0_thin[0], nudens_0_thin[1],
                                      nudens_0_thin[2], nudens_1_thin[0],
                                      nudens_1_thin[1], nudens_1_thin[2]);

        // NeutrinoDens assumes 3 species transport. Split heavy density if 4
        // species are used
        if (ng == 4) {
          nudens_0_thin[2] *= 0.5;
          nudens_1_thin[2] *= 0.5;
          nudens_0_thin[3] = nudens_0_thin[2];
          nudens_1_thin[3] = nudens_1_thin[2];
        }

        // Correct cross-sections for incoming neutrino energy
        for (int ig = 0; ig < ngroups * nspecies; ++ig) {
          int const i4D = layout2.linear(p.i, p.j, p.k, ig);

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
          nueave[i4D] = nudens_1 / nudens_0;

          // Correct absorption opacities for non-LTE effects
          // (kappa ~ E_nu^2)
          CCTK_REAL corr_fac = 1.0;
          corr_fac = (rJ[i4D] / rnnu[i4D]) * (nudens_0 / nudens_1);
          if (!isfinite(corr_fac)) {
            corr_fac = 1.0;
          }
          corr_fac *= corr_fac;
          corr_fac = max(1.0 / opacity_corr_fac_max,
                         min(corr_fac, opacity_corr_fac_max));

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
        }
      });
}

} // namespace nuX_M1
