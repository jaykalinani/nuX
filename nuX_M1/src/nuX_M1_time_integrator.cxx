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

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

extern "C" void nuX_M1_InitTimeIntegrator(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_M1_InitTimeIntegrator;
  DECLARE_CCTK_PARAMETERS

  if (verbose) {
    CCTK_INFO("nuX_M1_InitTimeIntegrator");
  }

  *TimeIntegratorStage = 2;
  *M1_OriginalTime = cctkGH->cctk_time;
  cctkGH->cctk_time -= cctkGH->cctk_delta_time / cctkGH->cctk_timefac;
}

extern "C" void nuX_M1_UpdateTime(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_M1_UpdateTime;
  DECLARE_CCTK_PARAMETERS

  if (verbose) {
    CCTK_INFO("nuX_M1_UpdateTime");
  }

  if (*TimeIntegratorStage == 1) {
    CCTK_REAL dt = cctkGH->cctk_delta_time / cctkGH->cctk_timefac;
    cctkGH->cctk_time = *M1_OriginalTime - (dt / 2);
  } else if (*TimeIntegratorStage == 0) {
    cctkGH->cctk_time = *M1_OriginalTime;
  }

  if (verbose) {
    CCTK_VINFO("Integrated to time %e", cctkGH->cctk_time);
  }
}

extern "C" void nuX_M1_FinalizeTimeIntegrator(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_M1_FinalizeTimeIntegrator;
  DECLARE_CCTK_PARAMETERS

  if (!QueryProlongating() || *TimeIntegratorStage != 0) {
    CCTK_VERROR("Unexpected prolongation state %d or time integrator stage %d. "
                "Expected 'true' and 0.",
                (int)QueryProlongating(), (int)*TimeIntegratorStage);
  }

  if (verbose) {
    CCTK_INFO("nuX_M1_FinalizeTimeIntegrator");
  }
}
