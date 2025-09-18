#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

namespace nuX_M1 {

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

  if (*TimeIntegratorStage != 0) {
    CCTK_VERROR("Unexpected time integrator stage %d. "
                "Expected 'true' and 0.",
                (int)*TimeIntegratorStage);
  }

  //  if (!QueryProlongating() || *TimeIntegratorStage != 0) {
  //    CCTK_VERROR("Unexpected prolongation state %d or time integrator stage
  //    %d. "
  //                "Expected 'true' and 0.",
  //                (int)QueryProlongating(), (int)*TimeIntegratorStage);
  //
  //  }

  if (verbose) {
    CCTK_INFO("nuX_M1_FinalizeTimeIntegrator");
  }
}

} // namespace nuX_M1
