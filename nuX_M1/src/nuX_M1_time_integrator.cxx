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
}

extern "C" void nuX_M1_FinalizeTimeIntegrator(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_M1_FinalizeTimeIntegrator;
  DECLARE_CCTK_PARAMETERS

  if (*TimeIntegratorStage != 0) {
    CCTK_VERROR("Unexpected time integrator stage %d. "
                "Expected 'true' and 0.",
                (int)*TimeIntegratorStage);
  }

  if (verbose) {
    CCTK_INFO("nuX_M1_FinalizeTimeIntegrator");
  }
}

} // namespace nuX_M1
