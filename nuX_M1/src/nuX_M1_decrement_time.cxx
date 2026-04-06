#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

namespace nuX_M1 {

extern "C" void nuX_M1_DecrementTimeIntegratorStage(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_nuX_M1_DecrementTimeIntegratorStage;
  DECLARE_CCTK_PARAMETERS;

  if (verbose && CCTK_MyProc(cctkGH) == 0) {
    CCTK_VINFO("nuX_M1_DecrementStage ENTER: current TimeIntegratorStage=%d "
               "cctk_time=%e",
               (int)*TimeIntegratorStage, (double)cctkGH->cctk_time);
  }

  /* Defensive: only decrement if positive; log if unexpected */
  if (*TimeIntegratorStage > 0) {
    --(*TimeIntegratorStage);
    if (verbose && CCTK_MyProc(cctkGH) == 0) {
      CCTK_VINFO("nuX_M1_DecrementStage: decremented -> %d",
                 (int)*TimeIntegratorStage);
    }
  } else {
    /* If this prints repeatedly, there's an ordering/race bug earlier */
    if (verbose && CCTK_MyProc(cctkGH) == 0) {
      CCTK_VINFO("nuX_M1_DecrementStage: nothing to do (TimeIntegratorStage=%d)",
                 (int)*TimeIntegratorStage);
    }
  }
}

} // namespace nuX_M1
