#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

namespace nuX_M1 {

extern "C" void nuX_M1_ParamCheck(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_M1_ParamCheck;
  DECLARE_CCTK_PARAMETERS

  if (optimize_prolongation) {
    if (cctk_nghostzones[0] < 4 || cctk_nghostzones[1] < 4 ||
        cctk_nghostzones[2] < 4) {
      CCTK_PARAMWARN("nuX_M1::optimize_prolongation requires at least "
                     "four ghost points");
    }
  }
}

} // namespace nuX_M1
