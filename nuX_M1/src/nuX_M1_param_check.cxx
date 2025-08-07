#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

namespace nuX_M1 {

extern "C" void nuX_M1_ParamCheck(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_M1_ParamCheck;
  DECLARE_CCTK_PARAMETERS

  if (CCTK_Equals(nuX_m1_test, "const")) {
    if (ngroups != 1 || nspecies != 3) {
      CCTK_PARAMWARN(
          "If \"nuX_M1::nuX_m1_test\" is set to \"const\", then"
          "\"nuX_M1::ngroups\" must be 1 and \"nuX_M1::nspecies\" must be 3");
    }
  }

  if (CCTK_Equals(nuX_m1_test, "shadow") ||
      CCTK_Equals(nuX_m1_test, "sphere")) {
    if (!CCTK_Equals(initial_hydro, "nuX_M1")) {
      CCTK_PARAMWARN("This test requires HydroBase::initial_data "
                     "to be \"nuX_M1\"");
    }
  }

  if (optimize_prolongation) {
    if (cctk_nghostzones[0] < 4 || cctk_nghostzones[1] < 4 ||
        cctk_nghostzones[2] < 4) {
      CCTK_PARAMWARN("nuX_M1::optimize_prolongation requires at least "
                     "four ghost points");
    }
  }
}

} // namespace nuX_M1
