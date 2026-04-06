#include <cmath>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

namespace nuX_M0 {

extern "C" void nuX_M0_ParamCheck(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_M0_ParamCheck;
  DECLARE_CCTK_PARAMETERS;

  int const ntheta = static_cast<int>(std::round(std::sqrt(nray / 2.0)));
  if (2 * ntheta * ntheta != nray) {
    CCTK_PARAMWARN("nray must equal 2*ntheta*ntheta for an integer ntheta");
  }

  if (bns_sep_threshold > 0.0 && !CCTK_IsThornActive("BNSTrackerGen")) {
    CCTK_PARAMWARN(
        "bns_sep_threshold requires BNSTrackerGen in this initial scaffold");
  }

  if (excision) {
    CCTK_PARAMWARN("The THC excision path is not ported yet");
  }
}

} // namespace nuX_M0
