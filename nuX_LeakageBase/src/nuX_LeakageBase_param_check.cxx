#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

namespace nuX_LeakageBase {

extern "C" void nuX_LeakageBase_ParamCheck(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_LeakageBase_ParamCheck;
  DECLARE_CCTK_PARAMETERS;

  if (atmo_rho <= 0.0) {
    CCTK_PARAMWARN("atmo_rho must be positive");
  }

  if (use_fakerates && !CCTK_IsThornActive("nuX_FakeRates")) {
    CCTK_PARAMWARN(
        "nuX_LeakageBase::use_fakerates = yes requires the nuX_FakeRates "
        "thorn to be active");
  }

  if (CCTK_Equals(init_method, "spherical")) {
    CCTK_PARAMWARN(
        "The THC-style spherical optical-depth initialization is not ported "
        "yet in this initial scaffold");
  }
}

} // namespace nuX_LeakageBase
