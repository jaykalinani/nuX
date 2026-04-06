#include "cctk.h"
#include "cctk_Arguments.h"

namespace nuX_LeakageBase {

extern "C" void nuX_LeakageBase_SelectBC(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_nuX_LeakageBase_SelectBC;

  // Do nothing. CarpetX performs the actual sync requested in schedule.ccl.
}

} // namespace nuX_LeakageBase
