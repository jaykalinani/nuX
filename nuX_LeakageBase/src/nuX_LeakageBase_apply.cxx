#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

namespace nuX_LeakageBase {

extern "C" void nuX_LeakageBase_Apply(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  static bool warned = false;
  if (!warned && CCTK_MyProc(cctkGH) == 0) {
    CCTK_WARN(CCTK_WARN_ALERT, "nuX_LeakageBase_Apply has not been ported yet");
    warned = true;
  }
}

} // namespace nuX_LeakageBase
