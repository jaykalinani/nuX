#include "cctk.h"

namespace nuX_LeakageBase {

extern "C" int nuX_LeakageBase_Init(void) {
  CCTK_RegisterBanner("nuX_LeakageBase: THC leakage infrastructure port");
  return 0;
}

} // namespace nuX_LeakageBase
