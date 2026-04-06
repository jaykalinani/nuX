#include "cctk.h"

namespace nuX_M0 {

extern "C" int nuX_M0_Init(void) {
  CCTK_RegisterBanner("nuX_M0: THC ray-by-ray M0 transport port");
  return 0;
}

} // namespace nuX_M0
