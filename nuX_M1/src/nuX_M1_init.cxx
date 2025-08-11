#include "cctk.h"

namespace nuX_M1 {

extern "C" int nuX_M1_Init(void) {
  CCTK_RegisterBanner(
      "nuX_M1: the best M1 neutrino transport solver for CarpetX");
  return 0;
}

} // namespace nuX_M1
