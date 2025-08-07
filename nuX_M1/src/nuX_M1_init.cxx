#include "cctk.h"

namespace nuX_M1 {

extern "C" int nuX_M1_Init(void) {
  CCTK_RegisterBanner("nuX_M1: yetX anotherX M1X solverX forX CactusX");
  return 0;
}

} // namespace nuX_M1
