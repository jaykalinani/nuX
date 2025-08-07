#include <cstring>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "nuX_utils.hxx"

namespace nuX_M1 {

extern "C" void nuX_M1_InitialCopy(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_M1_InitialCopy;
  DECLARE_CCTK_PARAMETERS

  if (verbose) {
    CCTK_INFO("nuX_M1_InitialCopy");
  }

  size_t siz = UTILS_GFSIZE(cctkGH) * nspecies * sizeof(CCTK_REAL);

  std::memcpy(rN, rN_p, siz);
  std::memcpy(rE, rE_p, siz);
  std::memcpy(rFx, rFx_p, siz);
  std::memcpy(rFy, rFy_p, siz);
  std::memcpy(rFz, rFz_p, siz);
}

} // namespace nuX_M1
