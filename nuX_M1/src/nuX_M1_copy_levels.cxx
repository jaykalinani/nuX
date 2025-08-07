#include <cstring>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "nuX_utils.hxx"

namespace nuX_M1 {

extern "C" void nuX_M1_CopyLevels(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_M1_CopyLevels;
  DECLARE_CCTK_PARAMETERS

  if (verbose) {
    CCTK_INFO("nuX_M1_CopyLevels");
  }

  size_t siz = UTILS_GFSIZE(cctkGH) * nspecies * sizeof(CCTK_REAL);

  std::memcpy(rN_p, rN, siz);
  std::memcpy(rE_p, rE, siz);
  std::memcpy(rFx_p, rFx, siz);
  std::memcpy(rFy_p, rFy, siz);
  std::memcpy(rFz_p, rFz, siz);

  if (timelevels > 2) {
    std::memcpy(rN_p_p, rN, siz);
    std::memcpy(rE_p_p, rE, siz);
    std::memcpy(rFx_p_p, rFx, siz);
    std::memcpy(rFy_p_p, rFy, siz);
    std::memcpy(rFz_p_p, rFz, siz);
  }
}

} // namespace nuX_M1
