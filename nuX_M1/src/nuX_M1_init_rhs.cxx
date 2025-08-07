#include <cstring>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "nuX_utils.hxx"

namespace nuX_M1 {

extern "C" void nuX_M1_InitRHS(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_M1_InitRHS;
  DECLARE_CCTK_PARAMETERS

  if (verbose) {
    CCTK_INFO("nuX_M1_InitRHS");
  }

  size_t siz = UTILS_GFSIZE(cctkGH) * nspecies * sizeof(CCTK_REAL);

  std::memset(rN_rhs, 0, siz);
  std::memset(rE_rhs, 0, siz);
  std::memset(rFx_rhs, 0, siz);
  std::memset(rFy_rhs, 0, siz);
  std::memset(rFz_rhs, 0, siz);
}

} // namespace nuX_M1
