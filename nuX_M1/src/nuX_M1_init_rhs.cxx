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

  const GF3D2layout layout2(cctkGH, {1, 1, 1});

  grid.loop_int_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        const int groupspec = ngroups * nspecies; // total components

        for (int ig = 0; ig < groupspec; ++ig) {
          const int i4D = layout2.linear(p.i, p.j, p.k, ig);

          // Zero the RHS buffers
          rN_rhs[i4D] = 0.0;
          rE_rhs[i4D] = 0.0;
          rFx_rhs[i4D] = 0.0;
          rFy_rhs[i4D] = 0.0;
          rFz_rhs[i4D] = 0.0;
        }
      });
}

} // namespace nuX_M1
