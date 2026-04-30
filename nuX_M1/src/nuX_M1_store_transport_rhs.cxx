#include <loop_device.hxx>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

namespace nuX_M1 {

using namespace Loop;

extern "C" void nuX_M1_StoreTransportRHS(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_M1_StoreTransportRHS;
  DECLARE_CCTK_PARAMETERS;

  if (verbose) {
    CCTK_INFO("nuX_M1_StoreTransportRHS");
  }

  const GridDescBaseDevice grid(cctkGH);
  const GF3D2layout layout_cc(cctkGH, {1, 1, 1});
  const int groupspec = ngroups * nspecies;

  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        for (int ig = 0; ig < groupspec; ++ig) {
          const int i4D = layout_cc.linear(p.i, p.j, p.k, ig);
          rN_rhs_transport[i4D] = rN_rhs[i4D];
          rE_rhs_transport[i4D] = rE_rhs[i4D];
          rFx_rhs_transport[i4D] = rFx_rhs[i4D];
          rFy_rhs_transport[i4D] = rFy_rhs[i4D];
          rFz_rhs_transport[i4D] = rFz_rhs[i4D];
        }
      });
}

} // namespace nuX_M1
