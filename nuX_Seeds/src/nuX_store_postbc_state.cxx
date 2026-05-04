#include <loop_device.hxx>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

namespace nuX_Seeds {

using namespace Loop;

extern "C" void nuX_Seeds_StorePostBCState(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_Seeds_StorePostBCState;
  DECLARE_CCTK_PARAMETERS;

  const GridDescBaseDevice grid(cctkGH);
  const GF3D2layout layout_cc(cctkGH, {1, 1, 1});
  const int groupspec = ngroups * nspecies;

  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        for (int ig = 0; ig < groupspec; ++ig) {
          const int i4D = layout_cc.linear(p.i, p.j, p.k, ig);
          rN_bc_snapshot[i4D] = rN[i4D];
          rE_bc_snapshot[i4D] = rE[i4D];
          rFx_bc_snapshot[i4D] = rFx[i4D];
          rFy_bc_snapshot[i4D] = rFy[i4D];
          rFz_bc_snapshot[i4D] = rFz[i4D];
        }
      });
}

} // namespace nuX_Seeds
