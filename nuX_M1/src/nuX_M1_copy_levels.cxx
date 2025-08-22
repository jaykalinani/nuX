#include <cstring>
#include <loop_device.hxx>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

namespace nuX_M1 {

using namespace Loop;

extern "C" void nuX_M1_CopyLevels(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_M1_CopyLevels;
  DECLARE_CCTK_PARAMETERS

  if (verbose) {
    CCTK_INFO("nuX_M1_CopyLevels");
  }

  const GridDescBaseDevice grid(cctkGH);
  const GF3D2layout layout2(cctkGH, {1, 1, 1});

  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        int const ijk = layout2.linear(p.i, p.j, p.k);

        // Loop over groups/species.
        // int groupspec = ngroups * nspecies;

        for (int ig = 0; ig < 3; ++ig) {
          int const i4D = layout2.linear(p.i, p.j, p.k, ig);

          rN_p[i4D] = rN[i4D];
          rE_p[i4D] = rE[i4D];
          rFx_p[i4D] = rFx[i4D];
          rFy_p[i4D] = rFy[i4D];
          rFz_p[i4D] = rFz[i4D];
          //printf("Setting prev timelevel to %e, from initial %e\n", rN_p[i4D], rN[i4D]);
        }
      });
}

} // namespace nuX_M1
