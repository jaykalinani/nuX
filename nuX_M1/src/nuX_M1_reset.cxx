#include <cstring>
#include <loop_device.hxx>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "nuX_utils.hxx"

namespace nuX_M1 {

using namespace Loop;

extern "C" void nuX_M1_Reset(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_M1_Reset;
  DECLARE_CCTK_PARAMETERS

  if (verbose) {
    CCTK_INFO("nuX_M1_Reset");
  }

  const GridDescBaseDevice grid(cctkGH);
  const GF3D2layout layout2(cctkGH, {1, 1, 1});

  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        const int ijk = layout2.linear(p.i, p.j, p.k);
        nuX_m1_mask[ijk] = 0.0;

        for (int ig = 0; ig < ngroups * nspecies; ++ig) {
          int const i4D = layout2.linear(p.i, p.j, p.k, ig);
          rE[i4D] = rad_E_floor;
          rN[i4D] = rad_N_floor;
          rFx[i4D] = 0.0;
          rFy[i4D] = 0.0;
          rFz[i4D] = 0.0;
        }
      });
}

} // namespace nuX_M1
