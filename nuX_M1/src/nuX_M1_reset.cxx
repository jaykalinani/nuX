#include <cstring>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "nuX_utils.hxx"
#include <loop_device.hxx>

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

  // size_t siz = UTILS_GFSIZE(cctkGH) * nspecies * sizeof(CCTK_REAL);
  // size_t isiz = UTILS_GFSIZE(cctkGH) * sizeof(CCTK_INT);

  // UTILS_LOOP3(nuX_m1_analysis, k, 0, cctk_lsh[2], j, 0, cctk_lsh[1], i, 0,
  //            cctk_lsh[0]) {
  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        for (int ig = 0; ig < ngroups * nspecies; ++ig) {
          int const i4D = layout2.linear(p.i, p.j, p.k, ig);
          rE[i4D] = rad_E_floor;
          rN[i4D] = rad_N_floor;
          rFx[i4D] = 0.0;
          rFy[i4D] = 0.0;
          rFz[i4D] = 0.0;
          nuX_m1_mask[i4D] = 0.0;
        }
      });
}

} // namespace nuX_M1
