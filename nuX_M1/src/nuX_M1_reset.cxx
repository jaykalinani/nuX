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

  if (verbose && CCTK_MyProc(cctkGH) == 0) {
    CCTK_INFO("nuX_M1_Reset");
  }

  const GridDescBaseDevice grid(cctkGH);
  const GF3D2layout layout_cc(cctkGH, {1, 1, 1});

  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        const int ijk = layout_cc.linear(p.i, p.j, p.k);
        nuX_m1_mask[ijk] = 0.0;

        for (int ig = 0; ig < ngroups * nspecies; ++ig) {
          int const i4D = layout_cc.linear(p.i, p.j, p.k, ig);
          rE[i4D] = rad_E_floor;
          rN[i4D] = rad_N_floor;
          rFx[i4D] = 0.0;
          rFy[i4D] = 0.0;
          rFz[i4D] = 0.0;

          rE_p[i4D] = rad_E_floor;
          rN_p[i4D] = rad_N_floor;
          rFx_p[i4D] = 0.0;
          rFy_p[i4D] = 0.0;
          rFz_p[i4D] = 0.0;

          abs_0[i4D] = 0.0;
          abs_1[i4D] = 0.0;
          eta_0[i4D] = 0.0;
          eta_1[i4D] = 0.0;
          scat_1[i4D] = 0.0;
          nueave[i4D] = 0.0;
        }
      });
}

} // namespace nuX_M1
