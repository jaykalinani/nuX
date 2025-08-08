#include <cassert>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

namespace nuX_M1 {

extern "C" void nuX_M1_InitVolform(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_nuX_M1_InitVolform;
  DECLARE_CCTK_PARAMETERS;

  if (verbose) {
    CCTK_INFO("nuX_M1_InitVolform");
  }

  {
    grid.loop_all_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          const smat<GF3D2<const CCTK_REAL>, dim> gf_g{gxx, gxy, gxz,
                                                       gyy, gyz, gzz};
          const smat<CCTK_REAL, 3> g_avg([&](int i, int j) ARITH_INLINE {
            return calc_avg_v2c(gf_g(i, j), p, dir);
          });

          /* determinant of spatial metric */
          const CCTK_REAL detg_avg = calc_det(g_avg);
          volform[p.I] = sqrt(detg_avg);
          assert(isfinite(volform[p.I]));
        });
  }
}

} // namespace nuX_M1
