#include <cassert>
#include <cmath>

#include <loop_device.hxx>
#include <mat.hxx>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "aster_utils.hxx"

namespace nuX_Seeds {

using namespace AsterUtils;
using namespace Loop;

extern "C" void nuX_Seeds_SetupNeutTest_equil(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_Seeds_SetupNeutTest_equil;
  DECLARE_CCTK_PARAMETERS;

  if (verbose)
    CCTK_INFO("nuX_Seeds_SetupNeutTest_equil");

  const GridDescBaseDevice grid(cctkGH);
  const GF3D2layout layout_cc(cctkGH, {1, 1, 1});
  const GF3D2layout layout_vc(cctkGH, {0, 0, 0});
  const smat<GF3D2<const CCTK_REAL8>, 3> gf_g{
      GF3D2<const CCTK_REAL8>(layout_vc, gxx),
      GF3D2<const CCTK_REAL8>(layout_vc, gxy),
      GF3D2<const CCTK_REAL8>(layout_vc, gxz),
      GF3D2<const CCTK_REAL8>(layout_vc, gyy),
      GF3D2<const CCTK_REAL8>(layout_vc, gyz),
      GF3D2<const CCTK_REAL8>(layout_vc, gzz)};

  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        const int ijk = layout_cc.linear(p.i, p.j, p.k);
        assert(ngroups == 1);
        assert(nspecies == 3);

        const smat<CCTK_REAL, 3> g_avg([&](int i, int j) ARITH_INLINE {
          return calc_avg_v2c(gf_g(i, j), p);
        });

        vec<CCTK_REAL, 3> v_up;
        vec<CCTK_REAL, 3> v_low;
        v_up(0) = velx[ijk];
        v_up(1) = vely[ijk];
        v_up(2) = velz[ijk];
        v_low = calc_contraction(g_avg, v_up);

        CCTK_REAL const W = calc_wlorentz(v_low, v_up);
        for (int is = 0; is < nspecies; ++is) {
          int const i4D = layout_cc.linear(p.i, p.j, p.k, is);
          CCTK_REAL const Jnu = equil_nudens_1[is];
          rE[i4D] = (4.0 * W * W - 1.0) / 3.0 * Jnu;
          rFx[i4D] = 4.0 / 3.0 * W * W * velx[ijk] * Jnu;
          rFy[i4D] = 4.0 / 3.0 * W * W * vely[ijk] * Jnu;
          rFz[i4D] = 4.0 / 3.0 * W * W * velz[ijk] * Jnu;
          rN[i4D] = equil_nudens_0[is] * W;
        }
      });
}

} // namespace nuX_Seeds
