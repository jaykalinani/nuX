#include <cmath>
#include <loop_device.hxx>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "nuX_M1_closure.hxx"
#include "nuX_utils.hxx"

namespace nuX_M1 {

using namespace Loop;

extern "C" void nuX_M1_RepairState(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_M1_RepairState;
  DECLARE_CCTK_PARAMETERS;

  const GridDescBaseDevice grid(cctkGH);
  const GF3D2layout layout_cc(cctkGH, {1, 1, 1});
  const GF3D2layout layout_vc(cctkGH, {0, 0, 0});
  tensor::slicing_geometry_const geom(layout_vc, layout_cc, alp, betax, betay,
                                      betaz, gxx, gxy, gxz, gyy, gyz, gzz, kxx,
                                      kxy, kxz, kyy, kyz, kzz);

  // ODESolvers forms intermediate and final states through tableau linear
  // combinations. Repair the stored radiation state once after each such
  // update so closure and transport never operate on a different temporary
  // state. This is a numerical invariant projection, not a collision solve,
  // and it does not update matter.
  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        tensor::inv_metric<4> g_uu;
        tensor::generic<CCTK_REAL, 4, 1> beta_u;
        geom.get_inv_metric(p, &g_uu);
        geom.get_shift_vec(p, &beta_u);

        for (int ig = 0; ig < nspecies * ngroups; ++ig) {
          const int i4D = layout_cc.linear(p.i, p.j, p.k, ig);
          CCTK_REAL E = rE[i4D];
          tensor::generic<CCTK_REAL, 4, 1> F_d;
          pack_F_d(beta_u(1), beta_u(2), beta_u(3), rFx[i4D], rFy[i4D],
                   rFz[i4D], &F_d);
          repair_moments(g_uu, &E, &F_d, rad_E_floor, rad_eps);

          rN[i4D] = isfinite(rN[i4D]) ? fmax(rN[i4D], rad_N_floor)
                                      : rad_N_floor;
          rE[i4D] = E;
          unpack_F_d(F_d, &rFx[i4D], &rFy[i4D], &rFz[i4D]);
        }
      });
}

} // namespace nuX_M1
