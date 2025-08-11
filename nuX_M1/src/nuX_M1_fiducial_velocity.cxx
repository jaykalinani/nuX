#include <algorithm>
#include <cstring>
#include <loop_device.hxx>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "aster_utils.hxx"
#include "nuX_utils.hxx"

#define CGS_GCC (1.619100425158886e-18)

namespace nuX_M1 {

using namespace Arith;
using namespace Loop;
using namespace AsterUtils;

extern "C" void nuX_M1_FiducialVelocity(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_M1_FiducialVelocity;
  DECLARE_CCTK_PARAMETERS

  if (verbose) {
    CCTK_INFO("nuX_M1_FiducialVelocity");
  }

  const GridDescBaseDevice grid(cctkGH);
  const GF3D2layout layout2(cctkGH, {1, 1, 1});
  const smat<GF3D2<const CCTK_REAL8>, 3> gf_g{
      GF3D2<const CCTK_REAL8>(layout2, gxx),
      GF3D2<const CCTK_REAL8>(layout2, gxy),
      GF3D2<const CCTK_REAL8>(layout2, gxz),
      GF3D2<const CCTK_REAL8>(layout2, gyy),
      GF3D2<const CCTK_REAL8>(layout2, gyz),
      GF3D2<const CCTK_REAL8>(layout2, gzz)};

  if (CCTK_Equals(fiducial_velocity, "fluid")) {

    grid.loop_all_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          const int ijk = layout2.linear(p.i, p.j, p.k);

          const smat<CCTK_REAL, 3> g_avg([&](int i, int j) ARITH_INLINE {
            return calc_avg_v2c(gf_g(i, j), p);
          });

          vec<CCTK_REAL, 3> v_up;
          vec<CCTK_REAL, 3> v_low;
          CCTK_REAL w_lorentz;

          fidu_velx[ijk] = velx[ijk];
          fidu_vely[ijk] = vely[ijk];
          fidu_velz[ijk] = velz[ijk];

          v_up(0) = velx[ijk];
          v_up(1) = vely[ijk];
          v_up(2) = velz[ijk];

          v_low = calc_contraction(g_avg, v_up);
          w_lorentz = calc_wlorentz(v_low, v_up);

          fidu_w_lorentz[ijk] = w_lorentz;
        });

  } else if (CCTK_Equals(fiducial_velocity, "mixed")) {

    grid.loop_all_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          const int ijk = layout2.linear(p.i, p.j, p.k);

          // Weight between fluid velocity and zero depending on density
          const CCTK_REAL fac =
              1.0 / fmax(dens[ijk], fiducial_velocity_rho_fluid * CGS_GCC);

          fidu_velx[ijk] = velx[ijk] * dens[ijk] * fac;
          fidu_vely[ijk] = vely[ijk] * dens[ijk] * fac;
          fidu_velz[ijk] = velz[ijk] * dens[ijk] * fac;

          // Metric-lowered components
          const CCTK_REAL fidu_vel_x = gxx[ijk] * fidu_velx[ijk] +
                                       gxy[ijk] * fidu_vely[ijk] +
                                       gxz[ijk] * fidu_velz[ijk];
          const CCTK_REAL fidu_vel_y = gxy[ijk] * fidu_velx[ijk] +
                                       gyy[ijk] * fidu_vely[ijk] +
                                       gyz[ijk] * fidu_velz[ijk];
          const CCTK_REAL fidu_vel_z = gxz[ijk] * fidu_velx[ijk] +
                                       gyz[ijk] * fidu_vely[ijk] +
                                       gzz[ijk] * fidu_velz[ijk];

          // Lorentz factor
          const CCTK_REAL v2 = fidu_vel_x * fidu_velx[ijk] +
                               fidu_vel_y * fidu_vely[ijk] +
                               fidu_vel_z * fidu_velz[ijk];

          fidu_w_lorentz[ijk] = 1.0 / sqrt(1.0 - v2);
        });

  } else {

    grid.loop_all_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          const int ijk = layout2.linear(p.i, p.j, p.k);

          fidu_velx[ijk] = CCTK_REAL(0);
          fidu_vely[ijk] = CCTK_REAL(0);
          fidu_velz[ijk] = CCTK_REAL(0);
          fidu_w_lorentz[ijk] = CCTK_REAL(1);
        });
  }
}

} // namespace nuX_M1
