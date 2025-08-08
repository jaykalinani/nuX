#include <algorithm>   // harmless to keep
#include <cstring>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include <loop_device.hxx>
#include "nuX_utils.hxx"

#define CGS_GCC (1.619100425158886e-18)

namespace nuX_M1 {

using namespace Loop;

extern "C" void nuX_M1_FiducialVelocity(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_M1_FiducialVelocity;
  DECLARE_CCTK_PARAMETERS

  if (verbose) {
    CCTK_INFO("nuX_M1_FiducialVelocity");
  }

  const GridDescBaseDevice grid(cctkGH);
  const GF3D2layout        layout2(cctkGH, {1, 1, 1});

  if (CCTK_Equals(fiducial_velocity, "fluid")) {

    grid.loop_all_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          const int ijk = layout2.linear(p.i, p.j, p.k);

          fidu_velx[ijk]      = velx[ijk];
          fidu_vely[ijk]      = vely[ijk];
          fidu_velz[ijk]      = velz[ijk];
          fidu_w_lorentz[ijk] = w_lorentz[ijk];
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

          fidu_velx[ijk]      = CCTK_REAL(0);
          fidu_vely[ijk]      = CCTK_REAL(0);
          fidu_velz[ijk]      = CCTK_REAL(0);
          fidu_w_lorentz[ijk] = CCTK_REAL(1);
        });
  }
}

} // namespace nuX_M1

