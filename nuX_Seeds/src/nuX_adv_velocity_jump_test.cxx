#include <cmath>
#include <cassert>
#include <loop_device.hxx>
#include <mat.hxx>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "setup_eos.hxx"
#include "aster_utils.hxx"

namespace nuX_M1 {

using namespace Loop;
using namespace EOSX;
using namespace AsterUtils;

// -----------------------------------------------------------------------------
// Main setup routine
// -----------------------------------------------------------------------------
extern "C" void nuX_Seeds_SetupTest_adv_velocity_jump(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_Seeds_SetupTest_adv_velocity_jump;
  DECLARE_CCTK_PARAMETERS;

  if (verbose)
    CCTK_INFO("nuX_Seeds_SetupTest_adv_velocity_jump");

  auto eos_3p_ig = global_eos_3p_ig;
  if (not CCTK_EQUALS(evolution_eos, "IdealGas")) {
    CCTK_VERROR("Invalid evolution EOS type '%s'. Please, set "
                "EOSX::evolution_eos = \"IdealGas\" in your parameter file.",
                evolution_eos);
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

  CCTK_REAL nx = test_nvec[0];
  CCTK_REAL ny = test_nvec[1];
  CCTK_REAL nz = test_nvec[2];
  CCTK_REAL n2 = nx * nx + ny * ny + nz * nz;

  if (n2 > 0) {
     CCTK_REAL nn = sqrt(n2);
     nx /= nn;
     ny /= nn;
     nz /= nn;
  } else {
     nx = 0.0;
     ny = 0.0;
     nz = 1.0;
  }

  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones, [=] CCTK_DEVICE(const PointDesc &p) {
        const int ijk = layout2.linear(p.i, p.j, p.k);
        for (int ig = 0; ig < ngroups * nspecies; ++ig) {
          int const i4D = layout2.linear(p.i, p.j, p.k, ig);
          CCTK_REAL const dotp3d = nx*p.x + ny*p.y + nz*p.z;
          if (dotp3d < 0.0) {
            velx[ijk] = static_velx;
            vely[ijk] = static_vely;
            velz[ijk] = static_velz;
            rE[i4D] = static_E;
          } else if (dotp3d >= 0.0) {
            velx[ijk] = -static_velx;
            vely[ijk] = -static_vely;
            velz[ijk] = -static_velz;
            rE[i4D] = 0.0;
          }
          rho[ijk] = static_rho;
          eps[ijk]  = static_eps;
          Ye[ijk]   = static_ye;
         
          press[ijk] = eos_3p_ig->press_from_valid_rho_eps_ye(
              rho[ijk], eps[ijk], Ye[ijk]);


          const smat<CCTK_REAL, 3> g_avg([&](int i, int j) ARITH_INLINE {
            return calc_avg_v2c(gf_g(i, j), p);
          });

          vec<CCTK_REAL, 3> v_up;
          vec<CCTK_REAL, 3> v_low;

          v_up(0) = velx[ijk];
          v_up(1) = vely[ijk];
          v_up(2) = velz[ijk];

          v_low = calc_contraction(g_avg, v_up);

          rN[i4D] = rE[i4D];
          CCTK_REAL const W = calc_wlorentz(v_low, v_up);
          CCTK_REAL const Jo3 = rE[i4D] / (4 * W * W - 1);
          rFx[i4D] = 4 * W * W * velx[ijk] * Jo3;
          rFy[i4D] = 4 * W * W * vely[ijk] * Jo3;
          rFz[i4D] = 4 * W * W * velz[ijk] * Jo3;
        }
      });
  grid.loop_all_device<1, 0, 0>(
      grid.nghostzones, [=] CCTK_DEVICE(const PointDesc &p)
                            CCTK_ATTRIBUTE_ALWAYS_INLINE { Avec_x(p.I) = 0.; });

  grid.loop_all_device<0, 1, 0>(
      grid.nghostzones, [=] CCTK_DEVICE(const PointDesc &p)
                            CCTK_ATTRIBUTE_ALWAYS_INLINE { Avec_y(p.I) = 0.; });

  grid.loop_all_device<0, 0, 1>(
      grid.nghostzones, [=] CCTK_DEVICE(const PointDesc &p)
                            CCTK_ATTRIBUTE_ALWAYS_INLINE { Avec_z(p.I) = 0.; });
}

} // namespace nuX_M1
