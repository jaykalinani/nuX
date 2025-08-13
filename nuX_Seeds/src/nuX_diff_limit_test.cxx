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
extern "C" void nuX_Seeds_SetupTest_diff_limit_test(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_Seeds_SetupTest_diff_limit_test;
  DECLARE_CCTK_PARAMETERS;

  if (verbose)
    CCTK_INFO("nuX_Seeds_SetupTest_diff_limit_test");

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

  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones, [=] CCTK_DEVICE(const PointDesc &p) {
        const int ijk = layout2.linear(p.i, p.j, p.k);
        for (int ig = 0; ig < ngroups * nspecies; ++ig) {
          int const i4D = layout2.linear(p.i, p.j, p.k, ig);
          
          // set the velocity to zero in param file
          velx[ijk] = static_velx;
          vely[ijk] = static_vely;
          velz[ijk] = static_velz;
          // 
          rho[ijk] = static_rho;
          eps[ijk] = static_eps; 
          Ye[ijk]  = static_ye;
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

          if (CCTK_EQUALS(nuX_test_case,"diff_limit_gaussian")){
            rE[i4D] = exp(-9.0*p.z*p.z);
          } else if (CCTK_EQUALS(nuX_test_case,"diff_limit_square")) {
            rE[i4D] = (p.z < 0.5)*(p.z > -0.5)*static_E;
          }
          rN[i4D] = rE[i4D];

          CCTK_REAL const W = calc_wlorentz(v_low, v_up);
          CCTK_REAL const Jo3 = rE[i4D] / (4 * W * W - 1);
          rFx[i4D] = 4 * W * W * velx[ijk] * Jo3;
          rFy[i4D] = 4 * W * W * vely[ijk] * Jo3;
          rFz[i4D] = 4 * W * W * velz[ijk] * Jo3;
        }
      });
}

} // namespace nuX_M1
