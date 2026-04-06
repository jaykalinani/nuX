#include <loop_device.hxx>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "setup_eos.hxx"

namespace nuX_Seeds {

using namespace EOSX;
using namespace Loop;

extern "C" void nuX_Seeds_SetupHydroTest_m0_homogeneous_sphere(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_Seeds_SetupHydroTest_m0_homogeneous_sphere;
  DECLARE_CCTK_PARAMETERS;

  auto eos_3p_ig = global_eos_3p_ig;
  if (!CCTK_EQUALS(evolution_eos, "IdealGas")) {
    CCTK_VERROR("Invalid evolution EOS type '%s'. Please set "
                "EOSX::evolution_eos = \"IdealGas\" for the M0 seed test.",
                evolution_eos);
  }

  const GridDescBaseDevice grid(cctkGH);
  const GF3D2layout layout_cc(cctkGH, {1, 1, 1});
  const GF3D2layout layout_x(cctkGH, {1, 0, 0});
  const GF3D2layout layout_y(cctkGH, {0, 1, 0});
  const GF3D2layout layout_z(cctkGH, {0, 0, 1});
  CCTK_REAL const eps_atmo =
      eos_3p_ig->eps_from_valid_rho_press_ye(rho_atmosphere, press_atmosphere,
                                             m0_static_ye);

  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        int const ijk = layout_cc.linear(p.i, p.j, p.k);
        CCTK_REAL const rp = ::sqrt(p.x * p.x + p.y * p.y + p.z * p.z);
        bool const inside = rp <= m0_test_radius;
        rho[ijk] = inside ? m0_static_rho : rho_atmosphere;
        eps[ijk] = inside ? m0_static_eps : eps_atmo;
        velx[ijk] = 0.0;
        vely[ijk] = 0.0;
        velz[ijk] = 0.0;
        Ye[ijk] = m0_static_ye;
        press[ijk] = inside
                         ? eos_3p_ig->press_from_valid_rho_eps_ye(
                               rho[ijk], eps[ijk], Ye[ijk])
                         : press_atmosphere;
        temperature[ijk] = eos_3p_ig->temp_from_valid_rho_eps_ye(
            rho[ijk], eps[ijk], Ye[ijk]);
        entropy[ijk] =
            eos_3p_ig->kappa_from_valid_rho_eps_ye(rho[ijk], eps[ijk], Ye[ijk]);
      });

  grid.loop_all_device<1, 0, 0>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        Avec_x[layout_x.linear(p.i, p.j, p.k)] = 0.0;
      });
  grid.loop_all_device<0, 1, 0>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        Avec_y[layout_y.linear(p.i, p.j, p.k)] = 0.0;
      });
  grid.loop_all_device<0, 0, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        Avec_z[layout_z.linear(p.i, p.j, p.k)] = 0.0;
      });
}

} // namespace nuX_Seeds
