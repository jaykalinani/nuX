#include <cmath>
#include <cassert>
#include <loop_device.hxx>
#include <mat.hxx>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "setup_eos.hxx"
#include "aster_utils.hxx"
#include "nuX_volume.hxx"

namespace nuX_Seeds {

using namespace Loop;
using namespace EOSX;
using namespace AsterUtils;
using namespace nuX_Seeds_volume;

// -----------------------------------------------------------------------------
// Main setup routine
// -----------------------------------------------------------------------------

extern "C" void nuX_Seeds_SetupHydroTest_sphere(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_Seeds_SetupHydroTest_sphere;
  DECLARE_CCTK_PARAMETERS;

  if (verbose)
    CCTK_INFO("nuX_Seeds_SetupHydroTest_sphere");

  auto eos_3p_ig = global_eos_3p_ig;
  if (not CCTK_EQUALS(evolution_eos, "IdealGas")) {
    CCTK_VERROR("Invalid evolution EOS type '%s'. Please, set "
                "EOSX::evolution_eos = \"IdealGas\" in your parameter file.",
                evolution_eos);
  }

  const GridDescBaseDevice grid(cctkGH);
  const GF3D2layout layout_cc(cctkGH, {1, 1, 1});
  const GF3D2layout layout3(cctkGH, {1, 0, 0});
  const GF3D2layout layout4(cctkGH, {0, 1, 0});
  const GF3D2layout layout5(cctkGH, {0, 0, 1});
  grid.loop_all_device<1, 1, 1>(grid.nghostzones, [=] CCTK_DEVICE(
                                                      const PointDesc &p) {
    const int ijk = layout_cc.linear(p.i, p.j, p.k);
    for (int ig = 0; ig < ngroups * nspecies; ++ig) {
      int const i4D = layout_cc.linear(p.i, p.j, p.k, ig);
      // Match THC sphere hydro setup: unit-radius density sphere at rest.
      rho[ijk] = volume_f(1.0, p.x, p.y, p.z, p.dx, p.dy, p.dz);
      eps[ijk] = static_eps;
      velx[ijk] = 0.0;
      vely[ijk] = 0.0;
      velz[ijk] = 0.0;
      Ye[ijk] = static_ye;

      press[ijk] =
          eos_3p_ig->press_from_rho_eps_ye(rho[ijk], eps[ijk], Ye[ijk]);
    }
  });

  grid.loop_all_device<1, 0, 0>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        const int ijk = layout3.linear(p.i, p.j, p.k);
        Avec_x[ijk] = 0.;
      });

  grid.loop_all_device<0, 1, 0>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        const int ijk = layout4.linear(p.i, p.j, p.k);
        Avec_y[ijk] = 0.;
      });

  grid.loop_all_device<0, 0, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        const int ijk = layout5.linear(p.i, p.j, p.k);
        Avec_z[ijk] = 0.;
      });
}

extern "C" void nuX_Seeds_SetupNeutTest_sphere(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_Seeds_SetupNeutTest_sphere;
  DECLARE_CCTK_PARAMETERS;

  if (verbose)
    CCTK_INFO("nuX_Seeds_SetupNeutTest_sphere");

  const GridDescBaseDevice grid(cctkGH);
  const GF3D2layout layout_cc(cctkGH, {1, 1, 1});
  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones, [=] CCTK_DEVICE(const PointDesc &p) {
        const int ijk = layout_cc.linear(p.i, p.j, p.k);
        for (int ig = 0; ig < ngroups * nspecies; ++ig) {
          int const i4D = layout_cc.linear(p.i, p.j, p.k, ig);
          rE[i4D] = rN[i4D] = rFx[i4D] = rFy[i4D] = rFz[i4D] = 0.0;
        }
      });
}

} // namespace nuX_Seeds
