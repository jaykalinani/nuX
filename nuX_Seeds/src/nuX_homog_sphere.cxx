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

extern "C" void nuX_Seeds_SetupHydroTest_homog_sphere(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_Seeds_SetupHydroTest_homog_sphere;
  DECLARE_CCTK_PARAMETERS;

  if (verbose)
    CCTK_INFO("nuX_Seeds_SetupHydroTest_homog_sphere");

  auto eos_3p_ig = global_eos_3p_ig;
  if (not CCTK_EQUALS(evolution_eos, "IdealGas")) {
    CCTK_VERROR("Invalid evolution EOS type '%s'. Please, set "
                "EOSX::evolution_eos = \"IdealGas\" in your parameter file.",
                evolution_eos);
  }

  const GridDescBaseDevice grid(cctkGH);
  const GF3D2layout layout2(cctkGH, {1, 1, 1});
  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones, [=] CCTK_DEVICE(const PointDesc &p) {
        const int ijk = layout2.linear(p.i, p.j, p.k);
        for (int ig = 0; ig < ngroups * nspecies; ++ig) {
          int const i4D = layout2.linear(p.i, p.j, p.k, ig);
          rho[ijk] = static_rho*volume_f(roi_radius, p.x, p.y, p.z, p.dx, p.dy, p.dz);
          eps[ijk] = static_eps*volume_f(roi_radius, p.x, p.y, p.z, p.dx, p.dy, p.dz);
          velx[ijk] = static_velx;
          vely[ijk] = static_vely;
          velz[ijk] = static_velz;
          Ye[ijk]   = static_ye;

          press[ijk] = eos_3p_ig->press_from_valid_rho_eps_ye(
              rho[ijk], eps[ijk], Ye[ijk]);
        }
      });
}

extern "C" void nuX_Seeds_SetupNeutTest_homog_sphere(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_Seeds_SetupNeutTest_homog_sphere;
  DECLARE_CCTK_PARAMETERS;

  if (verbose)
    CCTK_INFO("nuX_Seeds_SetupNeutTest_homog_sphere");

  const GridDescBaseDevice grid(cctkGH);
  const GF3D2layout layout2(cctkGH, {1, 1, 1});
  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones, [=] CCTK_DEVICE(const PointDesc &p) {
        const int ijk = layout2.linear(p.i, p.j, p.k);
        for (int ig = 0; ig < ngroups * nspecies; ++ig) {
          int const i4D = layout2.linear(p.i, p.j, p.k, ig);
          rE[i4D] = rN[i4D] = rFx[i4D] = rFy[i4D] = rFz[i4D] = 0.0;       
        }
      });
}

} // namespace nuX_Seeds
