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

extern "C" void nuX_Seeds_SetupHydroTest_shadow(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_Seeds_SetupHydroTest_shadow;
  DECLARE_CCTK_PARAMETERS;

  if (verbose)
    CCTK_INFO("nuX_Seeds_SetupHydroTest_shadow");

  auto eos_3p_ig = global_eos_3p_ig;
  if (not CCTK_EQUALS(evolution_eos, "IdealGas")) {
    CCTK_VERROR("Invalid evolution EOS type '%s'. Please, set "
                "EOSX::evolution_eos = \"IdealGas\" in your parameter file.",
                evolution_eos);
  }

  const GridDescBaseDevice grid(cctkGH);
  const GF3D2layout layout2(cctkGH, {1, 1, 1});

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
     nx = 1.0;
     ny = 0.0;
     nz = 0.0;
  }

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

extern "C" void nuX_Seeds_SetupNeutTest_shadow(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_Seeds_SetupNeutTest_shadow;
  DECLARE_CCTK_PARAMETERS;

  if (verbose)
    CCTK_INFO("nuX_Seeds_SetupNeutTest_shadow");

  const GridDescBaseDevice grid(cctkGH);
  const GF3D2layout layout2(cctkGH, {1, 1, 1});

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
     nx = 1.0;
     ny = 0.0;
     nz = 0.0;
  }

  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones, [=] CCTK_DEVICE(const PointDesc &p) {
        const int ijk = layout2.linear(p.i, p.j, p.k);
        for (int ig = 0; ig < ngroups * nspecies; ++ig) {
          int const i4D = layout2.linear(p.i, p.j, p.k, ig);
          rE[i4D] = rN[i4D] = rFx[i4D] = rFy[i4D] = rFz[i4D] = 0.0;       
        }
      });
   grid.loop_int_device<1, 1, 1>(
      grid.nghostzones, [=] CCTK_DEVICE(const PointDesc &p) {
        const int ijk = layout2.linear(p.i, p.j, p.k);
        for (int ig = 0; ig < ngroups * nspecies; ++ig) {
          int const i4D = layout2.linear(p.i, p.j, p.k, ig);
          if ((p.BI[0]== -1.0) && (nx == 1.0) && (abs(p.y) < beam_radius)) {
           rFx[i4D] = 1.0; // If on -X boundary, flux in +X
            rE[i4D] = 1.0;
            rN[i4D] = 1.0;
          }
          if ((p.BI[0]== 1.0) && (nx == -1.0) && (abs(p.y) < beam_radius)) {
           rFx[i4D] = -1.0; // If on +X boundary, flux in -X
            rE[i4D] = 1.0;
            rN[i4D] = 1.0;
          }
          if ((p.BI[1]== -1.0) && (ny == 1.0)&& (abs(p.x) < beam_radius)) {
           rFy[i4D] = 1.0; // If on -Y boundary, flux in +Y
            rE[i4D] = 1.0;
            rN[i4D] = 1.0;
          }
          if ((p.BI[1]== 1.0) && (ny == -1.0)&& (abs(p.x) < beam_radius)) {
           rFy[i4D] = -1.0; // If on +Y boundary, flux in -Y
            rE[i4D] = 1.0;
            rN[i4D] = 1.0;
          }
        }
      });
}

} // namespace nuX_Seeds
