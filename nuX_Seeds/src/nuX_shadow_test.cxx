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

CCTK_DEVICE CCTK_REAL volume(CCTK_REAL R, CCTK_REAL xp, CCTK_REAL yp,
                             CCTK_REAL zp, CCTK_REAL dx, CCTK_REAL dy,
                             CCTK_REAL dz) {
  const int NPOINTS = 10;
  int inside = 0, count = 0;
  for (int ii = 0; ii < NPOINTS; ++ii) {
    CCTK_REAL const myx = (xp - dx / 2.) + (ii + 0.5) * (dx / NPOINTS);
    for (int jj = 0; jj < NPOINTS; ++jj) {
      CCTK_REAL const myy = (yp - dy / 2.) + (jj + 0.5) * (dy / NPOINTS);
      for (int kk = 0; kk < NPOINTS; ++kk) {
        CCTK_REAL const myz = (zp - dz / 2.) + (kk + 0.5) * (dz / NPOINTS);
        ++count;
        if (myx * myx + myy * myy + myz * myz < R * R)
          ++inside;
      }
    }
  }
  return static_cast<CCTK_REAL>(inside) / static_cast<CCTK_REAL>(count);
}


extern "C" void nuX_Seeds_SetupTest_shadow(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_Seeds_SetupTest_shadow;
  DECLARE_CCTK_PARAMETERS;

  if (verbose)
    CCTK_INFO("nuX_Seeds_SetupTest_shadow");

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
     nx = 0.0;
     ny = 0.0;
     nz = 1.0;
  }

  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones, [=] CCTK_DEVICE(const PointDesc &p) {
        const int ijk = layout2.linear(p.i, p.j, p.k);
        for (int ig = 0; ig < ngroups * nspecies; ++ig) {
          int const i4D = layout2.linear(p.i, p.j, p.k, ig);
          rho[ijk] = static_rho*volume(shadow_radius, p.x, p.y, p.z, p.dx, p.dy, p.dz);
          eps[ijk] = static_eps*volume(shadow_radius, p.x, p.y, p.z, p.dx, p.dy, p.dz);
          velx[ijk] = static_velx;
          vely[ijk] = static_vely;
          velz[ijk] = static_velz;
          Ye[ijk]   = static_ye;

          press[ijk] = eos_3p_ig->press_from_valid_rho_eps_ye(
              rho[ijk], eps[ijk], Ye[ijk]);
          rE[i4D] = rN[i4D] = rFx[i4D] = rFy[i4D] = rFz[i4D] = 0.0;       
        }
      });
      grid.loop_int_device<1, 1, 1>(
      grid.nghostzones, [=] CCTK_DEVICE(const PointDesc &p) {
        const int ijk = layout2.linear(p.i, p.j, p.k);
        for (int ig = 0; ig < ngroups * nspecies; ++ig) {
          int const i4D = layout2.linear(p.i, p.j, p.k, ig);
          if ((p.BI[0]== -1.0) && (test_nvec[0] == 1.0) && (abs(p.y) < shadow_radius/2)) {
           rFx[i4D] = 1.0; // If on -X boundary, flux in +X
            rE[i4D] = 1.0;
            rN[i4D] = 1.0;
          }
          if ((p.BI[0]== 1.0) && (test_nvec[0] == -1.0) && (abs(p.y) < shadow_radius/2)) {
           rFx[i4D] = -1.0; // If on +X boundary, flux in -X
            rE[i4D] = 1.0;
            rN[i4D] = 1.0;
          }
          if ((p.BI[1]== -1.0) && (test_nvec[1] == 1.0)&& (abs(p.x) < shadow_radius/2)) {
           rFy[i4D] = 1.0; // If on -Y boundary, flux in +Y
            rE[i4D] = 1.0;
            rN[i4D] = 1.0;
          }
          if ((p.BI[1]== 1.0) && (test_nvec[1] == -1.0)&& (abs(p.x) < shadow_radius/2)) {
           rFy[i4D] = -1.0; // If on +Y boundary, flux in -Y
            rE[i4D] = 1.0;
            rN[i4D] = 1.0;
          }
        }
      });
}

} // namespace nuX_M1
