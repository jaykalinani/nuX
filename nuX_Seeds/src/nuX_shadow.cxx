#include <cmath>

#include <loop_device.hxx>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "setup_eos.hxx"
#include "nuX_volume.hxx"

namespace nuX_Seeds {

using namespace Loop;
using namespace EOSX;
using namespace nuX_Seeds_volume;

// -----------------------------------------------------------------------------
// Main setup routine
// -----------------------------------------------------------------------------

extern "C" void nuX_Seeds_SetupHydroTest_shadow(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_Seeds_SetupHydroTest_shadow;
  DECLARE_CCTK_PARAMETERS;

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

  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones, [=] CCTK_DEVICE(const PointDesc &p) {
        const int ijk = layout_cc.linear(p.i, p.j, p.k);
        // Match THC shadow hydro setup: unit-radius density sphere at rest.
        rho[ijk] = volume_f(1.0, p.x, p.y, p.z, p.dx, p.dy, p.dz);
        eps[ijk] = static_eps;
        velx[ijk] = 0.0;
        vely[ijk] = 0.0;
        velz[ijk] = 0.0;
        Ye[ijk] = static_ye;

        press[ijk] =
            eos_3p_ig->press_from_rho_eps_ye(rho[ijk], eps[ijk], Ye[ijk]);
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

extern "C" void nuX_Seeds_SetupNeutTest_shadow(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_Seeds_SetupNeutTest_shadow;
  DECLARE_CCTK_PARAMETERS;

  const GridDescBaseDevice grid(cctkGH);
  const GF3D2layout layout_cc(cctkGH, {1, 1, 1});

  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones, [=] CCTK_DEVICE(const PointDesc &p) {
        for (int ig = 0; ig < ngroups * nspecies; ++ig) {
          int const i4D = layout_cc.linear(p.i, p.j, p.k, ig);
          rE[i4D] = rN[i4D] = rFx[i4D] = rFy[i4D] = rFz[i4D] = 0.0;
        }
      });
}

extern "C" void nuX_Seeds_ShadowBCs(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_Seeds_ShadowBCs;
  DECLARE_CCTK_PARAMETERS;

  const GridDescBaseDevice grid(cctkGH);
  const GF3D2layout layout_cc(cctkGH, {1, 1, 1});

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

  grid.loop_int_device<1, 1, 1>(
      grid.nghostzones, [=] CCTK_DEVICE(const PointDesc &p) {
        for (int ig = 0; ig < ngroups * nspecies; ++ig) {
          int const i4D = layout_cc.linear(p.i, p.j, p.k, ig);
          if ((p.BI[0] == -1.0) && (nx > 0.5) && (abs(p.y) < beam_radius)) {
            rFx[i4D] = 1.0; // If on -X boundary, flux in +X
            rE[i4D] = 1.0;
            rN[i4D] = 1.0;
            rFy[i4D] = 0.0;
            rFz[i4D] = 0.0;
          }
          if ((p.BI[0] == 1.0) && (nx < -0.5) && (abs(p.y) < beam_radius)) {
            rFx[i4D] = -1.0; // If on +X boundary, flux in -X
            rE[i4D] = 1.0;
            rN[i4D] = 1.0;
            rFy[i4D] = 0.0;
            rFz[i4D] = 0.0;
          }
          if ((p.BI[1] == -1.0) && (ny > 0.5) && (abs(p.x) < beam_radius)) {
            rFy[i4D] = 1.0; // If on -Y boundary, flux in +Y
            rE[i4D] = 1.0;
            rN[i4D] = 1.0;
            rFx[i4D] = 0.0;
            rFz[i4D] = 0.0;
          }
          if ((p.BI[1] == 1.0) && (ny < -0.5) && (abs(p.x) < beam_radius)) {
            rFy[i4D] = -1.0; // If on +Y boundary, flux in -Y
            rE[i4D] = 1.0;
            rN[i4D] = 1.0;
            rFx[i4D] = 0.0;
            rFz[i4D] = 0.0;
          }
        }
      });
}

} // namespace nuX_Seeds
