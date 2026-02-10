#include <cmath>
#include <loop_device.hxx>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "setup_eos.hxx"

namespace nuX_Seeds {

using namespace Loop;
using namespace EOSX;

// -----------------------------------------------------------------------------
// Main setup routine
// -----------------------------------------------------------------------------
extern "C" void nuX_Seeds_SetupHydroTest_adv_velocity_jump(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_Seeds_SetupHydroTest_adv_velocity_jump;
  DECLARE_CCTK_PARAMETERS;

  if (verbose)
    CCTK_INFO("nuX_Seeds_SetupHydroTest_adv_velocity_jump");

  auto eos_3p_ig = global_eos_3p_ig;
  if (not CCTK_EQUALS(evolution_eos, "IdealGas")) {
    CCTK_VERROR("Invalid evolution EOS type '%s'. Please, set "
                "EOSX::evolution_eos = \"IdealGas\" in your parameter file.",
                evolution_eos);
  }

  const GridDescBaseDevice grid(cctkGH);
  const GF3D2layout layout2(cctkGH, {1, 1, 1});
  const GF3D2layout layout3(cctkGH, {1, 0, 0});
  const GF3D2layout layout4(cctkGH, {0, 1, 0});
  const GF3D2layout layout5(cctkGH, {0, 0, 1});

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

  grid.loop_all_device<1, 1, 1>(grid.nghostzones, [=] CCTK_DEVICE(
                                                      const PointDesc &p) {
    const int ijk = layout2.linear(p.i, p.j, p.k);
    for (int ig = 0; ig < ngroups * nspecies; ++ig) {
      int const i4D = layout2.linear(p.i, p.j, p.k, ig);
      CCTK_REAL const dotp3d = nx * p.x + ny * p.y + nz * p.z;
      if (dotp3d < 0.0) {
        velx[ijk] = static_velx;
        vely[ijk] = static_vely;
        velz[ijk] = static_velz;
      } else if (dotp3d >= 0.0) {
        velx[ijk] = -static_velx;
        vely[ijk] = -static_vely;
        velz[ijk] = -static_velz;
      }
      rho[ijk] = static_rho;
      eps[ijk] = static_eps;
      Ye[ijk] = static_ye;

      press[ijk] =
          eos_3p_ig->press_from_valid_rho_eps_ye(rho[ijk], eps[ijk], Ye[ijk]);
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

extern "C" void nuX_Seeds_SetupNeutTest_adv_velocity_jump(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_Seeds_SetupNeutTest_adv_velocity_jump;
  DECLARE_CCTK_PARAMETERS;

  if (verbose)
    CCTK_INFO("nuX_Seeds_SetupNeutTest_adv_velocity_jump");

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
        for (int ig = 0; ig < ngroups * nspecies; ++ig) {
          int const i4D = layout2.linear(p.i, p.j, p.k, ig);
          CCTK_REAL const dotp3d = nx * p.x + ny * p.y + nz * p.z;
          if (dotp3d < 0.0) {
            rE[i4D] = static_E;
          } else if (dotp3d >= 0.0) {
            rE[i4D] = 0.0;
          }
          rN[i4D] = rE[i4D];
          rFx[i4D] = rE[i4D] * nx;
          rFy[i4D] = rE[i4D] * ny;
          rFz[i4D] = rE[i4D] * nz;
        }
      });
}

} // namespace nuX_Seeds
