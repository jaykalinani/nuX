#include <cmath>

#include <loop_device.hxx>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "nuX_seed_utils.hxx"
#include "setup_eos.hxx"

namespace nuX_M1_Seeds {

using namespace Loop;
using namespace EOSX;

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
normalize_test_nvec(CCTK_REAL *nx, CCTK_REAL *ny, CCTK_REAL *nz) {
  CCTK_REAL const n2 = (*nx) * (*nx) + (*ny) * (*ny) + (*nz) * (*nz);
  if (n2 > 0) {
    CCTK_REAL const nn = sqrt(n2);
    *nx /= nn;
    *ny /= nn;
    *nz /= nn;
  } else {
    *nx = 1.0;
    *ny = 0.0;
    *nz = 0.0;
  }
}

extern "C" void nuX_M1_Seeds_SetupHydroTest_beam(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_M1_Seeds_SetupHydroTest_beam;
  DECLARE_CCTK_PARAMETERS;

  if (verbose)
    CCTK_INFO("nuX_M1_Seeds_SetupHydroTest_beam");

  auto eos_3p_ig = global_eos_3p_ig;
  const bool tabulated = CCTK_EQUALS(evolution_eos, "Tabulated3d");
  if (!tabulated && !CCTK_EQUALS(evolution_eos, "IdealGas")) {
    CCTK_VERROR("Invalid evolution EOS type '%s'. Please, set "
                "EOSX::evolution_eos to IdealGas or Tabulated3d.",
                evolution_eos);
  }
  CCTK_REAL tab_eps = static_eps, tab_temp = 0, tab_press = 0;
  if (tabulated) {
    const auto eos = global_eos_3p_tab3d;
    if (!eos || !std::isfinite(static_rho) || !std::isfinite(static_eps) ||
        !std::isfinite(static_ye) ||
        static_rho < eos->rgrho.min || static_rho > eos->rgrho.max ||
        static_ye < eos->rgye.min || static_ye > eos->rgye.max)
      CCTK_ERROR("Uniform beam matter is outside the initialized EOSX table");
    const auto eps_range = eos->range_eps_from_rho_ye(static_rho, static_ye);
    if (static_eps < eps_range.min || static_eps > eps_range.max)
      CCTK_ERROR("Uniform beam static_eps is outside the EOSX range at static_rho/static_ye");
    // The matter state is uniform: perform the tabulated EOS inversion once.
    tab_temp = eos->temp_from_rho_eps_ye(static_rho, tab_eps, static_ye);
    tab_press = eos->press_from_rho_temp_ye(static_rho, tab_temp, static_ye);
    if (!std::isfinite(tab_temp) || !std::isfinite(tab_press))
      CCTK_ERROR("Invalid tabulated EOSX state for the uniform beam matter");
  }

  const GridDescBaseDevice grid(cctkGH);
  const GF3D2layout layout_cc(cctkGH, {1, 1, 1});

  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        const int ijk = layout_cc.linear(p.i, p.j, p.k);
        velx[ijk] = static_velx;
        vely[ijk] = static_vely;
        velz[ijk] = static_velz;
        rho[ijk] = static_rho;
        eps[ijk] = static_eps;
        Ye[ijk] = static_ye;
        if (tabulated) {
          eps[ijk] = tab_eps;
          temperature[ijk] = tab_temp;
          press[ijk] = tab_press;
        } else {
          press[ijk] =
              eos_3p_ig->press_from_rho_eps_ye(rho[ijk], eps[ijk], Ye[ijk]);
          temperature[ijk] =
              eos_3p_ig->temp_from_rho_eps_ye(rho[ijk], eps[ijk], Ye[ijk]);
        }
      });

}

extern "C" void nuX_M1_Seeds_SetupNeutTest_beam(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_M1_Seeds_SetupNeutTest_beam;
  DECLARE_CCTK_PARAMETERS;

  if (verbose)
    CCTK_INFO("nuX_M1_Seeds_SetupNeutTest_beam");

  const GridDescBaseDevice grid(cctkGH);
  const GF3D2layout layout_cc(cctkGH, {1, 1, 1});
  int const ncomponents = radiation_components();

  CCTK_REAL nx = test_nvec[0];
  CCTK_REAL ny = test_nvec[1];
  CCTK_REAL nz = test_nvec[2];
  normalize_test_nvec(&nx, &ny, &nz);

  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        CCTK_REAL const proj = nx * p.x + ny * p.y + nz * p.z;
        CCTK_REAL const offset2 = (p.x - nx * p.x) * (p.x - nx * p.x) +
                                  (p.y - ny * p.y) * (p.y - ny * p.y) +
                                  (p.z - nz * p.z) * (p.z - nz * p.z);
        for (int ig = 0; ig < ncomponents; ++ig) {
          int const i4D = layout_cc.linear(p.i, p.j, p.k, ig);
          if (proj < beam_position && offset2 < beam_radius * beam_radius) {
            rE[i4D] = 1.0;
            rN[i4D] = 1.0;
            rFx[i4D] = nx;
            rFy[i4D] = ny;
            rFz[i4D] = nz;
          } else {
            rE[i4D] = 0.0;
            rN[i4D] = 0.0;
            rFx[i4D] = 0.0;
            rFy[i4D] = 0.0;
            rFz[i4D] = 0.0;
          }
        }
      });
}

} // namespace nuX_M1_Seeds
