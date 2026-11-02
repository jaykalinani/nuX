#include <loop_device.hxx>

#include <algorithm>
#include <limits>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "nuX_M1_closure.hxx"
#include "nuX_M1_nurates_cell.hxx"
#include "nuX_utils.hxx"
#include "setup_eos.hxx"

namespace nuX_M1 {
using namespace std;
using namespace Loop;
using namespace nuX_NuRates;
using namespace nuX_Utils;
using namespace EOSX;

// Keep invalid physical states/rates fatal in release GPU builds as well.
CCTK_HOST CCTK_DEVICE inline void
require_nurates(bool const valid, char const *message) {
  if (valid)
    return;
#if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
  printf("%s\n", message);
  device_abort();
#else
  CCTK_ERROR(message);
#endif
}

void CalcOpacityNuRates(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_M1_CalcOpacityNuRates;
  DECLARE_CCTK_PARAMETERS;

  if (verbose)
    CCTK_INFO("nuX_M1_CalcOpacityNuRates");
  if (nspecies != 3 || ngroups != 1)
    CCTK_ERROR("NuRates requires one group and three evolved species");
  const auto eos = global_eos_3p_tab3d;
  if (!eos)
    CCTK_ERROR("nuX_M1_CalcOpacityNuRates requires a tabulated EOS");

  const GridDescBaseDevice grid(cctkGH);
  const GF3D2layout layout_cc(cctkGH, {1, 1, 1});
  const GF3D2layout layout_vc(cctkGH, {0, 0, 0});
  const GF3D2<const CCTK_REAL> gf_gxx(layout_vc, gxx);
  const GF3D2<const CCTK_REAL> gf_gxy(layout_vc, gxy);
  const GF3D2<const CCTK_REAL> gf_gxz(layout_vc, gxz);
  const GF3D2<const CCTK_REAL> gf_gyy(layout_vc, gyy);
  const GF3D2<const CCTK_REAL> gf_gyz(layout_vc, gyz);
  const GF3D2<const CCTK_REAL> gf_gzz(layout_vc, gzz);

  MyQuadrature quadrature{};
  quadrature.type = kGauleg;
  quadrature.dim = 1;
  quadrature.nx = 10;
  quadrature.ny = quadrature.nz = 1;
  quadrature.x1 = 0;
  quadrature.x2 = 1;
  GaussLegendre(&quadrature);
  const auto flags = global_opac_flags;
  const auto parameters = global_opac_params;

  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones, [=] CCTK_DEVICE(const PointDesc &p) {
        const int ijk = layout_cc.linear(p.i, p.j, p.k);
        if (nuX_m1_mask[ijk]) {
          for (int sp = 0; sp < 3; ++sp) {
            const int index = layout_cc.linear(p.i, p.j, p.k, sp);
            abs_0[index] = abs_1[index] = scat_1[index] = 0;
            eta_0[index] = eta_1[index] = nueave[index] = 0;
          }
          return;
        }

        CCTK_REAL rho_cell;
        MyEOSParams matter{};
        require_nurates(
            prepare_nurates_matter(eos, rho[ijk], temperature[ijk], Ye[ijk],
                                    particle_mass, rho_cell, matter),
            "nuX_M1 NuRates: invalid matter/EOS state");
        const CCTK_REAL det = nuX_Utils::metric::spatial_det(
            tensor::interp_v2c(gf_gxx, p), tensor::interp_v2c(gf_gxy, p),
            tensor::interp_v2c(gf_gxz, p), tensor::interp_v2c(gf_gyy, p),
            tensor::interp_v2c(gf_gyz, p), tensor::interp_v2c(gf_gzz, p));
        require_nurates(isfinite(det) && det > 0,
                        "nuX_M1 NuRates: invalid spatial metric");
        const CCTK_REAL volume = sqrt(det);

        // rnnu and rJ are densitized COMOVING moments. The heavy group is a
        // sum of four flavors; the native library expects one member.
        M1Quantities radiation{};
        constexpr CCTK_REAL chi_slack =
            64 * numeric_limits<CCTK_REAL>::epsilon();
        for (int sp = 0; sp < 3; ++sp) {
          const int index = layout_cc.linear(p.i, p.j, p.k, sp);
          const CCTK_REAL flavors = sp == 2 ? 4 : 1;
          radiation.n[sp] = rnnu[index] / (volume * flavors) * nuX_ndens_conv;
          radiation.J[sp] = rJ[index] / (volume * flavors) * nuX_edens_conv;
          require_nurates(
              isfinite(chi[index]) && chi[index] >= CCTK_REAL(1) / 3 - chi_slack &&
                  chi[index] <= CCTK_REAL(1) + chi_slack,
              "nuX_M1 NuRates: invalid comoving closure factor");
          radiation.chi[sp] = min(max(chi[index], CCTK_REAL(1) / 3), CCTK_REAL(1));
        }
        radiation.n[3] = radiation.n[2];
        radiation.J[3] = radiation.J[2];
        radiation.chi[3] = radiation.chi[2];

        // Exactly the same spectrum, conversions and physical rates as the
        // nonlinear trial/final verification path. No floor-based reference
        // substitution or additional LTE grey correction is applied.
        const nurates_cell_rates_t callback{
            &quadrature, matter, flags, parameters, opacity_rate_max};
        nurates_opacities_t rates[3];
        require_nurates(callback(radiation, rates),
                        "nuX_M1 NuRates: invalid spectrum or collision rates");
        for (int sp = 0; sp < 3; ++sp) {
          const int index = layout_cc.linear(p.i, p.j, p.k, sp);
          abs_0[index] = rates[sp].abs_0;
          abs_1[index] = rates[sp].abs_1;
          scat_1[index] = rates[sp].scat_1;
          eta_0[index] = rates[sp].eta_0;
          eta_1[index] = rates[sp].eta_1;
          nueave[index] = rnnu[index] > 0 ? rJ[index] / rnnu[index] : 0;
        }
      });
}

extern "C" void nuX_M1_CalcOpacityNuRates(CCTK_ARGUMENTS) {
  CalcOpacityNuRates(cctkGH);
}

} // namespace nuX_M1
