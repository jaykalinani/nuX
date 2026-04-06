#include <loop_device.hxx>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "m1_opacities.hpp"
#include "nuX_LeakageBase_rates.hxx"
#include "nuX_fakerates.hxx"
#include "setup_eos.hxx"

namespace nuX_LeakageBase {

using namespace Loop;
using namespace nuX_FakeRates;
using namespace nuX_Rates;
using namespace EOSX;

extern "C" void nuX_LeakageBase_CalcOpacity(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_nuX_LeakageBase_CalcOpacity;
  DECLARE_CCTK_PARAMETERS;

  if (verbose && CCTK_MyProc(cctkGH) == 0) {
    CCTK_INFO("nuX_LeakageBase_CalcOpacity");
  }

  if (use_fakerates) {
    FakeRatesDef *myfakerates = global_fakerates;
    if (!myfakerates) {
      CCTK_ERROR(
          "nuX_LeakageBase_CalcOpacity requires nuX_FakeRates when "
          "nuX_LeakageBase::use_fakerates = yes");
    }

    grid.loop_all_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          M1Opacities const coeffs =
              myfakerates->ComputeFakeOpacities(rho(p.I));

          kappa_0_nue(p.I) =
              coeffs.kappa_0_a[id_nue] + coeffs.kappa_s[id_nue];
          kappa_0_nua(p.I) =
              coeffs.kappa_0_a[id_anue] + coeffs.kappa_s[id_anue];
          kappa_0_nux(p.I) =
              coeffs.kappa_0_a[id_nux] + coeffs.kappa_s[id_nux];

          kappa_1_nue(p.I) = coeffs.kappa_a[id_nue] + coeffs.kappa_s[id_nue];
          kappa_1_nua(p.I) =
              coeffs.kappa_a[id_anue] + coeffs.kappa_s[id_anue];
          kappa_1_nux(p.I) = coeffs.kappa_a[id_nux] + coeffs.kappa_s[id_nux];
        });
    return;
  }

  auto eos_3p = global_eos_3p_tab3d;
  if (!eos_3p) {
    CCTK_ERROR("nuX_LeakageBase_CalcOpacity requires the EOSX tabulated EOS");
  }

  MyQuadrature my_quad = {.type = kGauleg,
                          .alpha = -42.,
                          .dim = 1,
                          .nx = 10,
                          .ny = 1,
                          .nz = 1,
                          .x1 = 0.,
                          .x2 = 1.,
                          .y1 = -42.,
                          .y2 = -42.,
                          .z1 = -42.,
                          .z2 = -42.,
                          .points = {0},
                          .w = {0}};
  GaussLegendre(&my_quad);

  const OpacityFlags opacity_flags = global_opac_flags;
  const OpacityParams opacity_pars = global_opac_params;

  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        GreyOpacityParams grey_opacity_params;
        setup_equilibrium_grey_opacity_params(
            grey_opacity_params, opacity_flags, opacity_pars, eos_3p, rho(p.I),
            temperature(p.I), Ye(p.I), particle_mass);

        MyQuadrature gpu_quad;
        copy_quadrature(gpu_quad, my_quad);

        const M1Opacities coeffs =
            ComputeM1Opacities(&gpu_quad, &gpu_quad, &grey_opacity_params);

        kappa_0_nue(p.I) = coeffs.kappa_0_a[id_nue] * nuX_length_conv;
        kappa_0_nua(p.I) = coeffs.kappa_0_a[id_anue] * nuX_length_conv;
        kappa_0_nux(p.I) = coeffs.kappa_0_a[id_nux] * nuX_length_conv;

        kappa_1_nue(p.I) =
            (coeffs.kappa_a[id_nue] + coeffs.kappa_s[id_nue]) * nuX_length_conv;
        kappa_1_nua(p.I) =
            (coeffs.kappa_a[id_anue] + coeffs.kappa_s[id_anue]) *
            nuX_length_conv;
        kappa_1_nux(p.I) =
            (coeffs.kappa_a[id_nux] + coeffs.kappa_s[id_nux]) * nuX_length_conv;
      });
}

} // namespace nuX_LeakageBase
