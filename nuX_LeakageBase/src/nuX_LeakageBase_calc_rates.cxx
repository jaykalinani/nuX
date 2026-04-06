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

extern "C" void nuX_LeakageBase_Rates(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_nuX_LeakageBase_Rates;
  DECLARE_CCTK_PARAMETERS;

  if (verbose && CCTK_MyProc(cctkGH) == 0) {
    CCTK_INFO("nuX_LeakageBase_Rates");
  }

  if (use_fakerates) {
    FakeRatesDef *myfakerates = global_fakerates;
    if (!myfakerates) {
      CCTK_ERROR("nuX_LeakageBase_Rates requires nuX_FakeRates when "
                 "nuX_LeakageBase::use_fakerates = yes");
    }

    grid.loop_all_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          M1Opacities const coeffs =
              myfakerates->ComputeFakeOpacities(rho(p.I));

          CCTK_REAL num_nue, num_nua, num_nux;
          CCTK_REAL ene_nue, ene_nua, ene_nux;
          myfakerates->FakeNeutrinoDens(rho(p.I), num_nue, num_nua, num_nux,
                                        ene_nue, ene_nua, ene_nux);

          CCTK_REAL const r_free_nue = coeffs.eta_0[id_nue];
          CCTK_REAL const r_free_nua = coeffs.eta_0[id_anue];
          CCTK_REAL const r_free_nux = coeffs.eta_0[id_nux];
          CCTK_REAL const q_free_nue = coeffs.eta[id_nue];
          CCTK_REAL const q_free_nua = coeffs.eta[id_anue];
          CCTK_REAL const q_free_nux = coeffs.eta[id_nux];

          R_free_nue(p.I) = store_free_rates ? r_free_nue : 0.0;
          R_free_nua(p.I) = store_free_rates ? r_free_nua : 0.0;
          R_free_nux(p.I) = store_free_rates ? r_free_nux : 0.0;
          Q_free_nue(p.I) = store_free_rates ? q_free_nue : 0.0;
          Q_free_nua(p.I) = store_free_rates ? q_free_nua : 0.0;
          Q_free_nux(p.I) = store_free_rates ? q_free_nux : 0.0;

          R_eff_nue(p.I) =
              calc_eff_rate(r_free_nue, num_nue, kappa_0_nue(p.I),
                            optd_0_nue(p.I), DiffFact);
          R_eff_nua(p.I) =
              calc_eff_rate(r_free_nua, num_nua, kappa_0_nua(p.I),
                            optd_0_nua(p.I), DiffFact);
          R_eff_nux(p.I) =
              calc_eff_rate(r_free_nux, num_nux, kappa_0_nux(p.I),
                            optd_0_nux(p.I), DiffFact);
          Q_eff_nue(p.I) =
              calc_eff_rate(q_free_nue, ene_nue, kappa_1_nue(p.I),
                            optd_1_nue(p.I), DiffFact);
          Q_eff_nua(p.I) =
              calc_eff_rate(q_free_nua, ene_nua, kappa_1_nua(p.I),
                            optd_1_nua(p.I), DiffFact);
          Q_eff_nux(p.I) =
              calc_eff_rate(q_free_nux, ene_nux, kappa_1_nux(p.I),
                            optd_1_nux(p.I), DiffFact);
        });
    return;
  }

  auto eos_3p = global_eos_3p_tab3d;
  if (!eos_3p) {
    CCTK_ERROR("nuX_LeakageBase_Rates requires the EOSX tabulated EOS");
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
        constexpr CCTK_REAL heavy_species_factor = 4.0;

        GreyOpacityParams grey_opacity_params;
        setup_equilibrium_grey_opacity_params(
            grey_opacity_params, opacity_flags, opacity_pars, eos_3p, rho(p.I),
            temperature(p.I), Ye(p.I), particle_mass);

        MyQuadrature gpu_quad;
        copy_quadrature(gpu_quad, my_quad);

        const M1Opacities coeffs =
            ComputeM1Opacities(&gpu_quad, &gpu_quad, &grey_opacity_params);

        const CCTK_REAL r_free_nue =
            convert_number_rate(coeffs.eta_0[id_nue], 1.0);
        const CCTK_REAL r_free_nua =
            convert_number_rate(coeffs.eta_0[id_anue], 1.0);
        const CCTK_REAL r_free_nux =
            convert_number_rate(coeffs.eta_0[id_nux], heavy_species_factor);
        const CCTK_REAL q_free_nue =
            convert_energy_rate(coeffs.eta[id_nue], 1.0);
        const CCTK_REAL q_free_nua =
            convert_energy_rate(coeffs.eta[id_anue], 1.0);
        const CCTK_REAL q_free_nux =
            convert_energy_rate(coeffs.eta[id_nux], heavy_species_factor);

        const CCTK_REAL num_nue =
            convert_number_density(grey_opacity_params.m1_pars.n[id_nue], 1.0);
        const CCTK_REAL num_nua =
            convert_number_density(grey_opacity_params.m1_pars.n[id_anue], 1.0);
        const CCTK_REAL num_nux = convert_number_density(
            grey_opacity_params.m1_pars.n[id_nux], heavy_species_factor);
        const CCTK_REAL ene_nue =
            convert_energy_density(grey_opacity_params.m1_pars.J[id_nue], 1.0);
        const CCTK_REAL ene_nua =
            convert_energy_density(grey_opacity_params.m1_pars.J[id_anue], 1.0);
        const CCTK_REAL ene_nux = convert_energy_density(
            grey_opacity_params.m1_pars.J[id_nux], heavy_species_factor);

        R_free_nue(p.I) = store_free_rates ? r_free_nue : 0.0;
        R_free_nua(p.I) = store_free_rates ? r_free_nua : 0.0;
        R_free_nux(p.I) = store_free_rates ? r_free_nux : 0.0;
        Q_free_nue(p.I) = store_free_rates ? q_free_nue : 0.0;
        Q_free_nua(p.I) = store_free_rates ? q_free_nua : 0.0;
        Q_free_nux(p.I) = store_free_rates ? q_free_nux : 0.0;

        R_eff_nue(p.I) =
            calc_eff_rate(r_free_nue, num_nue, kappa_0_nue(p.I),
                          optd_0_nue(p.I), DiffFact);
        R_eff_nua(p.I) =
            calc_eff_rate(r_free_nua, num_nua, kappa_0_nua(p.I),
                          optd_0_nua(p.I), DiffFact);
        R_eff_nux(p.I) =
            calc_eff_rate(r_free_nux, num_nux, kappa_0_nux(p.I),
                          optd_0_nux(p.I), DiffFact);
        Q_eff_nue(p.I) =
            calc_eff_rate(q_free_nue, ene_nue, kappa_1_nue(p.I),
                          optd_1_nue(p.I), DiffFact);
        Q_eff_nua(p.I) =
            calc_eff_rate(q_free_nua, ene_nua, kappa_1_nua(p.I),
                          optd_1_nua(p.I), DiffFact);
        Q_eff_nux(p.I) =
            calc_eff_rate(q_free_nux, ene_nux, kappa_1_nux(p.I),
                          optd_1_nux(p.I), DiffFact);
      });
}

extern "C" void nuX_LeakageBase_NoAbsorption(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_nuX_LeakageBase_NoAbsorption;
  DECLARE_CCTK_PARAMETERS;

  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        abs_number(p.I) = 0.0;
        abs_energy(p.I) = 0.0;
      });
}

} // namespace nuX_LeakageBase
