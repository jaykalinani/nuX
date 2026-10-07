#include <loop_device.hxx>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "m1_opacities.hpp"
#include "nuX_LeakageBase_rates.hxx"
#include "setup_eos.hxx"

namespace nuX_LeakageBase {

using namespace Loop;
using namespace nuX_NuRates;

extern "C" void nuX_LeakageBase_CalcOpacityNuRates(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_LeakageBase_CalcOpacityNuRates;
  DECLARE_CCTK_PARAMETERS;

  if (verbose && CCTK_MyProc(cctkGH) == 0) {
    CCTK_INFO("nuX_LeakageBase_CalcOpacityNuRates");
  }

  auto eos_3p = EOSX::global_eos_3p_tab3d;
  if (!eos_3p) {
    CCTK_ERROR("nuX_LeakageBase_CalcOpacityNuRates requires a tabulated EOS");
  }

  MyQuadrature quadrature = {.type = kGauleg,
                             .alpha = -42.0,
                             .dim = 1,
                             .nx = 10,
                             .ny = 1,
                             .nz = 1,
                             .x1 = 0.0,
                             .x2 = 1.0,
                             .y1 = -42.0,
                             .y2 = -42.0,
                             .z1 = -42.0,
                             .z2 = -42.0,
                             .points = {0},
                             .w = {0}};
  GaussLegendre(&quadrature);
  const OpacityFlags opacity_flags = global_opac_flags;
  const OpacityParams opacity_params = global_opac_params;

  const GridDescBaseDevice grid(cctkGH);
  const GF3D2layout layout_cc(cctkGH, {1, 1, 1});
  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        const int ijk = layout_cc.linear(p.i, p.j, p.k);
        GreyOpacityParams grey_params;
        if (!setup_equilibrium_grey_opacity_params(
                grey_params, opacity_flags, opacity_params, eos_3p, rho[ijk],
                temperature[ijk], Ye[ijk], particle_mass)) {
          kappa_0_nue[ijk] = kappa_0_nua[ijk] = kappa_0_nux[ijk] = 0.0;
          kappa_1_nue[ijk] = kappa_1_nua[ijk] = kappa_1_nux[ijk] = 0.0;
          abs_0_nue[ijk] = abs_0_nua[ijk] = abs_0_nux[ijk] = 0.0;
          return;
        }

        MyQuadrature device_quadrature;
        copy_quadrature(device_quadrature, quadrature);
        const auto coeffs = ComputeM1Opacities(
            &device_quadrature, &device_quadrature, &grey_params);

        kappa_0_nue[ijk] =
            (coeffs.kappa_0_a[id_nue] + coeffs.kappa_s[id_nue]) *
            nuX_length_conv;
        kappa_0_nua[ijk] =
            (coeffs.kappa_0_a[id_anue] + coeffs.kappa_s[id_anue]) *
            nuX_length_conv;
        kappa_0_nux[ijk] =
            (coeffs.kappa_0_a[id_nux] + coeffs.kappa_s[id_nux]) *
            nuX_length_conv;
        kappa_1_nue[ijk] =
            (coeffs.kappa_a[id_nue] + coeffs.kappa_s[id_nue]) *
            nuX_length_conv;
        kappa_1_nua[ijk] =
            (coeffs.kappa_a[id_anue] + coeffs.kappa_s[id_anue]) *
            nuX_length_conv;
        kappa_1_nux[ijk] =
            (coeffs.kappa_a[id_nux] + coeffs.kappa_s[id_nux]) *
            nuX_length_conv;
        abs_0_nue[ijk] = coeffs.kappa_0_a[id_nue] * nuX_length_conv;
        abs_0_nua[ijk] = coeffs.kappa_0_a[id_anue] * nuX_length_conv;
        abs_0_nux[ijk] = coeffs.kappa_0_a[id_nux] * nuX_length_conv;
      });
}

extern "C" void nuX_LeakageBase_RatesNuRates(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_LeakageBase_RatesNuRates;
  DECLARE_CCTK_PARAMETERS;

  if (verbose && CCTK_MyProc(cctkGH) == 0) {
    CCTK_INFO("nuX_LeakageBase_RatesNuRates");
  }

  auto eos_3p = EOSX::global_eos_3p_tab3d;
  if (!eos_3p) {
    CCTK_ERROR("nuX_LeakageBase_RatesNuRates requires a tabulated EOS");
  }

  MyQuadrature quadrature = {.type = kGauleg,
                             .alpha = -42.0,
                             .dim = 1,
                             .nx = 10,
                             .ny = 1,
                             .nz = 1,
                             .x1 = 0.0,
                             .x2 = 1.0,
                             .y1 = -42.0,
                             .y2 = -42.0,
                             .z1 = -42.0,
                             .z2 = -42.0,
                             .points = {0},
                             .w = {0}};
  GaussLegendre(&quadrature);
  const OpacityFlags opacity_flags = global_opac_flags;
  const OpacityParams opacity_params = global_opac_params;

  const GridDescBaseDevice grid(cctkGH);
  const GF3D2layout layout_cc(cctkGH, {1, 1, 1});
  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        const int ijk = layout_cc.linear(p.i, p.j, p.k);
        constexpr CCTK_REAL heavy_species_factor = 4.0;
        GreyOpacityParams grey_params;
        if (!setup_equilibrium_grey_opacity_params(
                grey_params, opacity_flags, opacity_params, eos_3p, rho[ijk],
                temperature[ijk], Ye[ijk], particle_mass)) {
          R_free_nue[ijk] = R_free_nua[ijk] = R_free_nux[ijk] = 0.0;
          Q_free_nue[ijk] = Q_free_nua[ijk] = Q_free_nux[ijk] = 0.0;
          R_eff_nue[ijk] = R_eff_nua[ijk] = R_eff_nux[ijk] = 0.0;
          Q_eff_nue[ijk] = Q_eff_nua[ijk] = Q_eff_nux[ijk] = 0.0;
          return;
        }

        MyQuadrature device_quadrature;
        copy_quadrature(device_quadrature, quadrature);
        const auto coeffs = ComputeM1Opacities(
            &device_quadrature, &device_quadrature, &grey_params);

        const CCTK_REAL r_free[3] = {
            convert_number_rate(coeffs.eta_0[id_nue], 1.0),
            convert_number_rate(coeffs.eta_0[id_anue], 1.0),
            convert_number_rate(coeffs.eta_0[id_nux], heavy_species_factor)};
        const CCTK_REAL q_free[3] = {
            convert_energy_rate(coeffs.eta[id_nue], 1.0),
            convert_energy_rate(coeffs.eta[id_anue], 1.0),
            convert_energy_rate(coeffs.eta[id_nux], heavy_species_factor)};
        const CCTK_REAL number[3] = {
            convert_number_density(grey_params.m1_pars.n[id_nue], 1.0),
            convert_number_density(grey_params.m1_pars.n[id_anue], 1.0),
            convert_number_density(grey_params.m1_pars.n[id_nux],
                                   heavy_species_factor)};
        const CCTK_REAL energy[3] = {
            convert_energy_density(grey_params.m1_pars.J[id_nue], 1.0),
            convert_energy_density(grey_params.m1_pars.J[id_anue], 1.0),
            convert_energy_density(grey_params.m1_pars.J[id_nux],
                                   heavy_species_factor)};
        const CCTK_REAL kappa_0[3] = {
            kappa_0_nue[ijk], kappa_0_nua[ijk], kappa_0_nux[ijk]};
        const CCTK_REAL kappa_1[3] = {
            kappa_1_nue[ijk], kappa_1_nua[ijk], kappa_1_nux[ijk]};
        const CCTK_REAL tau_0[3] = {
            optd_0_nue[ijk], optd_0_nua[ijk], optd_0_nux[ijk]};
        const CCTK_REAL tau_1[3] = {
            optd_1_nue[ijk], optd_1_nua[ijk], optd_1_nux[ijk]};

        R_free_nue[ijk] = store_free_rates ? r_free[0] : 0.0;
        R_free_nua[ijk] = store_free_rates ? r_free[1] : 0.0;
        R_free_nux[ijk] = store_free_rates ? r_free[2] : 0.0;
        Q_free_nue[ijk] = store_free_rates ? q_free[0] : 0.0;
        Q_free_nua[ijk] = store_free_rates ? q_free[1] : 0.0;
        Q_free_nux[ijk] = store_free_rates ? q_free[2] : 0.0;

        R_eff_nue[ijk] =
            calc_eff_rate(r_free[0], number[0], kappa_0[0], tau_0[0], DiffFact);
        R_eff_nua[ijk] =
            calc_eff_rate(r_free[1], number[1], kappa_0[1], tau_0[1], DiffFact);
        R_eff_nux[ijk] =
            calc_eff_rate(r_free[2], number[2], kappa_0[2], tau_0[2], DiffFact);
        Q_eff_nue[ijk] =
            calc_eff_rate(q_free[0], energy[0], kappa_1[0], tau_1[0], DiffFact);
        Q_eff_nua[ijk] =
            calc_eff_rate(q_free[1], energy[1], kappa_1[1], tau_1[1], DiffFact);
        Q_eff_nux[ijk] =
            calc_eff_rate(q_free[2], energy[2], kappa_1[2], tau_1[2], DiffFact);
      });
}

} // namespace nuX_LeakageBase
