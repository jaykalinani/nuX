#include <loop_device.hxx>

#include <cmath>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "nuX_LeakageBase_rates.hxx"
#include "nuX_rate_units.hxx"
#include "nuX_weakrates.hxx"
#include "setup_eos.hxx"

namespace nuX_LeakageBase {

using namespace Loop;
using namespace nuX_Utils;

CCTK_HOST
    CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline nuX_WeakRates::EOSState
    make_weakrates_eos(EOSX::eos_3p_tabulated3d *const eos_3p,
                       const CCTK_REAL rho, const CCTK_REAL temperature,
                       const CCTK_REAL ye, const CCTK_REAL particle_mass) {
  if (!isfinite(rho) || rho <= 0.0 || !isfinite(temperature) ||
      temperature <= 0.0 || !isfinite(ye) || ye < 0.0 || ye > 1.0 ||
      !isfinite(particle_mass) || particle_mass <= 0.0) {
    return {};
  }
  using tabulated_eos = EOSX::eos_3p_tabulated3d;
  const CCTK_REAL rho_eos =
      fmin(fmax(rho, eos_3p->rgrho.min), eos_3p->rgrho.max);
  const CCTK_REAL temp_eos =
      fmin(fmax(temperature, eos_3p->rgtemp.min), eos_3p->rgtemp.max);
  const CCTK_REAL ye_eos = fmin(fmax(ye, eos_3p->rgye.min), eos_3p->rgye.max);
  const CCTK_REAL lr = log(rho_eos);
  const CCTK_REAL lt = log(temp_eos);
  const auto eos_state = eos_3p->interptable->interpolate<
      tabulated_eos::EV::MU_E, tabulated_eos::EV::MU_P, tabulated_eos::EV::MU_N,
      tabulated_eos::EV::XA, tabulated_eos::EV::XH, tabulated_eos::EV::XN,
      tabulated_eos::EV::XP, tabulated_eos::EV::ABAR, tabulated_eos::EV::ZBAR>(
      lr, lt, ye_eos);

  return {rho_eos * rate_units::code_density_to_g_nm3 /
              rate_units::per_cm3_to_per_nm3,
          temp_eos,
          ye_eos,
          particle_mass,
          eos_state[0],
          eos_state[1],
          eos_state[2],
          eos_state[3],
          eos_state[4],
          eos_state[5],
          eos_state[6],
          eos_state[7],
          eos_state[8]};
}

extern "C" void nuX_LeakageBase_CalcOpacityWeakRates(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_LeakageBase_CalcOpacityWeakRates;
  DECLARE_CCTK_PARAMETERS;

  if (verbose && CCTK_MyProc(cctkGH) == 0) {
    CCTK_INFO("nuX_LeakageBase_CalcOpacityWeakRates");
  }

  auto eos_3p = EOSX::global_eos_3p_tab3d;
  auto weakrates = nuX_WeakRates::global_weakrates;
  if (!eos_3p) {
    CCTK_ERROR("nuX_LeakageBase_CalcOpacityWeakRates requires a tabulated EOS");
  }
  if (!weakrates) {
    CCTK_ERROR("nuX_LeakageBase_CalcOpacityWeakRates requires nuX_WeakRates");
  }

  const GridDescBaseDevice grid(cctkGH);
  const GF3D2layout layout_cc(cctkGH, {1, 1, 1});
  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        const int ijk = layout_cc.linear(p.i, p.j, p.k);
        const auto weak_eos = make_weakrates_eos(
            eos_3p, rho[ijk], temperature[ijk], Ye[ijk], particle_mass);
        const auto coeffs = weakrates->compute_rates(weak_eos);

        kappa_0_nue[ijk] = (coeffs.kappa_0_a[0] + coeffs.kappa_0_s[0]) *
                           rate_units::per_cm_to_per_nm *
                           rate_units::code_length_to_nm;
        kappa_0_nua[ijk] = (coeffs.kappa_0_a[1] + coeffs.kappa_0_s[1]) *
                           rate_units::per_cm_to_per_nm *
                           rate_units::code_length_to_nm;
        kappa_0_nux[ijk] = (coeffs.kappa_0_a[2] + coeffs.kappa_0_s[2]) *
                           rate_units::per_cm_to_per_nm *
                           rate_units::code_length_to_nm;
        kappa_1_nue[ijk] = (coeffs.kappa_a[0] + coeffs.kappa_s[0]) *
                           rate_units::per_cm_to_per_nm *
                           rate_units::code_length_to_nm;
        kappa_1_nua[ijk] = (coeffs.kappa_a[1] + coeffs.kappa_s[1]) *
                           rate_units::per_cm_to_per_nm *
                           rate_units::code_length_to_nm;
        kappa_1_nux[ijk] = (coeffs.kappa_a[2] + coeffs.kappa_s[2]) *
                           rate_units::per_cm_to_per_nm *
                           rate_units::code_length_to_nm;
        abs_0_nue[ijk] = coeffs.kappa_0_a[0] * rate_units::per_cm_to_per_nm *
                         rate_units::code_length_to_nm;
        abs_0_nua[ijk] = coeffs.kappa_0_a[1] * rate_units::per_cm_to_per_nm *
                         rate_units::code_length_to_nm;
        abs_0_nux[ijk] = coeffs.kappa_0_a[2] * rate_units::per_cm_to_per_nm *
                         rate_units::code_length_to_nm;
      });
}

extern "C" void nuX_LeakageBase_RatesWeakRates(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_LeakageBase_RatesWeakRates;
  DECLARE_CCTK_PARAMETERS;

  if (verbose && CCTK_MyProc(cctkGH) == 0) {
    CCTK_INFO("nuX_LeakageBase_RatesWeakRates");
  }

  auto eos_3p = EOSX::global_eos_3p_tab3d;
  auto weakrates = nuX_WeakRates::global_weakrates;
  if (!eos_3p) {
    CCTK_ERROR("nuX_LeakageBase_RatesWeakRates requires a tabulated EOS");
  }
  if (!weakrates) {
    CCTK_ERROR("nuX_LeakageBase_RatesWeakRates requires nuX_WeakRates");
  }

  const GridDescBaseDevice grid(cctkGH);
  const GF3D2layout layout_cc(cctkGH, {1, 1, 1});
  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        const int ijk = layout_cc.linear(p.i, p.j, p.k);
        const auto weak_eos = make_weakrates_eos(
            eos_3p, rho[ijk], temperature[ijk], Ye[ijk], particle_mass);
        const auto coeffs = weakrates->compute_rates(weak_eos);
        const auto equilibrium = weakrates->equilibrium_densities(weak_eos);

        CCTK_REAL r_free[3], q_free[3], number[3], energy[3];
        for (int isp = 0; isp < 3; ++isp) {
          r_free[isp] = coeffs.eta_0[isp] * rate_units::per_cm3_to_per_nm3 /
                        rate_units::transport_number_density_to_nm3 *
                        rate_units::code_time_to_s;
          q_free[isp] = coeffs.eta[isp] * rate_units::per_cm3_to_per_nm3 /
                        rate_units::code_energy_density_to_mev_nm3 *
                        rate_units::code_time_to_s;
          number[isp] = equilibrium.number[isp] *
                        rate_units::per_cm3_to_per_nm3 /
                        rate_units::transport_number_density_to_nm3;
          energy[isp] = equilibrium.energy[isp] *
                        rate_units::per_cm3_to_per_nm3 /
                        rate_units::code_energy_density_to_mev_nm3;
        }

        const CCTK_REAL kappa_0[3] = {kappa_0_nue[ijk], kappa_0_nua[ijk],
                                      kappa_0_nux[ijk]};
        const CCTK_REAL kappa_1[3] = {kappa_1_nue[ijk], kappa_1_nua[ijk],
                                      kappa_1_nux[ijk]};
        const CCTK_REAL tau_0[3] = {optd_0_nue[ijk], optd_0_nua[ijk],
                                    optd_0_nux[ijk]};
        const CCTK_REAL tau_1[3] = {optd_1_nue[ijk], optd_1_nua[ijk],
                                    optd_1_nux[ijk]};

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
