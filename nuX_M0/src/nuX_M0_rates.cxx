#include <loop_device.hxx>

#include <cassert>
#include <cmath>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

#include "nuX_LeakageBase_rates.hxx"
#include "nuX_M0_device.hxx"
#include "nuX_M0_schedule.hxx"
#include "nuX_fakerates.hxx"
#include "nuX_rate_units.hxx"
#include "nuX_weakrates.hxx"
#include "setup_eos.hxx"

namespace nuX_M0 {

namespace {

struct RawRates {
  CCTK_REAL abs_0[3];
  CCTK_REAL kappa_0[3];
  CCTK_REAL kappa_1[3];
  CCTK_REAL r_free[3];
  CCTK_REAL q_free[3];
  CCTK_REAL number[3];
  CCTK_REAL energy[3];
};

struct EffectiveRates {
  CCTK_REAL abs_0[3];
  CCTK_REAL r_eff[3];
  CCTK_REAL q_eff[3];
};

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline EffectiveRates
make_effective_rates(const RawRates &raw, const CCTK_REAL tau_0[3],
                     const CCTK_REAL tau_1[3], const CCTK_REAL diff_fact) {
  EffectiveRates rates{};
  for (int isp = 0; isp < 3; ++isp) {
    rates.abs_0[isp] = raw.abs_0[isp];
    rates.r_eff[isp] =
        nuX_LeakageBase::calc_eff_rate(raw.r_free[isp], raw.number[isp],
                                       raw.kappa_0[isp], tau_0[isp], diff_fact);
    rates.q_eff[isp] =
        nuX_LeakageBase::calc_eff_rate(raw.q_free[isp], raw.energy[isp],
                                       raw.kappa_1[isp], tau_1[isp], diff_fact);
    assert(isfinite(rates.abs_0[isp]));
    assert(isfinite(rates.r_eff[isp]));
    assert(isfinite(rates.q_eff[isp]));
  }
  return rates;
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void store_rates(
    const int idx, const EffectiveRates &rates, CCTK_REAL *const abs_0_nue,
    CCTK_REAL *const abs_0_nua, CCTK_REAL *const abs_0_nux,
    CCTK_REAL *const r_nue, CCTK_REAL *const r_nua, CCTK_REAL *const r_nux,
    CCTK_REAL *const q_nue, CCTK_REAL *const q_nua, CCTK_REAL *const q_nux) {
  abs_0_nue[idx] = rates.abs_0[0];
  abs_0_nua[idx] = rates.abs_0[1];
  abs_0_nux[idx] = rates.abs_0[2];
  r_nue[idx] = rates.r_eff[0];
  r_nua[idx] = rates.r_eff[1];
  r_nux[idx] = rates.r_eff[2];
  q_nue[idx] = rates.q_eff[0];
  q_nua[idx] = rates.q_eff[1];
  q_nux[idx] = rates.q_eff[2];
}

int ray_point_count(const cGH *const cctkGH) {
  const int group_id = CCTK_GroupIndex("nuX_M0::nuX_M0_grid_vars");
  cGroupDynamicData group_data;
  const int ierr = CCTK_GroupDynamicData(cctkGH, group_id, &group_data);
  assert(!ierr);
  return group_data.ash[0] * group_data.ash[1];
}

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

  return {rho_eos * nuX_Utils::rate_units::code_density_to_g_nm3 /
              nuX_Utils::rate_units::per_cm3_to_per_nm3,
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

} // namespace

extern "C" void nuX_M0_CalcRatesFakeRates(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_M0_CalcRatesFakeRates;
  DECLARE_CCTK_PARAMETERS;

  if (!CarpetX_AllLevelsSynchronized() ||
      !iteration_is_due(cctk_iteration, compute_every) ||
      load_device_scalar(nuX_M0_is_on) == 0)
    return;

  const auto fakerates = nuX_FakeRates::global_fakerates;
  if (!fakerates)
    CCTK_ERROR("nuX_M0_CalcRatesFakeRates requires nuX_FakeRates");

  const int npts = ray_point_count(cctkGH);
  amrex::ParallelFor(npts, [=] CCTK_DEVICE(const int idx) {
    const auto coeffs = fakerates->ComputeFakeOpacities(nuX_M0_rho[idx]);
    RawRates raw{};
    fakerates->FakeNeutrinoDens(nuX_M0_rho[idx], raw.number[0], raw.number[1],
                                raw.number[2], raw.energy[0], raw.energy[1],
                                raw.energy[2]);
    for (int isp = 0; isp < 3; ++isp) {
      raw.abs_0[isp] = coeffs.kappa_0_a[isp];
      raw.kappa_0[isp] = coeffs.kappa_0_a[isp] + coeffs.kappa_s[isp];
      raw.kappa_1[isp] = coeffs.kappa_a[isp] + coeffs.kappa_s[isp];
      raw.r_free[isp] = coeffs.eta_0[isp];
      raw.q_free[isp] = coeffs.eta[isp];
    }
    const CCTK_REAL tau_0[3] = {nuX_M0_optd_0_nue[idx], nuX_M0_optd_0_nua[idx],
                                nuX_M0_optd_0_nux[idx]};
    const CCTK_REAL tau_1[3] = {nuX_M0_optd_1_nue[idx], nuX_M0_optd_1_nua[idx],
                                nuX_M0_optd_1_nux[idx]};
    const auto rates = make_effective_rates(raw, tau_0, tau_1, DiffFact);
    store_rates(idx, rates, nuX_M0_abs_0_nue, nuX_M0_abs_0_nua,
                nuX_M0_abs_0_nux, nuX_M0_R_nue, nuX_M0_R_nua, nuX_M0_R_nux,
                nuX_M0_Q_nue, nuX_M0_Q_nua, nuX_M0_Q_nux);
  });
}

extern "C" void nuX_M0_CalcRatesWeakRates(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_M0_CalcRatesWeakRates;
  DECLARE_CCTK_PARAMETERS;

  if (!CarpetX_AllLevelsSynchronized() ||
      !iteration_is_due(cctk_iteration, compute_every) ||
      load_device_scalar(nuX_M0_is_on) == 0)
    return;

  const auto eos_3p = EOSX::global_eos_3p_tab3d;
  const auto weakrates = nuX_WeakRates::global_weakrates;
  if (!eos_3p)
    CCTK_ERROR("nuX_M0_CalcRatesWeakRates requires a tabulated EOS");
  if (!weakrates)
    CCTK_ERROR("nuX_M0_CalcRatesWeakRates requires nuX_WeakRates");

  using namespace nuX_Utils;
  const int npts = ray_point_count(cctkGH);
  amrex::ParallelFor(npts, [=] CCTK_DEVICE(const int idx) {
    const auto eos_state =
        make_weakrates_eos(eos_3p, nuX_M0_rho[idx], nuX_M0_temp[idx],
                           nuX_M0_Ye[idx], particle_mass);
    const auto coeffs = weakrates->compute_rates(eos_state);
    const auto equilibrium = weakrates->equilibrium_densities(eos_state);
    RawRates raw{};
    for (int isp = 0; isp < 3; ++isp) {
      raw.abs_0[isp] = coeffs.kappa_0_a[isp] * rate_units::per_cm_to_per_nm *
                       rate_units::code_length_to_nm;
      raw.kappa_0[isp] = (coeffs.kappa_0_a[isp] + coeffs.kappa_0_s[isp]) *
                         rate_units::per_cm_to_per_nm *
                         rate_units::code_length_to_nm;
      raw.kappa_1[isp] = (coeffs.kappa_a[isp] + coeffs.kappa_s[isp]) *
                         rate_units::per_cm_to_per_nm *
                         rate_units::code_length_to_nm;
      raw.r_free[isp] = coeffs.eta_0[isp] * rate_units::per_cm3_to_per_nm3 /
                        rate_units::transport_number_density_to_nm3 *
                        rate_units::code_time_to_s;
      raw.q_free[isp] = coeffs.eta[isp] * rate_units::per_cm3_to_per_nm3 /
                        rate_units::code_energy_density_to_mev_nm3 *
                        rate_units::code_time_to_s;
      raw.number[isp] = equilibrium.number[isp] *
                        rate_units::per_cm3_to_per_nm3 /
                        rate_units::transport_number_density_to_nm3;
      raw.energy[isp] = equilibrium.energy[isp] *
                        rate_units::per_cm3_to_per_nm3 /
                        rate_units::code_energy_density_to_mev_nm3;
    }
    const CCTK_REAL tau_0[3] = {nuX_M0_optd_0_nue[idx], nuX_M0_optd_0_nua[idx],
                                nuX_M0_optd_0_nux[idx]};
    const CCTK_REAL tau_1[3] = {nuX_M0_optd_1_nue[idx], nuX_M0_optd_1_nua[idx],
                                nuX_M0_optd_1_nux[idx]};
    const auto rates = make_effective_rates(raw, tau_0, tau_1, DiffFact);
    store_rates(idx, rates, nuX_M0_abs_0_nue, nuX_M0_abs_0_nua,
                nuX_M0_abs_0_nux, nuX_M0_R_nue, nuX_M0_R_nua, nuX_M0_R_nux,
                nuX_M0_Q_nue, nuX_M0_Q_nua, nuX_M0_Q_nux);
  });
}

extern "C" void nuX_M0_CalcRatesNuRates(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_M0_CalcRatesNuRates;
  DECLARE_CCTK_PARAMETERS;

  if (!CarpetX_AllLevelsSynchronized() ||
      !iteration_is_due(cctk_iteration, compute_every) ||
      load_device_scalar(nuX_M0_is_on) == 0)
    return;

  const auto eos_3p = EOSX::global_eos_3p_tab3d;
  if (!eos_3p)
    CCTK_ERROR("nuX_M0_CalcRatesNuRates requires a tabulated EOS");

  using namespace nuX_NuRates;
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

  const int npts = ray_point_count(cctkGH);
  amrex::ParallelFor(npts, [=] CCTK_DEVICE(const int idx) {
    GreyOpacityParams grey_params;
    if (!nuX_LeakageBase::setup_equilibrium_grey_opacity_params(
            grey_params, opacity_flags, opacity_params, eos_3p,
            nuX_M0_rho[idx], nuX_M0_temp[idx], nuX_M0_Ye[idx],
            particle_mass)) {
      store_rates(idx, EffectiveRates{}, nuX_M0_abs_0_nue, nuX_M0_abs_0_nua,
                  nuX_M0_abs_0_nux, nuX_M0_R_nue, nuX_M0_R_nua,
                  nuX_M0_R_nux, nuX_M0_Q_nue, nuX_M0_Q_nua,
                  nuX_M0_Q_nux);
      return;
    }
    MyQuadrature device_quadrature;
    nuX_LeakageBase::copy_quadrature(device_quadrature, quadrature);
    const auto coeffs = ComputeM1Opacities(&device_quadrature,
                                           &device_quadrature, &grey_params);

    constexpr CCTK_REAL heavy_species_factor = 4.0;
    RawRates raw{};
    for (int isp = 0; isp < 3; ++isp) {
      const CCTK_REAL factor = isp == id_nux ? heavy_species_factor : 1.0;
      raw.abs_0[isp] = coeffs.kappa_0_a[isp] * nuX_length_conv;
      raw.kappa_0[isp] =
          (coeffs.kappa_0_a[isp] + coeffs.kappa_s[isp]) * nuX_length_conv;
      raw.kappa_1[isp] =
          (coeffs.kappa_a[isp] + coeffs.kappa_s[isp]) * nuX_length_conv;
      raw.r_free[isp] =
          nuX_LeakageBase::convert_number_rate(coeffs.eta_0[isp], factor);
      raw.q_free[isp] =
          nuX_LeakageBase::convert_energy_rate(coeffs.eta[isp], factor);
      raw.number[isp] = nuX_LeakageBase::convert_number_density(
          grey_params.m1_pars.n[isp], factor);
      raw.energy[isp] = nuX_LeakageBase::convert_energy_density(
          grey_params.m1_pars.J[isp], factor);
    }
    const CCTK_REAL tau_0[3] = {nuX_M0_optd_0_nue[idx], nuX_M0_optd_0_nua[idx],
                                nuX_M0_optd_0_nux[idx]};
    const CCTK_REAL tau_1[3] = {nuX_M0_optd_1_nue[idx], nuX_M0_optd_1_nua[idx],
                                nuX_M0_optd_1_nux[idx]};
    const auto rates = make_effective_rates(raw, tau_0, tau_1, DiffFact);
    store_rates(idx, rates, nuX_M0_abs_0_nue, nuX_M0_abs_0_nua,
                nuX_M0_abs_0_nux, nuX_M0_R_nue, nuX_M0_R_nua, nuX_M0_R_nux,
                nuX_M0_Q_nue, nuX_M0_Q_nua, nuX_M0_Q_nux);
  });
}

} // namespace nuX_M0
