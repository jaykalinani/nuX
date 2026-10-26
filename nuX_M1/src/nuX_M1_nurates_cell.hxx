#ifndef NUX_M1_NURATES_CELL_HXX
#define NUX_M1_NURATES_CELL_HXX

#include "m1_opacities.hpp"
#include "nuX_M1_nurates_equilibrium.hxx"
#include "setup_eos.hxx"

#include <cmath>

namespace nuX_M1 {

CCTK_HOST CCTK_DEVICE inline bool
nurates_rate_valid(const CCTK_REAL value, const CCTK_REAL maximum) {
  return isfinite(value) && value >= 0 && (maximum < 0 || value <= maximum);
}

// EOS interpolation is performed once per cell, not inside the radiation
// iteration. Matter is a fixed background during a local collision solve.
CCTK_HOST CCTK_DEVICE inline bool prepare_nurates_matter(
    const EOSX::eos_3p_tabulated3d *const eos, const CCTK_REAL rho_input,
    const CCTK_REAL temperature_input, const CCTK_REAL ye_input,
    const CCTK_REAL particle_mass, CCTK_REAL &rho, MyEOSParams &matter) {
  if (!isfinite(rho_input) || rho_input <= 0 || !isfinite(temperature_input) ||
      temperature_input <= 0 || !isfinite(ye_input) || ye_input < 0 ||
      ye_input > 1 || !isfinite(particle_mass) || particle_mass <= 0)
    return false;
  rho = fmin(fmax(rho_input, eos->rgrho.min), eos->rgrho.max);
  matter.temp = fmin(fmax(temperature_input, eos->rgtemp.min), eos->rgtemp.max);
  matter.ye = fmin(fmax(ye_input, eos->rgye.min), eos->rgye.max);
  matter.nb = rho * nuX_dens_conv / (particle_mass * kBS_MeVtog);
  using tabulated_eos = EOSX::eos_3p_tabulated3d;
  const auto composition = eos->interptable->interpolate<
      tabulated_eos::EV::MU_P, tabulated_eos::EV::MU_N, tabulated_eos::EV::MU_E,
      tabulated_eos::EV::XN, tabulated_eos::EV::XP>(
      log(rho), log(matter.temp), matter.ye);
  matter.mu_p = composition[0];
  matter.mu_n = composition[1];
  matter.mu_e = composition[2];
  if (!isfinite(matter.nb) || matter.nb <= 0 || !isfinite(matter.mu_p) ||
      !isfinite(matter.mu_n) || !isfinite(matter.mu_e))
    return false;
  matter.yn = isfinite(composition[3])
                  ? fmin(fmax(composition[3], CCTK_REAL(0)), CCTK_REAL(1))
                  : 1 - matter.ye;
  matter.yp = isfinite(composition[4])
                  ? fmin(fmax(composition[4], CCTK_REAL(0)), CCTK_REAL(1))
                  : matter.ye;
  return true;
}

// A lightweight per-cell callable for the local nonlinear solve. The same
// conversion and equilibrium helpers are used by CalcOpacityNuRates. Invalid
// Newton trials are rejected, not repaired and not made fatal mid-line-search.
struct nurates_cell_rates_t {
  const EOSX::eos_3p_tabulated3d *eos;
  const MyQuadrature *quadrature;
  MyEOSParams matter;
  OpacityFlags flags;
  OpacityParams parameters;
  CCTK_REAL rho;
  CCTK_REAL particle_mass;
  CCTK_REAL macro_proper_dt;
  CCTK_REAL volume;
  CCTK_REAL tau_trap;
  CCTK_REAL tau_delta;
  CCTK_REAL correction_max;
  CCTK_REAL rate_max;
  CCTK_REAL number_floor;
  CCTK_REAL energy_floor;

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_NOINLINE bool operator()(
      const M1Quantities &radiation, nurates_opacities_t *const result) const {
    GreyOpacityParams trial{};
    trial.eos_pars = matter;
    trial.opacity_flags = flags;
    trial.opacity_pars = parameters;
    trial.m1_pars = radiation;
    for (int sp = 0; sp < total_num_species; ++sp) {
      if (!equilibrium_moments_are_valid(radiation.n[sp], radiation.J[sp]) ||
          !isfinite(radiation.chi[sp]) || radiation.chi[sp] < CCTK_REAL(1) / 3 ||
          radiation.chi[sp] > 1)
        return false;
    }
    bool floor_only[3];
    nurates_rate_reference(trial.m1_pars, volume, matter.temp,
                            number_floor, energy_floor, floor_only);
    trial.distr_pars = CalculateDistrParamsFromM1(&trial.m1_pars, &matter);
    for (int sp = 0; sp < total_num_species; ++sp) {
      const auto &distribution = trial.distr_pars;
      if (!isfinite(distribution.temp_t[sp]) || distribution.temp_t[sp] <= 0 ||
          !isfinite(distribution.eta_t[sp]) ||
          !isfinite(distribution.temp_f[sp]) || distribution.temp_f[sp] <= 0 ||
          !isfinite(distribution.beta_f[sp]) || distribution.beta_f[sp] < 0 ||
          (distribution.w_t[sp] > 0 && trial.m1_pars.n[sp] > 0 &&
           !ThickDistributionMatchesM1(trial.m1_pars.n[sp], trial.m1_pars.J[sp],
                                       distribution.temp_t[sp],
                                       distribution.eta_t[sp])))
        return false;
    }

    const auto integrated = ComputeM1OpacitiesNonThermalSeparated(
        quadrature, quadrature, &trial);
    nurates_raw_opacities_t raw[3];
    for (int sp = 0; sp < 3; ++sp) {
      raw[sp] = convert_nurates_opacities(integrated, sp, sp == 2 ? 4 : 1);
      const auto &rate = raw[sp];
      if (!nurates_rate_valid(rate.abs_0, rate_max) ||
          !nurates_rate_valid(rate.abs_1_thermal, rate_max) ||
          !nurates_rate_valid(rate.abs_1_nonthermal, rate_max) ||
          !nurates_rate_valid(rate.scat_1, rate_max) ||
          !nurates_rate_valid(rate.eta_0, rate_max) ||
          !nurates_rate_valid(rate.eta_1_thermal, rate_max) ||
          !nurates_rate_valid(rate.eta_1_nonthermal, rate_max))
        return false;
    }
    CCTK_REAL equilibrium_number[3], equilibrium_energy[3];
    nurates_equilibrium_targets(
        eos, rho, matter, radiation, raw, particle_mass, macro_proper_dt,
        tau_trap, tau_delta, equilibrium_number, equilibrium_energy);
    for (int sp = 0; sp < 3; ++sp) {
      const CCTK_REAL factor = sp == 2 ? 4 : 1;
      result[sp] = correct_nurates_opacities(
          raw[sp], sp == 2, floor_only[sp],
          factor * radiation.n[sp] / nuX_ndens_conv,
          factor * radiation.J[sp] / nuX_edens_conv,
          trial.m1_pars.n[sp] > THRESHOLD_N,
          trial.m1_pars.J[sp] > THRESHOLD_J,
          equilibrium_number[sp], equilibrium_energy[sp], correction_max,
          number_floor, energy_floor);
      const auto &rate = result[sp];
      if (!nurates_rate_valid(rate.abs_0, rate_max) ||
          !nurates_rate_valid(rate.abs_1, rate_max) ||
          !nurates_rate_valid(rate.scat_1, rate_max) ||
          !nurates_rate_valid(rate.eta_0, rate_max) ||
          !nurates_rate_valid(rate.eta_1, rate_max))
        return false;
    }
    return true;
  }
};

} // namespace nuX_M1

#endif // NUX_M1_NURATES_CELL_HXX
