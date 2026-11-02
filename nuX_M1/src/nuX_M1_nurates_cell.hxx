#ifndef NUX_M1_NURATES_CELL_HXX
#define NUX_M1_NURATES_CELL_HXX

#include "m1_opacities.hpp"
#include "nuX_M1_nurates.hxx"
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
// conversion helpers are used by CalcOpacityNuRates. Invalid
// Newton trials are rejected, not repaired and not made fatal mid-line-search.
struct nurates_cell_rates_t {
  const MyQuadrature *quadrature;
  MyEOSParams matter;
  OpacityFlags flags;
  OpacityParams parameters;
  CCTK_REAL rate_max;

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
    if (!CalculateFermiDistrParamsFromM1(&trial.m1_pars, &matter,
                                         &trial.distr_pars))
      return false;

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
    for (int sp = 0; sp < 3; ++sp) {
      // All coefficients were integrated against this same trial spectrum.
      // Applying an additional LTE/mean-energy correction here would change
      // the collision operator and spoil pair lepton/number identities.
      result[sp] = native_nurates_opacities(raw[sp]);
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
