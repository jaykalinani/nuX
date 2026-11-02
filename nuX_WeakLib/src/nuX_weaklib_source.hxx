#ifndef NUX_WEAKLIB_SOURCE_HXX
#define NUX_WEAKLIB_SOURCE_HXX

#include "nuX_weaklib_collision.hxx"
#include "nuX_M1_spectral_types.hxx"
#include "nuX_rate_units.hxx"

namespace nuX_WeakLib {

// Shared by scheduled opacities and the implicit spectral callback. Selecting
// an isotropic scalar FD spectrum does not modify the M1 angular pressure or
// the comoving chi used by the transport closure.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool radiation_from_moments(const M1Quantities &moments, const MatterState &matter,
                             RadiationState &radiation) noexcept {
  MyEOSParams fixed_matter{};
  fixed_matter.temp = matter.temperature;
  // The reusable inverse carries NuRates' published h*c constant. Convert
  // its moment normalization to this provider's h*c convention explicitly;
  // do not mix the two phase-space prefactors even at sub-ppb differences.
  auto fitted_moments = moments;
  const double normalization = kBS_FourPi_hc3 /
      (phase_space * nuX_Utils::rate_units::per_cm3_to_per_nm3);
  for (int s = 0; s < total_num_species; ++s) {
    fitted_moments.n[s] *= normalization;
    fitted_moments.J[s] *= normalization;
  }
  NuDistributionParams distribution{};
  if (!CalculateFermiDistrParamsFromM1(&fitted_moments, &fixed_matter, &distribution))
    return false;
  for (int s = 0; s < 3; ++s) {
    const double flavors = s == 2 ? 4 : 1;
    radiation.number[s] = flavors * moments.n[s] / nuX_Utils::rate_units::per_cm3_to_per_nm3;
    radiation.energy[s] = flavors * moments.J[s] / nuX_Utils::rate_units::per_cm3_to_per_nm3;
    if (radiation.number[s] == 0) {
      radiation.temperature[s] = matter.temperature;
      radiation.eta[s] = 0; // unused finite vacuum placeholder
    } else {
      radiation.temperature[s] = distribution.temp_t[s];
      radiation.eta[s] = distribution.eta_t[s];
      if (!(std::isfinite(radiation.temperature[s]) && radiation.temperature[s] > 0 &&
            std::isfinite(radiation.eta[s])) ||
          !ThickDistributionMatchesM1(fitted_moments.n[s], fitted_moments.J[s],
              radiation.temperature[s], radiation.eta[s])) return false;
    }
  }
  return true;
}

// Continuous FD moments for the SAME inverse used by the collision callback.
// The legacy equilibrium model keeps its old tabulated-band moments. For the
// monopole model, initializing from a coarse-band integral and then inverting
// continuous FD moments would spuriously move the radiation away from LTE.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool equilibrium_radiation(const MatterState &matter, RadiationState &radiation) noexcept {
  if (!(std::isfinite(matter.temperature) && matter.temperature > 0 &&
        std::isfinite(matter.mu_nu))) return false;
  const double log_temp = std::log(matter.temperature);
  for (int s = 0; s < 3; ++s) {
    const double eta = s == 0 ? matter.mu_nu / matter.temperature :
        s == 1 ? -matter.mu_nu / matter.temperature : 0;
    double log_f2, log_f3;
    if (eta < -20) {
      const double e = std::exp(eta);
      log_f2 = eta + std::log(2.0) + std::log1p(-e / 8 + e * e / 27);
      log_f3 = eta + std::log(6.0) + std::log1p(-e / 16 + e * e / 81);
    } else {
      log_f2 = std::log(FDI_p2(eta));
      log_f3 = std::log(FDI_p3(eta));
    }
    const double log_scale = std::log((s == 2 ? 4 : 1) * phase_space);
    radiation.number[s] = std::exp(log_scale + 3 * log_temp + log_f2);
    radiation.energy[s] = std::exp(log_scale + 4 * log_temp + log_f3);
    radiation.temperature[s] = matter.temperature;
    radiation.eta[s] = eta;
    if (!nonnegative_finite(radiation.number[s]) || !nonnegative_finite(radiation.energy[s]) ||
        ((radiation.number[s] == 0) != (radiation.energy[s] == 0))) return false;
  }
  return true;
}

struct cell_rates_t {
  // Both referents outlive this synchronous per-cell solve. Avoid a private
  // copy of all resident-table axes/pointers in every GPU thread's workspace.
  const TableView &tables;
  const MatterState &matter;

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_NOINLINE
  bool operator()(const M1Quantities &moments,
                  nuX_M1::spectral_opacities_t result[3]) const noexcept {
    RadiationState radiation;
    if (!radiation_from_moments(moments, matter, radiation)) return false;
    Rates rates;
    if (evaluate_collision(tables, matter, radiation, rates) != Status::success)
      return false;
    using namespace nuX_Utils::rate_units;
    constexpr double inverse_cm = per_cm_to_per_nm * code_length_to_nm;
    constexpr double number_units = per_cm3_to_per_nm3 / transport_number_density_to_nm3;
    constexpr double energy_units = per_cm3_to_per_nm3 / code_energy_density_to_mev_nm3;
    for (int s = 0; s < 3; ++s)
      result[s] = {inverse_cm * rates.absorption_number[s],
                   inverse_cm * rates.absorption_energy[s],
                   inverse_cm * rates.scattering[s],
                   inverse_cm * number_units * rates.emission_number[s],
                   inverse_cm * energy_units * rates.emission_energy[s]};
    return true;
  }
};
} // namespace nuX_WeakLib

#endif
