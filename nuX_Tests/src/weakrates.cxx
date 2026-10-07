#include <cctk.h>
#include <cctk_Arguments.h>

#include <algorithm>
#include <cmath>

#include "nuX_weakrates.hxx"

namespace nuX_WeakRates {
namespace {

bool close_rel(const CCTK_REAL x, const CCTK_REAL y,
               const CCTK_REAL tolerance = 2.0e-12) {
  const CCTK_REAL scale = std::max(std::abs(x), std::abs(y));
  return scale == 0.0 || std::abs(x - y) <= tolerance * scale;
}

int check_rates(const WeakRates &weakrates, const EOSState &eos) {
  const auto rates = weakrates.compute_rates(eos);
  const auto equilibrium = weakrates.equilibrium_densities(eos);
  int failures = 0;
  for (int species = 0; species < nspecies; ++species) {
    const CCTK_REAL values[] = {
        rates.eta_0[species],     rates.eta[species],
        rates.kappa_0_a[species], rates.kappa_a[species],
        rates.kappa_0_s[species], rates.kappa_s[species]};
    for (const CCTK_REAL value : values) {
      failures += !std::isfinite(value) || value < 0.0;
    }
    failures += !close_rel(rates.eta_0[species],
                           detail::clight * rates.kappa_0_a[species] *
                               equilibrium.number[species]);
    failures +=
        !close_rel(rates.eta[species], detail::clight * rates.kappa_a[species] *
                                           equilibrium.energy[species]);
  }
  return failures;
}

int check_legacy_nonbeta_behavior(const WeakRates &weakrates,
                                  const EOSState &eos) {
  const auto rates = weakrates.compute_rates(eos);
  int failures = 0;
  for (int species = 0; species < nspecies; ++species) {
    failures +=
        !std::isfinite(rates.eta_0[species]) || rates.eta_0[species] <= 0.0;
    failures += !std::isfinite(rates.eta[species]) || rates.eta[species] <= 0.0;
    // Pair/plasmon emissivities were present in the legacy library, but no
    // inverse non-beta absorption was supplied.
    failures += rates.kappa_0_a[species] != 0.0;
    failures += rates.kappa_a[species] != 0.0;
  }
  return failures;
}

int check_zero_rates(const WeakRates &weakrates, const EOSState &eos) {
  const auto rates = weakrates.compute_rates(eos);
  const auto equilibrium = weakrates.equilibrium_densities(eos);
  int failures = 0;
  for (int species = 0; species < nspecies; ++species) {
    failures += rates.eta_0[species] != 0.0;
    failures += rates.eta[species] != 0.0;
    failures += rates.kappa_0_a[species] != 0.0;
    failures += rates.kappa_a[species] != 0.0;
    failures += rates.kappa_0_s[species] != 0.0;
    failures += rates.kappa_s[species] != 0.0;
    failures += equilibrium.number[species] != 0.0;
    failures += equilibrium.energy[species] != 0.0;
  }
  return failures;
}

} // namespace

extern "C" void nuX_WeakRates_TestRates(CCTK_ARGUMENTS) {
  const EOSState eos = {1.0e12, // g cm^-3
                        10.0,   // MeV
                        0.1,    930.18478708, 20.0, 900.0, 910.0, 0.0,
                        0.0,    0.9,          0.1,  0.0,   0.0};

  // A locally constructed object must match the public legacy defaults.  The
  // production singleton obtains the same values from param.ccl in init().
  const WeakRates legacy_defaults{};
  int failures = 0;
  failures += !legacy_defaults.include_beta;
  failures += !legacy_defaults.include_pair;
  failures += !legacy_defaults.include_plasmon;
  failures += legacy_defaults.include_nonbeta_absorption;
  failures += legacy_defaults.include_bremsstrahlung;
  failures += !legacy_defaults.include_elastic;
  failures += legacy_defaults.beta_low_density_threshold != 2.0e11;

  WeakRates weakrates{};
  weakrates.include_beta = false;
  weakrates.include_pair = true;
  weakrates.include_plasmon = true;
  weakrates.include_nonbeta_absorption = false;
  weakrates.include_bremsstrahlung = false;
  weakrates.include_elastic = false;
  weakrates.beta_low_density_threshold = 2.0e11;

  // Default/legacy mode keeps the raw pair and plasmon emissivities but does
  // not synthesize inverse rates.
  failures += check_legacy_nonbeta_behavior(weakrates, eos);

  // The correction is independently testable when explicitly enabled.
  weakrates.include_nonbeta_absorption = true;
  failures += check_rates(weakrates, eos);
  weakrates.include_pair = false;
  weakrates.include_plasmon = false;
  weakrates.include_bremsstrahlung = true;
  failures += check_rates(weakrates, eos);

  // The legacy beta factors contain an indeterminate 0/0 form when the
  // nucleon degeneracy and composition differences both vanish.  The
  // non-degenerate limiting form must keep this valid state finite.
  WeakRates beta_rates{};
  beta_rates.include_beta = true;
  beta_rates.include_pair = false;
  beta_rates.include_plasmon = false;
  beta_rates.include_nonbeta_absorption = false;
  beta_rates.include_bremsstrahlung = false;
  beta_rates.include_elastic = true;
  beta_rates.beta_low_density_threshold = 0.0;
  EOSState symmetric = eos;
  symmetric.xn = 0.5;
  symmetric.xp = 0.5;
  symmetric.mu_n = symmetric.mu_p + detail::qnp;
  const auto symmetric_rates = beta_rates.compute_rates(symmetric);
  for (int species = 0; species < nspecies; ++species) {
    failures += !std::isfinite(symmetric_rates.eta_0[species]);
    failures += !std::isfinite(symmetric_rates.eta[species]);
    failures += !std::isfinite(symmetric_rates.kappa_0_a[species]);
    failures += !std::isfinite(symmetric_rates.kappa_a[species]);
    failures += !std::isfinite(symmetric_rates.kappa_0_s[species]);
    failures += !std::isfinite(symmetric_rates.kappa_s[species]);
  }

  // A near-zero eta_hat is not a 0/0 limit when the composition difference is
  // finite.  Preserve the resulting large high-density blocking factors rather
  // than silently substituting the low-density approximation.
  EOSState near_singular = eos;
  near_singular.temp = 1.0;
  near_singular.xn = 0.4;
  near_singular.xp = 0.6;
  near_singular.mu_p = 0.0;
  near_singular.mu_n = detail::qnp - 1.0e-15;
  const auto near_singular_state =
      detail::derive_state(near_singular, CCTK_REAL(0));
  failures += !(near_singular_state.eta_np >
                CCTK_REAL(1.0e3) * near_singular_state.nb);
  failures += !(near_singular_state.eta_pn >
                CCTK_REAL(1.0e3) * near_singular_state.nb);

  // Exercise the corrected option together with the legacy beta and elastic
  // channels.  The added non-beta contribution must be exactly additive;
  // this does not impose a stronger detailed-balance identity on the legacy
  // beta approximation itself.
  WeakRates mixed_rates = beta_rates;
  mixed_rates.include_pair = true;
  mixed_rates.include_plasmon = true;
  mixed_rates.include_nonbeta_absorption = true;
  mixed_rates.include_bremsstrahlung = true;
  WeakRates nonbeta_rates = mixed_rates;
  nonbeta_rates.include_beta = false;
  nonbeta_rates.include_elastic = false;
  const auto mixed_coeffs = mixed_rates.compute_rates(eos);
  const auto beta_coeffs = beta_rates.compute_rates(eos);
  const auto nonbeta_coeffs = nonbeta_rates.compute_rates(eos);
  for (int species = 0; species < nspecies; ++species) {
    failures +=
        !close_rel(mixed_coeffs.eta_0[species],
                   beta_coeffs.eta_0[species] + nonbeta_coeffs.eta_0[species]);
    failures +=
        !close_rel(mixed_coeffs.eta[species],
                   beta_coeffs.eta[species] + nonbeta_coeffs.eta[species]);
    failures += !close_rel(mixed_coeffs.kappa_0_a[species],
                           beta_coeffs.kappa_0_a[species] +
                               nonbeta_coeffs.kappa_0_a[species]);
    failures += !close_rel(mixed_coeffs.kappa_a[species],
                           beta_coeffs.kappa_a[species] +
                               nonbeta_coeffs.kappa_a[species]);
  }

  EOSState invalid = eos;
  invalid.temp = 0.0;
  failures += check_zero_rates(beta_rates, invalid);
  invalid = eos;
  invalid.ye = -0.01;
  failures += check_zero_rates(beta_rates, invalid);
  invalid.ye = 1.01;
  failures += check_zero_rates(beta_rates, invalid);

  if (failures != 0) {
    CCTK_VERROR("WeakRates pair/bremsstrahlung self-test failed %d checks",
                failures);
  }
  CCTK_INFO("WeakRates pair/bremsstrahlung self-test passed");
}

} // namespace nuX_WeakRates
