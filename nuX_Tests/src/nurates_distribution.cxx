#include <cmath>
#include <limits>

#include <AMReX.H>
#include <AMReX_GpuControl.H>
#include <AMReX_GpuLaunch.H>
#include <AMReX_GpuMemory.H>

#include "cctk.h"
#include <loop_device.hxx>
#include "nuX_M1_nurates_cell.hxx"

namespace nuX_NuRates {
namespace {

enum class spectrum_component { thick, thin, total };

struct spectrum_moments {
  BS_REAL number{}, energy{};
};

// Independent forward quadrature: no Gamma helper or inverse moment formula.
// E/T_f=u^5 smooths the fractional thin-spectrum power at the lower endpoint.
CCTK_HOST CCTK_DEVICE spectrum_moments integrate_spectrum(
    const NuDistributionParams &distribution, const int sp,
    const spectrum_component component, const int intervals) {
  const BS_REAL step = pow(BS_REAL(80), BS_REAL(0.2)) / intervals;
  spectrum_moments result{};
  // The phase-space measure vanishes at u=0 for each component.
  for (int i = 1; i <= intervals; ++i) {
    const BS_REAL u = i * step;
    const BS_REAL energy = distribution.temp_f[sp] * POW5(u);
    const BS_REAL occupation = component == spectrum_component::thin
        ? NuFThin(energy, &distribution, sp)
        : component == spectrum_component::thick
              ? NuFThick(energy, &distribution, sp)
              : TotalNuF(energy, &distribution, sp);
    const int weight = i == intervals ? 1 : i % 2 ? 4 : 2;
    const BS_REAL measure = weight * 5 * distribution.temp_f[sp] *
                            POW4(u) * kBS_FourPi_hc3 * POW2(energy);
    result.number += measure * occupation;
    result.energy += measure * energy * occupation;
  }
  result.number *= step / 3;
  result.energy *= step / 3;
  return result;
}

CCTK_HOST CCTK_DEVICE int check_distribution_normalization() {
  int failures = 0;
  const BS_REAL temperatures[] = {0.1, 1, 3};
  const BS_REAL degeneracies[] = {-20, 0, 8};
  const BS_REAL eddington_factors[] = {BS_REAL(1) / 3, 0.6, 1};
  const int resolutions[] = {64, 128};
  constexpr BS_REAL tolerance = 2e-12;
  MyEOSParams matter{};
  matter.temp = 1;
  for (const BS_REAL temperature : temperatures) {
    for (const BS_REAL eta : degeneracies) {
      M1Quantities radiation{};
      for (int sp = 0; sp < total_num_species; ++sp) {
        // Distinct per-member spectra, including both heavy representatives.
        const BS_REAL T = (sp + 1) * temperature;
        radiation.n[sp] = kBS_FourPi_hc3 * POW3(T) * FDI_p2(eta);
        radiation.J[sp] = kBS_FourPi_hc3 * POW4(T) * FDI_p3(eta);
      }
      for (const BS_REAL chi : eddington_factors) {
        for (int sp = 0; sp < total_num_species; ++sp)
          radiation.chi[sp] = chi;
        const auto distribution = CalculateDistrParamsFromM1(&radiation, &matter);
        for (int sp = 0; sp < total_num_species; ++sp) {
          failures += distribution.c_f[sp] != BS_REAL(0.6);
          failures += fabs(distribution.w_t[sp] + distribution.w_f[sp] - 1) > tolerance;
          for (const int intervals : resolutions) {
            const auto thin = integrate_spectrum(
                distribution, sp, spectrum_component::thin, intervals);
            failures += !(isfinite(thin.number) &&
                fabs(thin.number / radiation.n[sp] - 1) <= tolerance);
            failures += !(isfinite(thin.energy) &&
                fabs(thin.energy / radiation.J[sp] - 1) <= tolerance);
            const auto thick = integrate_spectrum(
                distribution, sp, spectrum_component::thick, intervals);
            const auto total = integrate_spectrum(
                distribution, sp, spectrum_component::total, intervals);
            // Isolate thin normalization from the retained FD inverse fit.
            // This does not assert an exact forward normalization of that fit
            // or Pauli realizability of the FD/MB mixture.
            const BS_REAL number = distribution.w_t[sp] * thick.number +
                                    distribution.w_f[sp] * radiation.n[sp];
            const BS_REAL energy = distribution.w_t[sp] * thick.energy +
                                    distribution.w_f[sp] * radiation.J[sp];
            failures += !(isfinite(total.number) &&
                fabs(total.number / number - 1) <= tolerance);
            failures += !(isfinite(total.energy) &&
                fabs(total.energy / energy - 1) <= tolerance);
          }
        }
      }
    }
  }
  return failures;
}

CCTK_HOST CCTK_DEVICE int check_distribution_inverse() {
  int failures = 0;
  const BS_REAL degeneracies[] = {-40, -20, -2, 0, 8, 13, 20, 21, 33, 100, 860, 1e5};
  const BS_REAL temperatures[] = {0.1, 1, 3};
  M1Quantities radiation{};
  for (const BS_REAL temperature : temperatures) {
    for (const BS_REAL eta : degeneracies) {
      const BS_REAL n = kBS_FourPi_hc3 * POW3(temperature) * FDI_p2(eta);
      const BS_REAL J = kBS_FourPi_hc3 * POW4(temperature) * FDI_p3(eta);
      for (int sp = 0; sp < total_num_species; ++sp) {
        radiation.n[sp] = n;
        radiation.J[sp] = J;
        radiation.chi[sp] = BS_REAL(1) / 3;
      }
      NuDistributionParams distribution{};
      CalculateThickParamsFromM1(&radiation, &distribution);
      const BS_REAL recovered_T = distribution.temp_t[0];
      const BS_REAL recovered_eta = distribution.eta_t[0];
      failures += !ThickDistributionMatchesM1(n, J, recovered_T, recovered_eta);
      if (eta >= 13) {
        // Near the filled-sphere limit, rounded input moments determine eta
        // less accurately than they determine n/J. Do not confuse conditioning
        // with the former several-fold temperature error at eta=860.
        const BS_REAL tolerance = eta > 1000 ? BS_REAL(1e-4) : BS_REAL(2e-7);
        failures += !isfinite(recovered_T) ||
                    fabs(recovered_T / temperature - 1) > tolerance;
        failures += !isfinite(recovered_eta) ||
                    fabs(recovered_eta / eta - 1) > tolerance;
        const BS_REAL log_n = log(kBS_FourPi_hc3) + 3 * log(recovered_T) +
                              log(FDI_p2(recovered_eta));
        const BS_REAL log_J = log(kBS_FourPi_hc3) + 4 * log(recovered_T) +
                              log(FDI_p3(recovered_eta));
        failures += fabs(log_n - log(n)) > 2e-12;
        failures += fabs(log_J - log(J)) > 2e-12;
      }
    }
  }
  // The retained ordinary fit has a small pre-existing mismatch at its
  // matching point. Check both sides against its stated normalization
  // allowance, without claiming a globally smooth inverse.
  const BS_REAL joins[] = {BS_REAL(0.7) * (1 - 1e-8),
                           BS_REAL(0.7) * (1 + 1e-8)};
  for (const BS_REAL y : joins) {
    const BS_REAL n = kBS_FourPi_hc3;
    const BS_REAL J = kBS_FourPi_hc3 * cbrt(1 / y);
    for (int sp = 0; sp < total_num_species; ++sp) {
      radiation.n[sp] = n;
      radiation.J[sp] = J;
      radiation.chi[sp] = BS_REAL(1) / 3;
    }
    NuDistributionParams recovered{};
    CalculateThickParamsFromM1(&radiation, &recovered);
    failures += !ThickDistributionMatchesM1(n, J, recovered.temp_t[0], recovered.eta_t[0]);
  }

  // This populated state violates the isotropic FD Pauli bound; the inverse
  // must fail, not silently clamp its ratio or choose a maximum eta.
  for (int sp = 0; sp < total_num_species; ++sp) {
    radiation.n[sp] = kBS_FourPi_hc3;
    radiation.J[sp] = kBS_FourPi_hc3;
  }
  NuDistributionParams distribution{};
  CalculateThickParamsFromM1(&radiation, &distribution);
  failures += isfinite(distribution.eta_t[0]) || isfinite(distribution.temp_t[0]);
  // The same n/J is valid input to an UNUSED thick component of a thin
  // spectrum. Its finite placeholder cannot contaminate TotalNuF with 0*NaN.
  for (int sp = 0; sp < total_num_species; ++sp) radiation.chi[sp] = 1;
  CalculateThickParamsFromM1(&radiation, &distribution);
  failures += !isfinite(distribution.eta_t[0]) ||
              !isfinite(distribution.temp_t[0]) || distribution.w_t[0] != 0;
  MyEOSParams matter{};
  matter.temp = 1;
  const auto thin = CalculateDistrParamsFromM1(&radiation, &matter);
  failures += !isfinite(TotalNuF(1, &thin, 0));

  BS_REAL temperature, eta;
  // Exactly filled sphere / an unresolved finite width has no finite-T
  // representation here. Returning failure is explicit, not an eta cap.
  const BS_REAL log_n = log(kBS_FourPi_hc3 / 3);
  const BS_REAL log_J = log(kBS_FourPi_hc3 / 4);
  failures += CalculateDegenerateThickParams(log_n, log_J, &temperature, &eta);
  failures += CalculateDegenerateThickParams(
      std::numeric_limits<BS_REAL>::quiet_NaN(), log_J, &temperature, &eta);

  // Directly represented log moments exercise the existing dilute inverse
  // even though exp(-1000) by itself underflows.
  const BS_REAL hot = 1e100;
  const BS_REAL dilute_n = exp(log(kBS_FourPi_hc3) + 3 * log(hot) + log(BS_REAL(2)) - 1000);
  const BS_REAL dilute_J = exp(log(kBS_FourPi_hc3) + 4 * log(hot) + log(BS_REAL(6)) - 1000);
  for (int sp = 0; sp < total_num_species; ++sp) {
    radiation.n[sp] = dilute_n;
    radiation.J[sp] = dilute_J;
    radiation.chi[sp] = BS_REAL(1) / 3;
  }
  CalculateThickParamsFromM1(&radiation, &distribution);
  failures += !ThickDistributionMatchesM1(dilute_n, dilute_J,
      distribution.temp_t[0], distribution.eta_t[0]);
  for (int sp = 0; sp < total_num_species; ++sp) {
    radiation.n[sp] = radiation.J[sp] = 0;
  }
  const auto vacuum = CalculateDistrParamsFromM1(&radiation, &matter);
  failures += TotalNuF(0, &vacuum, 0) != 0 || TotalNuF(1, &vacuum, 0) != 0;
  return failures;
}

CCTK_HOST CCTK_DEVICE int check_cell_reconstruction(const MyQuadrature &quadrature) {
  MyEOSParams matter{};
  matter.nb = 1e17;
  matter.temp = 1;
  matter.ye = matter.yp = 0.3;
  matter.yn = 0.7;
  // No interactions/trapping: isolate the actual production callback's
  // reconstruction gate, without requiring an EOS or changing any rates.
  const nuX_M1::nurates_cell_rates_t callback{
      &quadrature, matter, OpacityFlags{}, OpacityParams{}, -1};
  int failures = 0;
  const BS_REAL degeneracies[] = {33, 860};
  for (const BS_REAL eta : degeneracies) {
    M1Quantities radiation{};
    for (int sp = 0; sp < total_num_species; ++sp) {
      radiation.n[sp] = kBS_FourPi_hc3 * FDI_p2(eta);
      radiation.J[sp] = kBS_FourPi_hc3 * FDI_p3(eta);
      radiation.chi[sp] = BS_REAL(1) / 3;
    }
    nuX_M1::nurates_opacities_t rates[3];
    if (!callback(radiation, rates)) {
      ++failures;
      continue;
    }
    for (int sp = 0; sp < 3; ++sp)
      failures += rates[sp].abs_0 != 0 || rates[sp].abs_1 != 0 ||
                  rates[sp].eta_0 != 0 || rates[sp].eta_1 != 0 || rates[sp].scat_1 != 0;
  }
  // Positive subnormal populations must not make disabled channels NaN via
  // an overflowing prefactor multiplied by an exactly zero integral.
  M1Quantities tiny{};
  for (int sp = 0; sp < total_num_species; ++sp) {
    tiny.n[sp] = tiny.J[sp] = 1e-320;
    tiny.chi[sp] = BS_REAL(1) / 3;
  }
  nuX_M1::nurates_opacities_t rates[3];
  failures += !callback(tiny, rates);
  for (int sp = 0; sp < 3; ++sp)
    failures += rates[sp].abs_0 != 0 || rates[sp].abs_1 != 0 ||
                rates[sp].eta_0 != 0 || rates[sp].eta_1 != 0 || rates[sp].scat_1 != 0;
  return failures;
}

CCTK_HOST CCTK_DEVICE int check_collision_spectrum() {
  int failures = 0;
  MyEOSParams matter{};
  matter.temp = 0.10944117450587615;
  const BS_REAL degeneracies[] = {-40, -20, -2, 0, 8, 13, 20, 33, 860};
  const BS_REAL closure_factors[] = {BS_REAL(1) / 3, 0.33333550562260805, 1};
  for (const BS_REAL eta : degeneracies) {
    M1Quantities radiation{};
    for (int sp = 0; sp < total_num_species; ++sp) {
      const BS_REAL T = matter.temp * (sp + 1);
      radiation.n[sp] = kBS_FourPi_hc3 * POW3(T) * FDI_p2(eta);
      radiation.J[sp] = kBS_FourPi_hc3 * POW4(T) * FDI_p3(eta);
    }
    for (const BS_REAL chi : closure_factors) {
      for (int sp = 0; sp < total_num_species; ++sp)
        radiation.chi[sp] = chi;
      NuDistributionParams distribution{};
      if (!CalculateFermiDistrParamsFromM1(&radiation, &matter, &distribution)) {
        ++failures;
        continue;
      }
      for (int sp = 0; sp < total_num_species; ++sp) {
        const BS_REAL T = distribution.temp_t[sp];
        const BS_REAL e = distribution.eta_t[sp];
        const BS_REAL log_n = log(kBS_FourPi_hc3) + 3 * log(T) + log(FDI_p2(e));
        const BS_REAL log_J = log(kBS_FourPi_hc3) + 4 * log(T) + log(FDI_p3(e));
        failures += fabs(log_n - log(radiation.n[sp])) > 2e-12;
        failures += fabs(log_J - log(radiation.J[sp])) > 2e-12;
        failures += radiation.chi[sp] != chi; // transport input is unchanged
        failures += distribution.w_t[sp] != 1 || distribution.w_f[sp] != 0;
        // Includes the problematic nearly filled low-energy sea and both
        // sides of the FD edge. Blocking must remain nonnegative without a clamp.
        for (int node = 0; node <= 64; ++node) {
          const BS_REAL energy = T * (fabs(e) + 10) * node / 32;
          const BS_REAL occupation = TotalNuF(energy, &distribution, sp);
          failures += !isfinite(occupation) || occupation < 0 || occupation > 1;
        }
      }
    }
  }
  M1Quantities vacuum{};
  NuDistributionParams distribution{};
  failures += !CalculateFermiDistrParamsFromM1(&vacuum, &matter, &distribution);
  failures += TotalNuF(0, &distribution, 0) != 0 || TotalNuF(1, &distribution, 0) != 0;
  failures += OneMinusTotalNuF(0, &distribution, 0) != 1;
  vacuum.J[0] = 1;
  failures += CalculateFermiDistrParamsFromM1(&vacuum, &matter, &distribution);
  vacuum.n[0] = -1;
  failures += CalculateFermiDistrParamsFromM1(&vacuum, &matter, &distribution);
  distribution.w_t[0] = 1;
  distribution.w_f[0] = 0;
  distribution.temp_t[0] = 1;
  distribution.eta_t[0] = 100;
  failures += TotalNuF(0, &distribution, 0) != 1;
  const BS_REAL hole = OneMinusTotalNuF(0, &distribution, 0);
  failures += !(hole > 0) || fabs(hole / exp(BS_REAL(-100)) - 1) > 1e-14;
  failures += NormalizeM1Opacity(0, 1e-320, kBS_FourPi_hc3) != 0;
  failures += NormalizeM1Opacity(1e-320, 1e-320, kBS_Clight) != 1;
  const BS_REAL normalized = NormalizeM1Opacity(1e200, 1e-200, 1e-200);
  failures += !isfinite(normalized) ||
              fabs(normalized / (1e200 / kBS_Clight) - 1) > 2e-12;
  return failures;
}

} // namespace

int distribution_unit_tests() {
  // Match production: prepare quadrature on the host, capture only its POD
  // data, and use the same nodes in host and device callback checks.
  MyQuadrature quadrature{};
  quadrature.type = kGauleg;
  quadrature.dim = 1;
  quadrature.nx = 6;
  quadrature.ny = quadrature.nz = 1;
  quadrature.x1 = 0;
  quadrature.x2 = 1;
  GaussLegendre(&quadrature);
  int failures = check_distribution_inverse() + check_distribution_normalization() +
                 check_collision_spectrum() + check_cell_reconstruction(quadrature);
  // Scalar helpers do not need AMReX; backend allocation needs its arena.
  // The established unit_tests.par activates CarpetX and exercises both
  // the host checks above and the active-backend checks below.
  if (!amrex::Initialized())
    return failures;
  const amrex::Gpu::LaunchSafeGuard launch_region(true);
  amrex::Gpu::DeviceScalar<int> device_failures(0);
  int *const output = device_failures.dataPtr();
  amrex::ParallelFor(1, [=] CCTK_DEVICE(int) {
    *output = check_distribution_inverse() + check_distribution_normalization() +
              check_collision_spectrum() + check_cell_reconstruction(quadrature);
  });
  amrex::Gpu::synchronize();
  failures += device_failures.dataValue();
  return failures;
}

} // namespace nuX_NuRates
