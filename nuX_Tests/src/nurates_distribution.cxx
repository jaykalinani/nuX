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
      nullptr, &quadrature, matter, OpacityFlags{}, OpacityParams{},
      1, 930.18478708, 1, 1, -1, 1, 1, -1, 0, 0};
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
  int failures = check_distribution_inverse() + check_cell_reconstruction(quadrature);
  // The established host-only test parfile activates no CarpetX driver.
  // Keep its scalar checks usable without allocating from an uninitialized
  // AMReX device arena. The CarpetX test additionally exercises its backend.
  if (!amrex::Initialized())
    return failures;
  const amrex::Gpu::LaunchSafeGuard launch_region(true);
  amrex::Gpu::DeviceScalar<int> device_failures(0);
  int *const output = device_failures.dataPtr();
  amrex::ParallelFor(1, [=] CCTK_DEVICE(int) {
    *output = check_distribution_inverse() + check_cell_reconstruction(quadrature);
  });
  amrex::Gpu::synchronize();
  failures += device_failures.dataValue();
  return failures;
}

} // namespace nuX_NuRates
