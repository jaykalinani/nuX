#include <cmath>
#include <vector>

#include <AMReX_GpuContainers.H>
#include <AMReX_GpuControl.H>
#include <AMReX_GpuLaunch.H>
#include <AMReX_GpuMemory.H>

#include "cctk.h"
#include "nuX_M1_closure.hxx"
#include "nuX_weaklib_collision.hxx"
#include "nuX_weaklib_source.hxx"

namespace nuX_WeakLib {
namespace {
void upload(amrex::Gpu::DeviceVector<double> &dest, const std::vector<double> &src) {
  dest.resize(src.size());
  amrex::Gpu::htod_memcpy(dest.data(), src.data(), src.size() * sizeof(double));
}

CCTK_HOST CCTK_DEVICE inline double relative_source(double eta, double loss) {
  return std::abs(eta - loss) / std::fmax(std::abs(eta) + std::abs(loss), 1e-300);
}

CCTK_HOST CCTK_DEVICE RadiationState sampled_fd(const TableView &v, double t,
                                                double eta, bool all_equal = false) {
  RadiationState r;
  for (int s = 0; s < 3; ++s) {
    r.temperature[s] = t;
    r.eta[s] = all_equal || s == 0 ? eta : s == 1 ? -eta : 0;
    const double flavors = s == 2 ? 4 : 1;
    for (int i = 0; i < v.ne; ++i) {
      const double e = v.energies[i], q = v.weights[i] * e * e;
      const double population = flavors * phase_space * q *
          fermi(e / t - r.eta[s]);
      r.number[s] += population;
      r.energy[s] += e * population;
    }
  }
  return r;
}

CCTK_HOST CCTK_DEVICE int check_collision_models(TableView view) {
  int failures = 0;
  const MatterState matter{1e12, 1, 0.3, 2, 5, 0.3, 0.7};
  const auto lte = sampled_fd(view, matter.temperature, 2);
  Rates rates;
  // Number-current kinematics must not change merely because the radiation
  // amplitude is below a configured numerical floor. For lab-isotropic T,
  // v=.6,W=1.25 gives J=1.75E and Gamma=5/7, not the material W.
  nuX_Utils::tensor::generic<CCTK_REAL, 4, 1> v, flux, u, h, current;
  for (int a = 0; a < 4; ++a) v(a) = flux(a) = u(a) = h(a) = 0;
  v(1) = 0.6;
  u(0) = 1.25;
  u(1) = 0.75;
  h(0) = -0.9375e-40;
  h(1) = -1.5625e-40;
  CCTK_REAL gamma = 0;
  failures +=
      !nuX_M1::compute_source_Gamma(1.25, v, 1.75e-40, 1e-40, flux, &gamma);
  failures += relative_source(gamma, 5.0 / 7) > 2e-12;
  nuX_M1::assemble_fnu(u, 1.75e-40, h, &current, 0);
  failures += relative_source(current(0), gamma) > 2e-12;
  failures += !nuX_M1::compute_source_Gamma(1.25, v, 0, 0, flux, &gamma) ||
              gamma != 1.25;
  auto full = lte;
  full.eta[0] = 51;
  // The occupied value rounds to one, but the physical hole must survive.
  failures += collision_detail::occupation(full, 0, 1) != 1;
  failures += relative_source(collision_detail::vacancy(full, 0, 1),
                                std::exp(-50)) > 2e-12;
  collision_detail::Integral nu, anti;
  collision_detail::pair_event(1, std::exp(-2), 1,
      fermi(-79), fermi(81), fermi(79), fermi(-81), 1, 1, nu, anti);
  failures += !(anti.emission_number > 0);
  failures += relative_source(anti.emission_number, anti.loss_number) > 2e-12;
  // Continuous equilibrium initialization and the callback inverse must use
  // one phase-space normalization, not the coarse table-grid moment sum.
  for (const double eta : {2.0, 40.0, 80.0}) {
    auto state = matter;
    state.mu_nu = eta * state.temperature;
    RadiationState exact, recovered;
    failures += !equilibrium_radiation(state, exact);
    M1Quantities moments{};
    for (int s = 0; s < 3; ++s) {
      const double flavors = s == 2 ? 4 : 1;
      moments.n[s] = exact.number[s] * nuX_Utils::rate_units::per_cm3_to_per_nm3 / flavors;
      moments.J[s] = exact.energy[s] * nuX_Utils::rate_units::per_cm3_to_per_nm3 / flavors;
    }
    moments.n[3] = moments.n[2];
    moments.J[3] = moments.J[2];
    failures += !radiation_from_moments(moments, state, recovered);
    for (int s = 0; s < 3; ++s) {
      failures += relative_source(recovered.temperature[s], state.temperature) > 2e-7;
      failures += std::abs(recovered.eta[s] - exact.eta[s]) > 2e-7 * (1 + eta);
    }
  }
  // Pair/Brem must annihilate any pair of FD distributions with opposite
  // chemical potentials at matter T, separately for all four heavy members.
  for (int channel = 0; channel < 3; ++channel) {
    view.include_pair = channel == 0;
    view.include_brem = channel == 1;
    view.include_plasmon = channel == 2;
    for (const double eta : {2.0, 40.0, 80.0}) {
      const auto equilibrium = sampled_fd(view, matter.temperature, eta);
      failures += evaluate_collision(view, matter, equilibrium, rates) != Status::success;
      for (int s = 0; s < 3; ++s) {
        failures += relative_source(rates.emission_number[s],
            rates.absorption_number[s] * equilibrium.number[s]) > 2e-12;
        failures += relative_source(rates.emission_energy[s],
            rates.absorption_energy[s] * equilibrium.energy[s]) > 2e-12;
      }
    }
    const auto hot = sampled_fd(view, 2, -0.4, true);
    failures += evaluate_collision(view, matter, hot, rates) != Status::success;
    const double snu = rates.emission_number[0] - rates.absorption_number[0] * hot.number[0];
    const double santi = rates.emission_number[1] - rates.absorption_number[1] * hot.number[1];
    const double scale = rates.emission_number[0] + rates.emission_number[1] +
        rates.absorption_number[0] * hot.number[0] +
        rates.absorption_number[1] * hot.number[1];
    failures += std::abs(snu - santi) > 2e-12 * scale;
    if (channel == 1) {
      failures += relative_source(rates.emission_number[2], 4 * rates.emission_number[0]) > 2e-12;
      failures += relative_source(rates.absorption_number[2], rates.absorption_number[0]) > 2e-12;
    }
    auto vacuum = hot;
    for (int s = 0; s < 3; ++s) vacuum.number[s] = vacuum.energy[s] = 0;
    failures += evaluate_collision(view, matter, vacuum, rates) != Status::success;
    for (int s = 0; s < 3; ++s)
      failures += !(std::isfinite(rates.emission_number[s]) && rates.emission_number[s] > 0 &&
          std::isfinite(rates.emission_energy[s]) && rates.emission_energy[s] > 0);
  }
  view.include_pair = view.include_brem = view.include_plasmon = false;
  view.include_nes = true;
  for (const double eta : {-3.0, 0.0, 2.0}) {
    const auto r = sampled_fd(view, matter.temperature, eta, true);
    failures += evaluate_collision(view, matter, r, rates) != Status::success;
    for (int s = 0; s < 3; ++s) {
      failures += rates.emission_number[s] != 0 || rates.absorption_number[s] != 0;
      failures += relative_source(rates.emission_energy[s],
          rates.absorption_energy[s] * r.energy[s]) > 2e-12;
      failures += !(rates.absorption_energy[s] > 0); // includes diagonal flux damping
    }
  }
  for (const double t : {0.4, 2.0}) {
    const auto r = sampled_fd(view, t, 0);
    failures += evaluate_collision(view, matter, r, rates) != Status::success;
    for (int s = 0; s < 3; ++s) {
      const double energy_source = rates.emission_energy[s] - rates.absorption_energy[s] * r.energy[s];
      failures += t < matter.temperature ? !(energy_source > 0) : !(energy_source < 0);
      failures += rates.emission_number[s] != 0 || rates.absorption_number[s] != 0;
    }
  }
  auto outside = matter;
  outside.mu_e = -1;
  failures += evaluate_collision(view, outside, lte, rates) != Status::out_of_range;
  return failures;
}
} // namespace

void collision_self_test() {
  const amrex::Gpu::LaunchSafeGuard launch_region(true);
  constexpr int ne = 64;
  std::vector<double> energy(ne), weight(ne);
  for (int i = 0; i < ne; ++i) energy[i] = std::pow(10.0, -5 + 8.0 * i / (ne - 1));
  for (int i = 0; i < ne; ++i)
    weight[i] = (energy[i + (i + 1 < ne)] - energy[i - (i > 0)]) / 2;
  amrex::Gpu::DeviceVector<double> energies, weights, temperatures, etas, densities;
  upload(energies, energy);
  upload(weights, weight);
  upload(temperatures, {8, 12});
  upload(etas, {-3, 3});
  upload(densities, {6, 16});
  wli::TableMeta<5> meta{};
  meta.n = {ne, ne, 4, 2, 2};
  std::vector<double> kernels(meta.size());
  // Deliberately unequal synthetic I/II diagonals test shared-event
  // conservation, NOT physical-table diagonal reciprocity/provenance.
  for (std::size_t i = 0; i < kernels.size(); ++i)
    kernels[i] = std::log10((int(i / (ne * ne)) % 4 + 1) * 1e-30);
  wli::ResidentTable<5> pair, nes, brem;
  pair.upload(meta, kernels.data(), kernels.size());
  nes.upload(meta, kernels.data(), kernels.size());
  meta.n = {ne, ne, 1, 2, 2};
  kernels.assign(meta.size(), -3); // S_sigma=1e-3 MeV^-1
  brem.upload(meta, kernels.data(), kernels.size());
  TableView view;
  view.ne = ne;
  view.energies = energies.data();
  view.weights = weights.data();
  view.pair = pair.view();
  view.nes = nes.view();
  view.brem = brem.view();
  view.pair_temp = view.nes_temp = view.brem_temp = {temperatures.data(), 2, 1e8, 1e12};
  view.pair_eta = view.nes_eta = {etas.data(), 2, 1e-3, 1e3};
  view.brem_rho = {densities.data(), 2, 1e6, 1e16};
  view.tail_tolerance = 1e-6;
  amrex::Gpu::DeviceScalar<int> failures(0);
  auto *const out = failures.dataPtr();
  amrex::ParallelFor(1, [=] CCTK_DEVICE(int) { *out = check_collision_models(view); });
  amrex::Gpu::synchronize();
  if (failures.dataValue())
    CCTK_VERROR("WeakLib monopole collision self-test failed (%d checks)", failures.dataValue());
  CCTK_INFO("WeakLib monopole conservation/detailed-balance self-test passed");
}
} // namespace nuX_WeakLib
