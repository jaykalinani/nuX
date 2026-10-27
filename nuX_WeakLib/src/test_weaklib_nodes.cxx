// Opt-in source tests: no compiled/backend execution is implied by this file.
#include "nuX_weaklib.hxx"

#include <array>
#include <cmath>
#include <vector>

#include <AMReX_GpuContainers.H>
#include <AMReX_GpuControl.H>
#include <AMReX_GpuLaunch.H>
#include <AMReX_GpuMemory.H>
#include "cctk.h"

namespace nuX_WeakLib {
namespace {
void upload(amrex::Gpu::DeviceVector<double> &dest, const std::vector<double> &src) {
  dest.resize(src.size());
  amrex::Gpu::htod_memcpy(dest.data(), src.data(), src.size() * sizeof(double));
}

// Analytic LOG-STORED manufactured data: nonlinear/discrete energy dependence,
// exactly trilinear in the GENUINE matter coordinates. No WLI call defines the
// oracle. Distinct species/moment offsets keep signed phi1 recoverable.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
double analytic_log_value(bool iso, int s, int m, int ie, double energy,
                          double ld, double lt, double ye) noexcept {
  const double x = (ld - 6) / 8, y = (lt - 8) / 4, z = (ye - 0.01) / 0.59;
  const double matter = x + 0.5 * y - 0.25 * z + 0.2 * x * y +
                        0.15 * x * z - 0.1 * y * z + 0.05 * x * y * z;
  const double le = std::log10(energy);
  if (!iso) return -2 + 0.2 * s + 0.12 * le * le + 0.01 * (ie % 5) + 0.05 * matter;
  return m == 0 ? 0.3 + 0.08 * s + 0.03 * le * le + 0.02 * (ie % 7) + 0.05 * matter
                : -0.4 + 0.06 * s + 0.015 * le * le + 0.01 * (ie % 3) + 0.03 * matter;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
double analytic_value(bool iso, int s, int m, int ie, double energy,
                      double ld, double lt, double ye) noexcept {
  const double offset = !iso ? 0.001 * (s + 1) : m == 0 ? 0.1 * (s + 1) : 2 + s;
  return std::pow(10.0, analytic_log_value(iso, s, m, ie, energy, ld, lt, ye)) - offset;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool near(double value, double expected) noexcept {
  return std::isfinite(value) && std::abs(value - expected) <= 1e-11 * std::abs(expected) + 1e-30;
}

// Independent PHYSICAL-E trapezoidal sums, direct FD (moderate degeneracy),
// analytic node values. No provider weights, wrapper or FD helper is the oracle.
Rates reference_rates(const std::vector<double> &energy, double rho, double temp,
                       double ye, double mu_nu) {
  Rates out;
  const double ld = std::log10(rho), lt = std::log10(temp / kboltz_mev_kelvin);
  for (int s = 0; s < 3; ++s) {
    const double mu = s == 0 ? mu_nu : s == 1 ? -mu_nu : 0;
    double n = 0, j = 0, a0 = 0, a1 = 0, scattering = 0;
    for (std::size_t i = 0; i < energy.size(); ++i) {
      const double e = energy[i];
      const double dE = (i ? (e - energy[i - 1]) / 2 : 0) +
                        (i + 1 < energy.size() ? (energy[i + 1] - e) / 2 : 0);
      const double f = 1 / (1 + std::exp((e - mu) / temp));
      const double wn = dE * e * e * f, wj = wn * e;
      n += wn;
      j += wj;
      if (s < 2) {
        const double chi = analytic_value(false, s, 0, int(i), e, ld, lt, ye);
        a0 += wn * chi;
        a1 += wj * chi;
      }
      double transport = 0;
      for (int c = 0; c < 2; ++c) {
        if (s < 2 && c != s) continue;
        const double p0 = analytic_value(true, c, 0, int(i), e, ld, lt, ye);
        const double p1 = analytic_value(true, c, 1, int(i), e, ld, lt, ye);
        transport += (s == 2 ? 0.5 : 1) * 4 * pi * e * e * (p0 - p1 / 3);
      }
      scattering += wj * transport;
    }
    const double prefactor = (s == 2 ? 4 : 1) * phase_space;
    out.number[s] = prefactor * n;
    out.energy[s] = prefactor * j;
    out.mean_energy[s] = j / n;
    out.absorption_number[s] = a0 / n;
    out.absorption_energy[s] = a1 / j;
    out.emission_number[s] = prefactor * a0;
    out.emission_energy[s] = prefactor * a1;
    out.scattering[s] = scattering / j;
  }
  return out;
}

void test_grid(const std::vector<double> &energy, bool full_integral) {
  const int ne = int(energy.size()); // test fixtures are bounded by 512 nodes
  std::vector<double> coordinates(ne), weights(ne);
  for (int i = 0; i < ne; ++i) {
    coordinates[i] = double(i);
    weights[i] = (energy[i + (i + 1 < ne)] - energy[i - (i > 0)]) / 2;
  }
  std::array<amrex::Gpu::DeviceVector<double>, 6> axes;
  upload(axes[0], energy);
  upload(axes[1], coordinates); // one node-coordinate buffer, not an added log-E buffer
  upload(axes[2], weights);
  upload(axes[3], {6, 14});
  upload(axes[4], {8, 12});
  upload(axes[5], {0.01, 0.6});
  TableView v;
  v.ne = ne;
  v.energies = axes[0].data();
  v.energy_node_coordinates = axes[1].data();
  v.weights = axes[2].data();
  v.emab_thermo = v.iso_thermo = {
      {axes[3].data(), 2, 1e6, 1e14}, {axes[4].data(), 2, 1e8, 1e12},
      {axes[5].data(), 2, 0.01, 0.6}};
  v.include_emab = v.include_iso = true;
  v.tail_tolerance = 1e-6;
  wli::ResidentTable<4> emab[2];
  wli::ResidentTable<5> iso[2];
  std::array<std::vector<double>, 2> absorption;
  for (int s = 0; s < 2; ++s) {
    v.emab_offsets[s] = 0.001 * (s + 1);
    v.iso_offsets[s][0] = 0.1 * (s + 1);
    v.iso_offsets[s][1] = 2 + s;
    wli::TableMeta<4> a{};
    a.n = {ne, 2, 2, 2};
    absorption[s].reserve(std::size_t(a.size()));
    for (int iy = 0; iy < 2; ++iy)
      for (int it = 0; it < 2; ++it)
        for (int id = 0; id < 2; ++id)
          for (int ie = 0; ie < ne; ++ie)
            absorption[s].push_back(analytic_log_value(false, s, 0, ie, energy[ie],
                                                       6 + 8 * id, 8 + 4 * it, 0.01 + 0.59 * iy));
    emab[s].upload(a, absorption[s].data(), absorption[s].size());
    v.emab[s] = emab[s].view();
    wli::TableMeta<5> b{};
    b.n = {ne, 2, 2, 2, 2};
    std::vector<double> scattering;
    scattering.reserve(std::size_t(b.size()));
    for (int iy = 0; iy < 2; ++iy)
      for (int it = 0; it < 2; ++it)
        for (int id = 0; id < 2; ++id)
          for (int m = 0; m < 2; ++m)
            for (int ie = 0; ie < ne; ++ie)
              scattering.push_back(analytic_log_value(true, s, m, ie, energy[ie],
                                                       6 + 8 * id, 8 + 4 * it, 0.01 + 0.59 * iy));
    iso[s].upload(b, scattering.data(), scattering.size());
    v.iso[s] = iso[s].view();
  }
  int host_bad = !physical_energy_nodes_are_valid(energy.data(), energy.size());
  // Sensitivity control for the old WRONG nonuniform physical-log-E locator.
  // This deliberately bypasses the adapter on HOST buffers only; it is not a
  // supported physical-energy API and must disagree with the nonlinear oracle.
  if (ne == 5 && energy[1] == 0.11) {
    std::vector<double> loge(ne);
    for (int i = 0; i < ne; ++i) loge[i] = std::log10(energy[i]);
    const double ld[] = {6, 14}, lt[] = {8, 12}, ye[] = {0.01, 0.6};
    const double wrong = wli::EmAbInterpolateSingleVariable4DPoint(
        loge[2], 9.2, 10.8, 0.23, loge.data(), ne, ld, 2, lt, 2, ye, 2,
        v.emab_offsets[0], absorption[0].data());
    const double expected = analytic_value(false, 0, 0, 2, energy[2], 9.2, 10.8, 0.23);
    host_bad += !std::isfinite(wrong) || std::abs(wrong - expected) <= 1e-6 * std::abs(expected);
  }
  amrex::GpuArray<Rates, 2> reference{};
  if (full_integral) {
    reference[0] = reference_rates(energy, 1e10, 5, 0.2, 2);
    reference[1] = reference_rates(energy, 3e12, 10, 0.45, -3);
  }
  amrex::Gpu::DeviceScalar<int> failures(0);
  int *const failure = failures.dataPtr();
  amrex::ParallelFor(1, [=] AMREX_GPU_DEVICE(int) noexcept {
    int bad = host_bad;
    const double points[4][3] = {{6, 8, 0.01}, {14, 12, 0.6},
                                 {9.2, 10.8, 0.23}, {12.5, 9.5, 0.48}};
    for (int p = 0; p < 4; ++p)
      for (int ie = 0; ie < ne; ++ie)
        for (int s = 0; s < 2; ++s) {
          const double ld = points[p][0], lt = points[p][1], ye = points[p][2];
          bad += !near(emab_at_energy_node(v, s, ie, ld, lt, ye),
                       analytic_value(false, s, 0, ie, v.energies[ie], ld, lt, ye));
          for (int m = 0; m < 2; ++m)
            bad += !near(iso_at_energy_node(v, s, m, ie, ld, lt, ye),
                         analytic_value(true, s, m, ie, v.energies[ie], ld, lt, ye));
        }
    bad += !std::isnan(emab_at_energy_node(v, 0, -1, 9, 10, 0.2));
    bad += !std::isnan(emab_at_energy_node(v, 0, ne, 9, 10, 0.2));
    bad += !std::isnan(iso_at_energy_node(v, 0, 0, -1, 9, 10, 0.2));
    bad += !std::isnan(iso_at_energy_node(v, 0, 0, ne, 9, 10, 0.2));
    bad += !std::isnan(emab_at_energy_node(v, -1, 0, 9, 10, 0.2));
    bad += !std::isnan(emab_at_energy_node(v, 2, 0, 9, 10, 0.2));
    bad += !std::isnan(iso_at_energy_node(v, -1, 0, 0, 9, 10, 0.2));
    bad += !std::isnan(iso_at_energy_node(v, 2, 0, 0, 9, 10, 0.2));
    bad += !std::isnan(iso_at_energy_node(v, 0, -1, 0, 9, 10, 0.2));
    bad += !std::isnan(iso_at_energy_node(v, 0, 2, 0, 9, 10, 0.2));
    auto malformed = v;
    malformed.iso[0].n[1] = 3;
    // The dynamic table extent cannot authorize indexing beyond the fixed
    // two-moment offset storage, even for a malformed manufactured view.
    bad += !std::isnan(iso_at_energy_node(malformed, 0, 2, 0, 9, 10, 0.2));
    const double duplicate[] = {1, 1}, zero[] = {0, 1}, descending[] = {2, 1};
    const double single[] = {1};
    const double negative[] = {-1, 1}, nan[] = {1, std::numeric_limits<double>::quiet_NaN()};
    const double nonfinite[] = {1, std::numeric_limits<double>::infinity()};
    bad += physical_energy_nodes_are_valid(duplicate, 2);
    bad += physical_energy_nodes_are_valid(zero, 2);
    bad += physical_energy_nodes_are_valid(descending, 2);
    bad += physical_energy_nodes_are_valid(negative, 2);
    bad += physical_energy_nodes_are_valid(nan, 2);
    bad += physical_energy_nodes_are_valid(nonfinite, 2);
    bad += physical_energy_nodes_are_valid(single, 1);
    bad += physical_energy_nodes_are_valid(nullptr, 2);
    bad += physical_energy_nodes_are_valid(nullptr, std::size_t(std::numeric_limits<int>::max()) + 1);
    if (full_integral)
      for (int p = 0; p < 2; ++p) {
        Rates actual;
        bad += evaluate(v, p == 0 ? 1e10 : 3e12, p == 0 ? 5 : 10,
                        p == 0 ? 0.2 : 0.45, p == 0 ? 2 : -3, actual) != Status::success;
        const auto expected = reference[p];
        for (int s = 0; s < 3; ++s) {
          bad += !near(actual.number[s], expected.number[s]);
          bad += !near(actual.energy[s], expected.energy[s]);
          bad += !near(actual.mean_energy[s], expected.mean_energy[s]);
          bad += !near(actual.absorption_number[s], expected.absorption_number[s]);
          bad += !near(actual.absorption_energy[s], expected.absorption_energy[s]);
          bad += !near(actual.emission_number[s], expected.emission_number[s]);
          bad += !near(actual.emission_energy[s], expected.emission_energy[s]);
          bad += !near(actual.scattering[s], expected.scattering[s]);
          bad += !near(actual.emission_number[s], actual.absorption_number[s] * actual.number[s]);
          bad += !near(actual.emission_energy[s], actual.absorption_energy[s] * actual.energy[s]);
        }
      }
    *failure = bad;
  });
  amrex::Gpu::synchronize();
  if (failures.dataValue())
    CCTK_VERROR("WeakLib exact energy-node tests failed for ne=%d: %d checks", ne, failures.dataValue());
}
} // namespace

void energy_node_self_test() {
  const amrex::Gpu::LaunchSafeGuard launch_region(true);
  for (const auto &energy : std::vector<std::vector<double>>{
           {0.1, 100}, {0.1, 0.11, 100}, {0.1, 0.11, 0.12, 0.13, 100},
           {1, 2, 4, 8, 16, 128}, {0.1, 1, 10, 100}})
    test_grid(energy, false);
  std::vector<double> energy(512);
  for (std::size_t i = 0; i < energy.size(); ++i) {
    const double x = double(i) / (energy.size() - 1);
    energy[i] = std::pow(10.0, -4 + 7 * (0.75 * x + 0.25 * x * x));
  }
  test_grid(energy, true);
  CCTK_INFO("WeakLib exact energy-node analytic and physical-E quadrature source tests passed");
}
} // namespace nuX_WeakLib
