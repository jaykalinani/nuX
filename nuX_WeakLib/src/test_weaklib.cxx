#include <array>
#include <cmath>
#include <limits>
#include <vector>

#include <AMReX_GpuContainers.H>
#include <AMReX_GpuControl.H>
#include <AMReX_GpuMemory.H>
#include <AMReX_GpuLaunch.H>

#include "cctk.h"
#include "nuX_weaklib.hxx"
#include "nuX_weaklib_io.hxx"

namespace nuX_WeakLib {
namespace {
void upload(amrex::Gpu::DeviceVector<double> &dest, const std::vector<double> &src) {
  dest.resize(src.size());
  amrex::Gpu::htod_memcpy(dest.data(), src.data(), src.size() * sizeof(double));
}
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool near(double got, double want) {
  return std::isfinite(got) && std::abs(got - want) <= 1e-11 * std::abs(want) + 1e-30;
}
} // namespace

// Opt-in tests run on the SAME backend as CarpetX. They are source-only until
// a future authorized build/run; no production table is needed by the math.
void self_test() {
  host_io::io_self_test(); // host files/collectives only, never inside a device lambda
  energy_node_self_test();
  const amrex::Gpu::LaunchSafeGuard launch_region(true);
  constexpr int ne = 256;
  const std::array<int, 4> emab_shape{ne, 2, 2, 2};
  const std::array<int, 5> iso_shape{ne, 2, 2, 2, 2};
  std::vector<double> energy(ne), node_coordinates(ne), weight(ne);
  for (int i = 0; i < ne; ++i) {
    energy[i] = std::pow(10.0, -4 + 7.0 * i / (ne - 1));
    node_coordinates[i] = double(i);
  }
  for (int i = 0; i < ne; ++i)
    weight[i] = (energy[i + (i + 1 < ne)] - energy[i - (i > 0)]) / 2;
  std::array<amrex::Gpu::DeviceVector<double>, 6> axes;
  upload(axes[0], energy);
  upload(axes[1], node_coordinates);
  upload(axes[2], weight);
  upload(axes[3], {6, 14});
  upload(axes[4], {8, 12});
  upload(axes[5], {0.01, 0.6});
  TableView v;
  v.energies = axes[0].data();
  v.energy_node_coordinates = axes[1].data();
  v.weights = axes[2].data();
  v.ne = ne;
  v.tail_tolerance = 1e-6;
  v.include_emab = v.include_iso = true;
  v.emab_thermo = v.iso_thermo = {
      {axes[3].data(), 2, 1e6, 1e14}, {axes[4].data(), 2, 1e8, 1e12},
      {axes[5].data(), 2, 0.01, 0.6}};
  wli::ResidentTable<4> emab[2];
  wli::ResidentTable<5> iso[2];
  for (int s = 0; s < 2; ++s) {
    wli::TableMeta<4> a{};
    for (int d = 0; d < 4; ++d) a.n[d] = emab_shape[d];
    std::vector<double> absorption(std::size_t(a.size()));
    for (std::size_t i = 0; i < absorption.size(); ++i)
      absorption[i] = std::log10((s + 1) * 0.01 * energy[i % ne] * energy[i % ne]);
    emab[s].upload(a, absorption.data(), absorption.size());
    v.emab[s] = emab[s].view();
    wli::TableMeta<5> b{};
    for (int d = 0; d < 5; ++d) b.n[d] = iso_shape[d];
    std::vector<double> scatter(std::size_t(b.size()));
    const double phi0 = s == 0 ? 2 : 5, phi1 = s == 0 ? -0.5 : 1.5;
    for (int m = 0; m < 2; ++m) v.iso_offsets[s][m] = 1 + s + 2 * m;
    for (std::size_t i = 0; i < scatter.size(); ++i) {
      const int m = int(i / ne) % 2;
      scatter[i] = std::log10((m == 0 ? phi0 : phi1) + v.iso_offsets[s][m]);
    }
    iso[s].upload(b, scatter.data(), scatter.size());
    v.iso[s] = iso[s].view();
  }
  // A large finite opacity can leave representable emission even when its
  // equilibrium population is below double precision's absolute range.
  wli::ResidentTable<4> strong_emab;
  wli::TableMeta<4> strong_meta{};
  for (int d = 0; d < 4; ++d) strong_meta.n[d] = emab_shape[d];
  const std::vector<double> strong_values(std::size_t(strong_meta.size()), 130);
  strong_emab.upload(strong_meta, strong_values.data(), strong_values.size());
  const auto strong_view = strong_emab.view();

  // Independent quadrature reference. No production helper computes these
  // moment powers or their number-versus-energy opacity averages for the test.
  double m2 = 0, m3 = 0, m4 = 0, m5 = 0;
  for (int i = 0; i < ne; ++i) {
    const double e = energy[i];
    const double f = 1 / (1 + std::exp(e / 5));
    const double x = weight[i] * e * e * f;
    m2 += x; m3 += x * e; m4 += x * e * e; m5 += x * e * e * e;
  }
  amrex::Gpu::DeviceScalar<int> failures(0);
  int *const failure = failures.dataPtr();
  amrex::ParallelFor(1, [=] AMREX_GPU_DEVICE(int) noexcept {
    int bad = 0;
    Rates r;
    bad += evaluate(v, 1e10, 5, 0.2, 0, r) != Status::success;
    for (int s = 0; s < 3; ++s) {
      const double multiplicity = s == 2 ? 4 : 1;
      bad += !near(r.number[s], multiplicity * phase_space * m2);
      bad += !near(r.energy[s], multiplicity * phase_space * m3);
      const double coefficient = s < 2 ? 0.01 * (s + 1) : 0;
      bad += !near(r.absorption_number[s], coefficient * m4 / m2);
      bad += !near(r.absorption_energy[s], coefficient * m5 / m3);
      // Independent angular integral of (1-mu)*(phi0+phi1*mu):
      // (1/2) integral[-1,1] = phi0-phi1/3. These are expansion
      // coefficients, not the already-normalized angular moments.
      const double transport = s == 0 ? 13.0 / 6 : s == 1 ? 9.0 / 2 : 10.0 / 3;
      bad += !near(r.scattering[s], 4 * pi * transport * m5 / m3);
      // Both collision components vanish at equilibrium. Emission remains
      // finite and positive at zero evolved population (which is not an input).
      bad += !near(r.emission_number[s], r.absorption_number[s] * r.number[s]);
      bad += !near(r.emission_energy[s], r.absorption_energy[s] * r.energy[s]);
      if (s < 2) bad += !(r.emission_number[s] > 0 && r.emission_energy[s] > 0);
    }
    Rates scratch;
    bad += evaluate(v, 1e6, 5, 0.01, 0, scratch) != Status::success;
    bad += evaluate(v, 1e14, 5, 0.6, 0, scratch) != Status::success;
    bad += evaluate(v, 1e10, 5, 0.2, 2, scratch) != Status::success;
    bad += !(scratch.number[0] > scratch.number[1]);
    bad += !near(scratch.number[2], r.number[2]);
    bad += !near(neutrino_chemical_potential(938.27208816 + 2, 939.56542052 + 6, 12, false),
                 neutrino_chemical_potential(2, 6, 12, true));
    bad += evaluate(v, 1e5, 5, 0.2, 0, scratch) != Status::out_of_range;
    bad += evaluate(v, 1e10, 5, 0.7, 0, scratch) != Status::out_of_range;
    bad += evaluate(v, 1e10, 5, 0.2, 1e6, scratch) != Status::insufficient_energy_coverage;
    bad += evaluate(v, 1e10, 0, 0.2, 0, scratch) != Status::invalid_state;
    bad += evaluate(v, 1e10, 5, 0.2, std::numeric_limits<double>::quiet_NaN(), scratch) != Status::invalid_state;
    auto elastic_only = v;
    elastic_only.include_emab = false;
    bad += evaluate(elastic_only, 1e10, 5, 0.2, 0, scratch) != Status::success;
    for (int s = 0; s < 3; ++s)
      bad += scratch.absorption_number[s] != 0 || scratch.emission_energy[s] != 0;
    auto absorption_only = v;
    absorption_only.include_iso = false;
    bad += evaluate(absorption_only, 1e10, 5, 0.2, 0, scratch) != Status::success;
    for (int s = 0; s < 3; ++s) bad += scratch.scattering[s] != 0;
    auto negative_kernel = v;
    negative_kernel.emab_offsets[0] = 1e10;
    bad += evaluate(negative_kernel, 1e10, 5, 0.2, 0, scratch) != Status::invalid_kernel;
    const double uniform[] = {-4, -2, 0, 2};
    const double nonuniform[] = {0, std::log10(2.0), std::log10(3.0), 2};
    const double narrow_nonuniform[] = {10, 10 + 1e-14, 10 + 2e-14, 10 + 1e-13};
    bad += !coordinate_axis_is_uniform(uniform, 4);
    bad += coordinate_axis_is_uniform(nonuniform, 4);
    bad += coordinate_axis_is_uniform(narrow_nonuniform, 4);
    // Large degeneracy suppresses one species below absolute representability,
    // but its normalized absorption and mean energy still have finite limits.
    double dilute2 = 0, dilute3 = 0, dilute4 = 0, dilute5 = 0;
    for (int ie = 0; ie < ne; ++ie) {
      const double e = v.energies[ie];
      const double w = v.weights[ie] * e * e * std::exp(-e / 0.01);
      dilute2 += w;
      dilute3 += w * e;
      dilute4 += w * e * e;
      dilute5 += w * e * e * e;
    }
    for (const double mu : {-10.0, 10.0}) {
      bad += evaluate(v, 1e10, 0.01, 0.2, mu, scratch) != Status::success;
      const int suppressed = mu > 0 ? 1 : 0;
      bad += scratch.number[suppressed] != 0 || scratch.energy[suppressed] != 0;
      bad += !near(scratch.mean_energy[suppressed], dilute3 / dilute2);
      const double coefficient = 0.01 * (suppressed + 1);
      bad += !near(scratch.absorption_number[suppressed], coefficient * dilute4 / dilute2);
      bad += !near(scratch.absorption_energy[suppressed], coefficient * dilute5 / dilute3);
    }
    auto strong_absorption = v;
    strong_absorption.emab[1] = strong_view;
    strong_absorption.emab_offsets[1] = 0;
    bad += evaluate(strong_absorption, 1e10, 0.01, 0.2, 10, scratch) != Status::success;
    bad += scratch.number[1] != 0 || scratch.energy[1] != 0;
    const double log_eta_number = std::log(phase_space) + std::log(dilute2) +
                                 130 * std::log(10.0) - 1000;
    const double log_eta_energy = std::log(phase_space) + std::log(dilute3) +
                                 130 * std::log(10.0) - 1000;
    // An absolute near() floor would miss loss of these tiny positive rates.
    bad += !(scratch.emission_number[1] > 0 &&
             near(std::log(scratch.emission_number[1]), log_eta_number));
    bad += !(scratch.emission_energy[1] > 0 &&
             near(std::log(scratch.emission_energy[1]), log_eta_energy));
    // exp(-750) alone is zero in double, but its product with a finite
    // phase-space prefactor need not underflow. Reconstruct in logarithms.
    const double restored = restore_fd_scale(1, 1e30, -750);
    bad += !(restored > 0 && near(std::log(restored), std::log(1e30) - 750));
    const double large_integral = restore_fd_scale(1e300, 1e30, -1000);
    bad += !(large_integral > 0 &&
             near(std::log(large_integral), std::log(1e300) + std::log(1e30) - 1000));
    *failure = bad;
  });
  amrex::Gpu::synchronize();
  if (failures.dataValue()) CCTK_VERROR("WeakLib synthetic self-tests failed: %d checks", failures.dataValue());
  CCTK_INFO("WeakLib synthetic grey integrals, multiplicities, balance, scattering normalization and range tests passed");
}
} // namespace nuX_WeakLib
