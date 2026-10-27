#ifndef NUX_WEAKLIB_HXX
#define NUX_WEAKLIB_HXX

#include <cmath>
#include <cstddef>
#include <limits>
#include <type_traits>

#include <AMReX_Array.H>
#include "wli_opacity_emab_iso.H"

namespace nuX_WeakLib {

// WeakLib tables use Kelvin, not EOSX's MeV temperature. h*c is in MeV cm.
inline constexpr double kboltz_mev_kelvin = 8.617333262145e-11;
inline constexpr double hc_mev_cm = 1.2398419843320026e-10;
inline constexpr double pi = 3.141592653589793238462643383279502884;
inline constexpr double neutron_proton_mass_difference = 1.29333236;
inline constexpr double phase_space = 4 * pi / (hc_mev_cm * hc_mev_cm * hc_mev_cm);

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
double neutrino_chemical_potential(double mup, double mun, double mue,
                                  bool nucleon_rest_mass_subtracted) noexcept {
  return mup + mue - mun - (nucleon_rest_mass_subtracted ? neutron_proton_mass_difference : 0);
}

struct AxisView {
  const double *log_nodes = nullptr;
  int n = 0;
  double min = 0, max = 0;
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  bool contains(double x) const noexcept {
    return std::isfinite(x) && x >= min && x <= max;
  }
};

struct ThermoView {
  AxisView rho, temp, ye; // ye.log_nodes actually contains RAW linear Ye
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  bool contains(double d, double t, double y) const noexcept {
    return rho.contains(d) && temp.contains(t) && ye.contains(y);
  }
};

// Only POD views enter grid kernels. The owners live on the host until shutdown.
struct TableView {
  ThermoView emab_thermo, iso_thermo;
  const double *energies = nullptr, *energy_node_coordinates = nullptr, *weights = nullptr;
  int ne = 0;
  wli::TableView<4> emab[2]{};
  wli::TableView<5> iso[2]{};
  double emab_offsets[2]{}, iso_offsets[2][2]{};
  bool include_emab = false, include_iso = false;
  double tail_tolerance = 0;
};
static_assert(std::is_trivially_copyable_v<TableView>);

struct Rates {
  double number[3]{}, energy[3]{}; // cm^-3, MeV cm^-3, nux is ALL FOUR species
  double mean_energy[3]{}; // MeV; remains defined when absolute LTE density underflows
  double absorption_number[3]{}, absorption_energy[3]{}, scattering[3]{}; // cm^-1
  double emission_number[3]{}, emission_energy[3]{}; // cm^-4, MeV cm^-4; multiply by c
};

enum class Status { success, invalid_state, out_of_range, insufficient_energy_coverage,
                    invalid_kernel, invalid_moment };

// VALUE ONLY at integer energy node ie. The resident coordinate buffer is
// {0,...,ne-1}, NOT log(E). WLI's existing value leaf then selects the exact
// stored energy plane while retaining genuine rho/T/Ye interpolation. Never
// use this coordinate for off-node physical energy queries or d/dE routines.
// Invalid public node/species indices return NaN; evaluate rejects invalid_kernel.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
double emab_at_energy_node(const TableView &v, int species, int ie,
                          double log_rho, double log_kelvin, double ye) noexcept {
  if (v.ne < 2 || ie < 0 || ie >= v.ne || species < 0 || species >= 2 ||
      !v.energy_node_coordinates)
    return std::numeric_limits<double>::quiet_NaN();
  const auto t = v.emab[species];
  const auto a = v.emab_thermo;
  return wli::EmAbInterpolateSingleVariable4DPoint(
      double(ie), log_rho, log_kelvin, ye, v.energy_node_coordinates, v.ne,
      a.rho.log_nodes, a.rho.n, a.temp.log_nodes, a.temp.n,
      a.ye.log_nodes, a.ye.n, v.emab_offsets[species], t.data);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
double iso_at_energy_node(const TableView &v, int species, int moment, int ie,
                          double log_rho, double log_kelvin, double ye) noexcept {
  if (v.ne < 2 || ie < 0 || ie >= v.ne || species < 0 || species >= 2 ||
      !v.energy_node_coordinates || moment < 0 || moment >= 2 ||
      moment >= v.iso[species].n[1])
    return std::numeric_limits<double>::quiet_NaN();
  const auto t = v.iso[species];
  const auto a = v.iso_thermo;
  return wli::IsoInterpolateSingleVariable5DPoint(
      double(ie), log_rho, log_kelvin, ye, v.energy_node_coordinates, v.ne,
      a.rho.log_nodes, a.rho.n, a.temp.log_nodes, a.temp.n,
      a.ye.log_nodes, a.ye.n, moment, t.n[1], v.iso_offsets[species][moment], t.data);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool physical_energy_nodes_are_valid(const double *nodes, std::size_t n) noexcept {
  if (!nodes || n < 2 || n > std::size_t(std::numeric_limits<int>::max())) return false;
  for (std::size_t i = 0; i < n; ++i)
    if (!std::isfinite(nodes[i]) || nodes[i] <= 0 || (i && nodes[i] <= nodes[i - 1]))
      return false;
  return true;
}

// Stable FD evaluation; no arbitrary exponent/chemical-potential clamp.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
double fermi(double x) noexcept {
  if (x >= 0) {
    const double e = std::exp(-x);
    return e / (1 + e);
  }
  return 1 / (1 + std::exp(x));
}

// WLI's coordinate-to-index formula requires uniform interpolation
// coordinates, not merely a monotone grid. Validate before uploading axes.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool coordinate_axis_is_uniform(const double *nodes, int n) noexcept {
  if (n < 2 || !std::isfinite(nodes[0]) || !std::isfinite(nodes[n - 1]))
    return false;
  const double step = (nodes[n - 1] - nodes[0]) / (n - 1);
  if (!std::isfinite(step) || step <= 0)
    return false;
  // Translation of the coordinate origin must not relax the spacing check.
  const double tolerance = 64 * std::numeric_limits<double>::epsilon() *
                           (nodes[n - 1] - nodes[0]);
  for (int i = 1; i < n; ++i)
    if (!std::isfinite(nodes[i]) || nodes[i] <= nodes[i - 1] ||
        std::abs((nodes[i] - nodes[0]) - i * step) > tolerance)
      return false;
  return true;
}

// Factor exp(mu/T) only for negative mu. Ratios of spectral moments then
// stay finite even when the absolute antineutrino population underflows.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
double scaled_fermi(double energy, double temp, double mu) noexcept {
  return mu < 0 ? std::exp(-energy / temp) /
                      (1 + std::exp((mu - energy) / temp))
                : fermi((energy - mu) / temp);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
double restore_fd_scale(double integral, double prefactor, double log_scale) noexcept {
  // Neither exp(log_scale) nor prefactor*integral need be representable on
  // its own when the final physical density or emissivity is finite.
  if (integral == 0)
    return 0;
  return log_scale == 0 ? prefactor * integral
                       : std::exp(std::log(integral) + std::log(prefactor) + log_scale);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool nonnegative_finite(double x) noexcept {
  return std::isfinite(x) && x >= 0;
}

// Upper bound on the omitted high-energy FD moment. f <= exp((mu-E)/T).
// The polynomial is the exact integer incomplete-gamma tail for E^q, q=2,3.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
double high_tail(int q, double e, double temp, double mu) noexcept {
  const double t2 = temp * temp;
  const double p = q == 2 ? e * e + 2 * e * temp + 2 * t2
                          : e * e * e + 3 * e * e * temp + 6 * e * t2 + 6 * t2 * temp;
  return std::exp((mu - e) / temp) * temp * p;
}

// Grey equilibrium-spectrum averages, not a grey energy-squared correction.
// Integrate over the SAME in-range energy nodes for N, J, emissivities and
// absorption averages. This makes Kirchhoff balance exact for the discrete
// spectral integral. No division by evolved N/J: vacuum emission stays finite.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Status evaluate(const TableView &v, double rho, double temp, double ye,
                double mu_nu, Rates &out) noexcept {
  out = {};
  if (!std::isfinite(rho) || rho <= 0 || !std::isfinite(temp) || temp <= 0 ||
      !std::isfinite(ye) || !std::isfinite(mu_nu) || v.ne < 2)
    return Status::invalid_state;
  const double kelvin = temp / kboltz_mev_kelvin;
  if ((v.include_emab && !v.emab_thermo.contains(rho, kelvin, ye)) ||
      (v.include_iso && !v.iso_thermo.contains(rho, kelvin, ye)))
    return Status::out_of_range;
  const double ld = std::log10(rho), lt = std::log10(kelvin);
  for (int s = 0; s < 3; ++s) {
    const double mu = s == 0 ? mu_nu : s == 1 ? -mu_nu : 0;
    const double log_scale = mu < 0 ? mu / temp : 0;
    double n = 0, j = 0, a0 = 0, a1 = 0, scat = 0;
    for (int ie = 0; ie < v.ne; ++ie) {
      const double e = v.energies[ie];
      const double w = v.weights[ie] * e * e * scaled_fermi(e, temp, mu);
      n += w;
      j += w * e;
      if (v.include_emab && s < 2) {
        const double chi = emab_at_energy_node(v, s, ie, ld, lt, ye);
        if (!nonnegative_finite(chi)) return Status::invalid_kernel;
        a0 += w * chi;
        a1 += w * e * chi;
      }
      if (v.include_iso) {
        // Neutral-current scattering is flavor independent, but neutrinos and
        // antineutrinos have distinct weak-magnetism tables. nux = 2 nu + 2 anti.
        const int first = s < 2 ? s : 0, last = s < 2 ? s : 1;
        double transport = 0;
        for (int c = first; c <= last; ++c) {
          double phi[2];
          for (int m = 0; m < 2; ++m)
            phi[m] = iso_at_energy_node(v, c, m, ie, ld, lt, ye);
          // WeakLib stores Legendre EXPANSION COEFFICIENTS: the creator
          // multiplies the l-th angular moment by (2*l+1)/2. Thus the first
          // coefficient contributes phi_1/3 to the transport moment, not
          // phi_1. It can be signed; do not clamp it.
          const double kappa = 4 * pi * e * e * (phi[0] - phi[1] / 3);
          if (!nonnegative_finite(phi[0]) || !std::isfinite(phi[1]) ||
              !nonnegative_finite(kappa)) return Status::invalid_kernel;
          transport += (s == 2 ? 0.5 : 1.0) * kappa;
        }
        scat += w * e * transport;
      }
    }
    if (!nonnegative_finite(n) || !nonnegative_finite(j))
      return Status::invalid_moment;
    if (n == 0 || j == 0)
      return Status::insufficient_energy_coverage;
    const double emin = v.energies[0], emax = v.energies[v.ne - 1];
    const double f0 = scaled_fermi(0, temp, mu);
    const double low2 = f0 * emin * emin * emin / 3;
    const double low3 = f0 * emin * emin * emin * emin / 4;
    const double tail_mu = mu < 0 ? 0 : mu;
    if (!((low2 + high_tail(2, emax, temp, tail_mu)) <= v.tail_tolerance * n &&
          (low3 + high_tail(3, emax, temp, tail_mu)) <= v.tail_tolerance * j))
      return Status::insufficient_energy_coverage;
    const double multiplicity = s == 2 ? 4 : 1;
    out.number[s] = restore_fd_scale(n, multiplicity * phase_space, log_scale);
    out.energy[s] = restore_fd_scale(j, multiplicity * phase_space, log_scale);
    out.mean_energy[s] = j / n;
    out.absorption_number[s] = a0 / n;
    out.absorption_energy[s] = a1 / j;
    out.scattering[s] = scat / j;
    out.emission_number[s] = restore_fd_scale(a0, multiplicity * phase_space, log_scale);
    out.emission_energy[s] = restore_fd_scale(a1, multiplicity * phase_space, log_scale);
    if (!nonnegative_finite(out.number[s]) || !nonnegative_finite(out.energy[s]) ||
        !nonnegative_finite(out.mean_energy[s]) ||
        !nonnegative_finite(out.absorption_number[s]) ||
        !nonnegative_finite(out.absorption_energy[s]) ||
        !nonnegative_finite(out.scattering[s]) ||
        !nonnegative_finite(out.emission_number[s]) ||
        !nonnegative_finite(out.emission_energy[s])) return Status::invalid_moment;
  }
  return Status::success;
}

void self_test();
void energy_node_self_test();
} // namespace nuX_WeakLib
#endif
