#ifndef NUX_WEAKLIB_COLLISION_HXX
#define NUX_WEAKLIB_COLLISION_HXX

#include "nuX_weaklib.hxx"
#include "cctk.h"
#include "nuX_weakrates.hxx"

namespace nuX_WeakLib {

// Fixed EOSX matter; cgs density, MeV temperature/chemical potentials.
struct MatterState {
  double rho, temperature, ye, mu_nu, mu_e, xp, xn;
};

// Number/energy are the actual comoving grey moments in cm^-3/MeV cm^-3.
// Heavy moments sum four species, but their scalar FD parameters are per
// member. The angular transport closure remains independent of this ansatz.
struct RadiationState {
  double number[3]{}, energy[3]{}, temperature[3]{}, eta[3]{};
};

namespace collision_detail {

struct Integral {
  double emission_number = 0, emission_energy = 0;
  double loss_number = 0, loss_energy = 0, scattering = 0;
};

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
double occupation(const RadiationState &r, int species, double energy) noexcept {
  return r.number[species] == 0 ? 0 :
      fermi(energy / r.temperature[species] - r.eta[species]);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
double vacancy(const RadiationState &r, int species, double energy) noexcept {
  // Do not subtract a nearly filled FD occupation from one: that loses a
  // representable hole long before exp(-|E/T-eta|) underflows.
  return r.number[species] == 0 ? 1 :
      fermi(r.eta[species] - energy / r.temperature[species]);
}

// Native Pair/NES index order is (Eprime,E); the upper Pair triangle also
// exchanges coupling component I/II and its offset. The zeroth Legendre
// coefficient has units cm^-1 MeV^-3, not an already integrated opacity.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
double pair_nes_kernel(const TableView &v, bool pair, bool heavy, bool anti,
                       int output, int partner, double temp, double eta_e) noexcept {
  const auto t = pair ? v.pair : v.nes;
  const auto ts = pair ? v.pair_temp : v.nes_temp;
  const auto es = pair ? v.pair_eta : v.nes_eta;
  const auto os = pair ? v.pair_offsets : v.nes_offsets;
  constexpr double cv = 0.96, ca = 0.5;
  const double c1 = heavy ? (cv + ca - 2) * (cv + ca - 2) : (cv + ca) * (cv + ca);
  const double c2 = (cv - ca) * (cv - ca);
  const double lt = std::log10(temp / kboltz_mev_kelvin);
  const double le = std::log10(eta_e);
  double value[2];
  for (int c = 0; c < 2; ++c) {
    if (partner <= output) {
      value[c] = wli::NESPairInterpolateSingleVariable2D2DAlignedPoint(
          lt, le, ts.log_nodes, ts.n, es.log_nodes, es.n,
          partner, output, v.ne, v.ne, c, t.n[2], os[c], t.data);
    } else if (pair) {
      value[c] = wli::PairCrossingSymmetryFillPoint(
          lt, le, ts.log_nodes, ts.n, es.log_nodes, es.n,
          partner, output, v.ne, v.ne, 1 - c, t.n[2], os[1 - c], t.data);
    } else {
      value[c] = wli::NESDetailedBalanceFillPoint(
          lt, le, ts.log_nodes, ts.n, es.log_nodes, es.n,
          partner, output, v.ne, v.ne, c, t.n[2], os[c], t.data,
          v.energies, temp);
    }
  }
  return 4 * pi * (anti ? c2 * value[0] + c1 * value[1]
                                 : c1 * value[0] + c2 * value[1]);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
double brem_kernel(const TableView &v, int i, int j, const MatterState &m) noexcept {
  constexpr double ca = -1.26 / 2, gf = 1.166e-11;
  // Match the native thornado HR98 response normalization, including its
  // declared Avogadro/hbar constants; this is not the EOSX baryon mass.
  constexpr double nref = 6.02214199e23 * 1e14;
  constexpr double hbar_c = 6.582119569e-22 * 2.99792458e10;
  constexpr double normalization = ca * ca * gf * gf * nref * hbar_c * hbar_c /
                                   ((2 * pi) * (2 * pi) * (2 * pi));
  const double fractions[3] = {m.xp, m.xn, std::sqrt(m.xp * m.xn)};
  const double weights[3] = {1, 1, 28.0 / 3};
  double response = 0;
  for (int a = 0; a < 3; ++a) {
    if (fractions[a] == 0) continue; // genuinely absent constituent
    const double d = m.rho * fractions[a];
    if (!v.brem_rho.contains(d))
      return std::numeric_limits<double>::quiet_NaN();
    const double component = wli::BremInterpolateSingleDensity2DAlignedPoint(
        std::log10(d), std::log10(m.temperature / kboltz_mev_kelvin),
        v.brem_rho.log_nodes, v.brem_rho.n, v.brem_temp.log_nodes, v.brem_temp.n,
        j, i, v.ne, v.ne, 0, v.brem.n[2], v.brem_offset, v.brem.data);
    if (!nonnegative_finite(component))
      return std::numeric_limits<double>::quiet_NaN();
    response += weights[a] * component;
  }
  // S_sigma already contains rho_14. Do not multiply by local nb again.
  return 4 * pi * 3 * normalization * response;
}

// A pair event is visited once, adding identical particle-number creation
// and destruction to its neutrino and antineutrino. Loss includes stimulated
// absorption, A+B, not just B; this also gives monopole angular damping.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void pair_event(double weight, double production, double annihilation,
                double gi, double gj, double bi, double bj, double ei, double ej,
                Integral &nu, Integral &anti, double multiplicity = 1) noexcept {
  const double a = multiplicity * weight * production;
  const double b = multiplicity * weight * annihilation;
  const double eta_i = a * bj, eta_j = a * bi;
  const double loss_i = (eta_i + b * gj) * gi;
  const double loss_j = (eta_j + b * gi) * gj;
  nu.emission_number += eta_i;
  anti.emission_number += eta_j;
  nu.loss_number += loss_i;
  anti.loss_number += loss_j;
  nu.emission_energy += ei * eta_i;
  anti.emission_energy += ej * eta_j;
  nu.loss_energy += ei * loss_i;
  anti.loss_energy += ej * loss_j;
}

// Analytic, monoenergetic plasmon approximation using ONLY the WeakRates
// plasmon formula. Its inverse is constructed from the same photon bath and
// actual partner FD occupations, not independent grey Kirchhoff corrections.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void plasmon(const MatterState &m, const RadiationState &r, Integral out[3]) noexcept {
  constexpr double cv = nuX_WeakRates::detail::cv;
  const double eta_e = m.mu_e / m.temperature;
  const auto kernel = nuX_WeakRates::detail::plasmon_kernel(m.temperature, eta_e);
  const double gamma = kernel.gamma, energy = kernel.energy;
  const double prefactor = kernel.prefactor / nuX_WeakRates::detail::clight;
  const double production = prefactor * std::exp(-gamma) * (1 + gamma);
  // 2*energy/T-gamma = 1+1/(1+gamma), evaluated without huge exponentials.
  const double annihilation = prefactor * std::exp(1 + 1 / (1 + gamma)) * (1 + gamma);
  const double nue = occupation(r, 0, energy), anti = occupation(r, 1, energy);
  pair_event(cv * cv, production, annihilation, nue, anti,
             vacancy(r, 0, energy), vacancy(r, 1, energy), energy, energy,
             out[0], out[1]);
  Integral heavy_nu, heavy_anti;
  const double heavy = occupation(r, 2, energy);
  pair_event((cv - 1) * (cv - 1), production, annihilation, heavy, heavy,
             vacancy(r, 2, energy), vacancy(r, 2, energy),
             energy, energy, heavy_nu, heavy_anti, 2);
  out[2].emission_number += heavy_nu.emission_number + heavy_anti.emission_number;
  out[2].loss_number += heavy_nu.loss_number + heavy_anti.loss_number;
  out[2].emission_energy += heavy_nu.emission_energy + heavy_anti.emission_energy;
  out[2].loss_energy += heavy_nu.loss_energy + heavy_anti.loss_energy;
}
} // namespace collision_detail

// Opt-in monopole collision model. No per-cell energy-pair matrix is stored.
// The radiation is reconstructed anew for every nonlinear trial by the caller.
// For h(E)=g(E)*H/J, the energy-weighted A+B gives the grey flux damping.
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_NOINLINE inline
Status evaluate_collision(const TableView &v, const MatterState &m,
                          const RadiationState &r, Rates &out) noexcept {
  using namespace collision_detail;
  out = {};
  if (!(std::isfinite(m.rho) && m.rho > 0 && std::isfinite(m.temperature) &&
        m.temperature > 0 && std::isfinite(m.ye) && std::isfinite(m.mu_nu) &&
        std::isfinite(m.mu_e) && v.ne >= 2) ||
      (v.include_brem && !(nonnegative_finite(m.xp) && m.xp <= 1 &&
        nonnegative_finite(m.xn) && m.xn <= 1)))
    return Status::invalid_state;
  const double kelvin = m.temperature / kboltz_mev_kelvin;
  const double eta_e = m.mu_e / m.temperature;
  if ((v.include_emab && !v.emab_thermo.contains(m.rho, kelvin, m.ye)) ||
      (v.include_iso && !v.iso_thermo.contains(m.rho, kelvin, m.ye)) ||
      (v.include_pair && !(v.pair_temp.contains(kelvin) && v.pair_eta.contains(eta_e))) ||
      (v.include_nes && !(v.nes_temp.contains(kelvin) && v.nes_eta.contains(eta_e))) ||
      (v.include_brem && !v.brem_temp.contains(kelvin)))
    return Status::out_of_range;
  for (int s = 0; s < 3; ++s) {
    if (!nonnegative_finite(r.number[s]) || !nonnegative_finite(r.energy[s]) ||
        (r.number[s] == 0) != (r.energy[s] == 0) ||
        !(std::isfinite(r.temperature[s]) && r.temperature[s] > 0 && std::isfinite(r.eta[s])))
      return Status::invalid_state;
    out.number[s] = r.number[s];
    out.energy[s] = r.energy[s];
    out.mean_energy[s] = r.number[s] > 0 ? r.energy[s] / r.number[s] : 0;
    if (!nonnegative_finite(out.mean_energy[s])) return Status::invalid_moment;
    if (r.number[s] > 0) {
      const double t = r.temperature[s], mu = r.eta[s] * t;
      const double log_scale = r.eta[s] < 0 ? r.eta[s] : 0;
      const double emin = v.energies[0], emax = v.energies[v.ne - 1];
      const double f0 = scaled_fermi(0, t, mu);
      const double low2 = f0 * emin * emin * emin / 3;
      const double low3 = f0 * emin * emin * emin * emin / 4;
      const double tail_mu = mu < 0 ? 0 : mu;
      const double flavors = s == 2 ? 4 : 1;
      const double n_tail = restore_fd_scale(low2 + high_tail(2, emax, t, tail_mu),
                                              flavors * phase_space, log_scale);
      const double j_tail = restore_fd_scale(low3 + high_tail(3, emax, t, tail_mu),
                                              flavors * phase_space, log_scale);
      if (!(n_tail <= v.tail_tolerance * r.number[s] &&
            j_tail <= v.tail_tolerance * r.energy[s]))
        return Status::insufficient_energy_coverage;
    }
  }
  Integral integral[3];
  const double ld = std::log10(m.rho), lt = std::log10(kelvin);
  for (int i = 0; i < v.ne; ++i) {
    const double e = v.energies[i], q = v.weights[i] * e * e;
    for (int s = 0; s < 3; ++s) {
      const double g = occupation(r, s, e), multiplicity = s == 2 ? 4 : 1;
      const double weight = multiplicity * phase_space * q;
      if (v.include_emab && s < 2) {
        const double chi = emab_at_energy_node(v, s, i, ld, lt, m.ye);
        if (!nonnegative_finite(chi)) return Status::invalid_kernel;
        const double mu = s == 0 ? m.mu_nu : -m.mu_nu;
        const double eta = restore_fd_scale(q * chi * scaled_fermi(e, m.temperature, mu),
            multiplicity * phase_space, mu < 0 ? mu / m.temperature : 0);
        integral[s].emission_number += eta;
        integral[s].emission_energy += e * eta;
        const double loss = r.number[s] == 0 ? 0 : restore_fd_scale(
            q * chi * scaled_fermi(e, r.temperature[s], r.eta[s] * r.temperature[s]),
            multiplicity * phase_space, r.eta[s] < 0 ? r.eta[s] : 0);
        integral[s].loss_number += loss;
        integral[s].loss_energy += e * loss;
      }
      if (v.include_iso) {
        double transport = 0;
        const int first = s == 2 ? 0 : s, last = s == 2 ? 1 : s;
        for (int a = first; a <= last; ++a) {
          const double p0 = iso_at_energy_node(v, a, 0, i, ld, lt, m.ye);
          const double p1 = iso_at_energy_node(v, a, 1, i, ld, lt, m.ye);
          const double rate = 4 * pi * e * e * (p0 - p1 / 3);
          if (!nonnegative_finite(p0) || !std::isfinite(p1) || !nonnegative_finite(rate))
            return Status::invalid_kernel;
          transport += (s == 2 ? 0.5 : 1) * rate;
        }
        integral[s].scattering += weight * e * g * transport;
      }
    }
    if (!(v.include_pair || v.include_brem || v.include_nes)) continue;
    for (int j = 0; j < v.ne; ++j) {
      const double ep = v.energies[j], qp = v.weights[j] * ep * ep;
      const double weight = phase_space * q * qp;
      const double brem = v.include_brem ? brem_kernel(v, i, j, m) : 0;
      if (!nonnegative_finite(brem)) return Status::invalid_kernel;
      if (v.include_pair || v.include_brem) {
        for (int heavy = 0; heavy < 2; ++heavy) {
          const double pair = v.include_pair ?
              pair_nes_kernel(v, true, heavy, false, i, j, m.temperature, eta_e) : 0;
          if (!nonnegative_finite(pair)) return Status::invalid_kernel;
          const double ann = brem + pair;
          if (!nonnegative_finite(ann)) return Status::invalid_kernel;
          const double pro = ann * std::exp(-(e + ep) / m.temperature);
          if (!heavy)
            pair_event(weight, pro, ann, occupation(r, 0, e), occupation(r, 1, ep),
                       vacancy(r, 0, e), vacancy(r, 1, ep),
                       e, ep, integral[0], integral[1]);
          else {
            Integral nu, anti;
            pair_event(weight, pro, ann, occupation(r, 2, e), occupation(r, 2, ep),
                       vacancy(r, 2, e), vacancy(r, 2, ep), e, ep, nu, anti, 2);
            integral[2].emission_number += nu.emission_number + anti.emission_number;
            integral[2].emission_energy += nu.emission_energy + anti.emission_energy;
            integral[2].loss_number += nu.loss_number + anti.loss_number;
            integral[2].loss_energy += nu.loss_energy + anti.loss_energy;
          }
        }
      }
      if (v.include_nes) {
        for (int s = 0; s < 3; ++s) {
          const double gi = occupation(r, s, e), gj = occupation(r, s, ep);
          const int first = s == 2 ? 0 : s, last = s == 2 ? 1 : s;
          double down = 0, up = 0;
          for (int anti = first; anti <= last; ++anti) {
            const double factor = s == 2 ? 2 : 1;
            // Native Phi(Eprime,E) is OUTGOING from E into Eprime.
            // Query its transpose for scattering into the current E bin.
            // The upper destination-energy triangle is Boltzmann suppressed.
            const double incoming = pair_nes_kernel(v, false, s == 2, anti,
                j, i, m.temperature, eta_e);
            const double outgoing = pair_nes_kernel(v, false, s == 2, anti,
                i, j, m.temperature, eta_e);
            if (!nonnegative_finite(incoming) || !nonnegative_finite(outgoing))
              return Status::invalid_kernel;
            down += factor * incoming;
            up += factor * outgoing;
          }
          const double emission = weight * down * gj;
          const double loss = weight * (down * gj + up * vacancy(r, s, ep)) * gi;
          // NES does not alter the integrated neutrino number, exactly.
          integral[s].emission_energy += e * emission;
          integral[s].loss_energy += e * loss;
        }
      }
    }
  }
  if (v.include_plasmon) plasmon(m, r, integral);
  for (int s = 0; s < 3; ++s) {
    const auto &a = integral[s];
    if (!nonnegative_finite(a.emission_number) || !nonnegative_finite(a.emission_energy) ||
        !nonnegative_finite(a.loss_number) || !nonnegative_finite(a.loss_energy) ||
        !nonnegative_finite(a.scattering)) return Status::invalid_kernel;
    out.emission_number[s] = a.emission_number;
    out.emission_energy[s] = a.emission_energy;
    out.absorption_number[s] = r.number[s] > 0 ? a.loss_number / r.number[s] : 0;
    out.absorption_energy[s] = r.energy[s] > 0 ? a.loss_energy / r.energy[s] : 0;
    out.scattering[s] = r.energy[s] > 0 ? a.scattering / r.energy[s] : 0;
    if (!nonnegative_finite(out.absorption_number[s]) ||
        !nonnegative_finite(out.absorption_energy[s]) ||
        !nonnegative_finite(out.scattering[s])) return Status::invalid_kernel;
  }
  return Status::success;
}

} // namespace nuX_WeakLib

#endif
