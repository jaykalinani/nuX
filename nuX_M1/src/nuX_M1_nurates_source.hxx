#ifndef NUX_M1_NURATES_SOURCE_HXX
#define NUX_M1_NURATES_SOURCE_HXX

// Reusable local, matter-frozen NuRates collision solve. The including source
// must enable NUX_M1_CLOSURE_IMPLEMENTATION and NUX_M1_SOURCES_IMPLEMENTATION
// before its first inclusion of nuX_M1_sources.hxx.
#include "nuX_M1_sources.hxx"
#include "nuX_M1_nurates.hxx"

namespace nuX_M1 {

// E, F_d and N are DENSITIZED evolved-group code-unit moments. Species 2 is
// the sum of four heavy flavors. chi is the selected closure's Eddington
// factor, not an independently evolved collision variable.
struct nurates_source_state_t {
  CCTK_REAL N, E;
  tensor::generic<CCTK_REAL, 4, 1> F_d;
  CCTK_REAL chi;
};

struct nurates_source_options_t {
  closure_t closure = CLOSURE_MINERBO;
  CCTK_INT maxiter = 64;
  CCTK_INT inner_maxiter = 100;
  CCTK_INT angular_maxiter = 32;
  CCTK_REAL epsrel = 1.0e-9;
  // Retained as a caller convenience, but deliberately NOT used to loosen
  // the radiation-scale residual: legacy source_epsabs can exceed E itself.
  CCTK_REAL epsabs = 0;
  // Relative step tolerance in the six (log Tnu, eta) coordinates. A
  // nonpositive value selects epsrel. Residual AND step tests are required.
  CCTK_REAL step_tolerance = 0;
  // False is a caller-proven invariant vacuum (e.g. heavy flavors in a
  // beta-only problem). Such a group is exactly zero, not an FD floor.
  // Every residual/final verification still requires its emission to vanish.
  bool active_species[3] = {true, true, true};
  // Optional LTE-like guess, ONLY a starting point. A failed attempt retries
  // from the positive baseline spectrum and mixed starting points.
  bool use_initial_spectrum = false;
  CCTK_REAL initial_log_temperature[3] = {};
  CCTK_REAL initial_eta[3] = {};
};

namespace nurates_source_detail {

using spectral_vector = Arith::vec<CCTK_REAL, 6>;
using spectral_matrix = Arith::mat<CCTK_REAL, 6>;

CCTK_HOST CCTK_DEVICE inline CCTK_REAL norm(const spectral_vector &v) {
  CCTK_REAL result = 0;
  for (int a = 0; a < 6; ++a)
    result = fmax(result, fabs(v(a)));
  return result;
}

CCTK_HOST CCTK_DEVICE inline bool valid(const nurates_opacities_t &r) {
  return isfinite(r.abs_0) && r.abs_0 >= 0 &&
         isfinite(r.abs_1) && r.abs_1 >= 0 &&
         isfinite(r.scat_1) && r.scat_1 >= 0 &&
         isfinite(r.eta_0) && r.eta_0 >= 0 &&
         isfinite(r.eta_1) && r.eta_1 >= 0;
}

// The small populations are represented in log coordinates; no lower
// population, temperature or degeneracy clamp is imposed. Failure/underflow
// rejects a trial and lets the outer line search choose another step.
CCTK_HOST CCTK_DEVICE inline bool spectral_moments(
    const spectral_vector &x, const nurates_source_options_t &options,
    M1Quantities &radiation) {
  for (int sp = 0; sp < 3; ++sp) {
    if (!options.active_species[sp]) {
      radiation.n[sp] = radiation.J[sp] = 0;
      continue;
    }
    const CCTK_REAL f2 = FDI_p2(x(2 * sp + 1));
    const CCTK_REAL f3 = FDI_p3(x(2 * sp + 1));
    if (!(isfinite(x(2 * sp)) && isfinite(x(2 * sp + 1)) &&
          isfinite(f2) && f2 > 0 && isfinite(f3) && f3 > 0))
      return false;
    radiation.n[sp] = exp(log(kBS_FourPi_hc3) + 3 * x(2 * sp) + log(f2));
    radiation.J[sp] = exp(log(kBS_FourPi_hc3) + 4 * x(2 * sp) + log(f3));
    if (!(isfinite(radiation.n[sp]) && radiation.n[sp] > 0 &&
          isfinite(radiation.J[sp]) && radiation.J[sp] > 0))
      return false;
  }
  radiation.n[3] = radiation.n[2];
  radiation.J[3] = radiation.J[2];
  return true;
}

CCTK_HOST CCTK_DEVICE inline bool vacuum(const nurates_source_state_t &s) {
  if (s.N != 0 || s.E != 0)
    return false;
  for (int a = 0; a < 4; ++a)
    if (s.F_d(a) != 0)
      return false;
  return true;
}

CCTK_HOST CCTK_DEVICE inline CCTK_REAL invariant(
    const nurates_source_state_t &s, const SourceUpdateContext &g) {
  return s.E - tensor::dot(s.F_d, g.v_u);
}

// This check cancels elastic scattering ANALYTICALLY. A dt*ks*E roundoff
// allowance here would erase the invariant in precisely the stiffest cells.
CCTK_HOST CCTK_DEVICE inline bool verify_invariant(
    const SourceUpdateContext &g, const nurates_source_state_t &old,
    const nurates_source_state_t &next, const CCTK_REAL J,
    const CCTK_REAL densitized_eta, const CCTK_REAL kabs,
    const CCTK_REAL epsrel) {
  const CCTK_REAL before = invariant(old, g), after = invariant(next, g);
  const CCTK_REAL h = g.alp * g.cdt / g.W;
  const CCTK_REAL change = h * (densitized_eta - kabs * J);
  CCTK_REAL subtraction_scale = fabs(old.E) + fabs(next.E);
  for (int a = 1; a < 4; ++a)
    subtraction_scale += fabs(g.v_u(a) * old.F_d(a)) +
                         fabs(g.v_u(a) * next.F_d(a));
  const CCTK_REAL tolerance = epsrel * fmax(fabs(before), fabs(after)) +
      256 * std::numeric_limits<CCTK_REAL>::epsilon() *
          (subtraction_scale + fabs(h * densitized_eta) + fabs(h * kabs * J));
  return isfinite(change) && isfinite(tolerance) &&
         fabs(after - before - change) <= tolerance;
}

// Strict closure evaluation, including for trial states. There is no
// alternate-closure fallback, atmosphere projection or repaired source state.
CCTK_HOST CCTK_DEVICE inline bool project(
    const SourceUpdateContext &g, const nurates_source_options_t &o,
    nurates_source_state_t &s, CCTK_REAL &J, CCTK_REAL &Gamma,
    tensor::generic<CCTK_REAL, 4, 1> *H_out = nullptr) {
  if (!(isfinite(s.N) && s.N >= 0 &&
        source_state_is_realizable(g.g_uu, s.E, s.F_d)))
    return false;
  tensor::symmetric2<CCTK_REAL, 4, 2> P_dd, T_dd;
  bool closure_valid = true;
  const CCTK_REAL closure_tol = fmin(g.closure_epsilon, 0.05 * o.epsrel);
  calc_closure(g.cctkGH, g.i, g.j, g.k, g.ig, o.closure, g.g_dd, g.g_uu,
               g.n_d, g.W, g.u_u, g.v_d, g.proj_ud, s.E, s.F_d, &s.chi,
               &P_dd, closure_tol, g.closure_maxiter, false, &closure_valid);
  if (!closure_valid || !isfinite(s.chi))
    return false;
  assemble_rT(g.n_d, s.E, s.F_d, P_dd, &T_dd);
  J = calc_J_from_rT(T_dd, g.u_u);
  tensor::generic<CCTK_REAL, 4, 1> H_d;
  calc_H_from_rT(T_dd, g.u_u, g.proj_ud, &H_d);
  if (!comoving_state_is_acceptable(o.closure, g.g_uu, J, H_d))
    return false;
  // Match the scheduled NuRates source/current map, without transport floors.
  // A vanishing energy population cannot carry a positive number population.
  if ((J == 0 && s.N != 0) ||
      !compute_source_Gamma(g.W, g.v_u, J, s.E, s.F_d, &Gamma))
    return false;
  if (H_out)
    *H_out = H_d;
  return true;
}

CCTK_HOST CCTK_DEVICE inline bool physical_moments(
    const SourceUpdateContext &g, const CCTK_REAL volume,
    const nurates_source_options_t &o, nurates_source_state_t states[3],
    M1Quantities &radiation) {
  for (int sp = 0; sp < 3; ++sp) {
    CCTK_REAL J, Gamma;
    if (!project(g, o, states[sp], J, Gamma))
      return false;
    const CCTK_REAL flavors = sp == 2 ? 4 : 1;
    radiation.n[sp] = states[sp].N / (Gamma * volume) *
                       (nuX_ndens_conv / flavors);
    radiation.J[sp] = J / volume * (nuX_edens_conv / flavors);
    radiation.chi[sp] = states[sp].chi;
    if (!(isfinite(radiation.n[sp]) && radiation.n[sp] >= 0 &&
          isfinite(radiation.J[sp]) && radiation.J[sp] >= 0))
      return false;
  }
  radiation.n[3] = radiation.n[2];
  radiation.J[3] = radiation.J[2];
  radiation.chi[3] = radiation.chi[2];
  return true;
}

template <class Rates> struct spectral_source_t {
  const SourceUpdateContext &geometry;
  const CCTK_REAL volume;
  const nurates_source_state_t *base;
  const Rates &rates;
  const nurates_source_options_t &options;

  // All species' coefficients are recomputed simultaneously for each trial
  // spectrum and each angular fixed-point iteration. Nothing is cached
  // between independent residual/Jacobian evaluations.
  CCTK_HOST CCTK_DEVICE bool evaluate(
      const spectral_vector &x, spectral_vector &residual,
      nurates_source_state_t *accepted = nullptr) const {
    M1Quantities radiation = {};
    if (!spectral_moments(x, options, radiation))
      return false;
    CCTK_REAL chi[3] = {base[0].chi, base[1].chi, base[2].chi};
    const CCTK_REAL lapse_dt = geometry.alp * geometry.cdt;
    const CCTK_REAL eps = std::numeric_limits<CCTK_REAL>::epsilon();
    for (int iteration = 0; iteration < options.angular_maxiter; ++iteration) {
      for (int sp = 0; sp < 4; ++sp)
        radiation.chi[sp] = chi[sp < 3 ? sp : 2];
      nurates_opacities_t coefficients[3];
      if (!rates(radiation, coefficients))
        return false;
      nurates_source_state_t states[3];
      CCTK_REAL angular_error = 0;
      for (int sp = 0; sp < 3; ++sp) {
        const auto &r = coefficients[sp];
        if (!valid(r))
          return false;
        const auto &g = geometry;
        const auto &old = base[sp];
        if (!options.active_species[sp]) {
          if (!vacuum(old) || r.eta_0 != 0 || r.eta_1 != 0)
            return false;
          states[sp] = old;
          // Anchor unused coordinates so the six-by-six Jacobian is
          // nonsingular while their physical populations stay exactly zero.
          residual(2 * sp) = x(2 * sp);
          residual(2 * sp + 1) = x(2 * sp + 1);
          continue;
        }
        const CCTK_REAL eta = volume * r.eta_1;
        auto &s = states[sp];
        // A realizable LTE-like INITIAL GUESS only. source_update retries
        // the identical selected-closure problem from the baseline if useful.
        // The amplitude is seeded by the scattering-free invariant, NOT by
        // an energy equation whose enormous scattering terms cancel.
        const CCTK_REAL proper_dt = lapse_dt / g.W;
        const CCTK_REAL Jguess = (invariant(old, g) + proper_dt * eta) /
                                 (1 + proper_dt * r.abs_1);
        s.E = (4 * g.W * g.W - 1) / 3 * Jguess;
        const CCTK_REAL factor = 4 * g.W * g.W / 3 * Jguess;
        pack_F_d(-g.alp * g.n_u(1), -g.alp * g.n_u(2), -g.alp * g.n_u(3),
                 factor * g.v_d(1), factor * g.v_d(2),
                 factor * g.v_d(3), &s.F_d);
        s.chi = chi[sp];
        const CCTK_REAL energy_scale = fmax(fabs(old.E), fabs(s.E));
        CCTK_REAL component_abs[4];
        component_abs[0] = 64 * eps *
            (fabs(old.E) + fabs(s.E) + fabs(proper_dt * eta) +
             fabs(proper_dt * r.abs_1 * Jguess));
        for (int a = 1; a < 4; ++a) {
          // At rest this depends on Fi, NEVER on E. In a moving frame only
          // the velocity-coupled energy can feed cancellation into this Fi.
          component_abs[a] = 64 * eps *
              (fabs(old.F_d(a)) + fabs(s.F_d(a)) +
               fabs(g.W * g.W * g.W * g.v_d(a)) * energy_scale);
        }
        for (int a = 0; a < 4; ++a)
          component_abs[a] = fmax(component_abs[a],
                                  std::numeric_limits<CCTK_REAL>::min());
        if (!(isfinite(s.E) && s.E >= 0 && isfinite(component_abs[0])))
          return false;
        const SourceUpdateContext ctx(
            g.cctkGH, g.i, g.j, g.k, sp,
            fmin(g.closure_epsilon, 0.05 * options.epsrel),
            g.closure_maxiter, false, g.cdt, g.alp, g.g_dd, g.g_uu, g.n_d,
            g.n_u, g.gamma_ud, g.u_d, g.u_u, g.v_d, g.v_u, g.proj_ud, g.W,
            old.E, old.F_d, old.E, old.F_d, eta, r.abs_1, r.scat_1,
            true, component_abs, true);
        const int status = source_update(
            ctx, options.closure, &s.chi, &s.E, &s.F_d, 1, -1, -1,
            options.inner_maxiter, 0, 0.05 * options.epsrel);
        if (status != NUX_M1_SOURCE_OK)
          return false;
        CCTK_REAL J, Gamma;
        // N is assigned after the energy/momentum solve. project only needs
        // a nonnegative N here; its exact value does not affect J or Gamma.
        s.N = old.N;
        if (!project(g, options, s, J, Gamma))
          return false;
        if (!verify_invariant(g, old, s, J, eta, r.abs_1, options.epsrel))
          return false;
        s.N = (old.N + lapse_dt * volume * r.eta_0) /
              (1 + lapse_dt * r.abs_0 / Gamma);
        const CCTK_REAL flavors = sp == 2 ? 4 : 1;
        const CCTK_REAL n = s.N / (Gamma * volume) *
                            (nuX_ndens_conv / flavors);
        const CCTK_REAL j = J / volume * (nuX_edens_conv / flavors);
        if (!(isfinite(s.N) && s.N >= 0 && isfinite(n) && n >= 0 &&
              isfinite(j) && j >= 0))
          return false;
        residual(2 * sp) =
            (radiation.n[sp] - n) / fmax(radiation.n[sp], n);
        residual(2 * sp + 1) =
            (radiation.J[sp] - j) / fmax(radiation.J[sp], j);
        angular_error = fmax(angular_error, fabs(s.chi - chi[sp]));
        chi[sp] = s.chi;
      }
      if (angular_error <= fmax(0.05 * options.epsrel, 64 * eps)) {
        if (!nuX_Utils::roots::finite(residual))
          return false;
        if (accepted)
          for (int sp = 0; sp < 3; ++sp)
            accepted[sp] = states[sp];
        return true;
      }
    }
    return false;
  }

  CCTK_HOST CCTK_DEVICE bool jacobian(
      const spectral_vector &x, const spectral_vector &f,
      spectral_matrix &matrix) const {
    for (int column = 0; column < 6; ++column) {
      if (!options.active_species[column / 2]) {
        for (int row = 0; row < 6; ++row)
          matrix(row, column) = row == column ? 1 : 0;
        continue;
      }
      const CCTK_REAL h = 1e-5 * fmax(CCTK_REAL(1), fabs(x(column)));
      spectral_vector xp = x, xm = x, fp, fm;
      xp(column) += h;
      xm(column) -= h;
      const bool plus = evaluate(xp, fp), minus = evaluate(xm, fm);
      if (!plus && !minus)
        return false;
      for (int row = 0; row < 6; ++row)
        matrix(row, column) = plus && minus ? (fp(row) - fm(row)) / (2 * h)
                               : plus ? (fp(row) - f(row)) / h
                                      : (f(row) - fm(row)) / h;
    }
    return nuX_Utils::roots::finite(matrix);
  }

  // Recompute coefficients from accepted physical moments, not an out-of-date
  // spectral trial. Independently check the scattering-free invariant with
  // no kscat allowance. Unscaled E/Fi equations are checked componentwise
  // with their arithmetic backward-error bounds: a stiff cancellation can
  // limit evaluation of those rows, but NEVER loosens the invariant check.
  CCTK_HOST CCTK_DEVICE bool verify(nurates_source_state_t states[3]) const {
    M1Quantities radiation = {};
    if (!physical_moments(geometry, volume, options, states, radiation))
      return false;
    nurates_opacities_t coefficients[3];
    if (!rates(radiation, coefficients))
      return false;
    const auto &g = geometry;
    const CCTK_REAL dt = g.cdt, lapse_dt = dt * g.alp;
    const CCTK_REAL eps = std::numeric_limits<CCTK_REAL>::epsilon();
    for (int sp = 0; sp < 3; ++sp) {
      const auto &r = coefficients[sp];
      if (!valid(r))
        return false;
      auto &s = states[sp];
      if (!options.active_species[sp]) {
        if (!vacuum(base[sp]) || !vacuum(s) || r.eta_0 != 0 || r.eta_1 != 0)
          return false;
        continue;
      }
      CCTK_REAL J, Gamma;
      tensor::generic<CCTK_REAL, 4, 1> H_d, S_d, Fsource;
      if (!project(g, options, s, J, Gamma, &H_d))
        return false;
      if (!verify_invariant(g, base[sp], s, J, volume * r.eta_1,
                             r.abs_1, options.epsrel))
        return false;
      calc_rad_sources(volume * r.eta_1, r.abs_1, r.scat_1, g.u_d, J, H_d, &S_d);
      const CCTK_REAL Esource = calc_rE_source(g.alp, g.n_u, S_d);
      calc_rF_source(g.alp, g.gamma_ud, S_d, &Fsource);
      const CCTK_REAL number_change =
          lapse_dt * (volume * r.eta_0 - r.abs_0 * s.N / Gamma);
      const CCTK_REAL number_scale = fmax(fabs(base[sp].N), fabs(s.N));
      const CCTK_REAL number_tol = options.epsrel * number_scale +
          256 * eps * (fabs(s.N) + fabs(base[sp].N) +
                       fabs(lapse_dt * volume * r.eta_0) +
                       fabs(lapse_dt * r.abs_0 * s.N / Gamma));
      if (!(isfinite(number_change) &&
            fabs(s.N - base[sp].N - number_change) <= number_tol))
        return false;
      for (int a = 0; a < 4; ++a) {
        const CCTK_REAL next = a == 0 ? s.E : s.F_d(a);
        const CCTK_REAL old = a == 0 ? base[sp].E : base[sp].F_d(a);
        const CCTK_REAL increment = dt * (a == 0 ? Esource : Fsource(a));
        const CCTK_REAL component_scale = fmax(fabs(old), fabs(next));
        const CCTK_REAL energy_scale = fmax(fabs(base[sp].E), fabs(s.E));
        const CCTK_REAL coupling = a == 0 ? g.W :
            fabs(g.W * g.W * g.W * g.v_d(a));
        const CCTK_REAL flux_scale = a == 0 ? energy_scale : fabs(s.F_d(a));
        const CCTK_REAL tolerance = options.epsrel * component_scale +
            256 * eps * (fabs(next) + fabs(old) + fabs(increment) +
                         fabs(lapse_dt * coupling * volume * r.eta_1) +
                         fabs(lapse_dt * (r.abs_1 + r.scat_1) *
                              (coupling * energy_scale + g.W * flux_scale)));
        if (!(isfinite(increment) &&
              fabs(next - old - increment) <= tolerance))
          return false;
      }
    }
    return true;
  }

  CCTK_HOST CCTK_DEVICE bool solve(spectral_vector x,
                                   nurates_source_state_t result[3]) const {
    spectral_vector f;
    if (!evaluate(x, f))
      return false;
    const CCTK_REAL step_tol = options.step_tolerance > 0
                                 ? options.step_tolerance : options.epsrel;
    for (int iteration = 0; iteration < options.maxiter; ++iteration) {
      spectral_matrix J, Q, R;
      if (!jacobian(x, f, J) ||
          nuX_Utils::roots::detail::qr_decomp_unpack(J, &Q, &R) != ROOTS_SUCCESS)
        return false;
      spectral_vector rhs, step;
      nuX_Utils::roots::detail::compute_qtf(Q, f, &rhs);
      for (int a = 0; a < 6; ++a)
        rhs(a) = -rhs(a);
      if (nuX_Utils::roots::detail::r_solve(R, rhs, &step) != ROOTS_SUCCESS ||
          !nuX_Utils::roots::finite(step))
        return false;
      CCTK_REAL relative_step = 0;
      for (int a = 0; a < 6; ++a)
        relative_step = fmax(relative_step,
                            fabs(step(a)) / fmax(CCTK_REAL(1), fabs(x(a))));
      if (norm(f) <= options.epsrel && relative_step <= step_tol) {
        nurates_source_state_t candidate[3];
        if (!evaluate(x, f, candidate) || norm(f) > options.epsrel)
          return false;
        if (verify(candidate)) {
          for (int sp = 0; sp < 3; ++sp)
            result[sp] = candidate[sp];
          return true;
        }
        // A very stiff coefficient can amplify a small spectral residual.
        // Continue Newton refinement until the original physical source
        // residual also passes; never accept the spectral test alone.
      }
      const CCTK_REAL old_norm = norm(f);
      CCTK_REAL damping = 1;
      bool accepted = false;
      for (int backtrack = 0; backtrack < 24; ++backtrack) {
        spectral_vector trial, trial_f;
        for (int a = 0; a < 6; ++a)
          trial(a) = x(a) + damping * step(a);
        if (evaluate(trial, trial_f) &&
            norm(trial_f) < (1 - 1e-4 * damping) * old_norm) {
          x = trial;
          f = trial_f;
          accepted = true;
          break;
        }
        damping *= 0.5;
      }
      if (!accepted)
        return false;
    }
    return false;
  }
};

} // namespace nurates_source_detail

// Rates must be a device-callable, non-aborting function object:
//   bool operator()(M1Quantities const &physical_per_flavor,
//                   nurates_opacities_t *evolved_group_coefficients) const;
// n is in nm^-3 and J in MeV nm^-3; species 3 mirrors species 2. The callback
// owns the complete SAME microphysics/distribution reconstruction used by the
// scheduled opacity path, including all-species coupling and equilibrium
// corrections. Its code-unit eta coefficients include heavy-group factor 4;
// they are NOT densitized. Absorption/scattering are group-independent inverse
// code lengths. False means invalid trial, never an accepted no-source update.
// This solves that supplied collision MODEL; simultaneous nonlinear coupling
// does not itself impose pair-number conservation. In particular, a provider
// that independently replaces electron-neutrino/antineutrino thermal rates by
// legacy grey Kirchhoff equilibria need not retain identical pair-only net
// number sources away from equilibrium. No full pair/lepton-channel balance
// is claimed here; that requires a separately conservative pair closure in
// the provider (and channel-specific conservation tests).
// geometry.cdt is the coordinate diagonal-stage interval: lapse enters once.
// Returns true only for a verified source solution. result is untouched on
// failure. No floors, clipping or repairs are applied to accepted moments.
template <class Rates>
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_NOINLINE bool coupled_nurates_source(
    const SourceUpdateContext &geometry, const CCTK_REAL volume,
    const nurates_source_state_t base[3], const Rates &rates,
    const nurates_source_options_t &options, nurates_source_state_t result[3]) {
  using namespace nurates_source_detail;
  if (!(isfinite(volume) && volume > 0 && isfinite(geometry.cdt) &&
        geometry.cdt >= 0 && isfinite(geometry.alp) && geometry.alp > 0 &&
        isfinite(options.epsrel) && options.epsrel > 0 &&
        options.maxiter > 0 && options.inner_maxiter > 0 &&
        options.angular_maxiter > 0))
    return false;
  nurates_source_state_t baseline[3];
  for (int sp = 0; sp < 3; ++sp)
    baseline[sp] = base[sp];
  for (int sp = 0; sp < 3; ++sp)
    if (!options.active_species[sp] && !vacuum(baseline[sp]))
      return false;
  M1Quantities baseline_physical = {};
  if (!physical_moments(geometry, volume, options, baseline, baseline_physical))
    return false;
  if (geometry.cdt == 0) {
    for (int sp = 0; sp < 3; ++sp)
      result[sp] = baseline[sp];
    return true;
  }
  // This is the EXACT zero-collision operator, not an asymptotic shortcut.
  // It also supports a true vacuum without inventing a positive population.
  nurates_opacities_t initial_rates[3];
  if (rates(baseline_physical, initial_rates)) {
    bool zero = true;
    for (int sp = 0; sp < 3; ++sp)
      zero = zero && valid(initial_rates[sp]) &&
             source_rates_are_zero(initial_rates[sp].eta_0,
                                    initial_rates[sp].eta_1,
                                    initial_rates[sp].abs_0,
                                    initial_rates[sp].abs_1,
                                    initial_rates[sp].scat_1);
    if (zero) {
      for (int sp = 0; sp < 3; ++sp)
        result[sp] = baseline[sp];
      return true;
    }
  }
  spectral_vector baseline_guess, supplied_guess;
  bool baseline_guess_valid = true, supplied_guess_valid = options.use_initial_spectrum;
  for (int sp = 0; sp < 3; ++sp)
    baseline_guess_valid = baseline_guess_valid &&
        (!options.active_species[sp] ||
         (baseline_physical.n[sp] > 0 && baseline_physical.J[sp] > 0));
  if (baseline_guess_valid) {
    NuDistributionParams distribution = {};
    CalculateThickParamsFromM1(&baseline_physical, &distribution);
    for (int sp = 0; sp < 3; ++sp) {
      baseline_guess(2 * sp) = options.active_species[sp]
          ? log(distribution.temp_t[sp]) : CCTK_REAL(0);
      baseline_guess(2 * sp + 1) = options.active_species[sp]
          ? distribution.eta_t[sp] : CCTK_REAL(0);
    }
    baseline_guess_valid = nuX_Utils::roots::finite(baseline_guess);
  }
  for (int sp = 0; sp < 3; ++sp) {
    supplied_guess(2 * sp) = options.active_species[sp]
        ? options.initial_log_temperature[sp] : CCTK_REAL(0);
    supplied_guess(2 * sp + 1) = options.active_species[sp]
        ? options.initial_eta[sp] : CCTK_REAL(0);
  }
  supplied_guess_valid = supplied_guess_valid &&
                         nuX_Utils::roots::finite(supplied_guess);
  const spectral_source_t<Rates> source{geometry, volume, baseline, rates, options};
  if (supplied_guess_valid && source.solve(supplied_guess, result))
    return true;
  if (baseline_guess_valid && source.solve(baseline_guess, result))
    return true;
  // Mixed spectral INITIAL GUESSES help number-conserving scattering when an
  // LTE population is a poor guess, without freezing any species' coefficient.
  if (baseline_guess_valid && supplied_guess_valid) {
    spectral_vector mixed = baseline_guess;
    for (int sp = 0; sp < 3; ++sp)
      mixed(2 * sp) = supplied_guess(2 * sp);
    if (source.solve(mixed, result))
      return true;
    mixed = supplied_guess;
    for (int sp = 0; sp < 3; ++sp)
      mixed(2 * sp) = baseline_guess(2 * sp);
    if (source.solve(mixed, result))
      return true;
  }
  return false;
}

} // namespace nuX_M1

#endif // NUX_M1_NURATES_SOURCE_HXX
