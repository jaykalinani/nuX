// Focused tests for the reusable local collision solve.
#define NUX_M1_CLOSURE_IMPLEMENTATION
#define NUX_M1_SOURCES_IMPLEMENTATION
#include "nuX_M1_nurates_source.hxx"

#include <cmath>
#include <limits>

namespace nuX_M1 {
namespace {

struct geometry_t {
  tensor::metric<4> gdd{};
  tensor::inv_metric<4> guu{};
  tensor::generic<CCTK_REAL, 4, 1> nd{}, nu{}, ud{}, uu{}, vd{}, vu{};
  tensor::generic<CCTK_REAL, 4, 2> projector{}, gamma{};
  CCTK_REAL alpha = 1, W;

  explicit geometry_t(const CCTK_REAL vx) {
    gdd(0, 0) = guu(0, 0) = -1;
    nd(0) = -1;
    nu(0) = 1;
    for (int a = 1; a < 4; ++a) {
      gdd(a, a) = guu(a, a) = gamma(a, a) = 1;
    }
    vu(1) = vd(1) = vx;
    W = 1 / std::sqrt(1 - vx * vx);
    uu(0) = W;
    uu(1) = W * vx;
    tensor::contract(gdd, uu, &ud);
    calc_proj(ud, uu, &projector);
  }
};

struct rates_t {
  nurates_opacities_t values[3]{};
  bool reject_vacuum = false;
  CCTK_HOST CCTK_DEVICE bool operator()(
      const M1Quantities &physical, nurates_opacities_t *out) const {
    for (int sp = 0; sp < 3; ++sp) {
      if (reject_vacuum && (physical.n[sp] == 0 || physical.J[sp] == 0))
        return false;
      out[sp] = values[sp];
    }
    return true;
  }
};

SourceUpdateContext context(cGH const *cctkGH, const geometry_t &g,
                             const nurates_source_state_t &old) {
  return SourceUpdateContext(
      cctkGH, 0, 0, 0, 0, 1e-12, 128, false, 1, g.alpha, g.gdd, g.guu,
      g.nd, g.nu, g.gamma, g.ud, g.uu, g.vd, g.vu, g.projector, g.W,
      old.E, old.F_d, old.E, old.F_d, 0, 0, 0, true);
}

nurates_source_state_t state(const CCTK_REAL N, const CCTK_REAL E,
                              const CCTK_REAL Fx) {
  nurates_source_state_t out{};
  out.N = N;
  out.E = E;
  pack_F_d(0, 0, 0, Fx, 0, 0, &out.F_d);
  out.chi = 1.0 / 3;
  return out;
}

int check_small_flux(cGH const *cctkGH, const CCTK_REAL scattering) {
  const geometry_t g(0);
  nurates_source_state_t old[3], next[3];
  rates_t rates;
  for (int sp = 0; sp < 3; ++sp) {
    old[sp] = state(1, 1, 1e-10);
    rates.values[sp].scat_1 = scattering;
  }
  auto options = nurates_source_options_t{};
  options.closure = CLOSURE_EDDINGTON;
  options.epsrel = 1e-10;
  const auto ctx = context(cctkGH, g, old[0]);
  int failures = !coupled_nurates_source(ctx, 1, old, rates, options, next);
  if (failures)
    return failures;
  for (int sp = 0; sp < 3; ++sp) {
    const CCTK_REAL exact_flux = 1e-10 / (1 + scattering);
    failures += std::abs(next[sp].F_d(1) - exact_flux) >
                2e-9 * exact_flux;
    failures += std::abs(next[sp].E - 1) > 2e-10;
    failures += std::abs(next[sp].N - 1) > 2e-10;
  }
  // The final check must reject the historical falsely accepted Fnew=0.
  nurates_source_state_t wrong[3];
  for (int sp = 0; sp < 3; ++sp)
    wrong[sp] = state(1, 1, 0);
  const nurates_source_detail::spectral_source_t<rates_t> source{
      ctx, 1, old, rates, options};
  failures += source.verify(wrong);
  return failures;
}

int check_stiff_invariant(cGH const *cctkGH) {
  const geometry_t g(0.5);
  nurates_source_state_t old[3], next[3], wrong[3];
  rates_t rates;
  for (int sp = 0; sp < 3; ++sp) {
    old[sp] = state(1, 1, 0);
    rates.values[sp].scat_1 = 1e30;
    const CCTK_REAL J = 0.5; // Wrong LTE amplitude, half the invariant.
    wrong[sp] = state(1, (4 * g.W * g.W - 1) / 3 * J,
                      4 * g.W * g.W / 3 * J * 0.5);
  }
  auto options = nurates_source_options_t{};
  options.closure = CLOSURE_EDDINGTON;
  options.epsrel = 1e-9;
  const auto ctx = context(cctkGH, g, old[0]);
  const nurates_source_detail::spectral_source_t<rates_t> source{
      ctx, 1, old, rates, options};
  int failures = source.verify(wrong);
  failures += nurates_source_detail::verify_invariant(
      ctx, old[0], wrong[0], 0.5, 0, 0, options.epsrel);
  if (!coupled_nurates_source(ctx, 1, old, rates, options, next))
    return failures + 1;
  for (int sp = 0; sp < 3; ++sp) {
    const CCTK_REAL I = next[sp].E - 0.5 * next[sp].F_d(1);
    failures += std::abs(I - 1) > 2e-9;
    failures += std::abs(next[sp].E - (4 * g.W * g.W - 1) / 3) > 2e-9;
  }
  return failures;
}

int check_invariant_jacobian(cGH const *cctkGH) {
  const geometry_t g(0.5);
  const auto old = state(1, 1.3, 0.2);
  const CCTK_REAL absolute[4] = {1e-14, 1e-14, 1e-14, 1e-14};
  const SourceUpdateContext ctx(
      cctkGH, 0, 0, 0, 0, 1e-12, 128, false, 0.3, g.alpha, g.gdd, g.guu,
      g.nd, g.nu, g.gamma, g.ud, g.uu, g.vd, g.vu, g.projector, g.W,
      old.E, old.F_d, old.E, old.F_d, 0.4, 0.7, 1e12,
      true, absolute, true);
  Params parameters(ctx, CLOSURE_EDDINGTON, 1.0 / 3);
  const arith_vector q{1.1, 0.3, 0.01, -0.02};
  arith_vector residual;
  arith_matrix analytic;
  if (impl_func_val_jac(q, &parameters, residual, analytic) != ROOTS_SUCCESS)
    return 1;
  int failures = 0;
  for (int col = 0; col < 4; ++col) {
    auto plus = q, minus = q;
    plus(col) += 1e-6;
    minus(col) -= 1e-6;
    arith_vector fp, fm;
    failures += impl_func_val(plus, &parameters, fp) != ROOTS_SUCCESS;
    failures += impl_func_val(minus, &parameters, fm) != ROOTS_SUCCESS;
    for (int row = 0; row < 4; ++row) {
      const CCTK_REAL derivative = (fp(row) - fm(row)) / 2e-6;
      failures += std::abs(derivative - analytic(row, col)) >
                  2e-7 * std::fmax(CCTK_REAL(1), std::abs(derivative));
    }
  }
  return failures;
}

int check_inactive_species(cGH const *cctkGH) {
  const geometry_t g(0);
  nurates_source_state_t old[3], next[3];
  rates_t rates;
  for (int sp = 0; sp < 3; ++sp) {
    old[sp] = state(0, 0, 0);
    if (sp < 2) {
      rates.values[sp].eta_0 = rates.values[sp].eta_1 = 1;
      rates.values[sp].abs_0 = rates.values[sp].abs_1 = 1;
    }
  }
  auto options = nurates_source_options_t{};
  options.closure = CLOSURE_EDDINGTON;
  options.epsrel = 1e-9;
  options.active_species[2] = false;
  options.use_initial_spectrum = true;
  M1Quantities target{};
  for (int sp = 0; sp < 4; ++sp) {
    // A positive analytic solution for the two emitting electron species.
    target.n[sp] = 0.5 * nuX_ndens_conv;
    target.J[sp] = 0.5 * nuX_edens_conv;
    target.chi[sp] = 1.0 / 3;
  }
  NuDistributionParams distribution{};
  CalculateThickParamsFromM1(&target, &distribution);
  for (int sp = 0; sp < 3; ++sp) {
    options.initial_log_temperature[sp] = std::log(distribution.temp_t[sp]);
    options.initial_eta[sp] = distribution.eta_t[sp];
  }
  const auto ctx = context(cctkGH, g, old[0]);
  int failures = !coupled_nurates_source(ctx, 1, old, rates, options, next);
  if (failures)
    return failures;
  failures += !nurates_source_detail::vacuum(next[2]);
  for (int sp = 0; sp < 2; ++sp) {
    failures += std::abs(next[sp].N - 0.5) > 2e-9;
    failures += std::abs(next[sp].E - 0.5) > 2e-9;
  }
  // A false proof of inactivity must fail, not discard spontaneous emission.
  rates.values[2].eta_0 = rates.values[2].eta_1 = 1;
  failures += coupled_nurates_source(ctx, 1, old, rates, options, next);
  return failures;
}

int check_log_residual() {
  using nurates_source_detail::log_moment_residual;
  const CCTK_REAL eps = std::numeric_limits<CCTK_REAL>::epsilon();
  int failures = log_moment_residual(1, 1) != 0;
  const CCTK_REAL near = std::nextafter(CCTK_REAL(1), CCTK_REAL(2));
  failures += std::abs(log_moment_residual(near, 1) - std::log1p(near - 1)) >
              eps * eps;
  const CCTK_REAL distant = log_moment_residual(1e300, 1e-300);
  failures += !std::isfinite(distant) ||
              std::abs(distant - 600 * std::log(CCTK_REAL(10))) >
                  4 * eps * std::abs(distant);
  failures += log_moment_residual(1e-300, 1e300) != -distant;
  // A logarithmic perturbation retains unit derivative even where the old
  // bounded relative residual rounded to exactly one.
  constexpr CCTK_REAL h = 1e-5;
  const CCTK_REAL derivative =
      (log_moment_residual(1e100 * std::exp(h), 1) -
       log_moment_residual(1e100 * std::exp(-h), 1)) / (2 * h);
  failures += std::abs(derivative - 1) > 1e-8;
  failures += std::isfinite(log_moment_residual(0, 1));
  failures += std::isfinite(log_moment_residual(1, 0));
  failures += std::isfinite(log_moment_residual(-1, 1));
  failures += std::isfinite(log_moment_residual(
      std::numeric_limits<CCTK_REAL>::infinity(), 1));
  failures += std::isfinite(log_moment_residual(
      1, std::numeric_limits<CCTK_REAL>::quiet_NaN()));
  return failures;
}

int check_remote_initial_spectrum(cGH const *cctkGH) {
  const geometry_t g(0);
  auto options = nurates_source_options_t{};
  options.closure = CLOSURE_EDDINGTON;
  options.epsrel = 1e-9;
  nurates_source_detail::spectral_vector equilibrium;
  for (int sp = 0; sp < 3; ++sp) {
    equilibrium(2 * sp) = std::log(CCTK_REAL(3));
    equilibrium(2 * sp + 1) = -2;
  }
  M1Quantities target{};
  if (!nurates_source_detail::spectral_moments(equilibrium, options, target))
    return 1;
  nurates_source_state_t old[3], next[3];
  rates_t rates;
  for (int sp = 0; sp < 3; ++sp) {
    const CCTK_REAL flavors = sp == 2 ? 4 : 1;
    old[sp] = state(flavors * target.n[sp] / nuX_ndens_conv,
                    flavors * target.J[sp] / nuX_edens_conv, 0);
    rates.values[sp].abs_0 = rates.values[sp].abs_1 = 1;
    rates.values[sp].eta_0 = old[sp].N;
    rates.values[sp].eta_1 = old[sp].E;
  }
  const auto ctx = context(cctkGH, g, old[0]);
  const nurates_source_detail::spectral_source_t<rates_t> source{
      ctx, 1, old, rates, options};
  int failures = 0;
  for (const CCTK_REAL offset : {-20.0, 20.0}) {
    auto trial = equilibrium;
    for (int sp = 0; sp < 3; ++sp)
      trial(2 * sp) += offset;
    // Test the remote start itself, not a successful fallback to the baseline.
    if (!source.solve(trial, next)) {
      ++failures;
      continue;
    }
    for (int sp = 0; sp < 3; ++sp) {
      failures += std::abs(next[sp].N / old[sp].N - 1) > 2e-8;
      failures += std::abs(next[sp].E / old[sp].E - 1) > 2e-8;
    }
  }
  return failures;
}

} // namespace

int nurates_adapter_unit_tests() {
  M1OpacitiesNonThermalSeparated raw{};
  raw.kappa_0_a[2] = 2;
  raw.kappa_0_a[3] = 6;
  raw.eta_0[2] = 3;
  raw.eta_0[3] = 7;
  const auto grouped = convert_nurates_opacities(raw, 2, 4);
  int failures = std::abs(grouped.abs_0 / nuX_length_conv - 4) > 1e-14;
  failures += std::abs(grouped.eta_0 * nuX_ndens_conv / nuX_time_conv - 20) >
              1e-13;

  M1Quantities moments{};
  constexpr CCTK_REAL volume = 2, temperature = 3;
  constexpr CCTK_REAL Nfloor = 1e-20, Efloor = 1e-20;
  for (int sp = 0; sp < 3; ++sp) {
    const CCTK_REAL flavors = sp == 2 ? 4 : 1;
    moments.n[sp] = Nfloor / (volume * flavors) * nuX_ndens_conv;
    moments.J[sp] = Efloor / (volume * flavors) * nuX_edens_conv;
  }
  const auto original = moments;
  bool floor_only[3];
  nurates_rate_reference(moments, volume, temperature, Nfloor, Efloor,
                         floor_only);
  for (int sp = 0; sp < 3; ++sp) {
    failures += !floor_only[sp];
    failures += moments.n[sp] > original.n[sp];
    failures += moments.J[sp] > original.J[sp];
    failures += std::abs(moments.J[sp] / moments.n[sp] - 3 * temperature) > 1e-13;
  }
  failures += moments.n[3] != moments.n[2] || moments.J[3] != moments.J[2];
  return failures;
}

int nurates_source_unit_tests(cGH const *cctkGH) {
  return check_small_flux(cctkGH, 1) + check_small_flux(cctkGH, 1e8) +
         check_stiff_invariant(cctkGH) +
         check_invariant_jacobian(cctkGH) + check_inactive_species(cctkGH) +
         check_log_residual() + check_remote_initial_spectrum(cctkGH) +
         nurates_adapter_unit_tests();
}

} // namespace nuX_M1
