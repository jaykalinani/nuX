#include "m1_unit_tests.hxx"
#define NUX_M1_CLOSURE_IMPLEMENTATION
#define NUX_M1_SOURCES_IMPLEMENTATION
#include "nuX_M1_sources.hxx"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>

namespace nuX_M1 {
namespace {

struct source_test_geometry_t {
  tensor::metric<4> g_dd{};
  tensor::inv_metric<4> g_uu{};
  tensor::generic<CCTK_REAL, 4, 1> n_d{};
  tensor::generic<CCTK_REAL, 4, 1> n_u{};
  tensor::generic<CCTK_REAL, 4, 2> gamma_ud{};
  tensor::generic<CCTK_REAL, 4, 1> u_d{};
  tensor::generic<CCTK_REAL, 4, 1> u_u{};
  tensor::generic<CCTK_REAL, 4, 1> v_d{};
  tensor::generic<CCTK_REAL, 4, 1> v_u{};
  tensor::generic<CCTK_REAL, 4, 2> proj_ud{};
  CCTK_REAL alpha{};
  CCTK_REAL W{};

  source_test_geometry_t(const CCTK_REAL lapse,
                         const std::array<CCTK_REAL, 3> &velocity)
      : alpha(lapse) {
    for (int a = 0; a < 4; ++a) {
      n_d(a) = n_u(a) = u_d(a) = u_u(a) = v_d(a) = v_u(a) = 0.0;
      for (int b = 0; b < 4; ++b)
        gamma_ud(a, b) = proj_ud(a, b) = 0.0;
    }
    g_dd(0, 0) = -alpha * alpha;
    g_uu(0, 0) = -1.0 / (alpha * alpha);
    for (int a = 1; a < 4; ++a) {
      g_dd(a, a) = 1.0;
      g_uu(a, a) = 1.0;
      gamma_ud(a, a) = 1.0;
      v_u(a) = velocity[a - 1];
      v_d(a) = velocity[a - 1];
    }
    n_d(0) = -alpha;
    n_u(0) = 1.0 / alpha;
    const CCTK_REAL v2 = velocity[0] * velocity[0] + velocity[1] * velocity[1] +
                         velocity[2] * velocity[2];
    W = 1.0 / std::sqrt(1.0 - v2);
    u_u(0) = W / alpha;
    for (int a = 1; a < 4; ++a)
      u_u(a) = W * velocity[a - 1];
    tensor::contract(g_dd, u_u, &u_d);
    calc_proj(u_d, u_u, &proj_ud);
  }
};

int fixed_chi_residual(const SourceUpdateContext &ctx, const arith_vector &q,
                       const CCTK_REAL chi, arith_vector *residual) {
  const CCTK_REAL E = q(0);
  if (!std::isfinite(E) || E < 0.0 || !nuX_Utils::roots::finite(q))
    return ROOTS_EBADFUNC;

  tensor::generic<CCTK_REAL, 4, 1> F_d;
  pack_F_d(-ctx.alp * ctx.n_u(1), -ctx.alp * ctx.n_u(2), -ctx.alp * ctx.n_u(3),
           q(1), q(2), q(3), &F_d);
  tensor::symmetric2<CCTK_REAL, 4, 2> P_dd;
  apply_closure(ctx.g_dd, ctx.g_uu, ctx.n_d, ctx.W, ctx.u_u, ctx.v_d,
                ctx.proj_ud, E, F_d, chi, &P_dd);
  tensor::symmetric2<CCTK_REAL, 4, 2> T_dd;
  assemble_rT(ctx.n_d, E, F_d, P_dd, &T_dd);
  const CCTK_REAL J = calc_J_from_rT(T_dd, ctx.u_u);
  tensor::generic<CCTK_REAL, 4, 1> H_d;
  calc_H_from_rT(T_dd, ctx.u_u, ctx.proj_ud, &H_d);
  tensor::generic<CCTK_REAL, 4, 1> S_d, tS_d;
  calc_rad_sources(ctx.eta, ctx.kabs, ctx.kscat, ctx.u_d, J, H_d, &S_d);
  const CCTK_REAL Edot = calc_rE_source(ctx.alp, ctx.n_u, S_d);
  calc_rF_source(ctx.alp, ctx.gamma_ud, S_d, &tS_d);
  (*residual)(0) = q(0) - ctx.Estar - ctx.cdt * Edot;
  for (int a = 1; a < 4; ++a)
    (*residual)(a) = q(a) - ctx.Fstar_d(a) - ctx.cdt * tS_d(a);
  return nuX_Utils::roots::finite(*residual) ? ROOTS_SUCCESS : ROOTS_EBADFUNC;
}

int check_source_jacobian(cGH const *cctkGH, const CCTK_REAL chi,
                          const arith_vector &q) {
  const source_test_geometry_t geom(0.7, {0.13, -0.21, 0.08});
  tensor::generic<CCTK_REAL, 4, 1> Fstar_d;
  pack_F_d(0.0, 0.0, 0.0, q(1), q(2), q(3), &Fstar_d);
  const SourceUpdateContext ctx(cctkGH, 0, 0, 0, 0, 1.0e-12, 64, true, 0.17,
                                geom.alpha, geom.g_dd, geom.g_uu, geom.n_d,
                                geom.n_u, geom.gamma_ud, geom.u_d, geom.u_u,
                                geom.v_d, geom.v_u, geom.proj_ud, geom.W, q(0),
                                Fstar_d, q(0), Fstar_d, 0.31, 1.7, 2.4);

  tensor::generic<CCTK_REAL, 4, 1> F_d, F_u;
  pack_F_d(0.0, 0.0, 0.0, q(1), q(2), q(3), &F_d);
  tensor::contract(geom.g_uu, F_d, &F_u);
  double q_data[4] = {q(0), q(1), q(2), q(3)};
  double Fup[4] = {F_u(0), F_u(1), F_u(2), F_u(3)};
  double vup[4] = {geom.v_u(0), geom.v_u(1), geom.v_u(2), geom.v_u(3)};
  double vdown[4] = {geom.v_d(0), geom.v_d(1), geom.v_d(2), geom.v_d(3)};
  arith_matrix analytic;
  source_jacobian_fixed_chi(
      q_data, Fup, tensor::dot(F_u, F_d), chi, ctx.kabs, ctx.kscat, vup, vdown,
      tensor::dot(geom.v_u, geom.v_d), geom.W, geom.alpha, ctx.cdt, analytic);

  int failures = 0;
  for (int column = 0; column < 4; ++column) {
    const CCTK_REAL step = 1.0e-6 * std::max(CCTK_REAL(1), std::abs(q(column)));
    arith_vector plus = q;
    arith_vector minus = q;
    plus(column) += step;
    minus(column) -= step;
    arith_vector fplus, fminus;
    if (fixed_chi_residual(ctx, plus, chi, &fplus) != ROOTS_SUCCESS ||
        fixed_chi_residual(ctx, minus, chi, &fminus) != ROOTS_SUCCESS) {
      ++failures;
      continue;
    }
    for (int row = 0; row < 4; ++row) {
      const CCTK_REAL numerical = (fplus(row) - fminus(row)) / (2.0 * step);
      const CCTK_REAL scale = std::max(
          {CCTK_REAL(1), std::abs(numerical), std::abs(analytic(row, column))});
      failures += std::abs(numerical - analytic(row, column)) > 5.0e-5 * scale;
    }
  }
  return failures;
}

int check_moment_repair() {
  source_test_geometry_t geom(1.0, {0.0, 0.0, 0.0});
  int failures = 0;

  CCTK_REAL E = 1.0;
  tensor::generic<CCTK_REAL, 4, 1> F_d;
  pack_F_d(0.0, 0.0, 0.0, 2.0, 0.0, 0.0, &F_d);
  moment_repair_t repair;
  repair_moments(geom.g_uu, &E, &F_d, 1.0e-12, 1.0e-4, &repair);
  const CCTK_REAL F2 = tensor::dot(geom.g_uu, F_d, F_d);
  failures += !repair.repaired_flux || repair.repaired_energy;
  failures += std::abs(F2 - (1.0 - 1.0e-4) * E * E) > 5.0e-14;
  failures += !source_state_is_realizable(geom.g_uu, E, F_d);

  const CCTK_REAL E_valid = E;
  const auto F_valid = F_d;
  repair_moments(geom.g_uu, &E, &F_d, 1.0e-12, 1.0e-4, &repair);
  failures += repair.repaired_energy || repair.repaired_flux || E != E_valid;
  for (int a = 0; a < 4; ++a)
    failures += F_d(a) != F_valid(a);

  E = -1.0;
  pack_F_d(0.0, 0.0, 0.0, 1.0, 0.0, 0.0, &F_d);
  repair_moments(geom.g_uu, &E, &F_d, 1.0e-6, 1.0e-4, &repair);
  failures += !repair.repaired_energy || !repair.repaired_flux ||
              E != CCTK_REAL(1.0e-6) ||
              tensor::dot(geom.g_uu, F_d, F_d) > (1.0 - 1.0e-4) * E * E;

  pack_F_d(0.0, 0.0, 0.0, 2.0, 0.0, 0.0, &F_d);
  failures += source_state_is_realizable(geom.g_uu, 1.0, F_d);
  failures += source_state_is_realizable(geom.g_uu, -1.0, F_d);
  return failures;
}

int check_flux_factor() {
  const source_test_geometry_t geom(1.0, {0.0, 0.0, 0.0});
  tensor::generic<CCTK_REAL, 4, 1> H_d;
  for (int a = 0; a < 4; ++a)
    H_d(a) = 0.0;
  H_d(1) = 0.6;
  H_d(2) = 0.8;

  int failures = 0;
  failures +=
      std::abs(flux_factor(geom.g_uu, 2.0, H_d, 1.0e-12) - 0.5) > 2.0e-15;
  failures += flux_factor(geom.g_uu, 0.0, H_d, 1.0e-12) != 0.0;
  H_d(1) = 4.0;
  failures += flux_factor(geom.g_uu, 2.0, H_d, 1.0e-12) != 1.0;
  return failures;
}

int check_eddington_comoving_domain() {
  // A fixed-Eddington tensor is not a realizability-preserving closure under a
  // boost.  This lab-realizable streaming state has positive J but |H| > J in
  // the fluid frame.  It must remain admissible only for the explicitly chosen
  // Eddington approximation; variable M1 closures must reject it.
  const source_test_geometry_t geom(1.0, {0.3, 0.0, 0.0});
  constexpr CCTK_REAL E = 1.0;
  tensor::generic<CCTK_REAL, 4, 1> F_d;
  pack_F_d(0.0, 0.0, 0.0, 1.0, 0.0, 0.0, &F_d);
  tensor::symmetric2<CCTK_REAL, 4, 2> P_dd;
  apply_closure(geom.g_dd, geom.g_uu, geom.n_d, geom.W, geom.u_u,
                geom.v_d, geom.proj_ud, E, F_d, CCTK_REAL(1.0 / 3.0),
                &P_dd);
  tensor::symmetric2<CCTK_REAL, 4, 2> T_dd;
  assemble_rT(geom.n_d, E, F_d, P_dd, &T_dd);
  const CCTK_REAL J = calc_J_from_rT(T_dd, geom.u_u);
  tensor::generic<CCTK_REAL, 4, 1> H_d;
  calc_H_from_rT(T_dd, geom.u_u, geom.proj_ud, &H_d);

  int failures = 0;
  failures += !(J > CCTK_REAL(0));
  failures += comoving_state_is_realizable(geom.g_uu, J, H_d);
  failures +=
      !comoving_state_is_acceptable(CLOSURE_EDDINGTON, geom.g_uu, J, H_d);
  failures +=
      comoving_state_is_acceptable(CLOSURE_MINERBO, geom.g_uu, J, H_d);

  // Eddington does not excuse a negative or nonfinite comoving energy.
  failures += comoving_state_is_acceptable(CLOSURE_EDDINGTON, geom.g_uu,
                                            CCTK_REAL(-1), H_d);
  return failures;
}

int check_number_source() {
  constexpr CCTK_REAL alpha = 0.7;
  constexpr CCTK_REAL volform = 1.3;
  constexpr CCTK_REAL eta = 0.4;
  constexpr CCTK_REAL kabs = 2.1;
  constexpr CCTK_REAL N = 0.6;
  constexpr CCTK_REAL Gamma = 1.2;
  const CCTK_REAL expected = alpha * (volform * eta - kabs * N / Gamma);
  const CCTK_REAL source = calc_rN_source(alpha, volform, eta, kabs, N, Gamma);
  return std::abs(source - expected) >
         4.0 * std::numeric_limits<CCTK_REAL>::epsilon();
}

int check_imex_number_update() {
  int failures = 0;
  constexpr CCTK_REAL Nstar = 2.0;
  constexpr CCTK_REAL alpha = 0.63;
  constexpr CCTK_REAL dt = 0.4;
  constexpr CCTK_REAL volform = 1.3;
  constexpr CCTK_REAL eta = 0.3;
  constexpr CCTK_REAL Gamma = 1.2;
  for (const CCTK_REAL kabs : {CCTK_REAL(0), CCTK_REAL(0.1),
                               CCTK_REAL(10), CCTK_REAL(1.0e8)}) {
    const auto update =
        source_number_update(Nstar, dt, alpha, volform, eta, kabs, Gamma);
    const CCTK_REAL residual =
        update.state - Nstar -
        dt * calc_rN_source(alpha, volform, eta, kabs, update.state, Gamma);
    failures += !std::isfinite(update.state) || update.state < 0.0;
    failures += std::abs(residual) >
                32.0 * std::numeric_limits<CCTK_REAL>::epsilon() * Nstar;
    failures += std::abs(update.state - Nstar - update.increment) >
                8.0 * std::numeric_limits<CCTK_REAL>::epsilon() * Nstar;
  }

  // Strong absorption leaves a representable positive root that would be
  // erased by reconstructing it from its rounded difference from Nstar.
  const auto tiny =
      source_number_update(Nstar, 1.0, alpha, 1.0, 0.0, 1.0e30, 1.0);
  failures += !(tiny.state > 0.0 && tiny.state < 1.0e-29);
  failures += Nstar + tiny.increment != 0.0;
  const auto transparent =
      source_number_update(Nstar, dt, alpha, volform, 0.0, 0.0, Gamma);
  failures += transparent.state != Nstar || transparent.increment != 0.0;

  // Near equilibrium, form the exchange from the source, not state-Nstar.
  const CCTK_REAL near_eta =
      1.0 + 64.0 * std::numeric_limits<CCTK_REAL>::epsilon();
  const auto near = source_number_update(1.0, dt, alpha, 1.0, near_eta, 1.0,
                                          1.0);
  const long double h = static_cast<long double>(dt) * alpha;
  const long double reference =
      h * (static_cast<long double>(near_eta) - 1.0L) / (1.0L + h);
  failures += !(near.increment > 0.0);
  failures += std::abs(static_cast<long double>(near.increment) - reference) >
              8.0L * std::numeric_limits<CCTK_REAL>::epsilon() *
                  std::abs(reference);
  return failures;
}

int check_physical_number_current() {
  int failures = 0;
  for (const auto &velocity : {std::array<CCTK_REAL, 3>{0.3, -0.2, 0.1},
                               std::array<CCTK_REAL, 3>{0.9, 0.0, 0.0}}) {
    const source_test_geometry_t geom(0.63, velocity);
    constexpr CCTK_REAL J = 1.0e-30;
    tensor::generic<CCTK_REAL, 4, 1> H_d;
    tensor::symmetric2<CCTK_REAL, 4, 2> K_dd, T_dd;
    for (int a = 0; a < 4; ++a) {
      H_d(a) = 0.0;
      for (int b = a; b < 4; ++b)
        K_dd(a, b) = J * (geom.g_dd(a, b) + geom.u_d(a) * geom.u_d(b)) / 3.0;
    }
    assemble_rT(geom.u_d, J, H_d, K_dd, &T_dd);
    const CCTK_REAL E = calc_J_from_rT(T_dd, geom.n_u);
    tensor::generic<CCTK_REAL, 4, 1> F_d;
    calc_H_from_rT(T_dd, geom.n_u, geom.gamma_ud, &F_d);
    CCTK_REAL Gamma = std::numeric_limits<CCTK_REAL>::quiet_NaN();
    failures += !compute_source_Gamma(geom.W, geom.v_u, J, E, F_d, &Gamma);
    failures += std::abs(Gamma - geom.W) > 2.0e-13 * geom.W;
    // This physical state is far below the transport floor, where the legacy
    // helper returns 1 even though an isotropic moving current has Gamma=W.
    failures += compute_Gamma(geom.W, geom.v_u, J, E, F_d, 1.0e-15,
                               1.0e-5) != 1.0;
    constexpr CCTK_REAL Nstar = 2.0e-30;
    constexpr CCTK_REAL dt = 0.4;
    constexpr CCTK_REAL eta = 0.3e-30;
    constexpr CCTK_REAL kabs = 4.0;
    const auto update =
        source_number_update(Nstar, dt, geom.alpha, 1.0, eta, kabs, Gamma);
    const CCTK_REAL expected =
        (Nstar + dt * geom.alpha * eta) /
        (1.0 + dt * geom.alpha * kabs / geom.W);
    failures += std::abs(update.state - expected) > 2.0e-13 * expected;
    const auto floor_current =
        source_number_update(Nstar, dt, geom.alpha, 1.0, eta, kabs, 1.0);
    failures += std::abs(update.state - floor_current.state) < 1.0e-3 * expected;
    if (velocity[0] > 0.8) {
      // A large optional transport margin must likewise not clamp F.v/E in
      // the accepted physical collision current.
      const CCTK_REAL margin_Gamma =
          compute_Gamma(geom.W, geom.v_u, J, E, F_d, 0.0, 0.5);
      failures += std::abs(margin_Gamma - Gamma) < 0.01 * Gamma;
    }
    pack_F_d(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, &F_d);
    failures += !compute_source_Gamma(geom.W, geom.v_u, 0.0, 0.0, F_d,
                                       &Gamma);
    failures += Gamma != geom.W;
    failures += compute_source_Gamma(geom.W, geom.v_u, 0.0, E, F_d, &Gamma);
  }
  return failures;
}

int check_source_convergence_contract() {
  nuX_Utils::roots::hybridsj_solver<CCTK_REAL, 4> solver{};
  const arith_vector baseline{1.0, 1.0e-8, -2.0e-8, 0.0};
  solver.x = baseline;
  for (int n = 0; n < 4; ++n)
    solver.f(n) = solver.dx(n) = 0.0;
  CCTK_REAL residual, step;
  const auto converged = [&]() {
    return source_solver_converged(solver, baseline, 1.0e-14, 1.0e-8,
                                    &residual, &step);
  };
  int failures = !converged();
  // A small step alone must never accept a failed residual. In particular,
  // the energy tolerance must not hide errors in a much smaller flux field.
  solver.f(1) = 2.0e-14;
  failures += converged();
  solver.f(1) = 0.0;
  solver.dx(1) = 2.0e-14;
  failures += converged();
  solver.dx(1) = 0.0;
  solver.f(2) = std::numeric_limits<CCTK_REAL>::quiet_NaN();
  failures += converged();
  return failures;
}

int check_closure_trial_failure(cGH const *cctkGH) {
  const source_test_geometry_t geom(0.7, {0.13, -0.21, 0.08});
  tensor::generic<CCTK_REAL, 4, 1> Fstar_d;
  pack_F_d(0.0, 0.0, 0.0, 0.4, 0.2, 0.1, &Fstar_d);
  SourceUpdateContext ctx(
      cctkGH, 0, 0, 0, 0, 1.0e-14, 1, false, 0.1, geom.alpha, geom.g_dd,
      geom.g_uu, geom.n_d, geom.n_u, geom.gamma_ud, geom.u_d, geom.u_u,
      geom.v_d, geom.v_u, geom.proj_ud, geom.W, 2.0, Fstar_d, 2.0, Fstar_d,
      0.3, 1.4, 2.1, true);
  const arith_vector q{2.0, 0.4, 0.2, 0.1};
  Params p(ctx, CLOSURE_MINERBO, 0.5);
  PreparedState state;
  // Deliberately insufficient closure iterations must reject the callback
  // without a device abort or an endpoint-closure substitution.
  int failures = prepare_closure(q, &p, &state) != ROOTS_EBADFUNC;
  ctx.closure_maxiter = 128;
  failures += prepare_closure(q, &p, &state) != ROOTS_SUCCESS;
  failures += !(p.chi > 1.0 / 3.0 && p.chi < 1.0);
  return failures;
}

int check_imex_diagonal_limits(cGH const *cctkGH) {
  const source_test_geometry_t geom(0.63, {0.0, 0.0, 0.0});
  constexpr CCTK_REAL Estar = 2.0;
  constexpr CCTK_REAL eta = 0.3;
  constexpr CCTK_REAL kabs = 1.7;
  constexpr CCTK_REAL kscat = 2.4;
  tensor::generic<CCTK_REAL, 4, 1> Fstar_d;
  pack_F_d(0.0, 0.0, 0.0, 0.3, -0.1, 0.05, &Fstar_d);
  int failures = 0;
  for (const closure_t closure : {CLOSURE_EDDINGTON, CLOSURE_MINERBO}) {
    for (const CCTK_REAL dt : {CCTK_REAL(0.001), CCTK_REAL(0.1),
                               CCTK_REAL(1), CCTK_REAL(50)}) {
      const SourceUpdateContext ctx(
          cctkGH, 0, 0, 0, 0, 1.0e-12, 128, false, dt, geom.alpha, geom.g_dd,
          geom.g_uu, geom.n_d, geom.n_u, geom.gamma_ud, geom.u_d, geom.u_u,
          geom.v_d, geom.v_u, geom.proj_ud, geom.W, Estar, Fstar_d, Estar,
          Fstar_d, eta, kabs, kscat, true);
      CCTK_REAL Enew = Estar * 0.9;
      auto Fnew_d = Fstar_d;
      CCTK_REAL chi = 1.0 / 3.0;
      // Enable every legacy threshold, and request no optional forcing.
      // The production IMEX policy must still select the strict residual.
      const int status = source_update(
          ctx, closure, &chi, &Enew, &Fnew_d,
          source_requires_nonlinear(false, false), 0.01, 0.01, 128,
          1.0e-14, 1.0e-10);
      failures += status != NUX_M1_SOURCE_OK;
      const CCTK_REAL h = dt * geom.alpha;
      const CCTK_REAL expected_E = (Estar + h * eta) / (1.0 + h * kabs);
      failures += std::abs(Enew - expected_E) > 3.0e-10 * Estar;
      for (int a = 1; a < 4; ++a) {
        const CCTK_REAL expected_F = Fstar_d(a) / (1.0 + h * (kabs + kscat));
        failures += std::abs(Fnew_d(a) - expected_F) > 3.0e-11;
      }
      Params p(ctx, closure, chi);
      arith_vector residual;
      const arith_vector q{Enew, Fnew_d(1), Fnew_d(2), Fnew_d(3)};
      failures += impl_func_val(q, &p, residual) != ROOTS_SUCCESS;
      for (int n = 0; n < 4; ++n) {
        const CCTK_REAL base = n == 0 ? Estar : Fstar_d(n);
        const CCTK_REAL scale =
            1.0e-14 + 1.0e-10 * std::max(std::abs(q(n)), std::abs(base));
        failures += !std::isfinite(residual(n)) || std::abs(residual(n)) > scale;
      }
    }
  }
  return failures;
}

int check_tiny_source_root(cGH const *cctkGH) {
  const source_test_geometry_t geom(0.63, {0.0, 0.0, 0.0});
  constexpr CCTK_REAL Estar = 2.0;
  constexpr CCTK_REAL kabs = 1.0e30;
  tensor::generic<CCTK_REAL, 4, 1> Fstar_d;
  pack_F_d(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, &Fstar_d);
  const SourceUpdateContext ctx(
      cctkGH, 0, 0, 0, 0, 1.0e-12, 128, false, 1.0, geom.alpha, geom.g_dd,
      geom.g_uu, geom.n_d, geom.n_u, geom.gamma_ud, geom.u_d, geom.u_u,
      geom.v_d, geom.v_u, geom.proj_ud, geom.W, Estar, Fstar_d, Estar, Fstar_d,
      0.0, kabs, 0.0, true);
  const CCTK_REAL expected = Estar / (1.0 + geom.alpha * kabs);
  CCTK_REAL Enew = expected;
  auto Fnew_d = Fstar_d;
  CCTK_REAL chi = 1.0 / 3.0;
  const int status = source_update(
      ctx, CLOSURE_EDDINGTON, &chi, &Enew, &Fnew_d,
      source_requires_nonlinear(false, false), 0.01, 0.01, 128,
      1.0e-14, 1.0e-10);
  int failures = status != NUX_M1_SOURCE_OK;
  failures += !(Enew > 0.0 && Enew < 1.0e-29);
  failures += std::abs(Enew - expected) > 1.0e-12 * expected;
  // The production IMEX path stores Enew itself to avoid this roundtrip.
  failures += Estar + (Enew - Estar) != 0.0;
  return failures;
}

int check_source_solve(cGH const *cctkGH) {
  const source_test_geometry_t geom(0.63, {0.11, -0.07, 0.04});
  constexpr CCTK_REAL Estar = 1.7;
  tensor::generic<CCTK_REAL, 4, 1> Fstar_d;
  pack_F_d(0.0, 0.0, 0.0, 0.24, -0.16, 0.09, &Fstar_d);
  const SourceUpdateContext ctx(cctkGH, 0, 0, 0, 0, 1.0e-12, 64, true, 0.13,
                                geom.alpha, geom.g_dd, geom.g_uu, geom.n_d,
                                geom.n_u, geom.gamma_ud, geom.u_d, geom.u_u,
                                geom.v_d, geom.v_u, geom.proj_ud, geom.W, Estar,
                                Fstar_d, Estar, Fstar_d, 0.29, 1.4, 2.1);

  // The thick-limit state is only an initial guess.  Every nonzero diagonal
  // stage in the paper-faithful mode must still satisfy the nonlinear
  // residual solve.
  CCTK_REAL Enew = 1.2;
  tensor::generic<CCTK_REAL, 4, 1> Fnew_d;
  pack_F_d(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, &Fnew_d);
  CCTK_REAL chi = 1.0 / 3.0;
  const int status =
      source_update(ctx, CLOSURE_EDDINGTON, &chi, &Enew, &Fnew_d, true, 20.0,
                    -1.0, 100, 1.0e-13, 1.0e-10);

  int failures = 0;
  failures += status == NUX_M1_SOURCE_FAIL;
  failures += !std::isfinite(Enew) || Enew < 0.0;
  for (int a = 0; a < 4; ++a)
    failures += !std::isfinite(Fnew_d(a));

  // Validate the accepted state independently from the solver return code.
  const arith_vector solution{Enew, Fnew_d(1), Fnew_d(2), Fnew_d(3)};
  arith_vector residual;
  failures +=
      fixed_chi_residual(ctx, solution, chi, &residual) != ROOTS_SUCCESS;
  for (int n = 0; n < 4; ++n) {
    const CCTK_REAL baseline = n == 0 ? ctx.Estar : ctx.Fstar_d(n);
    const CCTK_REAL scale =
        1.0e-13 +
        1.0e-10 * std::max(std::abs(solution(n)), std::abs(baseline));
    failures += !std::isfinite(residual(n)) || std::abs(residual(n)) > scale;
  }

  // A deliberately invalid analytic initial guess must not be accepted as a
  // solution.  Strict mode retries from the valid provisional state, so
  // require that recovery to converge to the actual implicit root.
  Enew = -1.0;
  pack_F_d(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, &Fnew_d);
  const int invalid_status =
      source_update(ctx, CLOSURE_EDDINGTON, &chi, &Enew, &Fnew_d, true, -1.0,
                    -1.0, 2, 1.0e-15, 1.0e-8);
  failures += invalid_status == NUX_M1_SOURCE_FAIL;
  failures += !std::isfinite(Enew) || Enew < 0.0;
  const arith_vector recovered{Enew, Fnew_d(1), Fnew_d(2), Fnew_d(3)};
  failures +=
      fixed_chi_residual(ctx, recovered, chi, &residual) != ROOTS_SUCCESS;
  for (int n = 0; n < 4; ++n) {
    const CCTK_REAL baseline = n == 0 ? ctx.Estar : ctx.Fstar_d(n);
    const CCTK_REAL scale =
        1.0e-15 +
        1.0e-8 * std::max(std::abs(recovered(n)), std::abs(baseline));
    failures += !std::isfinite(residual(n)) || std::abs(residual(n)) > scale;
  }
  return failures;
}

int check_repaired_source_baseline(cGH const *cctkGH) {
  const source_test_geometry_t geom(0.8, {0.0, 0.0, 0.0});

  // Mimic a provisional tableau state with negative energy and excessive
  // flux. The repaired state—not the invalid original—must define both the
  // closure and the implicit residual baseline.
  CCTK_REAL Estar = -2.0;
  tensor::generic<CCTK_REAL, 4, 1> Fstar_d;
  pack_F_d(0.0, 0.0, 0.0, 3.0, -1.0, 0.5, &Fstar_d);
  moment_repair_t repair;
  repair_moments(geom.g_uu, &Estar, &Fstar_d, 1.0e-6, 1.0e-4, &repair);

  const SourceUpdateContext ctx(cctkGH, 0, 0, 0, 0, 1.0e-12, 64, true, 0.3,
                                geom.alpha, geom.g_dd, geom.g_uu, geom.n_d,
                                geom.n_u, geom.gamma_ud, geom.u_d, geom.u_u,
                                geom.v_d, geom.v_u, geom.proj_ud, geom.W, Estar,
                                Fstar_d, Estar, Fstar_d, 0.0, 0.0, 0.0);

  CCTK_REAL Enew = Estar;
  tensor::generic<CCTK_REAL, 4, 1> Fnew_d = Fstar_d;
  CCTK_REAL chi = 1.0 / 3.0;
  const int status =
      source_update(ctx, CLOSURE_EDDINGTON, &chi, &Enew, &Fnew_d, true, -1.0,
                    -1.0, 16, 1.0e-15, 1.0e-12);

  int failures = 0;
  failures += !repair.repaired_energy || !repair.repaired_flux;
  failures += status == NUX_M1_SOURCE_FAIL;
  failures += Enew != Estar;
  for (int a = 0; a < 4; ++a)
    failures += Fnew_d(a) != Fstar_d(a);
  failures += !source_state_is_realizable(geom.g_uu, Enew, Fnew_d);
  return failures;
}

int check_vertex_velocity_repair() {
  const source_test_geometry_t geom(0.8, {0.0, 0.0, 0.0});
  int failures = 0;

  CCTK_REAL vx = 0.2;
  CCTK_REAL vy = -0.1;
  CCTK_REAL vz = 0.05;
  const CCTK_REAL vx_old = vx;
  const CCTK_REAL vy_old = vy;
  const CCTK_REAL vz_old = vz;
  CCTK_REAL W = repair_velocity_and_compute_W(geom.g_dd, &vx, &vy, &vz);
  failures += vx != vx_old || vy != vy_old || vz != vz_old;
  const CCTK_REAL v2 = vx * vx + vy * vy + vz * vz;
  failures += std::abs(W * W * (1.0 - v2) - 1.0) > 2.0e-15;

  vx = 2.0;
  vy = 0.0;
  vz = 0.0;
  W = repair_velocity_and_compute_W(geom.g_dd, &vx, &vy, &vz);
  const CCTK_REAL repaired_v2 = vx * vx + vy * vy + vz * vz;
  failures += !std::isfinite(W) || !(repaired_v2 < 1.0);
  failures += std::abs(W * W * (1.0 - repaired_v2) - 1.0) > 1.0e-4;

  vx = std::numeric_limits<CCTK_REAL>::quiet_NaN();
  vy = 1.0;
  vz = 1.0;
  W = repair_velocity_and_compute_W(geom.g_dd, &vx, &vy, &vz);
  failures += W != 1.0 || vx != 0.0 || vy != 0.0 || vz != 0.0;
  return failures;
}

int check_roundoff_closure(cGH const *cctkGH) {
  source_test_geometry_t geom(0.71, {0.0, 0.0, 0.0});
  for (int a = 1; a < 4; ++a) {
    geom.g_dd(a, a) = 1.72;
    geom.g_uu(a, a) = 1 / 1.72;
  }
  int failures = 0;
  for (const CCTK_REAL E : {CCTK_REAL(1), CCTK_REAL(1e-15)}) {
    for (const auto closure : {CLOSURE_MINERBO, CLOSURE_KERSHAW}) {
      tensor::generic<CCTK_REAL, 4, 1> F{};
      tensor::symmetric2<CCTK_REAL, 4, 2> P, T;
      CCTK_REAL chi = CCTK_REAL(1) / 3;
      bool valid = false;
      calc_closure(cctkGH, 0, 0, 0, 0, closure, geom.g_dd, geom.g_uu,
                   geom.n_d, geom.W, geom.u_u, geom.v_d, geom.proj_ud,
                   E, F, &chi, &P, 1e-12, 128, false, &valid);
      failures += !valid || chi != CCTK_REAL(1) / 3;

      // Force a one-ulp projector error so this regression also exercises
      // the cancellation on compilers that simplify the nominal geometry.
      auto rounded_projector = geom.proj_ud;
      rounded_projector(0, 0) = std::numeric_limits<CCTK_REAL>::epsilon();
      Parameters parameters(closure, geom.g_dd, geom.g_uu, geom.n_d, geom.W,
                             geom.u_u, geom.v_d, rounded_projector, E, F);
      apply_closure(geom.g_dd, geom.g_uu, geom.n_d, geom.W, geom.u_u,
                    geom.v_d, rounded_projector, E, F, CCTK_REAL(1) / 3, &P);
      assemble_rT(geom.n_d, E, F, P, &T);
      const CCTK_REAL J = calc_J_from_rT(T, geom.u_u);
      tensor::generic<CCTK_REAL, 4, 1> H;
      calc_H_from_rT(T, geom.u_u, rounded_projector, &H);
      failures += !(tensor::dot(geom.g_uu, H, H) < 0);
      failures += !comoving_state_is_realizable(geom.g_uu, J, H);
      failures += zFunction(0, &parameters) != 0;
      valid = false;
      calc_closure(cctkGH, 0, 0, 0, 0, closure, geom.g_dd, geom.g_uu,
                   geom.n_d, geom.W, geom.u_u, geom.v_d, rounded_projector,
                   E, F, &chi, &P, 1e-12, 128, false, &valid);
      failures += !valid || chi != CCTK_REAL(1) / 3;

      // Do not turn a small positive norm into an isotropic state.
      F(1) = 1e-10 * E;
      Parameters streaming(closure, geom.g_dd, geom.g_uu, geom.n_d, geom.W,
                            geom.u_u, geom.v_d, rounded_projector, E, F);
      failures += !(zFunction(0, &streaming) < 0);

      // A genuinely timelike flux from an inconsistent observer must fail.
      auto bad_u = geom.u_u;
      bad_u(0) *= 2;
      tensor::generic<CCTK_REAL, 4, 1> bad_ud;
      tensor::generic<CCTK_REAL, 4, 2> bad_projector;
      tensor::contract(geom.g_dd, bad_u, &bad_ud);
      calc_proj(bad_ud, bad_u, &bad_projector);
      F(1) = 0;
      Parameters bad(closure, geom.g_dd, geom.g_uu, geom.n_d, geom.W,
                       bad_u, geom.v_d, bad_projector, E, F);
      failures += std::isfinite(zFunction(0, &bad));
      valid = true;
      calc_closure(cctkGH, 0, 0, 0, 0, closure, geom.g_dd, geom.g_uu,
                   geom.n_d, geom.W, bad_u, geom.v_d, bad_projector,
                   E, F, &chi, &P, 1e-12, 128, false, &valid);
      failures += valid;
    }
  }
  return failures;
}

int check_lapse_relaxation(cGH const *cctkGH) {
  int failures = 0;
  constexpr CCTK_REAL Estar = 2.0;
  constexpr CCTK_REAL eta = 0.3;
  constexpr CCTK_REAL kabs = 4.0;
  constexpr CCTK_REAL dt = 0.4;

  for (const CCTK_REAL alpha : {CCTK_REAL(0.25), CCTK_REAL(0.8)}) {
    const source_test_geometry_t geom(alpha, {0.0, 0.0, 0.0});
    tensor::generic<CCTK_REAL, 4, 1> Fstar_d;
    pack_F_d(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, &Fstar_d);
    const SourceUpdateContext ctx(
        cctkGH, 0, 0, 0, 0, 1.0e-12, 64, true, dt, alpha, geom.g_dd, geom.g_uu,
        geom.n_d, geom.n_u, geom.gamma_ud, geom.u_d, geom.u_u, geom.v_d,
        geom.v_u, geom.proj_ud, geom.W, Estar, Fstar_d, Estar, Fstar_d, eta,
        kabs, 0.0);
    CCTK_REAL Enew = Estar;
    tensor::generic<CCTK_REAL, 4, 1> Fnew_d = Fstar_d;
    CCTK_REAL chi = 1.0 / 3.0;
    const int status =
        source_update(ctx, CLOSURE_EDDINGTON, &chi, &Enew, &Fnew_d, true, -1.0,
                      -1.0, 100, 1.0e-14, 1.0e-12);
    const CCTK_REAL expected =
        (Estar + dt * alpha * eta) / (1.0 + dt * alpha * kabs);
    failures += status == NUX_M1_SOURCE_FAIL;
    failures += std::abs(Enew - expected) > 2.0e-12;
    failures += std::abs(Fnew_d(1)) > 2.0e-14 ||
                std::abs(Fnew_d(2)) > 2.0e-14 || std::abs(Fnew_d(3)) > 2.0e-14;
  }
  return failures;
}

} // namespace

int source_math_unit_tests(cGH const *cctkGH) {
  int failures = 0;
  const auto run = [&failures](const char *const name, const int count) {
    if (count != 0)
      CCTK_VINFO("nuX M1 self-test '%s' failed %d checks", name, count);
    failures += count;
  };
  run("source Jacobian transition",
      check_source_jacobian(cctkGH, 0.72,
                            arith_vector{2.0, 0.31, -0.27, 0.19}));
  run("source Jacobian diffusion",
      check_source_jacobian(cctkGH, 0.36,
                            arith_vector{2.0, 0.025, -0.018, 0.011}));
  run("source Jacobian streaming",
      check_source_jacobian(cctkGH, 0.99,
                            arith_vector{2.0, 1.45, -0.62, 0.31}));
  run("source Jacobian isotropic",
      check_source_jacobian(cctkGH, 1.0 / 3.0,
                            arith_vector{2.0, 0.0, 0.0, 0.0}));
  run("moment repair", check_moment_repair());
  run("flux factor", check_flux_factor());
  run("Eddington comoving domain", check_eddington_comoving_domain());
  run("number source", check_number_source());
  run("IMEX source policy", !source_requires_nonlinear(false, false) +
                                !source_requires_nonlinear(false, true) +
                                source_requires_nonlinear(true, false) +
                                !source_requires_nonlinear(true, true));
  run("IMEX number update", check_imex_number_update());
  run("physical IMEX number current", check_physical_number_current());
  run("source convergence contract", check_source_convergence_contract());
  run("recoverable closure trial failure", check_closure_trial_failure(cctkGH));
  run("IMEX diagonal limits", check_imex_diagonal_limits(cctkGH));
  run("IMEX tiny absorbing root", check_tiny_source_root(cctkGH));
  run("source solve", check_source_solve(cctkGH));
  run("repaired source baseline", check_repaired_source_baseline(cctkGH));
  run("vertex velocity repair", check_vertex_velocity_repair());
  run("lapse relaxation", check_lapse_relaxation(cctkGH));
  run("roundoff-negative closure norm", check_roundoff_closure(cctkGH));

  // A low lapse and finite W reduce the collision interval in fluid proper
  // time. This case would be classified stiff if alpha/W were omitted.
  const CCTK_REAL proper_dt = 0.1 * 1.0 / 2.0;
  run("proper-time stiffness classification",
      !source_is_nonstiff(proper_dt, 2.0, 2.0) +
          source_is_nonstiff(1.0, 2.0, 2.0));
  run("proper-time number thermalization classification",
      source_uses_thermalized_number_limit(proper_dt, 10.0, 0.6) +
          !source_uses_thermalized_number_limit(proper_dt, 13.0, 0.6) +
          source_uses_thermalized_number_limit(proper_dt, 13.0, -1.0));
  return failures;
}

} // namespace nuX_M1
