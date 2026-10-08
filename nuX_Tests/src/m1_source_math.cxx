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
  run("source solve", check_source_solve(cctkGH));
  run("repaired source baseline", check_repaired_source_baseline(cctkGH));
  run("vertex velocity repair", check_vertex_velocity_repair());
  run("lapse relaxation", check_lapse_relaxation(cctkGH));

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
