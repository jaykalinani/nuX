#include "m1_unit_tests.hxx"
#define NUX_M1_CLOSURE_IMPLEMENTATION
#include "nuX_M1_hlle.hxx"

#include <algorithm>
#include <cmath>

namespace nuX_M1 {
namespace {

bool near(const CCTK_REAL a, const CCTK_REAL b,
          const CCTK_REAL tolerance = 2.0e-11) {
  return std::isfinite(a) && std::isfinite(b) &&
         std::abs(a - b) <= tolerance *
             std::max({CCTK_REAL(1), std::abs(a), std::abs(b)});
}

int check_hlle_algebra() {
  int failures = 0;
  for (const CCTK_REAL A : {CCTK_REAL(0), CCTK_REAL(0.03), CCTK_REAL(1)}) {
    for (const bool momentum : {false, true}) {
      // A constant state must retain its complete physical flux, not half
      // the flux (the normalization printed in paper Eq. 64 is incorrect).
      failures += !near(hlle_flux(2, 2, 3, 3, -0.4, 0.8, A, momentum), 3);
      failures += !near(hlle_flux(2, 2, 3, 3, 0, 0, A, momentum), 3);
    }
  }
  failures += !near(hlle_flux(1, 3, 5, 7, 0, 2, 1, false), 5);
  failures += !near(hlle_flux(1, 3, 5, 7, -2, 0, 1, false), 7);
  failures += !near(hlle_flux(1, 3, 5, 7, 0, 2, 1, true), 5);
  failures += !near(hlle_flux(1, 3, 5, 7, -2, 0, 1, true), 7);
  failures += !near(hlle_flux(1, 3, 5, 7, 0, 2, 0, true), 6);
  // Signed HLLE reference: (2*5-(-1)*7+2*(-1)*(3-1))/(2-(-1)).
  failures += !near(hlle_flux(1, 3, 5, 7, -1, 2, 1, false), 13.0 / 3.0);
  failures += !near(hlle_flux(1, 3, 5, 7, -1, 2, 0, false), 17.0 / 3.0);
  failures += !near(hlle_flux(1, 3, 5, 7, -1, 2, 0, true), 6);
  // O(A) and O(A^2) jump viscosity at high optical depth.
  failures += !near(hlle_flux(1, 3, 0, 0, -1, 1, 0.1, false), -0.1);
  failures += !near(hlle_flux(1, 3, 0, 0, -1, 1, 0.1, true), -0.01);
  failures += !near(hlle_diffusion_factor(0, 0, 0.5), 1);
  failures += !near(hlle_diffusion_factor(100, 300, 0.5), 0.01);
  failures += !near(hlle_slope(0, 1, 2, 2), 1);
  failures += !near(hlle_slope(2, 1, 0, 2), -1);
  failures += !near(hlle_slope(0, 2, 1, 2), 0);
  failures += !near(hlle_slope(0, 0.1, 2, 1), 0.1);
  failures += !near(hlle_slope(0, 0.1, 2, 2), 0.2);
  return failures;
}

int check_face_input_validation() {
  int failures = 0;
  const CCTK_REAL nan = std::numeric_limits<CCTK_REAL>::quiet_NaN();
  const CCTK_REAL inf = std::numeric_limits<CCTK_REAL>::infinity();
  failures += !hlle_velocity_is_finite(0.2, -0.1, 0.05);
  for (const CCTK_REAL invalid : {nan, inf, -inf}) {
    failures += hlle_velocity_is_finite(invalid, 0, 0);
    failures += hlle_velocity_is_finite(0, invalid, 0);
    failures += hlle_velocity_is_finite(0, 0, invalid);
  }
  const CCTK_REAL vacuum[5] = {};
  failures += !hlle_reconstruction_is_finite(0, vacuum);
  CCTK_REAL q[5] = {0.4, 0.15, -0.08, 0.06, 2};
  failures += !hlle_reconstruction_is_finite(q[HLLE_E], q);
  for (int a = 0; a < 5; ++a) {
    const CCTK_REAL saved = q[a];
    q[a] = nan;
    failures += hlle_reconstruction_is_finite(2, q);
    q[a] = saved;
  }
  // Finite primitive ratios can still overflow when multiplied by E. Reject
  // that arithmetic before cone repair can silently turn the flux into zero.
  q[HLLE_FX] = 2;
  failures += hlle_reconstruction_is_finite(
      std::numeric_limits<CCTK_REAL>::max(), q);
  return failures;
}

int check_radiation_speeds() {
  int failures = 0;
  CCTK_REAL lm, lp;
  // Exact 1D fixed-v thick roots are the relativistic addition of v and
  // +/-1/sqrt(3). In particular, alpha*v/W fails this moving-medium test.
  const CCTK_REAL cs = 1 / std::sqrt(CCTK_REAL(3));
  for (const CCTK_REAL v : {CCTK_REAL(-0.8), CCTK_REAL(0), CCTK_REAL(0.5)}) {
    const CCTK_REAL W = 1 / std::sqrt(1 - v * v);
    hlle_radiation_speeds(0.7, 0.2, 1, W, v, 0, 0, 1,
                          CLOSURE_EDDINGTON, v == 0, &lm, &lp);
    failures += !near(lm, -0.2 + 0.7 * (v - cs) / (1 - v * cs));
    failures += !near(lp, -0.2 + 0.7 * (v + cs) / (1 + v * cs));
  }
  // The paper's chi-interpolated estimate underbounds this exact static
  // Minerbo Jacobian: lambda+ = .7541494487418, not .6460308504.
  hlle_static_speeds(0.5, 1, CLOSURE_MINERBO, &lm, &lp);
  failures += lm > -0.2541494487418 || lp < 0.7541494487418;
  failures += !near(lm, -0.2541494487418, 2e-10);
  failures += !near(lp, 0.7541494487418, 2e-10);
  // Independent 4x4 pressure-Jacobian eigenvalues for oblique/transverse F.
  hlle_static_speeds(0.5, 0.3, CLOSURE_MINERBO, &lm, &lp);
  failures += !near(lm, -0.4933818225584, 2e-10);
  failures += !near(lp, 0.6315309560944, 2e-10);
  hlle_static_speeds(0.5, 0, CLOSURE_MINERBO, &lm, &lp);
  failures += !near(lm, -0.5682575707077, 2e-10);
  failures += !near(lp, 0.5682575707077, 2e-10);
  hlle_static_speeds(0.5, 0.3, CLOSURE_KERSHAW, &lm, &lp);
  failures += !near(lm, -0.4024937810560, 2e-10);
  failures += !near(lp, 0.6024937810560, 2e-10);
  hlle_static_speeds(1, 0.4, CLOSURE_KERSHAW, &lm, &lp);
  failures += !near(lm, 0.4 / 3);
  failures += !near(lp, 0.4);
  hlle_static_speeds(1, -0.4, CLOSURE_KERSHAW, &lm, &lp);
  failures += !near(lm, -0.4);
  failures += !near(lp, -0.4 / 3);
  // Shift sign, lapse scaling, inverse metric, and oblique streaming. With
  // gamma^{dd}=.25, F^d=-.2, |F|=.5, the unit-normal cosine is -.8.
  hlle_radiation_speeds(0.6, -0.3, 0.25, 1, 0, -0.2, 0.5, 0.5,
                        CLOSURE_MINERBO, true, &lm, &lp);
  failures += !near(lm, 0.3 - 0.24);
  failures += !near(lp, 0.3 - 0.24);
  // The isotropic exact limit does not divide by |F| or a cubic discriminant.
  hlle_radiation_speeds(0.6, 0.1, 0.25, 1, 0, 0, 0, 1,
                        CLOSURE_MINERBO, true, &lm, &lp);
  failures += !near(lm, -0.1 - 0.3 / std::sqrt(CCTK_REAL(3)));
  failures += !near(lp, -0.1 + 0.3 / std::sqrt(CCTK_REAL(3)));
  // Check the current moving estimate's coordinate mapping only. This does
  // NOT establish enclosure of the inherited moving pressure Jacobian.
  hlle_radiation_speeds(0.6, 0.1, 0.25, 1.1, 0.1, 0.1, 0.3, 1,
                        CLOSURE_MINERBO, false, &lm, &lp);
  failures += !near(lm, -0.4);
  failures += !near(lp, 0.2);
  for (const CCTK_REAL f : {CCTK_REAL(1e-12), CCTK_REAL(1 - 1e-12)}) {
    hlle_static_speeds(f, 0.4, CLOSURE_MINERBO, &lm, &lp);
    failures += !(std::isfinite(lm) && std::isfinite(lp) && lm <= lp &&
                  lm >= -1 && lp <= 1);
  }
  return failures;
}

int check_face_mapping(cGH const *cctkGH) {
  int failures = 0;
  tensor::metric<3> gamma;
  gamma(0, 0) = 1.4;
  gamma(0, 1) = 0.12;
  gamma(0, 2) = -0.07;
  gamma(1, 1) = 0.9;
  gamma(1, 2) = 0.04;
  gamma(2, 2) = 1.2;
  tensor::generic<CCTK_REAL, 4, 1> beta;
  beta(0) = 0;
  beta(1) = 0.2;
  beta(2) = -0.1;
  beta(3) = 0.03;
  hlle_geometry_t g;
  g.set(0.7, beta, gamma, 0, 0, 0);
  const CCTK_REAL q[5] = {0.4, 0.15, -0.08, 0.06, 2};
  for (int dir = 0; dir < 3; ++dir) {
    const auto state = hlle_state(cctkGH, 0, 0, 0, 0, dir, g, q,
                                  CLOSURE_EDDINGTON, 1e-15, 1e-15, 1e-12,
                                  1e-12, 128, true);
    CCTK_REAL Fup = 0;
    for (int a = 0; a < 3; ++a)
      Fup += g.gamma_uu(dir, a) * q[a + 1] * q[HLLE_E];
    failures += !near(state.u[HLLE_N], g.volume * q[HLLE_N] * q[HLLE_E]);
    failures += !near(state.flux[HLLE_N],
                       g.volume * q[HLLE_N] *
                           (g.alpha * Fup - beta(dir + 1) * q[HLLE_E]));
    failures += !near(state.flux[HLLE_E],
                       g.volume * (g.alpha * Fup - beta(dir + 1) * q[HLLE_E]));
    for (int a = 0; a < 3; ++a) {
      failures += !near(state.u[a + 1], g.volume * q[a + 1] * q[HLLE_E]);
      failures += !near(
          state.flux[a + 1],
          g.volume * (g.alpha * q[HLLE_E] * (a == dir) / 3 -
                      beta(dir + 1) * q[a + 1] * q[HLLE_E]));
    }
    // Equal face states remain consistent across opacity and component map.
    for (int a = 0; a < 5; ++a)
      failures += !near(
          hlle_flux(state.u[a], state.u[a], state.flux[a], state.flux[a],
                     std::min(CCTK_REAL(0), state.lm),
                     std::max(CCTK_REAL(0), state.lp), 0.01, a >= 1 && a <= 3),
          state.flux[a]);
  }
  // A componentwise reconstructed ratio may leave the cone; the face
  // projection must repair it using the non-diagonal inverse face metric.
  const CCTK_REAL unphysical[5] = {0.4, 2, -3, 4, 2};
  const auto repaired = hlle_state(cctkGH, 0, 0, 0, 0, 1, g, unphysical,
                                    CLOSURE_MINERBO, 1e-15, 1e-15, 1e-8,
                                    1e-12, 256, true);
  CCTK_REAL F2 = 0;
  for (int a = 0; a < 3; ++a)
    for (int b = 0; b < 3; ++b)
      F2 += g.gamma_uu(a, b) * repaired.u[a + 1] * repaired.u[b + 1];
  failures += F2 > repaired.u[HLLE_E] * repaired.u[HLLE_E] * (1 - 0.5e-8);
  // W comes from the face velocity and metric, not an independently
  // interpolated Lorentz factor. This also covers moving-medium projection.
  g.set(0.7, beta, gamma, 0.2, -0.1, 0.05);
  const CCTK_REAL v2 = 1.4 * 0.04 + 0.9 * 0.01 + 1.2 * 0.0025 +
                       2 * (0.12 * -0.02 - 0.07 * 0.01 + 0.04 * -0.005);
  failures += !near(g.W, 1 / std::sqrt(1 - v2));
  // Both zero floors are legal. A true zero-moment vacuum must not divide by
  // E or J, and has exactly zero numerical/physical flux in any direction.
  const CCTK_REAL vacuum[5] = {};
  for (const closure_t model : {CLOSURE_EDDINGTON, CLOSURE_MINERBO,
                                 CLOSURE_KERSHAW}) {
    for (int dir = 0; dir < 3; ++dir) {
      const auto empty = hlle_state(cctkGH, 0, 0, 0, 0, dir, g, vacuum,
                                     model, 0, 0, 0, 1e-12, 128, false);
      for (int a = 0; a < 5; ++a)
        failures += empty.u[a] != 0 || empty.flux[a] != 0;
      failures += !(std::isfinite(empty.lm) && std::isfinite(empty.lp) &&
                    empty.lm <= empty.lp);
    }
  }
  return failures;
}

} // namespace

int hlle_unit_tests(cGH const *cctkGH) {
  return check_hlle_algebra() + check_face_input_validation() +
         check_radiation_speeds() + check_face_mapping(cctkGH);
}

} // namespace nuX_M1
