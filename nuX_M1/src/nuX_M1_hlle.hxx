#ifndef NUX_M1_HLLE_HXX
#define NUX_M1_HLLE_HXX

// Private helpers for the optional radiation solver; the Lax-Friedrichs
// kernels do not include this header. Components follow nu_flux:
// N, Fx, Fy, Fz, E.
#include "nuX_M1_closure.hxx"

namespace nuX_M1 {

enum hlle_component : int { HLLE_N = 0, HLLE_FX, HLLE_FY, HLLE_FZ, HLLE_E };

CCTK_HOST CCTK_DEVICE inline bool
hlle_velocity_is_finite(const CCTK_REAL vx, const CCTK_REAL vy,
                        const CCTK_REAL vz) {
  return isfinite(vx) && isfinite(vy) && isfinite(vz);
}

CCTK_HOST CCTK_DEVICE inline bool
hlle_reconstruction_is_finite(const CCTK_REAL E, const CCTK_REAL q[5]) {
  if (!isfinite(E))
    return false;
  for (int a = 0; a < 5; ++a)
    if (!isfinite(q[a]) || (a != HLLE_E && !isfinite(E * q[a])))
      return false;
  return true;
}

CCTK_HOST CCTK_DEVICE inline CCTK_REAL
hlle_slope(const CCTK_REAL qm, const CCTK_REAL q, const CCTK_REAL qp,
           const CCTK_REAL theta) {
  const CCTK_REAL dl = q - qm, dr = qp - q;
  if ((dl > 0 && dr > 0) || (dl < 0 && dr < 0)) {
    const CCTK_REAL mag =
        fmin(theta * fmin(fabs(dl), fabs(dr)), fabs(0.5 * dl + 0.5 * dr));
    return copysign(mag, dl);
  }
  return 0;
}

CCTK_HOST CCTK_DEVICE inline CCTK_REAL
hlle_diffusion_factor(const CCTK_REAL kappa_l, const CCTK_REAL kappa_r,
                      const CCTK_REAL proper_dx) {
  // Opacities are inverse proper lengths. Do not use the Lax-Friedrichs
  // mindiss floor:
  // a nonzero opacity-independent floor destroys the diffusion limit.
  if (!(isfinite(kappa_l) && kappa_l >= 0 && isfinite(kappa_r) &&
        kappa_r >= 0 && isfinite(proper_dx) && proper_dx > 0))
    device_abort();
  const CCTK_REAL tau = (0.5 * kappa_l + 0.5 * kappa_r) * proper_dx;
  return tau > 1 ? 1 / tau : CCTK_REAL(1);
}

// Signed convention: sm <= 0 <= sp. This preserves F(U,U)=F(U) and
// becomes exactly upwind in the transparent limit when every characteristic
// has the same sign. The optically thick momentum correction also acts on
// one-sided fans, as required by the centered diffusion-limit blend.
CCTK_HOST CCTK_DEVICE inline CCTK_REAL
hlle_flux(const CCTK_REAL ul, const CCTK_REAL ur, const CCTK_REAL fl,
          const CCTK_REAL fr, const CCTK_REAL sm, const CCTK_REAL sp,
          const CCTK_REAL A, const bool momentum) {
  const CCTK_REAL span = sp - sm;
  if (!(span > 0))
    return 0.5 * fl + 0.5 * fr;
  // Scale before multiplying speeds to avoid both sp*sm and a tiny
  // denominator. The momentum blend has unit normalization at A=0 and 1.
  const CCTK_REAL wl = sp / span, wr = -sm / span;
  const CCTK_REAL weighted = wl * fl + wr * fr;
  const CCTK_REAL jump = (wl * sm) * (ur - ul);
  if (momentum) {
    const CCTK_REAL A2 = A * A;
    return A2 * (weighted + jump) + (1 - A2) * (0.5 * fl + 0.5 * fr);
  }
  return weighted + A * jump;
}

// Musolino & Rezzolla (2023), Sec. 3.1: blend optically thin and thick
// speed estimates with the SAME chi as the face pressure. These are not
// exact bounds for the intermediate variable-closure flux Jacobian.
CCTK_HOST CCTK_DEVICE inline void
hlle_radiation_speeds(const CCTK_REAL alpha, const CCTK_REAL beta_dir,
                      const CCTK_REAL gamma_uu_dir, const CCTK_REAL W,
                      const CCTK_REAL v_dir, const CCTK_REAL F_up_dir,
                      const CCTK_REAL Fnorm, const CCTK_REAL chi,
                      CCTK_REAL *lm, CCTK_REAL *lp) {
  if (!(isfinite(alpha) && alpha > 0 && isfinite(beta_dir) &&
        isfinite(gamma_uu_dir) && gamma_uu_dir > 0 && isfinite(W) && W >= 1 &&
        isfinite(v_dir) && isfinite(F_up_dir) && isfinite(Fnorm) && Fnorm >= 0 &&
        isfinite(chi) && chi >= CCTK_REAL(1) / 3 && chi <= 1))
    device_abort();

  // Fixed-v thick roots. p=alpha*v^d, not the printed alpha*v^d/W;
  // lower MIN/upper MAX includes the transverse advective characteristic.
  const CCTK_REAL p = alpha * v_dir;
  const CCTK_REAL W2 = W * W, den = 2 * W2 + 1;
  const CCTK_REAL radicand =
      alpha * alpha * gamma_uu_dir * den - 2 * W2 * p * p;
  const CCTK_REAL root = sqrt(fmax(CCTK_REAL(0), radicand));
  const CCTK_REAL adv = -beta_dir + p;
  const CCTK_REAL thick_m = fmin(adv, -beta_dir + (2 * W2 * p - root) / den);
  const CCTK_REAL thick_p = fmax(adv, -beta_dir + (2 * W2 * p + root) / den);
  const CCTK_REAL dthin = (3 * chi - 1) / 2, dthick = 1 - dthin;
  *lm = thick_m;
  *lp = thick_p;
  if (dthin > 0) {
    // Sort the signed thin endpoints for either sign of F^d. At F=0 the
    // direction is undefined; use the two coordinate light-cone endpoints.
    // Isotropic chi=1/3 needs no thin direction and takes no division here.
    const CCTK_REAL thin_offset = Fnorm > 0
        ? alpha * (fabs(F_up_dir) / Fnorm)
        : alpha * sqrt(gamma_uu_dir);
    *lm = dthick * thick_m + dthin * (-beta_dir - thin_offset);
    *lp = dthick * thick_p + dthin * (-beta_dir + thin_offset);
  }
  if (!(isfinite(*lm) && isfinite(*lp) && *lm <= *lp))
    device_abort();
}

struct hlle_geometry_t {
  tensor::metric<4> g_dd;
  tensor::inv_metric<4> g_uu;
  tensor::inv_metric<3> gamma_uu;
  tensor::generic<CCTK_REAL, 4, 1> beta_u, n_d, u_u, v_u, v_d;
  tensor::generic<CCTK_REAL, 4, 2> proj_ud;
  CCTK_REAL alpha, volume, W;

  CCTK_HOST CCTK_DEVICE inline void
  set(const CCTK_REAL lapse, const tensor::generic<CCTK_REAL, 4, 1> &shift,
      const tensor::metric<3> &gamma, CCTK_REAL vx, CCTK_REAL vy, CCTK_REAL vz) {
    alpha = lapse;
    beta_u = shift;
    const CCTK_REAL det = gamma.det();
    if (!(isfinite(alpha) && alpha > 0 && isfinite(det) && det > 0 &&
          gamma(0, 0) > 0 &&
          gamma(0, 0) * gamma(1, 1) > gamma(0, 1) * gamma(0, 1)))
      device_abort();
    for (int a = 1; a < 4; ++a)
      if (!isfinite(beta_u(a)))
        device_abort();
    // The shared helper repairs finite interpolation overshoots, but also
    // maps NaN/Inf to rest. Do not hide corrupt face inputs behind that policy.
    if (!hlle_velocity_is_finite(vx, vy, vz))
      device_abort();
    volume = sqrt(det);
    g_dd.from_adm(alpha, beta_u(1), beta_u(2), beta_u(3), gamma(0, 0),
                  gamma(0, 1), gamma(0, 2), gamma(1, 1), gamma(1, 2),
                  gamma(2, 2));
    g_uu.from_adm(alpha, beta_u(1), beta_u(2), beta_u(3), gamma(0, 0),
                  gamma(0, 1), gamma(0, 2), gamma(1, 1), gamma(1, 2),
                  gamma(2, 2));
    gamma_uu.from_metric(gamma);
    n_d(0) = -alpha;
    n_d(1) = n_d(2) = n_d(3) = 0;
    // An independently averaged W is inconsistent in curved geometry.
    W = repair_velocity_and_compute_W(g_dd, &vx, &vy, &vz);
    pack_v_u(vx, vy, vz, &v_u);
    tensor::contract(g_dd, v_u, &v_d);
    nuX_Utils::valencia::uvel(alpha, beta_u(1), beta_u(2), beta_u(3), W,
                            vx, vy, vz, u_u.data());
    tensor::generic<CCTK_REAL, 4, 1> u_d;
    tensor::contract(g_dd, u_u, &u_d);
    calc_proj(u_d, u_u, &proj_ud);
  }
};

struct hlle_state_t {
  CCTK_REAL u[5], flux[5];
  CCTK_REAL lm, lp;
};

// Kept out of the caller's inline kernel to bound live closure temporaries.
// q = (N/E, Fx/E, Fy/E, Fz/E, E) contains undensitized face primitives.
static CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_NOINLINE hlle_state_t
hlle_state(cGH const *cctkGH, const int i, const int j, const int k,
           const int ig, const int dir, const hlle_geometry_t &g,
           const CCTK_REAL q[5], const closure_t closure_fun,
           const CCTK_REAL E_floor, const CCTK_REAL N_floor,
           const CCTK_REAL eps, const CCTK_REAL closure_tol,
           const CCTK_INT closure_iters, const bool fallback) {
  hlle_state_t result;
  // Reject invalid reconstruction before floors, vacuum handling or cone
  // repair can turn non-finite input into a finite radiation state.
  // Floors are specified for densitized evolved moments elsewhere in M1.
  CCTK_REAL E = fmax(q[HLLE_E], E_floor / g.volume);
  if (!hlle_reconstruction_is_finite(E, q))
    device_abort();
  if (E == 0 && N_floor == 0) {
    // Exact vacuum: all physical moments and fluxes vanish, independently of
    // reconstructed direction ratios. Do not divide 0/0 or rootfind a closure
    // whose entire vacuum equation is identically zero. Use the isotropic
    // regular-limit speed estimate, without assigning a thin flux direction.
    for (int a = 0; a < 5; ++a)
      result.u[a] = result.flux[a] = 0;
    hlle_radiation_speeds(g.alpha, g.beta_u(dir + 1), g.gamma_uu(dir, dir),
                          g.W, g.v_u(dir + 1), 0, 0, CCTK_REAL(1) / 3,
                          &result.lm, &result.lp);
    return result;
  }
  tensor::generic<CCTK_REAL, 4, 1> F_d, F_u;
  pack_F_d(g.beta_u(1), g.beta_u(2), g.beta_u(3), E * q[HLLE_FX],
           E * q[HLLE_FY], E * q[HLLE_FZ], &F_d);
  for (int a = 0; a < 4; ++a)
    if (!isfinite(F_d(a)))
      device_abort();
  // Componentwise limiting of Fi/E alone is not a multidimensional cone
  // constraint. Project once in the COMMON face metric after reconstruction.
  repair_moments(g.g_uu, &E, &F_d, E_floor / g.volume, eps);
  tensor::contract(g.g_uu, F_d, &F_u);

  const CCTK_REAL Fnorm = sqrt(fmax(CCTK_REAL(0), tensor::dot(F_u, F_d)));
  const bool static_frame =
      g.v_u(1) == 0 && g.v_u(2) == 0 && g.v_u(3) == 0;
  CCTK_REAL chi;
  tensor::symmetric2<CCTK_REAL, 4, 2> P_dd;
  if (static_frame && (closure_fun == CLOSURE_MINERBO ||
                       closure_fun == CLOSURE_KERSHAW)) {
    // Here J=E, |H|=|F|, so no nonlinear closure inversion is required.
    // The face pressure and speed interpolation use this same chi(f).
    chi = eval_closure(closure_fun, fmin(CCTK_REAL(1), Fnorm / E));
    apply_closure(g.g_dd, g.g_uu, g.n_d, g.W, g.u_u, g.v_d, g.proj_ud,
                  E, F_d, chi, &P_dd);
  } else {
    calc_closure(cctkGH, i, j, k, ig, closure_fun, g.g_dd, g.g_uu, g.n_d,
                 g.W, g.u_u, g.v_d, g.proj_ud, E, F_d, &chi, &P_dd,
                 closure_tol, closure_iters, fallback);
  }
  tensor::generic<CCTK_REAL, 4, 2> P_ud;
  tensor::contract(g.g_uu, P_dd, &P_ud);
  tensor::symmetric2<CCTK_REAL, 4, 2> T_dd;
  assemble_rT(g.n_d, E, F_d, P_dd, &T_dd);
  const CCTK_REAL J = calc_J_from_rT(T_dd, g.u_u);
  tensor::generic<CCTK_REAL, 4, 1> H_d, H_u, fnu_u;
  calc_H_from_rT(T_dd, g.u_u, g.proj_ud, &H_d);
  if (!comoving_state_is_acceptable(closure_fun, g.g_uu, J, H_d))
    device_abort();
  tensor::contract(g.g_uu, H_d, &H_u);
  assemble_fnu(g.u_u, J, H_u, &fnu_u, E_floor / g.volume);
  // Normalize the CURRENT actually used above: N=sqrt(gamma)*alpha*n*fnu^0.
  // This also gives Gamma=W when the floor regularization removes H/J.
  // Independently falling back to Gamma=1 would mis-map moving floor states.
  const CCTK_REAL Gamma = g.alpha * fnu_u(0);
  if (!(isfinite(Gamma) && Gamma > 0))
    device_abort();
  const CCTK_REAL N = fmax(E * fmax(CCTK_REAL(0), q[HLLE_N]),
                           N_floor / g.volume);
  const CCTK_REAL number_speed = g.alpha * fnu_u(dir + 1) / Gamma;
  result.u[HLLE_N] = g.volume * N;
  result.flux[HLLE_N] = result.u[HLLE_N] * number_speed;
  result.u[HLLE_E] = g.volume * E;
  result.flux[HLLE_E] =
      g.volume * calc_E_flux(g.alpha, g.beta_u, E, F_u, dir + 1);
  for (int a = 1; a < 4; ++a) {
    result.u[a] = g.volume * F_d(a);
    result.flux[a] =
        g.volume * calc_F_flux(g.alpha, g.beta_u, F_d, P_ud, dir + 1, a);
  }
  hlle_radiation_speeds(g.alpha, g.beta_u(dir + 1), g.gamma_uu(dir, dir),
                        g.W, g.v_u(dir + 1), F_u(dir + 1), Fnorm, chi,
                        &result.lm, &result.lp);
  // N is advected by the same closure-derived number current. Include its
  // scalar eigenvalue, which need not be inside the interpolated estimates.
  result.lm = fmin(result.lm, number_speed);
  result.lp = fmax(result.lp, number_speed);
  for (int a = 0; a < 5; ++a)
    if (!(isfinite(result.u[a]) && isfinite(result.flux[a])))
      device_abort();
  return result;
}

} // namespace nuX_M1

#endif // NUX_M1_HLLE_HXX
