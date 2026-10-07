#ifndef NUX_M1_CLOSURE_HXX
#define NUX_M1_CLOSURE_HXX

#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>
#include <loop_device.hxx>

#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "nuX_utils.hxx"

namespace nuX_M1 {

using namespace nuX_Utils;
using namespace Loop;
using namespace std;

[[noreturn]] CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
device_abort() {
#if defined(__CUDA_ARCH__)
  asm volatile("trap;");
  while (true) {
  }
#else
  __builtin_trap();
#endif
}

enum closure_t : int {
  CLOSURE_EDDINGTON = 0,
  CLOSURE_KERSHAW = 1,
  CLOSURE_MINERBO = 2,
  CLOSURE_THIN = 3,
};

#ifdef NUX_M1_CLOSURE_IMPLEMENTATION

struct Parameters {
  CCTK_HOST CCTK_DEVICE Parameters(
      closure_t _closure, tensor::metric<4> const &_g_dd,
      tensor::inv_metric<4> const &_g_uu,
      tensor::generic<CCTK_REAL, 4, 1> const &_n_d, CCTK_REAL const _w_lorentz,
      tensor::generic<CCTK_REAL, 4, 1> const &_u_u,
      tensor::generic<CCTK_REAL, 4, 1> const &_v_d,
      tensor::generic<CCTK_REAL, 4, 2> const &_proj_ud, CCTK_REAL const _E,
      tensor::generic<CCTK_REAL, 4, 1> const &_F_d)
      : closure(_closure), g_dd(_g_dd), g_uu(_g_uu), n_d(_n_d),
        w_lorentz(_w_lorentz), u_u(_u_u), v_d(_v_d), proj_ud(_proj_ud), E(_E),
        F_d(_F_d) {}
  closure_t closure;
  tensor::metric<4> const &g_dd;
  tensor::inv_metric<4> const &g_uu;
  tensor::generic<CCTK_REAL, 4, 1> const &n_d;
  CCTK_REAL const w_lorentz;
  tensor::generic<CCTK_REAL, 4, 1> const &u_u;
  tensor::generic<CCTK_REAL, 4, 1> const &v_d;
  tensor::generic<CCTK_REAL, 4, 2> const &proj_ud;
  CCTK_REAL const E;
  tensor::generic<CCTK_REAL, 4, 1> const &F_d;
};

#endif // NUX_M1_CLOSURE_IMPLEMENTATION

enum ClosFlag : CCTK_INT {
  CLOS_OK = 0,   // closure success
  CLOS_I = 1,    // initial value closure NaN or inf
  ROOT_I = 2,    // initial root solver error
  CLOS_IT = 3,   // iteration closure NaN or inf
  ROOT_IT = 4,   // iterative root solver error
  ROOT_MAXIT = 5 // root solver reached max iterations
};

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
pack_F_d(CCTK_REAL const betax, CCTK_REAL const betay, CCTK_REAL const betaz,
         CCTK_REAL const Fx, CCTK_REAL const Fy, CCTK_REAL const Fz,
         tensor::generic<CCTK_REAL, 4, 1> *F_d) {
  // F_0 = g_0i F^i = beta_i F^i = beta^i F_i
  F_d->at(0) = betax * Fx + betay * Fy + betaz * Fz;
  F_d->at(1) = Fx;
  F_d->at(2) = Fy;
  F_d->at(3) = Fz;
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
unpack_F_d(tensor::generic<CCTK_REAL, 4, 1> const &F_d, CCTK_REAL *Fx,
           CCTK_REAL *Fy, CCTK_REAL *Fz) {
  *Fx = F_d(1);
  *Fy = F_d(2);
  *Fz = F_d(3);
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
pack_F_d(CCTK_REAL const Fx, CCTK_REAL const Fy, CCTK_REAL const Fz,
         tensor::generic<CCTK_REAL, 3, 1> *F_d) {
  F_d->at(0) = Fx;
  F_d->at(1) = Fy;
  F_d->at(2) = Fz;
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
unpack_F_d(tensor::generic<CCTK_REAL, 3, 1> const &F_d, CCTK_REAL *Fx,
           CCTK_REAL *Fy, CCTK_REAL *Fz) {
  *Fx = F_d(0);
  *Fy = F_d(1);
  *Fz = F_d(2);
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
pack_H_d(CCTK_REAL const Ht, CCTK_REAL const Hx, CCTK_REAL const Hy,
         CCTK_REAL const Hz, tensor::generic<CCTK_REAL, 4, 1> *H_d) {
  H_d->at(0) = Ht;
  H_d->at(1) = Hx;
  H_d->at(2) = Hy;
  H_d->at(3) = Hz;
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
unpack_H_d(tensor::generic<CCTK_REAL, 4, 1> const &H_d, CCTK_REAL *Ht,
           CCTK_REAL *Hx, CCTK_REAL *Hy, CCTK_REAL *Hz) {
  *Ht = H_d(0);
  *Hx = H_d(1);
  *Hy = H_d(2);
  *Hz = H_d(3);
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
pack_P_dd(CCTK_REAL const betax, CCTK_REAL const betay, CCTK_REAL const betaz,
          CCTK_REAL const Pxx, CCTK_REAL const Pxy, CCTK_REAL const Pxz,
          CCTK_REAL const Pyy, CCTK_REAL const Pyz, CCTK_REAL const Pzz,
          tensor::symmetric2<CCTK_REAL, 4, 2> *P_dd) {
  CCTK_REAL const Pbetax = Pxx * betax + Pxy * betay + Pxz * betaz;
  CCTK_REAL const Pbetay = Pxy * betax + Pyy * betay + Pyz * betaz;
  CCTK_REAL const Pbetaz = Pxz * betax + Pyz * betay + Pzz * betaz;

  // P_00 = g_0i g_k0 P^ik = beta^i beta^k P_ik
  P_dd->at(0, 0) = Pbetax * betax + Pbetay * betay + Pbetaz * betaz;

  // P_0i = g_0j g_ki P^jk = beta_j P_i^j = beta^j P_ij
  P_dd->at(0, 1) = Pbetax;
  P_dd->at(0, 2) = Pbetay;
  P_dd->at(0, 3) = Pbetaz;

  P_dd->at(1, 1) = Pxx;
  P_dd->at(1, 2) = Pxy;
  P_dd->at(1, 3) = Pxz;
  P_dd->at(2, 2) = Pyy;
  P_dd->at(2, 3) = Pyz;
  P_dd->at(3, 3) = Pzz;
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
unpack_P_dd(tensor::symmetric2<CCTK_REAL, 4, 2> const &P_dd, CCTK_REAL *Pxx,
            CCTK_REAL *Pxy, CCTK_REAL *Pxz, CCTK_REAL *Pyy, CCTK_REAL *Pyz,
            CCTK_REAL *Pzz) {
  *Pxx = P_dd(1, 1);
  *Pxy = P_dd(1, 2);
  *Pxz = P_dd(1, 3);
  *Pyy = P_dd(2, 2);
  *Pyz = P_dd(2, 3);
  *Pzz = P_dd(3, 3);
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
pack_P_dd(CCTK_REAL const Pxx, CCTK_REAL const Pxy, CCTK_REAL const Pxz,
          CCTK_REAL const Pyy, CCTK_REAL const Pyz, CCTK_REAL const Pzz,
          tensor::symmetric2<CCTK_REAL, 3, 2> *P_dd) {
  P_dd->at(0, 0) = Pxx;
  P_dd->at(0, 1) = Pxy;
  P_dd->at(0, 2) = Pxz;
  P_dd->at(1, 1) = Pyy;
  P_dd->at(1, 2) = Pyz;
  P_dd->at(2, 2) = Pzz;
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
unpack_P_dd(tensor::symmetric2<CCTK_REAL, 3, 2> const &P_dd, CCTK_REAL *Pxx,
            CCTK_REAL *Pxy, CCTK_REAL *Pxz, CCTK_REAL *Pyy, CCTK_REAL *Pyz,
            CCTK_REAL *Pzz) {
  *Pxx = P_dd(0, 0);
  *Pxy = P_dd(0, 1);
  *Pxz = P_dd(0, 2);
  *Pyy = P_dd(1, 1);
  *Pyz = P_dd(1, 2);
  *Pzz = P_dd(2, 2);
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
pack_v_u(CCTK_REAL const velx, CCTK_REAL const vely, CCTK_REAL const velz,
         tensor::generic<CCTK_REAL, 4, 1> *v_u) {
  v_u->at(0) = 0.0;
  v_u->at(1) = velx;
  v_u->at(2) = vely;
  v_u->at(3) = velz;
}

// Fluid projector: delta^a_b + u^a u_b
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
calc_proj(tensor::generic<CCTK_REAL, 4, 1> const &u_d,
          tensor::generic<CCTK_REAL, 4, 1> const &u_u,
          tensor::generic<CCTK_REAL, 4, 2> *proj_ud) {
  for (int a = 0; a < 4; ++a)
    for (int b = 0; b < 4; ++b) {
      proj_ud->at(a, b) = tensor::delta(a, b) + u_u(a) * u_d(b);
    }
}

// Compute the closure in the thin limit
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
calc_Pthin(tensor::inv_metric<4> const &g_uu, CCTK_REAL const E,
           tensor::generic<CCTK_REAL, 4, 1> const &F_d,
           tensor::symmetric2<CCTK_REAL, 4, 2> *P_dd) {
  CCTK_REAL const F2 = tensor::dot(g_uu, F_d, F_d);
  CCTK_REAL fac = 0.0;
  if (isfinite(E) && isfinite(F2) && F2 > 0.0) {
    fac = E / F2;
    if (!isfinite(fac)) {
      fac = 0.0;
    }
  }
  for (int a = 0; a < 4; ++a)
    for (int b = a; b < 4; ++b) {
      P_dd->at(a, b) = fac * F_d(a) * F_d(b);
    }
}

// Compute the closure in the thick limit
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
calc_Pthick(tensor::metric<4> const &g_dd, tensor::inv_metric<4> const &g_uu,
            tensor::generic<CCTK_REAL, 4, 1> const &n_d, CCTK_REAL const W,
            tensor::generic<CCTK_REAL, 4, 1> const &v_d, CCTK_REAL const E,
            tensor::generic<CCTK_REAL, 4, 1> const &F_d,
            tensor::symmetric2<CCTK_REAL, 4, 2> *P_dd) {
  CCTK_REAL const v_dot_F = tensor::dot(g_uu, v_d, F_d);

  CCTK_REAL const W2 = W * W;
  CCTK_REAL const coef = 1. / (2. * W2 + 1.);

  // J/3
  CCTK_REAL const Jo3 = coef * ((2. * W2 - 1.) * E - 2. * W2 * v_dot_F);

  // tH = gamma_ud H_d
  tensor::generic<CCTK_REAL, 4, 1> tH_d;
  for (int a = 0; a < 4; ++a) {
    tH_d(a) = F_d(a) / W +
              coef * W * v_d(a) * ((4. * W2 + 1.) * v_dot_F - 4. * W2 * E);
  }

  for (int a = 0; a < 4; ++a)
    for (int b = a; b < 4; ++b) {
      P_dd->at(a, b) =
          Jo3 * (4. * W2 * v_d(a) * v_d(b) + g_dd(a, b) + n_d(a) * n_d(b));
      P_dd->at(a, b) += W * (tH_d(a) * v_d(b) + tH_d(b) * v_d(a));
    }
}

// Computes the comoving flux factor xi = sqrt(H_a H^a) / J.
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
flux_factor(tensor::inv_metric<4> const &g_uu, CCTK_REAL const J,
            tensor::generic<CCTK_REAL, 4, 1> const &H_d,
            CCTK_REAL rad_E_floor) {
  if (!(J > rad_E_floor) || !isfinite(J))
    return CCTK_REAL(0);
  const CCTK_REAL H2 = tensor::dot(g_uu, H_d, H_d);
  if (!isfinite(H2) || H2 <= CCTK_REAL(0))
    return CCTK_REAL(0);
  const CCTK_REAL xi = sqrt(H2) / J;
  return max(CCTK_REAL(0), min(xi, CCTK_REAL(1)));
}

// Closures
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
eddington(CCTK_REAL const xi) {
  return 1.0 / 3.0;
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
kershaw(CCTK_REAL const xi) {
  return 1.0 / 3.0 + 2.0 / 3.0 * xi * xi;
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
minerbo(CCTK_REAL const xi) {
  return 1.0 / 3.0 + xi * xi * (6.0 - 2.0 * xi + 6.0 * xi * xi) / 15.0;
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
thin(CCTK_REAL const xi) {
  return 1.0;
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline bool
closure_is_eddington(closure_t const closure) {
  return closure == CLOSURE_EDDINGTON;
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline bool
closure_is_thin(closure_t const closure) {
  return closure == CLOSURE_THIN;
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
eval_closure(closure_t const closure, CCTK_REAL const xi) {
  switch (closure) {
  case CLOSURE_EDDINGTON:
    return eddington(xi);
  case CLOSURE_KERSHAW:
    return kershaw(xi);
  case CLOSURE_MINERBO:
    return minerbo(xi);
  case CLOSURE_THIN:
    return thin(xi);
  }
  return eddington(xi);
}

// Computes the closure in the lab frame given the Eddington factor chi
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
apply_closure(tensor::metric<4> const &g_dd, tensor::inv_metric<4> const &g_uu,
              tensor::generic<CCTK_REAL, 4, 1> const &n_d,
              CCTK_REAL const w_lorentz,
              tensor::generic<CCTK_REAL, 4, 1> const &u_u,
              tensor::generic<CCTK_REAL, 4, 1> const &v_d,
              tensor::generic<CCTK_REAL, 4, 2> const &proj_ud,
              CCTK_REAL const E, tensor::generic<CCTK_REAL, 4, 1> const &F_d,
              CCTK_REAL const chi, tensor::symmetric2<CCTK_REAL, 4, 2> *P_dd) {
  CCTK_REAL const chi_phys =
      max(CCTK_REAL(1.0 / 3.0), min(CCTK_REAL(1.0), chi));
  CCTK_REAL const dthick = 3. * (1 - chi_phys) / 2.;
  CCTK_REAL const dthin = 1. - dthick;
  CCTK_REAL const coeff_eps =
      CCTK_REAL(64.0) * std::numeric_limits<CCTK_REAL>::epsilon();
  bool const use_thick = isfinite(dthick) && abs(dthick) > coeff_eps;
  bool const use_thin = isfinite(dthin) && abs(dthin) > coeff_eps;

  tensor::symmetric2<CCTK_REAL, 4, 2> Pthin_dd;
  tensor::symmetric2<CCTK_REAL, 4, 2> Pthick_dd;

  if (use_thin) {
    calc_Pthin(g_uu, E, F_d, &Pthin_dd);
  }
  if (use_thick) {
    calc_Pthick(g_dd, g_uu, n_d, w_lorentz, v_d, E, F_d, &Pthick_dd);
  }

  for (int a = 0; a < 4; ++a)
    for (int b = a; b < 4; ++b) {
      CCTK_REAL const thick_term = use_thick ? dthick * Pthick_dd(a, b) : 0.0;
      CCTK_REAL const thin_term = use_thin ? dthin * Pthin_dd(a, b) : 0.0;
      P_dd->at(a, b) = thick_term + thin_term;
    }
}

// Assemble the unit-norm radiation number current
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
assemble_fnu(tensor::generic<CCTK_REAL, 4, 1> const &u_u, CCTK_REAL const J,
             tensor::generic<CCTK_REAL, 4, 1> const &H_u,
             tensor::generic<CCTK_REAL, 4, 1> *fnu_u, CCTK_REAL rad_E_floor) {
  for (int a = 0; a < 4; ++a) {
    fnu_u->at(a) = u_u(a) + (J > rad_E_floor ? H_u(a) / J : 0);
  }
}

// Compute the ratio of neutrino number densities in the lab and fluid frame
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
compute_Gamma(CCTK_REAL const W, tensor::generic<CCTK_REAL, 4, 1> const &v_u,
              CCTK_REAL const J, CCTK_REAL const E,
              tensor::generic<CCTK_REAL, 4, 1> const &F_d,
              CCTK_REAL rad_E_floor, CCTK_REAL rad_eps) {
  if (!(isfinite(W) && W >= CCTK_REAL(1) && isfinite(E) && E > rad_E_floor &&
        isfinite(J) && J > rad_E_floor))
    return CCTK_REAL(1);

  const CCTK_REAL raw_f_dot_v = tensor::dot(F_d, v_u) / E;
  if (!isfinite(raw_f_dot_v))
    return CCTK_REAL(1);
  // repair_moments enforces F_a F^a/E^2 <= 1-rad_eps, so the corresponding
  // bound on the (unsquared) flux factor is sqrt(1-rad_eps).
  const CCTK_REAL f_limit = sqrt(max(CCTK_REAL(0), CCTK_REAL(1) - rad_eps));
  const CCTK_REAL f_dot_v = max(-f_limit, min(raw_f_dot_v, f_limit));
  const CCTK_REAL Gamma = W * (E / J) * (CCTK_REAL(1) - f_dot_v);
  return isfinite(Gamma) && Gamma > CCTK_REAL(0) ? Gamma : CCTK_REAL(1);
}

// Assemble the radiation stress tensor in any frame
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
assemble_rT(tensor::generic<CCTK_REAL, 4, 1> const &u_d, CCTK_REAL const J,
            tensor::generic<CCTK_REAL, 4, 1> const &H_d,
            tensor::symmetric2<CCTK_REAL, 4, 2> const &K_dd,
            tensor::symmetric2<CCTK_REAL, 4, 2> *rT_dd) {
  for (int a = 0; a < 4; ++a)
    for (int b = a; b < 4; ++b) {
      rT_dd->at(a, b) =
          J * u_d(a) * u_d(b) + H_d(a) * u_d(b) + H_d(b) * u_d(a) + K_dd(a, b);
    }
}

// Project out the radiation energy (in any frame)
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
calc_J_from_rT(tensor::symmetric2<CCTK_REAL, 4, 2> const &rT_dd,
               tensor::generic<CCTK_REAL, 4, 1> const &u_u) {
  return tensor::dot(rT_dd, u_u, u_u);
}

// Project out the radiation fluxes (in any frame)
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
calc_H_from_rT(tensor::symmetric2<CCTK_REAL, 4, 2> const &rT_dd,
               tensor::generic<CCTK_REAL, 4, 1> const &u_u,
               tensor::generic<CCTK_REAL, 4, 2> const &proj_ud,
               tensor::generic<CCTK_REAL, 4, 1> *H_d) {
  for (int a = 0; a < 4; ++a) {
    H_d->at(a) = 0.0;
    for (int b = 0; b < 4; ++b)
      for (int c = 0; c < 4; ++c) {
        H_d->at(a) -= proj_ud(b, a) * u_u(c) * rT_dd(b, c);
      }
  }
}

// Project out the radiation pressure tensor (in any frame)
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
calc_K_from_rT(tensor::symmetric2<CCTK_REAL, 4, 2> const &rT_dd,
               tensor::generic<CCTK_REAL, 4, 2> const &proj_ud,
               tensor::symmetric2<CCTK_REAL, 4, 2> *K_dd) {
  for (int a = 0; a < 4; ++a)
    for (int b = a; b < 4; ++b) {
      K_dd->at(a, b) = 0.0;
      for (int c = 0; c < 4; ++c)
        for (int d = 0; d < 4; ++d) {
          K_dd->at(a, b) += proj_ud(c, a) * proj_ud(d, b) * rT_dd(c, d);
        }
    }
}

// Compute the radiation energy flux
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
calc_E_flux(CCTK_REAL const alp, tensor::generic<CCTK_REAL, 4, 1> const &beta_u,
            CCTK_REAL const E, tensor::generic<CCTK_REAL, 4, 1> const &F_u,
            int const dir) {
  return alp * F_u(dir) - beta_u(dir) * E;
}

// Compute the flux of neutrino energy flux
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
calc_F_flux(CCTK_REAL const alp, tensor::generic<CCTK_REAL, 4, 1> const &beta_u,
            tensor::generic<CCTK_REAL, 4, 1> const &F_d,
            tensor::generic<CCTK_REAL, 4, 2> const &P_ud, int const dir,
            int const comp) {
  return alp * P_ud(dir, comp) - beta_u(dir) * F_d(comp);
}

// Computes the sources S_a = [eta - k_abs J] u_a - [k_abs + k_scat] H_a
// WARNING: be consistent with the densitization of eta, J, and H_d
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
calc_rad_sources(CCTK_REAL const eta, CCTK_REAL const kabs,
                 CCTK_REAL const kscat,
                 tensor::generic<CCTK_REAL, 4, 1> const &u_d, CCTK_REAL const J,
                 tensor::generic<CCTK_REAL, 4, 1> const H_d,
                 tensor::generic<CCTK_REAL, 4, 1> *S_d) {
  for (int a = 0; a < 4; ++a) {
    S_d->at(a) = (eta - kabs * J) * u_d(a) - (kabs + kscat) * H_d(a);
  }
}

// Computes the source term for E: -alp n^a S_a
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
calc_rE_source(CCTK_REAL const alp, tensor::generic<CCTK_REAL, 4, 1> const &n_u,
               tensor::generic<CCTK_REAL, 4, 1> const &S_d) {
  return -alp * tensor::dot(n_u, S_d);
}

// Computes the physical source term for the densitized number density.
// N must be the accepted (and, when necessary, separately repaired) stage
// state.  Numerical floors belong to the state projection and must not be
// introduced again in this collision source.
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
calc_rN_source(CCTK_REAL const alp, CCTK_REAL const volform,
               CCTK_REAL const eta_0, CCTK_REAL const kabs_0, CCTK_REAL const N,
               CCTK_REAL const Gamma) {
  return alp * (volform * eta_0 - kabs_0 * N / Gamma);
}

// Computes the source term for F_a: alp gamma^b_a S_b
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
calc_rF_source(CCTK_REAL const alp,
               tensor::generic<CCTK_REAL, 4, 2> const gamma_ud,
               tensor::generic<CCTK_REAL, 4, 1> const &S_d,
               tensor::generic<CCTK_REAL, 4, 1> *tS_d) {
  for (int a = 0; a < 4; ++a) {
    tS_d->at(a) = 0.0;
    for (int b = 0; b < 4; ++b) {
      tS_d->at(a) += alp * gamma_ud(b, a) * S_d(b);
    }
  }
}

#ifdef NUX_M1_CLOSURE_IMPLEMENTATION

// Function to rootfind in order to determine the closure
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline double
zFunction(double xi, void *params) {
  Parameters *p = reinterpret_cast<Parameters *>(params);

  tensor::symmetric2<CCTK_REAL, 4, 2> P_dd;
  apply_closure(p->g_dd, p->g_uu, p->n_d, p->w_lorentz, p->u_u, p->v_d,
                p->proj_ud, p->E, p->F_d, eval_closure(p->closure, xi), &P_dd);

  tensor::symmetric2<CCTK_REAL, 4, 2> rT_dd;
  assemble_rT(p->n_d, p->E, p->F_d, P_dd, &rT_dd);

  CCTK_REAL const J = calc_J_from_rT(rT_dd, p->u_u);

  tensor::generic<CCTK_REAL, 4, 1> H_d;
  calc_H_from_rT(rT_dd, p->u_u, p->proj_ud, &H_d);

  CCTK_REAL const H2 = tensor::dot(p->g_uu, H_d, H_d);
  return (J * xi) * (J * xi) - H2;
}

// Computes the closure in the lab frame with a rootfinding procedure
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
closure_abort_if_no_fallback(bool use_fallback) {
  if (!use_fallback)
    device_abort();
}

// Computes the closure in the lab frame with a rootfinding procedure
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_NOINLINE void calc_closure(
    cGH const *cctkGH, int const i, int const j, int const k, int const ig,
    closure_t closure_fun, tensor::metric<4> const &g_dd,
    tensor::inv_metric<4> const &g_uu,
    tensor::generic<CCTK_REAL, 4, 1> const &n_d, CCTK_REAL const w_lorentz,
    tensor::generic<CCTK_REAL, 4, 1> const &u_u,
    tensor::generic<CCTK_REAL, 4, 1> const &v_d,
    tensor::generic<CCTK_REAL, 4, 2> const &proj_ud, CCTK_REAL const E,
    tensor::generic<CCTK_REAL, 4, 1> const &F_d, CCTK_REAL *chi,
    tensor::symmetric2<CCTK_REAL, 4, 2> *P_dd, CCTK_REAL closure_epsilon,
    CCTK_INT closure_maxiter, bool use_fallback) {
  // These are special cases for which no root finding is needed
  if (closure_is_eddington(closure_fun)) {
    *chi = 1. / 3.;
    apply_closure(g_dd, g_uu, n_d, w_lorentz, u_u, v_d, proj_ud, E, F_d, *chi,
                  P_dd);
    return;
  }
  if (closure_is_thin(closure_fun)) {
    *chi = 1.0;
    apply_closure(g_dd, g_uu, n_d, w_lorentz, u_u, v_d, proj_ud, E, F_d, *chi,
                  P_dd);
    return;
  }

  Parameters params(closure_fun, g_dd, g_uu, n_d, w_lorentz, u_u, v_d, proj_ud,
                    E, F_d);
  auto fn = [&params](auto x) { return zFunction(x, &params); };
  auto fallback_chi = [&]() {
    // Score the exact endpoint tensors that the fallback can return.  For the
    // nonlinear closure families xi=0 gives chi=1/3 (Eddington), whereas xi=1
    // gives chi=1 (thin).  Evaluating fn(1/3) here would score an intermediate
    // nonlinear pressure tensor and then return a different Eddington tensor.
    CCTK_REAL const z_ed = fn(CCTK_REAL(0.0));
    CCTK_REAL const z_th = fn(CCTK_REAL(1.0));
    if (isfinite(z_th) && isfinite(z_ed)) {
      return (abs(z_th) < abs(z_ed)) ? CCTK_REAL(1.0) : CCTK_REAL(1.0 / 3.0);
    }
    if (isfinite(z_th)) {
      return CCTK_REAL(1.0);
    }
    return CCTK_REAL(1.0 / 3.0);
  };

  double x_lo = 0.0;
  double x_hi = 1.0;
  CCTK_REAL const f_lo = fn(x_lo);
  CCTK_REAL const f_hi = fn(x_hi);

  if (isfinite(f_lo) && isfinite(f_hi) &&
      (f_lo == CCTK_REAL(0.0) || f_hi == CCTK_REAL(0.0))) {
    // Match GSL/THC endpoint-root handling. If both endpoints solve the
    // degenerate equation, THC returns the upper endpoint for this test.
    CCTK_REAL const xi =
        (f_hi == CCTK_REAL(0.0)) ? CCTK_REAL(x_hi) : CCTK_REAL(x_lo);
    *chi = eval_closure(closure_fun, xi);
    apply_closure(g_dd, g_uu, n_d, w_lorentz, u_u, v_d, proj_ud, E, F_d, *chi,
                  P_dd);
    return;
  }

  // No root, most likely because of high velocities in the fluid
  // We use very simple approximation in this case
  if (!isfinite(f_lo) || !isfinite(f_hi) || f_lo * f_hi > 0.0) {
    closure_abort_if_no_fallback(use_fallback);
    *chi = fallback_chi();
    apply_closure(g_dd, g_uu, n_d, w_lorentz, u_u, v_d, proj_ud, E, F_d, *chi,
                  P_dd);
    return;
  }

  nuX_Utils::roots::brent_solver<CCTK_REAL> solver{};
  int ierr = nuX_Utils::roots::solver_set(&solver, fn, x_lo, x_hi);

  if (ierr != nuX_Utils::roots::code(nuX_Utils::roots::status::success)) {
    closure_abort_if_no_fallback(use_fallback);
    *chi = fallback_chi();
    apply_closure(g_dd, g_uu, n_d, w_lorentz, u_u, v_d, proj_ud, E, F_d, *chi,
                  P_dd);
    return;
  }

  int iter = 0;
  int test_status =
      nuX_Utils::roots::code(nuX_Utils::roots::status::continue_iter);
  bool solver_failed = false;

  do {
    ++iter;
    ierr = nuX_Utils::roots::solver_iterate(&solver, fn);
    if (ierr != nuX_Utils::roots::code(nuX_Utils::roots::status::success)) {
      solver_failed = true;
      break;
    }

    test_status = nuX_Utils::roots::test_interval(
        nuX_Utils::roots::solver_x_lower(&solver),
        nuX_Utils::roots::solver_x_upper(&solver), closure_epsilon,
        CCTK_REAL(0.0));

    if (test_status !=
            nuX_Utils::roots::code(nuX_Utils::roots::status::success) &&
        test_status !=
            nuX_Utils::roots::code(nuX_Utils::roots::status::continue_iter)) {
      solver_failed = true;
      break;
    }
  } while (test_status == nuX_Utils::roots::code(
                              nuX_Utils::roots::status::continue_iter) &&
           iter < closure_maxiter);

  if (solver_failed || test_status != nuX_Utils::roots::code(
                                          nuX_Utils::roots::status::success)) {
    closure_abort_if_no_fallback(use_fallback);
    *chi = fallback_chi();
    apply_closure(g_dd, g_uu, n_d, w_lorentz, u_u, v_d, proj_ud, E, F_d, *chi,
                  P_dd);
    return;
  }

  CCTK_REAL const xi = nuX_Utils::roots::solver_root(&solver);
  if (!isfinite(xi)) {
    closure_abort_if_no_fallback(use_fallback);
    *chi = fallback_chi();
  } else {
    CCTK_REAL const chi_try = eval_closure(closure_fun, xi);
    if (!isfinite(chi_try)) {
      closure_abort_if_no_fallback(use_fallback);
      *chi = fallback_chi();
    } else {
      *chi = chi_try;
    }
  }
  // We are done, update the closure with the newly found chi
  apply_closure(g_dd, g_uu, n_d, w_lorentz, u_u, v_d, proj_ud, E, F_d, *chi,
                P_dd);
}

#endif // NUX_M1_CLOSURE_IMPLEMENTATION

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_NOINLINE void calc_closure(
    cGH const *cctkGH, int const i, int const j, int const k, int const ig,
    closure_t closure_fun, tensor::metric<4> const &g_dd,
    tensor::inv_metric<4> const &g_uu,
    tensor::generic<CCTK_REAL, 4, 1> const &n_d, CCTK_REAL const w_lorentz,
    tensor::generic<CCTK_REAL, 4, 1> const &u_u,
    tensor::generic<CCTK_REAL, 4, 1> const &v_d,
    tensor::generic<CCTK_REAL, 4, 2> const &proj_ud, CCTK_REAL const E,
    tensor::generic<CCTK_REAL, 4, 1> const &F_d, CCTK_REAL *chi,
    tensor::symmetric2<CCTK_REAL, 4, 2> *P_dd, CCTK_REAL closure_epsilon,
    CCTK_INT closure_maxiter, bool use_fallback);

struct moment_repair_t {
  CCTK_REAL delta_E;
  tensor::generic<CCTK_REAL, 4, 1> delta_F_d;
  bool repaired_energy;
  bool repaired_flux;
};

// F_a is spatial with respect to the Eulerian normal, so F_a F^a is
// nonnegative.  Its four-dimensional contraction can nevertheless acquire a
// tiny negative value from cancellation between lapse/shift terms.  Scale the
// tolerance with the absolute contraction terms, not only with E^2, so a
// harmless roundoff error is not mistaken for an invalid radiation state.
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
spacelike_norm_roundoff_tolerance(
    tensor::symmetric2<CCTK_REAL, 4, 2> const &g_uu,
    tensor::generic<CCTK_REAL, 4, 1> const &F_d,
    CCTK_REAL const E = CCTK_REAL(0)) {
  CCTK_REAL contraction_scale = abs(E * E);
  for (int a = 0; a < 4; ++a)
    for (int b = 0; b < 4; ++b)
      contraction_scale += abs(g_uu(a, b) * F_d(a) * F_d(b));
  contraction_scale =
      max(contraction_scale, std::numeric_limits<CCTK_REAL>::min());
  return CCTK_REAL(256) * std::numeric_limits<CCTK_REAL>::epsilon() *
         contraction_scale;
}

// Check that a comoving radiation energy and flux belong to the physical
// moment cone. This is a validation only: independently modifying J or H_a
// would make them inconsistent with the stress tensor from which they were
// projected.
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline bool
comoving_state_is_realizable(tensor::symmetric2<CCTK_REAL, 4, 2> const &g_uu,
                             CCTK_REAL const J,
                             tensor::generic<CCTK_REAL, 4, 1> const &H_d) {
  if (!isfinite(J) || J < CCTK_REAL(0))
    return false;
  for (int a = 0; a < 4; ++a)
    if (!isfinite(H_d(a)))
      return false;

  const CCTK_REAL H2 = tensor::dot(g_uu, H_d, H_d);
  const CCTK_REAL tol = spacelike_norm_roundoff_tolerance(g_uu, H_d, J);
  return isfinite(H2) && H2 >= -tol && H2 <= J * J + tol;
}

// A fixed Eddington tensor is a diffusion approximation rather than a
// realizability-preserving M1 closure.  A lab-frame state can therefore have a
// finite, positive comoving energy while its projected flux lies outside the
// M1 moment cone.  Such a state is still usable by the explicitly requested
// Eddington model.  Variable closures, on the other hand, must preserve the
// complete comoving cone.  Keep this model distinction identical in the
// closure and collision-source paths.
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline bool
comoving_state_is_acceptable(
    closure_t const closure,
    tensor::symmetric2<CCTK_REAL, 4, 2> const &g_uu, CCTK_REAL const J,
    tensor::generic<CCTK_REAL, 4, 1> const &H_d) {
  if (!closure_is_eddington(closure))
    return comoving_state_is_realizable(g_uu, J, H_d);

  if (!isfinite(J) || J < CCTK_REAL(0))
    return false;
  for (int a = 0; a < 4; ++a)
    if (!isfinite(H_d(a)))
      return false;
  const CCTK_REAL H2 = tensor::dot(g_uu, H_d, H_d);
  const CCTK_REAL tol = spacelike_norm_roundoff_tolerance(g_uu, H_d, J);
  return isfinite(H2) && H2 >= -tol;
}

// Reconstruct a Lorentz factor from a metric and an interpolated velocity.
// High-order interpolation can overshoot the timelike domain; repair that
// velocity before constructing u^a rather than interpolating W independently.
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
repair_velocity_and_compute_W(tensor::metric<4> const &g_dd,
                              CCTK_REAL *const vx, CCTK_REAL *const vy,
                              CCTK_REAL *const vz,
                              CCTK_REAL const v2_margin = CCTK_REAL(1.0e-12)) {
  if (!isfinite(*vx) || !isfinite(*vy) || !isfinite(*vz)) {
    *vx = *vy = *vz = 0.0;
  }

  CCTK_REAL v2 =
      g_dd(1, 1) * (*vx) * (*vx) + g_dd(2, 2) * (*vy) * (*vy) +
      g_dd(3, 3) * (*vz) * (*vz) +
      2.0 * (g_dd(1, 2) * (*vx) * (*vy) + g_dd(1, 3) * (*vx) * (*vz) +
             g_dd(2, 3) * (*vy) * (*vz));
  const CCTK_REAL v2_limit = 1.0 - v2_margin;
  if (!isfinite(v2) || v2 < 0.0 || !(v2_limit > 0.0)) {
    *vx = *vy = *vz = 0.0;
    v2 = 0.0;
  } else if (v2 >= v2_limit) {
    const CCTK_REAL fac = sqrt(v2_limit / v2);
    *vx *= fac;
    *vy *= fac;
    *vz *= fac;
    // Recompute from the rounded repaired components.  Reusing the ideal
    // target v2_limit can make W inconsistent with the velocity when W is
    // large, even though the scaling itself is correct.
    v2 = g_dd(1, 1) * (*vx) * (*vx) + g_dd(2, 2) * (*vy) * (*vy) +
         g_dd(3, 3) * (*vz) * (*vz) +
         2.0 * (g_dd(1, 2) * (*vx) * (*vy) +
                g_dd(1, 3) * (*vx) * (*vz) +
                g_dd(2, 3) * (*vy) * (*vz));
    if (!isfinite(v2) || v2 < CCTK_REAL(0) || v2 >= CCTK_REAL(1)) {
      *vx = *vy = *vz = CCTK_REAL(0);
      v2 = CCTK_REAL(0);
    }
  }
  return CCTK_REAL(1) / sqrt(CCTK_REAL(1) - v2);
}

// Enforce E >= rad_E_floor and F_a F^a <= (1-rad_eps) E^2.
//
// This routine is shared by cell-centred source states and vertex-interpolated
// stress-energy states.  It records the numerical change separately from any
// physical collision increment.  Already-realizable finite states are left
// bitwise unchanged.
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
repair_moments(tensor::symmetric2<CCTK_REAL, 4, 2> const &g_uu, CCTK_REAL *E,
               tensor::generic<CCTK_REAL, 4, 1> *F_d,
               CCTK_REAL const rad_E_floor, CCTK_REAL const rad_eps,
               moment_repair_t *const repair = nullptr) {

  const CCTK_REAL Eold = *E;
  const tensor::generic<CCTK_REAL, 4, 1> Fold_d = *F_d;
  bool repaired_energy = false;
  bool repaired_flux = false;

  if (!isfinite(*E) || *E < rad_E_floor) {
    *E = rad_E_floor;
    repaired_energy = true;
  }

  bool flux_finite = true;
  for (int a = 0; a < 4; ++a)
    flux_finite = flux_finite && isfinite(F_d->at(a));

  CCTK_REAL F2 = flux_finite ? tensor::dot(g_uu, *F_d, *F_d) : -1.0;
  const CCTK_REAL norm_tol =
      flux_finite && isfinite(F2)
          ? spacelike_norm_roundoff_tolerance(g_uu, *F_d, *E)
          : CCTK_REAL(0);
  if (!flux_finite || !isfinite(F2) || F2 < -norm_tol) {
    for (int a = 0; a < 4; ++a)
      F_d->at(a) = 0.0;
    F2 = 0.0;
    repaired_flux = true;
  } else if (F2 < 0.0) {
    F2 = 0.0;
  }

  const CCTK_REAL lim = (*E) * (*E) * (1.0 - rad_eps);
  if (F2 > lim) {
    const CCTK_REAL fac = sqrt(lim / F2);
    for (int a = 0; a < 4; ++a)
      F_d->at(a) *= fac;
    repaired_flux = true;
  }

  if (repair != nullptr) {
    repair->delta_E = *E - Eold;
    for (int a = 0; a < 4; ++a)
      repair->delta_F_d(a) = F_d->at(a) - Fold_d(a);
    repair->repaired_energy = repaired_energy;
    repair->repaired_flux = repaired_flux;
  }
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
apply_floor(tensor::symmetric2<CCTK_REAL, 4, 2> const &g_uu, CCTK_REAL *E,
            tensor::generic<CCTK_REAL, 4, 1> *F_d, CCTK_REAL const rad_E_floor,
            CCTK_REAL const rad_eps) {
  repair_moments(g_uu, E, F_d, rad_E_floor, rad_eps);
}

} // namespace nuX_M1

#endif
