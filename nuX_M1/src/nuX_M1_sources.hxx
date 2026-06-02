#ifndef NUX_M1_SOURCES_HXX
#define NUX_M1_SOURCES_HXX

#include <algorithm>
#include <cassert>
#include <cmath>
#include <loop_device.hxx>

#include "nuX_M1_closure.hxx"
#include "nuX_M1_macro.hxx"
#include "nuX_utils.hxx"

#include <cctk_Parameters.h>

#define NUX_M1_SOURCE_OK 0
#define NUX_M1_SOURCE_THIN 1
#define NUX_M1_SOURCE_EDDINGTON 2
#define NUX_M1_SOURCE_EQUIL 3
#define NUX_M1_SOURCE_SCAT 4
#define NUX_M1_SOURCE_FAIL 5

enum : int {
  ROOTS_SUCCESS = static_cast<int>(nuX_Utils::roots::status::success),
  ROOTS_CONTINUE = static_cast<int>(nuX_Utils::roots::status::continue_iter),
  ROOTS_EBADFUNC = static_cast<int>(nuX_Utils::roots::status::ebadfunc),
  ROOTS_ENOPROG = static_cast<int>(nuX_Utils::roots::status::enoprog),
  ROOTS_ENOPROGJ = static_cast<int>(nuX_Utils::roots::status::enoprogj),
};

// Fixed-size vector/matrix aliases for the implicit 4D root solve.
typedef Arith::vec<CCTK_REAL, 4> arith_vector;
typedef Arith::mat<CCTK_REAL, 4> arith_matrix;

namespace nuX_M1 {

using namespace nuX_Utils;
using namespace Loop;
using namespace std;

struct SourceUpdateContext {
  CCTK_HOST CCTK_DEVICE SourceUpdateContext(
      cGH const *_cctkGH, int const _i, int const _j, int const _k,
      int const _ig, CCTK_REAL _closure_epsilon,
      CCTK_INT _closure_maxiter, bool _closure_use_fallback,
      CCTK_REAL const _cdt, CCTK_REAL const _alp,
      tensor::metric<4> const &_g_dd, tensor::inv_metric<4> const &_g_uu,
      tensor::generic<CCTK_REAL, 4, 1> const &_n_d,
      tensor::generic<CCTK_REAL, 4, 1> const &_n_u,
      tensor::generic<CCTK_REAL, 4, 2> const &_gamma_ud,
      tensor::generic<CCTK_REAL, 4, 1> const &_u_d,
      tensor::generic<CCTK_REAL, 4, 1> const &_u_u,
      tensor::generic<CCTK_REAL, 4, 1> const &_v_d,
      tensor::generic<CCTK_REAL, 4, 1> const &_v_u,
      tensor::generic<CCTK_REAL, 4, 2> const &_proj_ud, CCTK_REAL const _W,
      CCTK_REAL const _Eold, tensor::generic<CCTK_REAL, 4, 1> const &_Fold_d,
      CCTK_REAL const _Estar,
      tensor::generic<CCTK_REAL, 4, 1> const &_Fstar_d,
      CCTK_REAL const _eta, CCTK_REAL const _kabs, CCTK_REAL const _kscat)
      : cctkGH(_cctkGH), i(_i), j(_j), k(_k), ig(_ig),
        closure_epsilon(_closure_epsilon), closure_maxiter(_closure_maxiter),
        closure_use_fallback(_closure_use_fallback), cdt(_cdt), alp(_alp),
        g_dd(_g_dd), g_uu(_g_uu), n_d(_n_d), n_u(_n_u), gamma_ud(_gamma_ud),
        u_d(_u_d), u_u(_u_u), v_d(_v_d), v_u(_v_u), proj_ud(_proj_ud), W(_W),
        Eold(_Eold), Fold_d(_Fold_d), Estar(_Estar), Fstar_d(_Fstar_d),
        eta(_eta), kabs(_kabs),
        kscat(_kscat) {}
  cGH const *cctkGH;
  int const i;
  int const j;
  int const k;
  int const ig;
  CCTK_REAL closure_epsilon;
  CCTK_INT closure_maxiter;
  bool closure_use_fallback;
  CCTK_REAL const cdt;
  CCTK_REAL const alp;
  tensor::metric<4> const &g_dd;
  tensor::inv_metric<4> const &g_uu;
  tensor::generic<CCTK_REAL, 4, 1> const &n_d;
  tensor::generic<CCTK_REAL, 4, 1> const &n_u;
  tensor::generic<CCTK_REAL, 4, 2> const &gamma_ud;
  tensor::generic<CCTK_REAL, 4, 1> const &u_d;
  tensor::generic<CCTK_REAL, 4, 1> const &u_u;
  tensor::generic<CCTK_REAL, 4, 1> const &v_d;
  tensor::generic<CCTK_REAL, 4, 1> const &v_u;
  tensor::generic<CCTK_REAL, 4, 2> const &proj_ud;
  CCTK_REAL const W;
  CCTK_REAL const Eold;
  tensor::generic<CCTK_REAL, 4, 1> const &Fold_d;
  CCTK_REAL const Estar;
  tensor::generic<CCTK_REAL, 4, 1> const &Fstar_d;
  CCTK_REAL const eta;
  CCTK_REAL const kabs;
  CCTK_REAL const kscat;
};

#ifdef NUX_M1_SOURCES_IMPLEMENTATION

struct Params {
  CCTK_HOST CCTK_DEVICE Params(SourceUpdateContext const &ctx,
                               closure_t const _closure, CCTK_REAL const _chi)
      : ctx(&ctx), closure(_closure), chi(_chi) {}
  SourceUpdateContext const *ctx;
  closure_t closure;
  CCTK_REAL chi;
};

struct PreparedState {
  CCTK_REAL E;
  tensor::generic<CCTK_REAL, 4, 1> F_d;
  tensor::generic<CCTK_REAL, 4, 1> F_u;
  tensor::symmetric2<CCTK_REAL, 4, 2> P_dd;
  tensor::symmetric2<CCTK_REAL, 4, 2> T_dd;
  CCTK_REAL J;
  tensor::generic<CCTK_REAL, 4, 1> H_d;
  tensor::generic<CCTK_REAL, 4, 1> S_d;
  CCTK_REAL Edot;
  tensor::generic<CCTK_REAL, 4, 1> tS_d;
};

/* return dthin */
CCTK_HOST CCTK_DEVICE inline double radM1_set_dthin(double chi) {
  return 1.5 * chi - 0.5;
}

/* return dthick */
CCTK_HOST CCTK_DEVICE inline double radM1_set_dthick(double chi) {
  return 1.5 * (1 - chi);
}

/*
double sign(double x) {
  return (0 < x) - (x < 0);
}
*/

// Low level kernel computing the Jacobian matrix
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
__source_jacobian_low_level(double *qpre, double Fup[4], double F2, double chi,
                            double kapa, double kaps, double vup[4],
                            double vdown[4], double v2, double W, double alpha,
                            double cdt, arith_matrix &J) {
  const double kapas = kapa + kaps;
  const double alpW = alpha * W;

  const double dthin = radM1_set_dthin(chi);
  const double dthick = radM1_set_dthick(chi);

  const double vx = vdown[1];
  const double vy = vdown[2];
  const double vz = vdown[3];
  const double W2 = SQ(W);
  const double W3 = W2 * W;

  const double vdotF =
      Fup[1] * vdown[1] + Fup[2] * vdown[2] + Fup[3] * vdown[3];
  const double normF = sqrt(F2);
  const double inormF = (normF > 0 ? 1 / normF : 0);
  const double vdothatf = vdotF * inormF;
  const double vdothatf2 = SQ(vdothatf);
  const double hatfx = qpre[1] * inormF; // hatf_i
  const double hatfy = qpre[2] * inormF;
  const double hatfz = qpre[3] * inormF;
  const double hatfupx = Fup[1] * inormF; // hatf^i
  const double hatfupy = Fup[2] * inormF;
  const double hatfupz = Fup[3] * inormF;
  const double e = qpre[0];
  const double eonormF = min(e * inormF, 1.0); // with factor dthin ...

  // drvts of J
  double JdE = W2 + dthin * vdothatf2 * W2 +
               (dthick * (3 - 2 * W2) * (-1 + W2)) / (1 + 2 * W2);

  double JdFv = 2 * W2 *
                (-1 + (dthin * eonormF * vdothatf) +
                 (2 * dthick * (-1 + W2)) / (1 + 2 * W2));
  double JdFf = (-2 * dthin * eonormF * vdothatf2 * W2);

  double JdFx = JdFv * vup[1] + JdFf * hatfupx;
  double JdFy = JdFv * vup[2] + JdFf * hatfupy;
  double JdFz = JdFv * vup[3] + JdFf * hatfupz;

  // drvts of Hi
  double HdEv =
      W3 * (-1 - dthin * vdothatf2 + (dthick * (-3 + 2 * W2)) / (1 + 2 * W2));
  double HdEf = -(dthin * vdothatf * W);

  double HxdE = HdEv * vx + HdEf * hatfx;
  double HydE = HdEv * vy + HdEf * hatfy;
  double HzdE = HdEv * vz + HdEf * hatfz;

  double HdFdelta = (1 - dthick * v2 - (dthin * eonormF * vdothatf)) * W;
  double HdFvv = (2 * (1 - dthin * eonormF * vdothatf) * W3) +
                 dthick * W * (2 - 2 * W2 + 1 / (-1 - 2 * W2));
  double HdFff = (2 * dthin * eonormF * vdothatf * W);
  double HdFvf = (2 * dthin * eonormF * vdothatf2 * W3);
  double HdFfv = -(dthin * eonormF * W);

  double HxdFx = HdFdelta + HdFvv * vx * vup[1] + HdFff * hatfx * hatfupx +
                 HdFvf * vx * hatfupx + HdFfv * hatfx * vup[1];
  double HydFx = HdFvv * vy * vup[1] + HdFff * hatfy * hatfupx +
                 HdFvf * vy * hatfupx + HdFfv * hatfy * vup[1];
  double HzdFx = HdFvv * vz * vup[1] + HdFff * hatfz * hatfupx +
                 HdFvf * vz * hatfupx + HdFfv * hatfz * vup[1];

  double HxdFy = HdFvv * vx * vup[2] + HdFff * hatfx * hatfupy +
                 HdFvf * vx * hatfupy + HdFfv * hatfx * vup[2];
  double HydFy = HdFdelta + HdFvv * vy * vup[2] + HdFff * hatfy * hatfupy +
                 HdFvf * vy * hatfupy + HdFfv * hatfy * vup[2];
  double HzdFy = HdFvv * vz * vup[2] + HdFff * hatfz * hatfupy +
                 HdFvf * vz * hatfupy + HdFfv * hatfz * vup[2];

  double HxdFz = HdFvv * vx * vup[3] + HdFff * hatfx * hatfupz +
                 HdFvf * vx * hatfupz + HdFfv * hatfx * vup[3];
  double HydFz = HdFvv * vy * vup[3] + HdFff * hatfy * hatfupz +
                 HdFvf * vy * hatfupz + HdFfv * hatfy * vup[3];
  double HzdFz = HdFdelta + HdFvv * vz * vup[3] + HdFff * hatfz * hatfupz +
                 HdFvf * vz * hatfupz + HdFfv * hatfz * vup[3];

  // Build the Jacobian
  double J00 = -alpW * (kapas - kaps * JdE);

  double J0x = +alpW * kaps * JdFx + alpW * kapas * vup[1];
  double J0y = +alpW * kaps * JdFy + alpW * kapas * vup[2];
  double J0z = +alpW * kaps * JdFz + alpW * kapas * vup[3];

  double Jx0 = -alpha * (kapas * HxdE + W * kapa * vx * JdE);
  double Jy0 = -alpha * (kapas * HydE + W * kapa * vy * JdE);
  double Jz0 = -alpha * (kapas * HzdE + W * kapa * vz * JdE);

  double Jxx = -alpha * (kapas * HxdFx + W * kapa * vx * JdFx);
  double Jxy = -alpha * (kapas * HxdFy + W * kapa * vx * JdFy);
  double Jxz = -alpha * (kapas * HxdFz + W * kapa * vx * JdFz);

  double Jyy = -alpha * (kapas * HydFx + W * kapa * vy * JdFx);
  double Jyx = -alpha * (kapas * HydFy + W * kapa * vy * JdFy);
  double Jyz = -alpha * (kapas * HydFz + W * kapa * vy * JdFz);

  double Jzx = -alpha * (kapas * HzdFx + W * kapa * vz * JdFx);
  double Jzy = -alpha * (kapas * HzdFy + W * kapa * vz * JdFy);
  double Jzz = -alpha * (kapas * HzdFz + W * kapa * vz * JdFz);

  // Store Jacobian into J
  double A_data[4][4] = {
      1 - cdt * J00, -cdt * J0x,    -cdt * J0y,    -cdt * J0z,
      -cdt * Jx0,    1 - cdt * Jxx, -cdt * Jxy,    -cdt * Jxz,
      -cdt * Jy0,    -cdt * Jyx,    1 - cdt * Jyy, -cdt * Jyz,
      -cdt * Jz0,    -cdt * Jzx,    -cdt * Jzy,    1 - cdt * Jzz,
  };
  for (int a = 0; a < 4; ++a)
    for (int b = 0; b < 4; ++b) {
      J(a, b) = A_data[a][b];
    }
}

CCTK_HOST CCTK_DEVICE int prepare_closure(const arith_vector &q, Params *p,
                                          PreparedState *s) {
  SourceUpdateContext const &c = *p->ctx;
  s->E = max(q(0), 0.0);
  if (s->E < 0) {
    return ROOTS_EBADFUNC;
  }
  pack_F_d(-c.alp * c.n_u(1), -c.alp * c.n_u(2), -c.alp * c.n_u(3), q(1),
           q(2), q(3), &s->F_d);
  tensor::contract(c.g_uu, s->F_d, &s->F_u);

  calc_closure(c.cctkGH, c.i, c.j, c.k, c.ig, p->closure, c.g_dd, c.g_uu,
               c.n_d, c.W, c.u_u, c.v_d, c.proj_ud, s->E, s->F_d, &p->chi,
               &s->P_dd, c.closure_epsilon, c.closure_maxiter,
               c.closure_use_fallback);

  return ROOTS_SUCCESS;
}

CCTK_HOST CCTK_DEVICE int prepare_sources(Params *p, PreparedState *s) {
  SourceUpdateContext const &c = *p->ctx;
  assemble_rT(c.n_d, s->E, s->F_d, s->P_dd, &s->T_dd);

  s->J = calc_J_from_rT(s->T_dd, c.u_u);
  calc_H_from_rT(s->T_dd, c.u_u, c.proj_ud, &s->H_d);

  calc_rad_sources(c.eta, c.kabs, c.kscat, c.u_d, s->J, s->H_d, &s->S_d);

  s->Edot = calc_rE_source(c.alp, c.n_u, s->S_d);
  calc_rF_source(c.alp, c.gamma_ud, s->S_d, &s->tS_d);

  return ROOTS_SUCCESS;
}

CCTK_HOST CCTK_DEVICE int prepare(const arith_vector &q, Params *p,
                                  PreparedState *s) {
  int ierr = prepare_closure(q, p, s);
  if (ierr != ROOTS_SUCCESS) {
    return ierr;
  }

  ierr = prepare_sources(p, s);
  if (ierr != ROOTS_SUCCESS) {
    return ierr;
  }

  return ROOTS_SUCCESS;
}

// Function to rootfind for
//    f(q) = q - q^* - dt S[q]
CCTK_HOST CCTK_DEVICE int impl_func_val(const arith_vector &q, Params *p,
                                        arith_vector &f) {
  PreparedState s;
  int ierr = prepare(q, p, &s);
  if (ierr != ROOTS_SUCCESS) {
    return ierr;
  }
  SourceUpdateContext const &c = *p->ctx;

#define EVALUATE_ZFUNC                                                         \
  f(0) = q(0) - c.Estar - c.cdt * s.Edot;                                      \
  f(1) = q(1) - c.Fstar_d(1) - c.cdt * s.tS_d(1);                              \
  f(2) = q(2) - c.Fstar_d(2) - c.cdt * s.tS_d(2);                              \
  f(3) = q(3) - c.Fstar_d(3) - c.cdt * s.tS_d(3);

  EVALUATE_ZFUNC

  return ROOTS_SUCCESS;
}

// Jacobian of the implicit function
CCTK_HOST CCTK_DEVICE int impl_func_jac(const arith_vector &q, Params *p,
                                        arith_matrix &J) {
  PreparedState s;
  int ierr = prepare(q, p, &s);
  if (ierr != ROOTS_SUCCESS) {
    return ierr;
  }
  SourceUpdateContext const &c = *p->ctx;

#define EVALUATE_ZJAC                                                          \
  double m_q[] = {s.E, s.F_d(1), s.F_d(2), s.F_d(3)};                          \
  double m_Fup[] = {s.F_u(0), s.F_u(1), s.F_u(2), s.F_u(3)};                   \
  double m_F2 = tensor::dot(s.F_u, s.F_d);                                     \
  double m_chi = p->chi;                                                       \
  double m_kscat = c.kscat;                                                    \
  double m_kabs = c.kabs;                                                      \
  double m_vup[] = {c.v_u(0), c.v_u(1), c.v_u(2), c.v_u(3)};                   \
  double m_vdw[] = {c.v_d(0), c.v_d(1), c.v_d(2), c.v_d(3)};                   \
  double m_v2 = tensor::dot(c.v_u, c.v_d);                                     \
  double m_W = c.W;                                                            \
  double m_alpha = c.alp;                                                      \
  double m_cdt = c.cdt;                                                        \
                                                                               \
  __source_jacobian_low_level(m_q, m_Fup, m_F2, m_chi, m_kabs, m_kscat, m_vup, \
                              m_vdw, m_v2, m_W, m_alpha, m_cdt, J);

  EVALUATE_ZJAC

  return ROOTS_SUCCESS;
}

// Function and Jacobian evaluation
CCTK_HOST CCTK_DEVICE int impl_func_val_jac(const arith_vector &q, Params *p,
                                            arith_vector &f, arith_matrix &J) {
  PreparedState s;
  int ierr = prepare(q, p, &s);
  if (ierr != ROOTS_SUCCESS) {
    return ierr;
  }
  SourceUpdateContext const &c = *p->ctx;

  EVALUATE_ZFUNC
  EVALUATE_ZJAC

  return ROOTS_SUCCESS;
}

#undef EVALUATE_ZFUNC
#undef EVALUATE_ZJAC

#if 0
void scattering_limit(
        Params * p,
        CCTK_REAL J,
        CCTK_REAL * Enew,
        tensor::generic<CCTK_REAL, 4, 1> * Fnew_d) {
    tensor::symmetric2<CCTK_REAL, 4, 2> T_dd;
    for (int a = 0; a < 4; ++a)
    for (int b = a; b < 4; ++b) {
        T_dd(a,b) = (4./3.) * J * p->u_d(a) * p->u_d(b) +
                    (1./3.) * J * p->g_dd(a,b);
    }
    *Enew = calc_J_from_rT(T_dd, p->n_u);
    calc_H_from_rT(T_dd, p->n_u, p->gamma_ud, Fnew_d);
}

void thermal_equilibrium(
        Params * p,
        CCTK_REAL * Enew,
        tensor::generic<CCTK_REAL, 4, 1> * Fnew_d) {
    scattering_limit(p, (p->kabs > 0 ? p->eta/p->kabs : 0), Enew, Fnew_d);
}
#endif

CCTK_HOST CCTK_DEVICE void
explicit_update(Params *p, PreparedState const &s, CCTK_REAL *Enew,
                tensor::generic<CCTK_REAL, 4, 1> *Fnew_d) {
  SourceUpdateContext const &c = *p->ctx;
  *Enew = c.Estar + c.cdt * s.Edot;
  Fnew_d->at(1) = c.Fstar_d(1) + c.cdt * s.tS_d(1);
  Fnew_d->at(2) = c.Fstar_d(2) + c.cdt * s.tS_d(2);
  Fnew_d->at(3) = c.Fstar_d(3) + c.cdt * s.tS_d(3);
  // F_0 = g_0i F^i = beta_i F^i = beta^i F_i
  Fnew_d->at(0) = -c.alp * c.n_u(1) * Fnew_d->at(1) -
                  c.alp * c.n_u(2) * Fnew_d->at(2) -
                  c.alp * c.n_u(3) * Fnew_d->at(3);
}

#endif // NUX_M1_SOURCES_IMPLEMENTATION

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline bool
source_rates_are_zero(CCTK_REAL const eta_0, CCTK_REAL const eta_1,
                      CCTK_REAL const kabs_0, CCTK_REAL const kabs_1,
                      CCTK_REAL const kscat_1) {
  return eta_0 == CCTK_REAL(0) && eta_1 == CCTK_REAL(0) &&
         kabs_0 == CCTK_REAL(0) && kabs_1 == CCTK_REAL(0) &&
         kscat_1 == CCTK_REAL(0);
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline bool
source_is_nonstiff(CCTK_REAL const cdt, CCTK_REAL const kabs,
                   CCTK_REAL const kscat) {
  return cdt * kabs < CCTK_REAL(1) && cdt * kscat < CCTK_REAL(1);
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline bool
source_uses_thick_limit(CCTK_REAL const cdt, CCTK_REAL const kabs,
                        CCTK_REAL const kscat,
                        CCTK_REAL const source_thick_limit) {
  return source_thick_limit > CCTK_REAL(0) &&
         SQ(cdt) * (kabs * (kabs + kscat)) > SQ(source_thick_limit);
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline bool
source_uses_scat_limit(CCTK_REAL const cdt, CCTK_REAL const kscat,
                       CCTK_REAL const source_scat_limit) {
  return source_scat_limit > CCTK_REAL(0) &&
         cdt * kscat > source_scat_limit;
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline bool
source_can_use_no_source_fallback(
    CCTK_REAL const cdt, CCTK_REAL const kabs, CCTK_REAL const kscat,
    CCTK_REAL const Estar, tensor::generic<CCTK_REAL, 4, 1> const &Fstar_d) {
  constexpr CCTK_REAL source_fail_tau_safe = CCTK_REAL(1.0e-3);
  const CCTK_REAL tau_tot = cdt * (kabs + kscat);
  const bool star_finite = isfinite(Estar) && isfinite(Fstar_d(1)) &&
                           isfinite(Fstar_d(2)) && isfinite(Fstar_d(3));
  return isfinite(tau_tot) && tau_tot <= source_fail_tau_safe && star_finite;
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
source_apply_no_source_fallback(
    CCTK_REAL const Estar, tensor::generic<CCTK_REAL, 4, 1> const &Fstar_d,
    CCTK_REAL *chi, CCTK_REAL *Enew,
    tensor::generic<CCTK_REAL, 4, 1> *Fnew_d) {
  *Enew = Estar;
  *Fnew_d = Fstar_d;
  *chi = CCTK_REAL(1.0 / 3.0);
}

CCTK_HOST CCTK_DEVICE int source_update_nonstiff(
    SourceUpdateContext const &ctx, closure_t closure_fun, CCTK_REAL *chi,
    CCTK_REAL *Enew,
    tensor::generic<CCTK_REAL, 4, 1> *Fnew_d);

CCTK_HOST CCTK_DEVICE int source_update(
    SourceUpdateContext const &ctx, closure_t closure_fun, CCTK_REAL *chi,
    CCTK_REAL *Enew, tensor::generic<CCTK_REAL, 4, 1> *Fnew_d,
    CCTK_REAL source_thick_limit, CCTK_REAL source_scat_limit,
    CCTK_INT source_maxiter, CCTK_REAL source_epsabs, CCTK_REAL source_epsrel);

#ifdef NUX_M1_SOURCES_IMPLEMENTATION

CCTK_HOST CCTK_DEVICE int source_update_nonstiff_attempt(
    SourceUpdateContext const &ctx, closure_t closure_fun, CCTK_REAL *chi,
    CCTK_REAL *Enew,
    tensor::generic<CCTK_REAL, 4, 1> *Fnew_d) {

  Params p(ctx, closure_fun, *chi);

  arith_vector qold{ctx.Eold, ctx.Fold_d(1), ctx.Fold_d(2), ctx.Fold_d(3)};

  PreparedState state;
  prepare(qold, &p, &state);
  explicit_update(&p, state, Enew, Fnew_d);

  const bool state_finite =
      isfinite(*Enew) && isfinite(Fnew_d->at(0)) && isfinite(Fnew_d->at(1)) &&
      isfinite(Fnew_d->at(2)) && isfinite(Fnew_d->at(3));
  if (!state_finite) {
    return NUX_M1_SOURCE_FAIL;
  }

  arith_vector q{*Enew, Fnew_d->at(1), Fnew_d->at(2), Fnew_d->at(3)};
  PreparedState closure_state;
  int ierr_cl = prepare_closure(q, &p, &closure_state);
  if (ierr_cl != ROOTS_SUCCESS || !isfinite(p.chi)) {
    return NUX_M1_SOURCE_FAIL;
  }
  *chi = p.chi;

  return NUX_M1_SOURCE_THIN;
}

CCTK_HOST CCTK_DEVICE int source_update_nonstiff(
    SourceUpdateContext const &ctx, closure_t closure_fun, CCTK_REAL *chi,
    CCTK_REAL *Enew,
    tensor::generic<CCTK_REAL, 4, 1> *Fnew_d) {

  int ierr =
      source_update_nonstiff_attempt(ctx, closure_fun, chi, Enew, Fnew_d);
  if (ierr != NUX_M1_SOURCE_FAIL) {
    return ierr;
  }

  if (!closure_is_eddington(closure_fun)) {
    ierr = source_update_nonstiff_attempt(ctx, CLOSURE_EDDINGTON, chi, Enew,
                                          Fnew_d);
    if (ierr != NUX_M1_SOURCE_FAIL) {
      return (ierr == NUX_M1_SOURCE_OK) ? NUX_M1_SOURCE_EDDINGTON : ierr;
    }
  }

  if (source_can_use_no_source_fallback(ctx.cdt, ctx.kabs, ctx.kscat, ctx.Estar,
                                        ctx.Fstar_d)) {
    source_apply_no_source_fallback(ctx.Estar, ctx.Fstar_d, chi, Enew, Fnew_d);
    return NUX_M1_SOURCE_EDDINGTON;
  }

  return NUX_M1_SOURCE_FAIL;
}

CCTK_HOST CCTK_DEVICE int source_update_implicit_attempt(
    SourceUpdateContext const &ctx, closure_t closure_fun, CCTK_REAL *chi,
    CCTK_REAL *Enew, tensor::generic<CCTK_REAL, 4, 1> *Fnew_d,
    CCTK_INT source_maxiter, CCTK_REAL source_epsabs,
    CCTK_REAL source_epsrel) {

  Params p(ctx, closure_fun, *chi);

  // Initial guess for the solution
  arith_vector q_initial_guess{*Enew, Fnew_d->at(1), Fnew_d->at(2),
                               Fnew_d->at(3)};

  auto fn_nd_val = [&p](arith_vector q, arith_vector &f) -> int {
    return impl_func_val(q, &p, f);
  };
  auto fn_nd_jac = [&p](arith_vector q, arith_matrix &J) -> int {
    return impl_func_jac(q, &p, J);
  };
  auto fn_nd_val_jac =
      [&p](arith_vector q, arith_vector &f, arith_matrix &J) -> int {
    return impl_func_val_jac(q, &p, f, J);
  };

  nuX_Utils::roots::hybridsj_solver<CCTK_REAL, 4> solver{};
  int ierr =
      nuX_Utils::roots::hybridsj_set(&solver, fn_nd_val_jac, q_initial_guess);
  int iter = 0;
  bool failed = false;

  if (ierr != ROOTS_SUCCESS) {
    failed = true;
  } else {
    do {
      if (iter < source_maxiter) {
        ierr = nuX_Utils::roots::hybridsj_iterate(&solver, fn_nd_val, fn_nd_jac);
        iter++;
      }

      if (ierr == ROOTS_ENOPROG || ierr == ROOTS_ENOPROGJ ||
          ierr == ROOTS_EBADFUNC || iter >= source_maxiter) {
        failed = true;
        break;
      } else if (ierr != ROOTS_SUCCESS) {
        assert(false && "Unexpected error in hybridsj_iterate");
        failed = true;
        break;
      }

      ierr = nuX_Utils::roots::test_delta(solver.dx, solver.x, source_epsabs,
                                          source_epsrel);
    } while (ierr == ROOTS_CONTINUE);

    if (!failed && ierr != ROOTS_SUCCESS) {
      assert(false && "Unexpected error in roots::test_delta");
      failed = true;
    }
  }

  if (failed) {
    return NUX_M1_SOURCE_FAIL;
  }

  // Assign the output of the multiroot solve to the updated state.
  *Enew = solver.x(0);
  Fnew_d->at(1) = solver.x(1);
  Fnew_d->at(2) = solver.x(2);
  Fnew_d->at(3) = solver.x(3);

  // F_0 = g_0i F^i = beta_i F^i = beta^i F_i
  Fnew_d->at(0) = -ctx.alp * ctx.n_u(1) * Fnew_d->at(1) -
                  ctx.alp * ctx.n_u(2) * Fnew_d->at(2) -
                  ctx.alp * ctx.n_u(3) * Fnew_d->at(3);

  const bool state_finite = isfinite(*Enew) && isfinite(Fnew_d->at(0)) &&
                            isfinite(Fnew_d->at(1)) &&
                            isfinite(Fnew_d->at(2)) && isfinite(Fnew_d->at(3));
  if (!state_finite) {
    return NUX_M1_SOURCE_FAIL;
  }

  PreparedState closure_state;
  int ierr_cl = prepare_closure(solver.x, &p, &closure_state);
  if (ierr_cl != ROOTS_SUCCESS || !isfinite(p.chi)) {
    return NUX_M1_SOURCE_FAIL;
  }
  *chi = p.chi;

  return NUX_M1_SOURCE_OK;
}

// Solves the implicit problem
// .  q^new = q^star + dt S[q^new]
// The source term is S^a = (eta - ka J) u^a - (ka + ks) H^a and includes
// also emission.

CCTK_HOST CCTK_DEVICE int source_update(
    SourceUpdateContext const &ctx, closure_t closure_fun, CCTK_REAL *chi,
    CCTK_REAL *Enew, tensor::generic<CCTK_REAL, 4, 1> *Fnew_d,
    CCTK_REAL source_thick_limit, CCTK_REAL source_scat_limit,
    CCTK_INT source_maxiter, CCTK_REAL source_epsabs, CCTK_REAL source_epsrel) {

  // Non stiff limit, use explicit update
  if (source_is_nonstiff(ctx.cdt, ctx.kabs, ctx.kscat)) {
    return source_update_nonstiff(ctx, closure_fun, chi, Enew, Fnew_d);
  }

  // Our scheme cannot capture this dynamics (tau << dt), so we go
  // directly to the equilibrium
  if (source_uses_thick_limit(ctx.cdt, ctx.kabs, ctx.kscat,
                              source_thick_limit)) {
    return NUX_M1_SOURCE_EQUIL;
  }

  // This handles the scattering dominated limit
  if (source_uses_scat_limit(ctx.cdt, ctx.kscat, source_scat_limit)) {
    return NUX_M1_SOURCE_SCAT;
  }

  int ierr = source_update_implicit_attempt(
      ctx, closure_fun, chi, Enew, Fnew_d, source_maxiter, source_epsabs,
      source_epsrel);
  if (ierr != NUX_M1_SOURCE_FAIL) {
    return ierr;
  }

#ifdef WARN_FOR_SRC_FIX
  printf("hybridsj failed in the implicit solve!\n");
#endif

  if (!closure_is_eddington(closure_fun)) {
#ifdef WARN_FOR_SRC_FIX
    printf("Eddington closure\n");
#endif
    ierr = source_update_implicit_attempt(ctx, CLOSURE_EDDINGTON, chi, Enew,
                                          Fnew_d, source_maxiter,
                                          source_epsabs, source_epsrel);
    if (ierr == NUX_M1_SOURCE_OK) {
      return NUX_M1_SOURCE_EDDINGTON;
    }
  }

#ifdef WARN_FOR_SRC_FIX
  printf("using initial guess\n");
#endif
  if (source_can_use_no_source_fallback(ctx.cdt, ctx.kabs, ctx.kscat,
                                        ctx.Estar, ctx.Fstar_d)) {
    source_apply_no_source_fallback(ctx.Estar, ctx.Fstar_d, chi, Enew, Fnew_d);
    return NUX_M1_SOURCE_EDDINGTON;
  }

  return NUX_M1_SOURCE_FAIL;
}

#endif // NUX_M1_SOURCES_IMPLEMENTATION

} // namespace nuX_M1

#endif
