#ifndef NUX_M1_SOURCES_HXX
#define NUX_M1_SOURCES_HXX

#include <algorithm>
#include <cassert>
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

namespace nuX_M1 {

using namespace nuX_Utils;
using namespace Loop;
using namespace std;

struct Params {
  Params(cGH const *_cctkGH, int const _i, int const _j, int const _k,
         int const _ig, closure_t _closure, gsl_root_fsolver *_gsl_solver_1d,
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
         CCTK_REAL const _Estar,
         tensor::generic<CCTK_REAL, 4, 1> const &_Fstar_d, CCTK_REAL const _chi,
         CCTK_REAL const _eta, CCTK_REAL const _kabs, CCTK_REAL const _kscat)
      : cctkGH(_cctkGH), i(_i), j(_j), k(_k), ig(_ig),
        closure(_closure), // gsl_solver_1d(_gsl_solver_1d),
        cdt(_cdt), alp(_alp), g_dd(_g_dd), g_uu(_g_uu), n_d(_n_d), n_u(_n_u),
        gamma_ud(_gamma_ud), u_d(_u_d), u_u(_u_u), v_d(_v_d), v_u(_v_u),
        proj_ud(_proj_ud), W(_W), Estar(_Estar), Fstar_d(_Fstar_d), chi(_chi),
        eta(_eta), kabs(_kabs), kscat(_kscat) {}
  cGH const *cctkGH;
  int const i;
  int const j;
  int const k;
  int const ig;
  closure_t closure;
  gsl_root_fsolver *gsl_solver_1d;
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
  CCTK_REAL const Estar;
  tensor::generic<CCTK_REAL, 4, 1> const &Fstar_d;
  CCTK_REAL chi;
  CCTK_REAL const eta;
  CCTK_REAL const kabs;
  CCTK_REAL const kscat;

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
CCTK_HOST CCTK_DEVICE void
__source_jacobian_low_level(double *qpre, double Fup[4], double F2, double chi,
                            double kapa, double kaps, double vup[4],
                            double vdown[4], double v2, double W, double alpha,
                            double cdt, double *qstar),
//        gsl_matrix * J)
{
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
      // gsl_matrix_set(J, a, b, A_data[a][b]);
    }
}

CCTK_HOST CCTK_DEVICE int prepare_closure(gsl_vector const *q, Params *p) {
  p->E = max(gsl_vector_get(q, 0), 0.0);
  if (p->E < 0) {
    return GSL_EBADFUNC;
  }
  pack_F_d(-p->alp * p->n_u(1), -p->alp * p->n_u(2), -p->alp * p->n_u(3),
           gsl_vector_get(q, 1), gsl_vector_get(q, 2), gsl_vector_get(q, 3),
           &p->F_d);
  tensor::contract(p->g_uu, p->F_d, &p->F_u);

  calc_closure(p->cctkGH, p->i, p->j, p->k, p->ig, p->closure, p->gsl_solver_1d,
               p->g_dd, p->g_uu, p->n_d, p->W, p->u_u, p->v_d, p->proj_ud, p->E,
               p->F_d, &p->chi, &p->P_dd);

  return GSL_SUCCESS;
}

CCTK_HOST CCTK_DEVICE int prepare_sources(gsl_vector const *q, Params *p) {
  assemble_rT(p->n_d, p->E, p->F_d, p->P_dd, &p->T_dd);

  p->J = calc_J_from_rT(p->T_dd, p->u_u);
  calc_H_from_rT(p->T_dd, p->u_u, p->proj_ud, &p->H_d);

  calc_rad_sources(p->eta, p->kabs, p->kscat, p->u_d, p->J, p->H_d, &p->S_d);

  p->Edot = calc_rE_source(p->alp, p->n_u, p->S_d);
  calc_rF_source(p->alp, p->gamma_ud, p->S_d, &p->tS_d);

  return GSL_SUCCESS;
}

CCTK_HOST CCTK_DEVICE int prepare(gsl_vector const *q, Params *p) {
  int ierr = prepare_closure(q, p);
  if (ierr != GSL_SUCCESS) {
    return ierr;
  }

  ierr = prepare_sources(q, p);
  if (ierr != GSL_SUCCESS) {
    return ierr;
  }

  return GSL_SUCCESS;
}

// Function to rootfind for
//    f(q) = q - q^* - dt S[q]
CCTK_HOST CCTK_DEVICE int impl_func_val(gsl_vector const *q, void *params,
                                        gsl_vector *f) {
  Params *p = reinterpret_cast<Params *>(params);
  int ierr = prepare(q, p);
  if (ierr != GSL_SUCCESS) {
    return ierr;
  }

#define EVALUATE_ZFUNC                                                         \
  gsl_vector_set(f, 0, gsl_vector_get(q, 0) - p->Estar - p->cdt * p->Edot);    \
  gsl_vector_set(f, 1,                                                         \
                 gsl_vector_get(q, 1) - p->Fstar_d(1) - p->cdt * p->tS_d(1));  \
  gsl_vector_set(f, 2,                                                         \
                 gsl_vector_get(q, 2) - p->Fstar_d(2) - p->cdt * p->tS_d(2));  \
  gsl_vector_set(f, 3,                                                         \
                 gsl_vector_get(q, 3) - p->Fstar_d(3) - p->cdt * p->tS_d(3));

  EVALUATE_ZFUNC

  return GSL_SUCCESS;
}

// Jacobian of the implicit function
CCTK_HOST CCTK_DEVICE int impl_func_jac(gsl_vector const *q, void *params,
                                        gsl_matrix *J) {
  Params *p = reinterpret_cast<Params *>(params);
  int ierr = prepare(q, p);
  if (ierr != GSL_SUCCESS) {
    return ierr;
  }

#define EVALUATE_ZJAC                                                          \
  double m_q[] = {p->E, p->F_d(1), p->F_d(2), p->F_d(3)};                      \
  double m_Fup[] = {p->F_u(0), p->F_u(1), p->F_u(2), p->F_u(3)};               \
  double m_F2 = tensor::dot(p->F_u, p->F_d);                                   \
  double m_chi = p->chi;                                                       \
  double m_kscat = p->kscat;                                                   \
  double m_kabs = p->kabs;                                                     \
  double m_vup[] = {p->v_u(0), p->v_u(1), p->v_u(2), p->v_u(3)};               \
  double m_vdw[] = {p->v_d(0), p->v_d(1), p->v_d(2), p->v_d(3)};               \
  double m_v2 = tensor::dot(p->v_u, p->v_d);                                   \
  double m_W = p->W;                                                           \
  double m_alpha = p->alp;                                                     \
  double m_cdt = p->cdt;                                                       \
  double m_qstar[] = {p->Estar, p->Fstar_d(1), p->Fstar_d(2), p->Fstar_d(3)};  \
                                                                               \
  __source_jacobian_low_level(m_q, m_Fup, m_F2, m_chi, m_kscat, m_kabs, m_vup, \
                              m_vdw, m_v2, m_W, m_alpha, m_cdt, m_qstar, J);

  EVALUATE_ZJAC

  return GSL_SUCCESS;
}

// Function and Jacobian evaluation
CCTK_HOST CCTK_DEVICE int impl_func_val_jac(gsl_vector const *q, void *params,
                                            gsl_vector *f, gsl_matrix *J) {
  Params *p = reinterpret_cast<Params *>(params);
  int ierr = prepare(q, p);
  if (ierr != 0) {
    return ierr;
  }

  EVALUATE_ZFUNC
  EVALUATE_ZJAC

  return GSL_SUCCESS;
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
explicit_update(Params *p, CCTK_REAL *Enew,
                tensor::generic<CCTK_REAL, 4, 1> *Fnew_d) {
  *Enew = p->Estar + p->cdt * p->Edot;
  Fnew_d->at(1) = p->Fstar_d(1) + p->cdt * p->tS_d(1);
  Fnew_d->at(2) = p->Fstar_d(2) + p->cdt * p->tS_d(2);
  Fnew_d->at(3) = p->Fstar_d(3) + p->cdt * p->tS_d(3);
  // F_0 = g_0i F^i = beta_i F^i = beta^i F_i
  Fnew_d->at(0) = -p->alp * p->n_u(1) * Fnew_d->at(1) -
                  p->alp * p->n_u(2) * Fnew_d->at(2) -
                  p->alp * p->n_u(3) * Fnew_d->at(3);
}

// Solves the implicit problem
// .  q^new = q^star + dt S[q^new]
// The source term is S^a = (eta - ka J) u^a - (ka + ks) H^a and includes
// also emission.

CCTK_HOST CCTK_DEVICE inline int source_update(
    cGH const *cctkGH, int const i, int const j, int const k, int const ig,
    closure_t closure_fun,
    //        gsl_root_fsolver * gsl_solver_1d,
    //        gsl_multiroot_fdfsolver * gsl_solver_nd,
    CCTK_REAL const cdt, CCTK_REAL const alp, tensor::metric<4> const &g_dd,
    tensor::inv_metric<4> const &g_uu,
    tensor::generic<CCTK_REAL, 4, 1> const &n_d,
    tensor::generic<CCTK_REAL, 4, 1> const &n_u,
    tensor::generic<CCTK_REAL, 4, 2> const &gamma_ud,
    tensor::generic<CCTK_REAL, 4, 1> const &u_d,
    tensor::generic<CCTK_REAL, 4, 1> const &u_u,
    tensor::generic<CCTK_REAL, 4, 1> const &v_d,
    tensor::generic<CCTK_REAL, 4, 1> const &v_u,
    tensor::generic<CCTK_REAL, 4, 2> const &proj_ud, CCTK_REAL const W,
    CCTK_REAL const Eold, tensor::generic<CCTK_REAL, 4, 1> const &Fold_d,
    CCTK_REAL const Estar, tensor::generic<CCTK_REAL, 4, 1> const &Fstar_d,
    CCTK_REAL const eta, CCTK_REAL const kabs, CCTK_REAL const kscat,
    CCTK_REAL *chi, CCTK_REAL *Enew, tensor::generic<CCTK_REAL, 4, 1> *Fnew_d) {

  DECLARE_CCTK_PARAMETERS;
  Params p(cctkGH, i, j, k, ig, closure_fun, cdt, alp, g_dd, g_uu, n_d, n_u,
           gamma_ud, u_d, u_u, v_d, v_u, proj_ud, W, Estar, Fstar_d, *chi, eta,
           kabs, kscat);

  //    gsl_multiroot_function_fdf zfunc = {
  //        impl_func_val,
  //        impl_func_jac,
  //        impl_func_val_jac,
  //        4, &p};

  // Old solution
  CCTK_REAL qold[] = {Eold, Fold_d(1), Fold_d(2), Fold_d(3)};
  //    gsl_vector_view xold = gsl_vector_view_array(qold, 4);

  // Non stiff limit, use explicit update
  if (cdt * kabs < 1 && cdt * kscat < 1) {
    prepare(&xold.vector, &p);
    explicit_update(&p, Enew, Fnew_d);

    CCTK_REAL q[4] = {*Enew, Fnew_d->at(1), Fnew_d->at(2), Fnew_d->at(3)};
    //        gsl_vector_view x = gsl_vector_view_array(q, 4);
    prepare_closure(&x.vector, &p);
    *chi = p.chi;

    return NUX_M1_SOURCE_THIN;
  }

  // Our scheme cannot capture this dynamics (tau << dt), so we go
  // directly to the equilibrium
  if (source_thick_limit > 0 &&
      SQ(cdt) * (kabs * (kabs + kscat)) > SQ(source_thick_limit)) {
    return NUX_M1_SOURCE_EQUIL;
  }

  // This handles the scattering dominated limit
  if (source_scat_limit > 0 && cdt * kscat > source_scat_limit) {
    return NUX_M1_SOURCE_SCAT;
  }

  // Initial guess for the solution
  CCTK_REAL q[4] = {*Enew, Fnew_d->at(1), Fnew_d->at(2), Fnew_d->at(3)};
  //    gsl_vector_view x = gsl_vector_view_array(q, 4);

  //    int ierr = gsl_multiroot_fdfsolver_set(gsl_solver_nd, &zfunc,
  //    &x.vector);
  int iter = 0;
  do {
    if (iter < source_maxiter) {
      //           ierr = gsl_multiroot_fdfsolver_iterate(gsl_solver_nd);
      iter++;
    }
    // The nonlinear solver is stuck.
    if (ierr == GSL_ENOPROG || ierr == GSL_ENOPROGJ || ierr == GSL_EBADFUNC ||
        iter >= source_maxiter) {

      // If we are here, then we are in trouble
#ifdef WARN_FOR_SRC_FIX
      ostringstream ss;
      if (ierr == GSL_EBADFUNC) {
        ss << "NaNs or Infs found in the implicit solve.\n";
      } else if (iter > source_maxiter) {
        ss << "Source solver exceeded the maximum number of iterations\n";
      } else {
        ss << "Stuck nonlinear solver.\n";
      }
      ss << "Trying to save the day... ";
#endif

      // We are optically thick, suggest to retry with Eddington closure
      if (closure_fun != eddington) {
#ifdef WARN_FOR_SRC_FIX
        ss << "Eddington closure\n";
//                print_stuff(cctkGH, i, j, k, ig, &p, ss);
//                Printer::print_warn(ss.str());
#endif
        ierr = source_update(
            cctkGH, i, j, k, ig,
            //                    eddington, gsl_solver_1d, gsl_solver_nd, cdt,
            eddington, cdt, alp, g_dd, g_uu, n_d, n_u, gamma_ud, u_d, u_u, v_d,
            v_u, proj_ud, W, Eold, Fold_d, Estar, Fstar_d, eta, kabs, kscat,
            chi, Enew, Fnew_d);
        if (ierr == NUX_M1_SOURCE_OK) {
          return NUX_M1_SOURCE_EDDINGTON;
        } else {
          return ierr;
        }
      } else {
#ifdef WARN_FOR_SRC_FIX
        ss << "using initial guess\n";
//                print_stuff(cctkGH, i, j, k, ig, &p, ss);
//                Printer::print_warn(ss.str());
#endif
        return NUX_M1_SOURCE_FAIL;
      }
    } else if (ierr != GSL_SUCCESS) {
      char msg[BUFSIZ];
      snprintf(msg, BUFSIZ,
               "Unexpected error in "
               "gsl_multirootroot_fdfsolver_iterate, error code \"%d\"",
               ierr);
#pragma omp critical
      CCTK_ERROR(msg);
    }
    //        ierr = gsl_multiroot_test_delta(gsl_solver_nd->dx,
    //        gsl_solver_nd->x,
    //                source_epsabs, source_epsrel);
  } while (ierr == GSL_CONTINUE);

  //    *Enew = gsl_vector_get(gsl_solver_nd->x, 0);
  //    Fnew_d->at(1) = gsl_vector_get(gsl_solver_nd->x, 1);
  //    Fnew_d->at(2) = gsl_vector_get(gsl_solver_nd->x, 2);
  //    Fnew_d->at(3) = gsl_vector_get(gsl_solver_nd->x, 3);

  //    // F_0 = g_0i F^i = beta_i F^i = beta^i F_i
  Fnew_d->at(0) = -alp * n_u(1) * Fnew_d->at(1) - alp * n_u(2) * Fnew_d->at(2) -
                  alp * n_u(3) * Fnew_d->at(3);

  //    prepare_closure(gsl_solver_nd->x, &p);
  *chi = p.chi;

  return NUX_M1_SOURCE_OK;
}

} // namespace nuX_M1

#endif
