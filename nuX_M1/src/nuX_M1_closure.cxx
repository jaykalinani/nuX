//  Templated Hydrodynamics Code: an hydro code built on top of HRSCCore
//  Copyright (C) 2020, David Radice <david.radice@psu.edu>
//
//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.


#include <cassert>
#include <iostream>
#include <sstream>
#include <string>

// TODO: Dependence on GSL

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>

#include "thc_printer.hh"

#include "thc_M1_closure.hh"
#include "thc_M1_macro.h"

#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

using namespace std;

// TODO: Rename namespaces

namespace thc {
namespace m1 {

namespace {

// TODO: tensor:: objects

struct Parameters {
    Parameters(
        closure_t _closure,
        tensor::metric<4> const & _g_dd,
        tensor::inv_metric<4> const & _g_uu,
        tensor::generic<CCTK_REAL, 4, 1> const & _n_d,
        CCTK_REAL const _w_lorentz,
        tensor::generic<CCTK_REAL, 4, 1> const & _u_u,
        tensor::generic<CCTK_REAL, 4, 1> const & _v_d,
        tensor::generic<CCTK_REAL, 4, 2> const & _proj_ud,
        CCTK_REAL const _E,
        tensor::generic<CCTK_REAL, 4, 1> const & _F_d):
            closure(_closure), g_dd(_g_dd), g_uu(_g_uu), n_d(_n_d),
            w_lorentz(_w_lorentz), u_u(_u_u), v_d(_v_d), proj_ud(_proj_ud),
            E(_E), F_d(_F_d) {}
    closure_t closure;
    tensor::metric<4> const & g_dd;
    tensor::inv_metric<4> const & g_uu;
    tensor::generic<CCTK_REAL, 4, 1> const & n_d;
    CCTK_REAL const w_lorentz;
    tensor::generic<CCTK_REAL, 4, 1> const & u_u;
    tensor::generic<CCTK_REAL, 4, 1> const & v_d;
    tensor::generic<CCTK_REAL, 4, 2> const & proj_ud;
    CCTK_REAL const E;
    tensor::generic<CCTK_REAL, 4, 1> const & F_d;
};

// TODO: *cctkGH pointer might not work
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void print_stuff(
        cGH const * cctkGH,
        const PointDesc &p,
        int const ig,
        Parameters const * p,
        ostream & ss) {
    DECLARE_CCTK_ARGUMENTS;
    // TODO: layout2 definition
    int const ijk = layout2.linear(p.i, p.j, p.k);
    // TODO: might not work
    int const i4D = CCTK_VectGFIndex3D(cctkGH, i, j, k, ig);

    ss << "Iteration = " << cctkGH->cctk_iteration << endl;
    ss << "Reflevel = " << ilogb(cctkGH->cctk_levfac[0]) << endl;
    ss << "(i, j, k) = (" << i << ", " << j << ", " << k << ")\n";
    ss << "ig = " << ig << endl;
    ss << "(x, y, z) = (" << x[ijk] << ", " << y[ijk] << ", " << z[ijk] << ")\n";
    ss << "rho = " << rho[ijk] << endl;
    ss << "temperature = " << temperature[ijk] << endl;
    ss << "Y_e = " << Y_e[ijk] << endl;
    ss << "alp = " << alp[ijk] << endl;
    ss << "abs_0 = " << abs_0[i4D] << endl;
    ss << "abs_1 = " << abs_1[i4D] << endl;
    ss << "eta_0 = " << eta_0[i4D] << endl;
    ss << "eta_1 = " << eta_1[i4D] << endl;
    ss << "scat_1 = " << scat_1[i4D] << endl;
    ss << "alp = " << alp[ijk] << endl;
    ss << "beta = (" << betax[ijk] << ", "
                     << betay[ijk] << ", "
                     << betaz[ijk] << ")\n";
    ss << "g_uu = (";
    for (int a = 0; a < 4; ++a)
    for (int b = 0; b < 4; ++b) {
        ss << p->g_uu(a,b) << ", ";
    }
    ss << "\b\b)\n";
    ss << "g_dd = (";
    for (int a = 0; a < 4; ++a)
    for (int b = 0; b < 4; ++b) {
        ss << p->g_dd(a,b) << ", ";
    }
    ss << "\b\b)\n";
    ss << "w_lorentz = " << p->w_lorentz << endl;
    ss << "n_d = (";
    for (int a = 0; a < 4; ++a) {
        ss << p->n_d(a) << ", ";
    }
    ss << "\b\b)\n";
    ss << "u_u = (";
    for (int a = 0; a < 4; ++a) {
        ss << p->u_u(a) << ", ";
    }
    ss << "\b\b)\n";
    ss << "v_d = (";
    for (int a = 0; a < 4; ++a) {
        ss << p->v_d(a) << ", ";
    }
    ss << "\b\b)\n";
    ss << "E = " << p->E << endl;
    ss << "F_d = (";
    for (int a = 0; a < 4; ++a) {
        ss << p->F_d(a) << ", ";
    }
    ss << "\b\b)\n";
}

// Function to rootfind in order to determine the closure
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline double zFunction(double xi, void * params) {
    Parameters *p = reinterpret_cast<Parameters *>(params);

    tensor::symmetric2<CCTK_REAL, 4, 2> P_dd;
    apply_closure(p->g_dd, p->g_uu, p->n_d, p->w_lorentz, p->u_u, p->v_d,
            p->proj_ud, p->E, p->F_d, p->closure(xi), &P_dd);

    tensor::symmetric2<CCTK_REAL, 4, 2> rT_dd;
    assemble_rT(p->n_d, p->E, p->F_d, P_dd, &rT_dd);

    CCTK_REAL const J = calc_J_from_rT(rT_dd, p->u_u);

    tensor::generic<CCTK_REAL, 4, 1> H_d;
    calc_H_from_rT(rT_dd, p->u_u, p->proj_ud, &H_d);

    CCTK_REAL const H2 = tensor::dot(p->g_uu, H_d, H_d);
    return SQ(J*xi) - H2;
}

} // anonymous namespace

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void pack_F_d(
        CCTK_REAL const betax, CCTK_REAL const betay, CCTK_REAL const betaz,
        CCTK_REAL const Fx, CCTK_REAL const Fy, CCTK_REAL const Fz,
        tensor::generic<CCTK_REAL, 4, 1> * F_d) {
    // F_0 = g_0i F^i = beta_i F^i = beta^i F_i
    F_d->at(0) = betax*Fx + betay*Fy + betaz*Fz;
    F_d->at(1) = Fx;
    F_d->at(2) = Fy;
    F_d->at(3) = Fz;
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void unpack_F_d(
        tensor::generic<CCTK_REAL, 4, 1> const & F_d,
        CCTK_REAL * Fx, CCTK_REAL * Fy, CCTK_REAL * Fz) {
    *Fx = F_d(1);
    *Fy = F_d(2);
    *Fz = F_d(3);
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void pack_F_d(CCTK_REAL const Fx, CCTK_REAL const Fy, CCTK_REAL const Fz,
        tensor::generic<CCTK_REAL, 3, 1> * F_d) {
    F_d->at(0) = Fx;
    F_d->at(1) = Fy;
    F_d->at(2) = Fz;
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void unpack_F_d(
        tensor::generic<CCTK_REAL, 3, 1> const & F_d,
        CCTK_REAL * Fx, CCTK_REAL * Fy, CCTK_REAL * Fz) {
    *Fx = F_d(0);
    *Fy = F_d(1);
    *Fz = F_d(2);
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void pack_H_d(
        CCTK_REAL const Ht, CCTK_REAL const Hx, CCTK_REAL const Hy, CCTK_REAL const Hz,
        tensor::generic<CCTK_REAL, 4, 1> * H_d) {
    H_d->at(0) = Ht;
    H_d->at(1) = Hx;
    H_d->at(2) = Hy;
    H_d->at(3) = Hz;
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void unpack_H_d(
        tensor::generic<CCTK_REAL, 4, 1> const & H_d,
        CCTK_REAL * Ht, CCTK_REAL * Hx, CCTK_REAL * Hy, CCTK_REAL * Hz) {
    *Ht = H_d(0);
    *Hx = H_d(1);
    *Hy = H_d(2);
    *Hz = H_d(3);
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void pack_P_dd(
        CCTK_REAL const betax, CCTK_REAL const betay, CCTK_REAL const betaz,
        CCTK_REAL const Pxx, CCTK_REAL const Pxy, CCTK_REAL const Pxz,
        CCTK_REAL const Pyy, CCTK_REAL const Pyz, CCTK_REAL const Pzz,
        tensor::symmetric2<CCTK_REAL, 4, 2> * P_dd) {
    CCTK_REAL const Pbetax = Pxx*betax + Pxy*betay + Pxz*betaz;
    CCTK_REAL const Pbetay = Pxy*betax + Pyy*betay + Pyz*betaz;
    CCTK_REAL const Pbetaz = Pxz*betax + Pyz*betay + Pzz*betaz;

    // P_00 = g_0i g_k0 P^ik = beta^i beta^k P_ik
    P_dd->at(0,0) = Pbetax*betax + Pbetay*betay + Pbetaz*betaz;

    // P_0i = g_0j g_ki P^jk = beta_j P_i^j = beta^j P_ij
    P_dd->at(0,1) = Pbetax;
    P_dd->at(0,2) = Pbetay;
    P_dd->at(0,3) = Pbetaz;

    P_dd->at(1,1) = Pxx;
    P_dd->at(1,2) = Pxy;
    P_dd->at(1,3) = Pxz;
    P_dd->at(2,2) = Pyy;
    P_dd->at(2,3) = Pyz;
    P_dd->at(3,3) = Pzz;
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void unpack_P_dd(
        tensor::symmetric2<CCTK_REAL, 4, 2> const & P_dd,
        CCTK_REAL * Pxx, CCTK_REAL * Pxy, CCTK_REAL * Pxz,
        CCTK_REAL * Pyy, CCTK_REAL * Pyz, CCTK_REAL * Pzz) {
    *Pxx = P_dd(1,1);
    *Pxy = P_dd(1,2);
    *Pxz = P_dd(1,3);
    *Pyy = P_dd(2,2);
    *Pyz = P_dd(2,3);
    *Pzz = P_dd(3,3);
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void pack_P_dd(
        CCTK_REAL const Pxx, CCTK_REAL const Pxy, CCTK_REAL const Pxz,
        CCTK_REAL const Pyy, CCTK_REAL const Pyz, CCTK_REAL const Pzz,
        tensor::symmetric2<CCTK_REAL, 3, 2> * P_dd) {
    P_dd->at(0,0) = Pxx;
    P_dd->at(0,1) = Pxy;
    P_dd->at(0,2) = Pxz;
    P_dd->at(1,1) = Pyy;
    P_dd->at(1,2) = Pyz;
    P_dd->at(2,2) = Pzz;
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void unpack_P_dd(
        tensor::symmetric2<CCTK_REAL, 3, 2> const & P_dd,
        CCTK_REAL * Pxx, CCTK_REAL * Pxy, CCTK_REAL * Pxz,
        CCTK_REAL * Pyy, CCTK_REAL * Pyz, CCTK_REAL * Pzz) {
    *Pxx = P_dd(0,0);
    *Pxy = P_dd(0,1);
    *Pxz = P_dd(0,2);
    *Pyy = P_dd(1,1);
    *Pyz = P_dd(1,2);
    *Pzz = P_dd(2,2);
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void pack_v_u(
        CCTK_REAL const velx, CCTK_REAL const vely, CCTK_REAL const velz,
        tensor::generic<CCTK_REAL, 4, 1> * v_u) {
    v_u->at(0) = 0.0;
    v_u->at(1) = velx;
    v_u->at(2) = vely;
    v_u->at(3) = velz;
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void calc_proj(
        tensor::generic<CCTK_REAL, 4, 1> const & u_d,
        tensor::generic<CCTK_REAL, 4, 1> const & u_u,
        tensor::generic<CCTK_REAL, 4, 2> * proj_ud) {
    for (int a = 0; a < 4; ++a)
    for (int b = 0; b < 4; ++b) {
        proj_ud->at(a,b) = tensor::delta(a,b) + u_u(a)*u_d(b);
    }
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void calc_Pthin(
        tensor::inv_metric<4> const & g_uu,
        CCTK_REAL const E,
        tensor::generic<CCTK_REAL, 4, 1> const & F_d,
        tensor::symmetric2<CCTK_REAL, 4, 2> * P_dd) {
    CCTK_REAL const F2 = tensor::dot(g_uu, F_d, F_d);
    CCTK_REAL fac = (F2 > 0 ? E/F2 : 0);
    for (int a = 0; a < 4; ++a)
    for (int b = a; b < 4; ++b) {
        P_dd->at(a,b) = fac * F_d(a) * F_d(b);
    }
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void calc_Pthick(
        tensor::metric<4> const & g_dd,
        tensor::inv_metric<4> const & g_uu,
        tensor::generic<CCTK_REAL, 4, 1> const & n_d,
        CCTK_REAL const W,
        tensor::generic<CCTK_REAL, 4, 1> const & v_d,
        CCTK_REAL const E,
        tensor::generic<CCTK_REAL, 4, 1> const & F_d,
        tensor::symmetric2<CCTK_REAL, 4, 2> * P_dd) {
    CCTK_REAL const v_dot_F = tensor::dot(g_uu, v_d, F_d);

    CCTK_REAL const W2 = W*W;
    CCTK_REAL const coef = 1./(2.*W2 + 1.);

    // J/3
    CCTK_REAL const Jo3 = coef*((2.*W2 - 1.)*E - 2.*W2*v_dot_F);

    // tH = gamma_ud H_d
    tensor::generic<CCTK_REAL, 4, 1> tH_d;
    for (int a = 0; a < 4; ++a) {
        tH_d(a) = F_d(a)/W + coef*W*v_d(a)*((4.*W2 + 1.)*v_dot_F - 4.*W2*E);
    }

    for (int a = 0; a < 4; ++a)
    for (int b = a; b < 4; ++b) {
        P_dd->at(a,b)  = Jo3*(4.*W2*v_d(a)*v_d(b) + g_dd(a,b) + n_d(a)*n_d(b));
        P_dd->at(a,b) += W*(tH_d(a)*v_d(b) + tH_d(b)*v_d(a));
    }
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL flux_factor(
        tensor::inv_metric<4> const & g_uu,
        CCTK_REAL const J,
        tensor::generic<CCTK_REAL, 4, 1> const & H_d) {
    DECLARE_CCTK_PARAMETERS;
    CCTK_REAL xi = (J > rad_E_floor ? tensor::dot(g_uu, H_d, H_d)/SQ(J) : 0);
    return max(0.0, min(xi, 1.0));
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL eddington(CCTK_REAL const xi) {
    return 1.0/3.0;
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL kershaw(CCTK_REAL const xi) {
    return 1.0/3.0 + 2.0/3.0*xi*xi;
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL minerbo(CCTK_REAL const xi) {
    return 1.0/3.0 + xi*xi*(6.0 - 2.0*xi + 6.0*xi*xi)/15.0;
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL thin(CCTK_REAL const xi) {
    return 1.0;
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void apply_closure(
        tensor::metric<4> const & g_dd,
        tensor::inv_metric<4> const & g_uu,
        tensor::generic<CCTK_REAL, 4, 1> const & n_d,
        CCTK_REAL const w_lorentz,
        tensor::generic<CCTK_REAL, 4, 1> const & u_u,
        tensor::generic<CCTK_REAL, 4, 1> const & v_d,
        tensor::generic<CCTK_REAL, 4, 2> const & proj_ud,
        CCTK_REAL const E,
        tensor::generic<CCTK_REAL, 4, 1> const & F_d,
        CCTK_REAL const chi,
        tensor::symmetric2<CCTK_REAL, 4, 2> * P_dd) {
    tensor::symmetric2<CCTK_REAL, 4, 2> Pthin_dd;
    tensor::symmetric2<CCTK_REAL, 4, 2> Pthick_dd;

    calc_Pthin(g_uu, E, F_d, &Pthin_dd);
    calc_Pthick(g_dd, g_uu, n_d, w_lorentz, v_d, E, F_d, &Pthick_dd);

    CCTK_REAL const dthick = 3.*(1 - chi)/2.;
    CCTK_REAL const dthin = 1. - dthick;

    for (int a = 0; a < 4; ++a)
    for (int b = a; b < 4; ++b) {
        P_dd->at(a,b) = dthick*Pthick_dd(a,b) + dthin*Pthin_dd(a,b);
    }
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void calc_closure(
        cGH const * cctkGH,
        int const i, int const j, int const k,
        int const ig,
        closure_t closure_fun,
        gsl_root_fsolver * fsolver,
        tensor::metric<4> const & g_dd,
        tensor::inv_metric<4> const & g_uu,
        tensor::generic<CCTK_REAL, 4, 1> const & n_d,
        CCTK_REAL const w_lorentz,
        tensor::generic<CCTK_REAL, 4, 1> const & u_u,
        tensor::generic<CCTK_REAL, 4, 1> const & v_d,
        tensor::generic<CCTK_REAL, 4, 2> const & proj_ud,
        CCTK_REAL const E,
        tensor::generic<CCTK_REAL, 4, 1> const & F_d,
        CCTK_REAL * chi,
        tensor::symmetric2<CCTK_REAL, 4, 2> * P_dd) {
    DECLARE_CCTK_PARAMETERS;

    // These are special cases for which no root finding is needed
    if (closure_fun == eddington) {
        *chi = 1./3.;
        apply_closure(g_dd, g_uu, n_d, w_lorentz, u_u, v_d, proj_ud,
                E, F_d, *chi, P_dd);
        return;
    }
    if (closure_fun == thin) {
        *chi = 1.0;
        apply_closure(g_dd, g_uu, n_d, w_lorentz, u_u, v_d, proj_ud,
                E, F_d, *chi, P_dd);
        return;
    }

    Parameters params(closure_fun, g_dd, g_uu, n_d,
            w_lorentz, u_u, v_d, proj_ud, E, F_d);
    gsl_function F;
    F.function = &zFunction;
    F.params = reinterpret_cast<void *>(&params);

    double x_lo = 0.0;
    double x_hi = 1.0;

    int ierr = gsl_root_fsolver_set(fsolver, &F, x_lo, x_hi);
    // No root, most likely because of high velocities in the fluid
    // We use very simple approximation in this case
    if (ierr == GSL_EINVAL) {
        double const z_ed = zFunction(1./3., F.params);
        double const z_th = zFunction(1., F.params);
        if (abs(z_th) < abs(z_ed)) {
            *chi = 1.0;
        }
        else {
            *chi = 1.0/3.0;
        }
        apply_closure(g_dd, g_uu, n_d, w_lorentz, u_u, v_d, proj_ud,
                E, F_d, *chi, P_dd);
        return;
    }
    else if (ierr == GSL_EBADFUNC) {
        ostringstream ss;
        ss << "NaN or Inf found in closure!\n";
        print_stuff(cctkGH, i, j, k, ig, &params, ss);
        Printer::print_err(ss.str());
        CCTK_ERROR("NaN or Inf found in closure!");
    }
    else if (ierr != 0) {
        ostringstream ss;
        ss << "Unexpected error in gsl_root_fsolver_set, error code \""
           << ierr << "\"\n";
        print_stuff(cctkGH, i, j, k, ig, &params, ss);
        Printer::print_err(ss.str());
        CCTK_ERROR("Unexpected error in gsl_root_fsolver_set");
    }

    // Rootfinding
    int iter = 0;
    do {
        ++iter;
        ierr = gsl_root_fsolver_iterate(fsolver);
        // Some nans in the evaluation. This should not happen.
        if (ierr == GSL_EBADFUNC) {
            ostringstream ss;
            ss << "NaNs or Infs found when computing the closure!" << endl;
            print_stuff(cctkGH, i, j, k, ig, &params, ss);
            Printer::print_err(ss.str());
            CCTK_ERROR("NaNs or Infs found when computing the closure!");
        }
        else if (ierr != 0) {
            ostringstream ss;
            ss << "Unexpected error in gsl_root_fsolver_iterate,  error code \""
               << ierr << "\"\n";
            print_stuff(cctkGH, i, j, k, ig, &params, ss);
            Printer::print_err(ss.str());
            CCTK_ERROR("Unexpected error in gsl_root_fsolver_iterate");
        }
        *chi = closure_fun(gsl_root_fsolver_root(fsolver));
        x_lo = gsl_root_fsolver_x_lower(fsolver);
        x_hi = gsl_root_fsolver_x_upper(fsolver);
        ierr = gsl_root_test_interval(x_lo, x_hi, closure_epsilon, 0);
    } while (ierr == GSL_CONTINUE && iter < closure_maxiter);

#ifdef WARN_FOR_MAXITER
    if (ierr != GSL_SUCCESS) {
        ostringstream ss;
        ss << "Maximum number of iterations exceeded "
                "when computing the M1 closure";
        print_stuff(cctkGH, i, j, k, ig, &params, ss);
        Printer::print_warn(ss.str());
    }
#endif

    // We are done, update the closure with the newly found chi
    apply_closure(g_dd, g_uu, n_d, w_lorentz, u_u, v_d, proj_ud,
            E, F_d, *chi, P_dd);
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void assemble_fnu(
        tensor::generic<CCTK_REAL, 4, 1> const & u_u,
        CCTK_REAL const J,
        tensor::generic<CCTK_REAL, 4, 1> const & H_u,
        tensor::generic<CCTK_REAL, 4, 1> * fnu_u) {
    DECLARE_CCTK_PARAMETERS;
    for (int a = 0; a < 4; ++a) {
        fnu_u->at(a) = u_u(a) + (J > rad_E_floor ? H_u(a)/J : 0);
    }
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL compute_Gamma(
        CCTK_REAL const W,
        tensor::generic<CCTK_REAL, 4, 1> const & v_u,
        CCTK_REAL const J,
        CCTK_REAL const E,
        tensor::generic<CCTK_REAL, 4, 1> const & F_d) {
    DECLARE_CCTK_PARAMETERS;
    if (E > rad_E_floor && J > rad_E_floor) {
        CCTK_REAL f_dot_v = min(tensor::dot(F_d, v_u)/E, 1 - rad_eps);
        return W*(E/J)*(1 - f_dot_v);
    }
    else {
        return 1;
    }
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void assemble_rT(
        tensor::generic<CCTK_REAL, 4, 1> const & u_d,
        CCTK_REAL const J,
        tensor::generic<CCTK_REAL, 4, 1> const & H_d,
        tensor::symmetric2<CCTK_REAL, 4, 2> const & K_dd,
        tensor::symmetric2<CCTK_REAL, 4, 2> * rT_dd) {
    for (int a = 0; a < 4; ++a)
    for (int b = a; b < 4; ++b) {
        rT_dd->at(a,b) = J*u_d(a)*u_d(b) + H_d(a)*u_d(b) + H_d(b)*u_d(a) + K_dd(a,b);
    }
}

// Project out the radiation energy (in any frame)
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL calc_J_from_rT(
        tensor::symmetric2<CCTK_REAL, 4, 2> const & rT_dd,
        tensor::generic<CCTK_REAL, 4, 1> const & u_u) {
    return tensor::dot(rT_dd, u_u, u_u);
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void calc_H_from_rT(
        tensor::symmetric2<CCTK_REAL, 4, 2> const & rT_dd,
        tensor::generic<CCTK_REAL, 4, 1> const & u_u,
        tensor::generic<CCTK_REAL, 4, 2> const & proj_ud,
        tensor::generic<CCTK_REAL, 4, 1> * H_d) {
    for (int a = 0; a < 4; ++a) {
        H_d->at(a) = 0.0;
        for (int b = 0; b < 4; ++b)
        for (int c = 0; c < 4; ++c) {
            H_d->at(a) -= proj_ud(b,a)*u_u(c)*rT_dd(b,c);
        }
    }
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void calc_K_from_rT(
        tensor::symmetric2<CCTK_REAL, 4, 2> const & rT_dd,
        tensor::generic<CCTK_REAL, 4, 2> const & proj_ud,
        tensor::symmetric2<CCTK_REAL, 4, 2> * K_dd) {
    for (int a = 0; a < 4; ++a)
    for (int b = a; b < 4; ++b) {
        K_dd->at(a,b) = 0.0;
        for (int c = 0; c < 4; ++c)
        for (int d = 0; d < 4; ++d) {
            K_dd->at(a,b) += proj_ud(c,a)*proj_ud(d,b)*rT_dd(c,d);
        }
    }
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL calc_E_flux(
        CCTK_REAL const alp,
        tensor::generic<CCTK_REAL, 4, 1> const & beta_u,
        CCTK_REAL const E,
        tensor::generic<CCTK_REAL, 4, 1> const & F_u,
        int const dir) {
    return alp*F_u(dir) - beta_u(dir)*E;
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL calc_F_flux(
        CCTK_REAL const alp,
        tensor::generic<CCTK_REAL, 4, 1> const & beta_u,
        tensor::generic<CCTK_REAL, 4, 1> const & F_d,
        tensor::generic<CCTK_REAL, 4, 2> const & P_ud,
        int const dir,
        int const comp) {
    return alp*P_ud(dir,comp) - beta_u(dir)*F_d(comp);
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void calc_rad_sources(
        CCTK_REAL const eta,
        CCTK_REAL const kabs,
        CCTK_REAL const kscat,
        tensor::generic<CCTK_REAL, 4, 1> const & u_d,
        CCTK_REAL const J,
        tensor::generic<CCTK_REAL, 4, 1> const H_d,
        tensor::generic<CCTK_REAL, 4, 1> * S_d) {
    for (int a = 0; a < 4; ++a) {
        S_d->at(a) = (eta - kabs*J)*u_d(a) - (kabs + kscat)*H_d(a);
    }
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL calc_rE_source(
        CCTK_REAL const alp,
        tensor::generic<CCTK_REAL, 4, 1> const & n_u,
        tensor::generic<CCTK_REAL, 4, 1> const & S_d) {
    return - alp*tensor::dot(n_u, S_d);
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void calc_rF_source(
        CCTK_REAL const alp,
        tensor::generic<CCTK_REAL, 4, 2> const gamma_ud,
        tensor::generic<CCTK_REAL, 4, 1> const & S_d,
        tensor::generic<CCTK_REAL, 4, 1> * tS_d) {
    for (int a = 0; a < 4; ++a) {
        tS_d->at(a) = 0.0;
        for (int b = 0; b < 4; ++b) {
            tS_d->at(a) += alp*gamma_ud(b,a)*S_d(b);
        }
    }
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void apply_floor(
        tensor::symmetric2<CCTK_REAL, 4, 2> const & g_uu,
        CCTK_REAL * E,
        tensor::generic<CCTK_REAL, 4, 1> * F_d) {
    DECLARE_CCTK_PARAMETERS;

    *E = max(rad_E_floor, *E);

    CCTK_REAL const F2 = tensor::dot(g_uu, *F_d, *F_d);
    CCTK_REAL const lim = (*E)*(*E)*(1 - rad_eps);
    if (F2 > lim) {
        CCTK_REAL fac = lim/F2;
        for (int a = 0; a < 4; ++a) {
            F_d->at(a) *= fac;
        }
    }
}

} // namespace m1
} // namespace thc
