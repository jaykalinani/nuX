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


#ifndef NUX_M1_CLOSURE_HH
#define NUX_M1_CLOSURE_HH

// TODO: GSL dependence

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>

#include "cctk.h"

#include "utils_tensor.hh"

// TODO: Rename namespaces

namespace nuX_M1 {

using namespace utils;

// TODO: Dependence on tensor::

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void pack_F_d(
        CCTK_REAL const betax, CCTK_REAL const betay, CCTK_REAL const betaz,
        CCTK_REAL const Fx, CCTK_REAL const Fy, CCTK_REAL const Fz,
        tensor::generic<CCTK_REAL, 4, 1> * F_d);

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void unpack_F_d(
        tensor::generic<CCTK_REAL, 4, 1> const & F_d,
        CCTK_REAL * Fx, CCTK_REAL * Fy, CCTK_REAL * Fz);

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void pack_F_d(
        CCTK_REAL const Fx, CCTK_REAL const Fy, CCTK_REAL const Fz,
        tensor::generic<CCTK_REAL, 3, 1> * F_d);

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void unpack_F_d(
        tensor::generic<CCTK_REAL, 3, 1> const & F_d,
        CCTK_REAL * Fx, CCTK_REAL * Fy, CCTK_REAL * Fz);

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void pack_H_d(
        CCTK_REAL const Ht, CCTK_REAL const Hx, CCTK_REAL const Hy, CCTK_REAL const Hz,
        tensor::generic<CCTK_REAL, 4, 1> * H_d);

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void unpack_H_d(
        tensor::generic<CCTK_REAL, 4, 1> const & H_d,
        CCTK_REAL * Ht, CCTK_REAL * Hx, CCTK_REAL * Hy, CCTK_REAL * Hz);

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void pack_P_dd(
        CCTK_REAL const betax, CCTK_REAL const betay, CCTK_REAL const betaz,
        CCTK_REAL const Pxx, CCTK_REAL const Pxy, CCTK_REAL const Pxz,
        CCTK_REAL const Pyy, CCTK_REAL const Pyz, CCTK_REAL const Pzz,
        tensor::symmetric2<CCTK_REAL, 4, 2> * P_dd);

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void unpack_P_dd(
        tensor::symmetric2<CCTK_REAL, 4, 2> const & P_dd,
        CCTK_REAL * Pxx, CCTK_REAL * Pxy, CCTK_REAL * Pxz,
        CCTK_REAL * Pyy, CCTK_REAL * Pyz, CCTK_REAL * Pzz);

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void pack_P_dd(
        CCTK_REAL const Pxx, CCTK_REAL const Pxy, CCTK_REAL const Pxz,
        CCTK_REAL const Pyy, CCTK_REAL const Pyz, CCTK_REAL const Pzz,
        tensor::symmetric2<CCTK_REAL, 3, 2> * P_dd);

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void unpack_P_dd(
        tensor::symmetric2<CCTK_REAL, 3, 2> const & P_dd,
        CCTK_REAL * Pxx, CCTK_REAL * Pxy, CCTK_REAL * Pxz,
        CCTK_REAL * Pyy, CCTK_REAL * Pyz, CCTK_REAL * Pzz);

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void pack_v_u(
        CCTK_REAL const velx, CCTK_REAL const vely, CCTK_REAL const velz,
        tensor::generic<CCTK_REAL, 4, 1> * v_u);

// Fluid projector: delta^a_b + u^a u_b
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void calc_proj(
        tensor::generic<CCTK_REAL, 4, 1> const & u_d,
        tensor::generic<CCTK_REAL, 4, 1> const & u_u,
        tensor::generic<CCTK_REAL, 4, 2> * proj_ud);

// Compute the closure in the thin limit
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void calc_Pthin(
        tensor::inv_metric<4> const & g_uu,
        CCTK_REAL const E,
        tensor::generic<CCTK_REAL, 4, 1> const & F_d,
        tensor::symmetric2<CCTK_REAL, 4, 2> * P_dd);

// Compute the closure in the thick limit
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void calc_Pthick(
        tensor::metric<4> const & g_dd,
        tensor::inv_metric<4> const & g_uu,
        tensor::generic<CCTK_REAL, 4, 1> const & n_d,
        CCTK_REAL const w_lorentz,
        tensor::generic<CCTK_REAL, 4, 1> const & v_d,
        CCTK_REAL const E,
        tensor::generic<CCTK_REAL, 4, 1> const & F_d,
        tensor::symmetric2<CCTK_REAL, 4, 2> * P_dd);

// Computes the flux factor xi = H_a H^a / J^2
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL flux_factor(
        tensor::inv_metric<4> const & g_uu,
        CCTK_REAL const J,
        tensor::generic<CCTK_REAL, 4, 1> const & H_d);

// Closures
typedef CCTK_REAL (*closure_t)(CCTK_REAL const);
CCTK_REAL eddington(CCTK_REAL const xi);
CCTK_REAL kershaw(CCTK_REAL const xi);
CCTK_REAL minerbo(CCTK_REAL const xi);
CCTK_REAL thin(CCTK_REAL const xi);

// Computes the closure in the lab frame given the Eddington factor chi
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
        tensor::symmetric2<CCTK_REAL, 4, 2> * P_dd);

// Computes the closure in the lab frame with a rootfinding procedure
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void calc_closure(
        cGH const * cctkGH,
        int const i, int const j, int const k,
        int const ig,
        closure_t closure,
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
        tensor::symmetric2<CCTK_REAL, 4, 2> * P_dd);

// Assemble the unit-norm radiation number current
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void assemble_fnu(
        tensor::generic<CCTK_REAL, 4, 1> const & u_u,
        CCTK_REAL const J,
        tensor::generic<CCTK_REAL, 4, 1> const & H_u,
        tensor::generic<CCTK_REAL, 4, 1> * fnu_u);

// Compute the ratio of neutrino number densities in the lab and fluid frame
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL compute_Gamma(
        CCTK_REAL const W,
        tensor::generic<CCTK_REAL, 4, 1> const & v_u,
        CCTK_REAL const J,
        CCTK_REAL const E,
        tensor::generic<CCTK_REAL, 4, 1> const & F_d);

// Assemble the radiation stress tensor in any frame
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void assemble_rT(
        tensor::generic<CCTK_REAL, 4, 1> const & u_d,
        CCTK_REAL const J,
        tensor::generic<CCTK_REAL, 4, 1> const & H_d,
        tensor::symmetric2<CCTK_REAL, 4, 2> const & K_dd,
        tensor::symmetric2<CCTK_REAL, 4, 2> * rT_dd);

// Project out the radiation energy (in any frame)
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL calc_J_from_rT(
        tensor::symmetric2<CCTK_REAL, 4, 2> const & rT_dd,
        tensor::generic<CCTK_REAL, 4, 1> const & u_u);

// Project out the radiation fluxes (in any frame)
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void calc_H_from_rT(
        tensor::symmetric2<CCTK_REAL, 4, 2> const & rT_dd,
        tensor::generic<CCTK_REAL, 4, 1> const & u_u,
        tensor::generic<CCTK_REAL, 4, 2> const & proj_ud,
        tensor::generic<CCTK_REAL, 4, 1> * H_d);

// Project out the radiation pressure tensor (in any frame)
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void calc_K_from_rT(
        tensor::symmetric2<CCTK_REAL, 4, 2> const & rT_dd,
        tensor::generic<CCTK_REAL, 4, 2> const & proj_ud,
        tensor::symmetric2<CCTK_REAL, 4, 2> * K_dd);

// Compute the radiation energy flux
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL calc_E_flux(
        CCTK_REAL const alp,
        tensor::generic<CCTK_REAL, 4, 1> const & beta_u,
        CCTK_REAL const E,
        tensor::generic<CCTK_REAL, 4, 1> const & F_u,
        int const direction);                           // 1:x 2:y 3:z

// Compute the flux of neutrino energy flux
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL calc_F_flux(
        CCTK_REAL const alp,
        tensor::generic<CCTK_REAL, 4, 1> const & beta_u,
        tensor::generic<CCTK_REAL, 4, 1> const & F_d,
        tensor::generic<CCTK_REAL, 4, 2> const & P_ud,
        int const direction,                             // 1:x 2:y 3:z
        int const component);                            // 1:x 2:y 3:z

// Computes the sources S_a = [eta - k_abs J] u_a - [k_abs + k_scat] H_a
// WARNING: be consistent with the densitization of eta, J, and H_d
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void calc_rad_sources(
        CCTK_REAL const eta,
        CCTK_REAL const kabs,
        CCTK_REAL const kscat,
        tensor::generic<CCTK_REAL, 4, 1> const & u_d,
        CCTK_REAL const J,
        tensor::generic<CCTK_REAL, 4, 1> const H_d,
        tensor::generic<CCTK_REAL, 4, 1> * S_d);

// Computes the source term for E: -alp n^a S_a
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL calc_rE_source(
        CCTK_REAL const alp,
        tensor::generic<CCTK_REAL, 4, 1> const & n_u,
        tensor::generic<CCTK_REAL, 4, 1> const & S_d);

// Computes the source term for F_a: alp gamma^b_a S_b
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void calc_rF_source(
        CCTK_REAL const alp,
        tensor::generic<CCTK_REAL, 4, 2> const gamma_ud,
        tensor::generic<CCTK_REAL, 4, 1> const & S_d,
        tensor::generic<CCTK_REAL, 4, 1> * tS_d);

// Enforce that E > rad_E_floor and F_a F^a < (1 - rad_eps) E^2
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void apply_floor(
        tensor::symmetric2<CCTK_REAL, 4, 2> const & g_uu,
        CCTK_REAL * E,
        tensor::generic<CCTK_REAL, 4, 1> * F_d);

} // namespace nuX_M1

#endif
