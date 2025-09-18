#include <cassert>
#include <cstring>
#include <cmath>

#include <loop_device.hxx>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include <mat.hxx>
#include <vec.hxx>
#include <sum.hxx>
#include <simd.hxx>

#include "nuX_M1_closure.hxx"
#include "nuX_M1_macro.hxx"
#include "nuX_utils.hxx"
#include "aster_utils.hxx"

#define GFINDEX1D(__k, ig, iv)                                                 \
  ((iv) + (ig) * 5 + (__k) * (5 * ngroups * nspecies))

#define PINDEX1D(ig, iv) ((iv) + (ig) * 5)

namespace nuX_M1 {

using namespace Loop;
using namespace Arith;
using namespace AsterUtils;
using namespace nuX_Utils;
using namespace std;

CCTK_HOST CCTK_DEVICE inline CCTK_REAL minmod2(CCTK_REAL rl, CCTK_REAL rp,
                                               CCTK_REAL th) {
  CCTK_REAL val = 0.0;
  if (rl * rp > 0.0) {
    val = copysign(min(th * fabs(rl), th * fabs(rp)), rl);
  }
  return min(1.0, val);
}

CCTK_DEVICE CCTK_HOST inline int face_stride(int dir, const int *lsh) {
  // number of grid points in one face-centered component
  // vcc: (Nx+1, Ny, Nz), cvc: (Nx, Ny+1, Nz), ccv: (Nx, Ny, Nz+1)
  return (dir == 0)   ? (lsh[0] + 1) * lsh[1] * lsh[2]
         : (dir == 1) ? lsh[0] * (lsh[1] + 1) * lsh[2]
                      : lsh[0] * lsh[1] * (lsh[2] + 1);
}

CCTK_DEVICE CCTK_HOST inline vect<int, 3> face_centering(int dir) {
  // CI,CJ,CK template flags used by Loop; 0 means "face/vertex" in that axis.
  return {dir == 0 ? 0 : 1, dir == 1 ? 0 : 1, dir == 2 ? 0 : 1};
}

// Read a “conserved” scalar U^{(iv)} for limiter purposes from the 5-tuple.
// iv: 0->rN, 1->rFx, 2->rFy, 3->rFz, 4->rE
CCTK_HOST CCTK_DEVICE inline CCTK_REAL read_cons_by_iv(
    int iv, const CCTK_REAL *__restrict__ rN, const CCTK_REAL *__restrict__ rFx,
    const CCTK_REAL *__restrict__ rFy, const CCTK_REAL *__restrict__ rFz,
    const CCTK_REAL *__restrict__ rE, int idx) {
  return (iv == 0)   ? rN[idx]
         : (iv == 1) ? rFx[idx]
         : (iv == 2) ? rFy[idx]
         : (iv == 3) ? rFz[idx]
                     : rE[idx];
}

//---------------------------------------------------------------------
// (A) Flux calculation kernel (face-centred write-out) — with L/R recon,
//     Lax–Friedrichs blending, minmod limiter φ, and opacity factor A_{i+1/2}
//     per the paper’s §3.1 transport scheme.
//     Assumes ghost zones have been filled and cctk_nghostzones[dir] >= 2.
//---------------------------------------------------------------------
template <int dir> void M1_CalcFlux(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_M1_CalcFluxes;
  DECLARE_CCTK_PARAMETERS;

  // We need LL,L | R,RR -> 2 ghost layers in the sweep direction
  assert(cctk_nghostzones[dir] >= 2);

  const GridDescBaseDevice grid(cctkGH);

  // cell-centred and face-centred layouts
  const GF3D2layout layout_cc(cctkGH, {1, 1, 1});
  const GF3D2layout layout_fc(cctkGH, face_centering(dir));

  // packed face buffer for this direction
  CCTK_REAL *nu_flux_dir = (dir == 0   ? nu_flux_x
                            : dir == 1 ? nu_flux_y
                                       : nu_flux_z);

  const int STRIDE = face_stride(dir, cctk_lsh);

  // geometry helpers (currently sampling cell-centred; replace with proper
  // vertex→face interpolation once available)
  tensor::slicing_geometry_const geom(alp, betax, betay, betaz, gxx, gxy, gxz,
                                      gyy, gyz, gzz, kxx, kxy, kxz, kyy, kyz,
                                      kzz, volform);

  // fiducial velocity (cell-centred)
  tensor::fluid_velocity_field_const fidu(alp, betax, betay, betaz,
                                          fidu_w_lorentz, fidu_velx, fidu_vely,
                                          fidu_velz);

  // grid spacings
  const CCTK_REAL delta[3] = {CCTK_DELTA_SPACE(0), CCTK_DELTA_SPACE(1),
                              CCTK_DELTA_SPACE(2)};

  // Loop over interior faces for this direction
  grid.loop_int_device<(dir == 0 ? 0 : 1), (dir == 1 ? 0 : 1),
                       (dir == 2 ? 0 : 1)>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        // Face linear index where we STORE flux
        const int ijk_fc = layout_fc.linear(p.i, p.j, p.k);

        // Adjacent cell centers for this face: L = face − ê_dir, R = face
        const int iL = p.i - (dir == 0), jL = p.j - (dir == 1),
                  kL = p.k - (dir == 2);
        const int iR = p.i, jR = p.j, kR = p.k;

        const int idxL = layout_cc.linear(iL, jL, kL);
        const int idxR = layout_cc.linear(iR, jR, kR);

        // Two more along the ray for the limiter: LL and RR (ghosts are valid)
        const int iLL = iL - (dir == 0), jLL = jL - (dir == 1),
                  kLL = kL - (dir == 2);
        const int iRR = iR + (dir == 0), jRR = jR + (dir == 1),
                  kRR = kR + (dir == 2);

        const int idxLL = layout_cc.linear(iLL, jLL, kLL);
        const int idxRR = layout_cc.linear(iRR, jRR, kRR);

        // Geometry at L and R (cell-centred); cheap face averages for α,β,γ^-1
        tensor::inv_metric<3> gamma_uu_L, gamma_uu_R;
        tensor::inv_metric<4> g_uu_L, g_uu_R;
        tensor::generic<CCTK_REAL, 4, 1> beta_u_L, beta_u_R;

        geom.get_inv_metric(idxL, &gamma_uu_L);
        geom.get_inv_metric(idxL, &g_uu_L);
        geom.get_shift_vec(idxL, &beta_u_L);

        geom.get_inv_metric(idxR, &gamma_uu_R);
        geom.get_inv_metric(idxR, &g_uu_R);
        geom.get_shift_vec(idxR, &beta_u_R);

        const CCTK_REAL alpha_L = alp[idxL];
        const CCTK_REAL alpha_R = alp[idxR];
        const CCTK_REAL alphaF = 0.5 * (alpha_L + alpha_R);

        tensor::generic<CCTK_REAL, 4, 1> beta_uF;
        beta_uF(1) = 0.5 * (beta_u_L(1) + beta_u_R(1));
        beta_uF(2) = 0.5 * (beta_u_L(2) + beta_u_R(2));
        beta_uF(3) = 0.5 * (beta_u_L(3) + beta_u_R(3));

        tensor::inv_metric<3> gamma_uuF;
        for (int a = 1; a <= 3; ++a)
          for (int b = 1; b <= 3; ++b)
            gamma_uuF(a, b) = 0.5 * (gamma_uu_L(a, b) + gamma_uu_R(a, b));

        // Face wavespeed c_face = max(|β_dir|+α sqrt(γ^{dir,dir})) on L/R
        auto face_wavespeed = [](const tensor::inv_metric<3> &gamma_uu,
                                 const tensor::generic<CCTK_REAL, 4, 1> &beta_u,
                                 CCTK_REAL alpha,
                                 int idir) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          const CCTK_REAL gdd = gamma_uu(idir + 1, idir + 1);
          const CCTK_REAL cgeom = alpha * sqrt(fabs(gdd));
          const CCTK_REAL bdir = fabs(beta_u(idir + 1));
          return cgeom + bdir;
        };

        const CCTK_REAL cL = face_wavespeed(gamma_uu_L, beta_u_L, alpha_L, dir);
        const CCTK_REAL cR = face_wavespeed(gamma_uu_R, beta_u_R, alpha_R, dir);
        const CCTK_REAL c_face = fmax(cL, cR);

        // Opacity average for A (abs+scat)
        const CCTK_REAL kappa_L = abs_1[idxL] + scat_1[idxL];
        const CCTK_REAL kappa_R = abs_1[idxR] + scat_1[idxR];
        const CCTK_REAL kappa_ave = 0.5 * (kappa_L + kappa_R);

        CCTK_REAL A = 1.0;
        if (kappa_ave * delta[dir] > 1.0) {
          A = fmin(1.0, 1.0 / (delta[dir] * kappa_ave));
          A = fmax(A, mindiss);
        }

        // Fiducial 4-velocity and 3-velocity at L and R
        tensor::generic<CCTK_REAL, 4, 1> u_u_L, u_u_R, v_u_L, v_u_R;
        fidu.get(idxL, &u_u_L);
        fidu.get(idxR, &u_u_R);
        pack_v_u(fidu_velx[idxL], fidu_vely[idxL], fidu_velz[idxL], &v_u_L);
        pack_v_u(fidu_velx[idxR], fidu_vely[idxR], fidu_velz[idxR], &v_u_R);

        const int groupspec = ngroups * nspecies;

        for (int ig = 0; ig < groupspec; ++ig) {

          // Assemble LEFT tensors and physical fluxes
          const int i4D_L = layout_cc.linear(iL, jL, kL, ig);

          tensor::generic<CCTK_REAL, 4, 1> H_d_L, H_u_L, F_u_L, F_d_L, fnu_u_L;
          tensor::symmetric2<CCTK_REAL, 4, 2> P_dd_L;
          tensor::generic<CCTK_REAL, 4, 2> P_ud_L;

          pack_F_d(betax[idxL], betay[idxL], betaz[idxL], rFx[i4D_L],
                   rFy[i4D_L], rFz[i4D_L], &F_d_L);
          pack_H_d(rHt[i4D_L], rHx[i4D_L], rHy[i4D_L], rHz[i4D_L], &H_d_L);
          pack_P_dd(betax[idxL], betay[idxL], betaz[idxL], rPxx[i4D_L],
                    rPxy[i4D_L], rPxz[i4D_L], rPyy[i4D_L], rPyz[i4D_L],
                    rPzz[i4D_L], &P_dd_L);

          tensor::contract(g_uu_L, H_d_L, &H_u_L);
          tensor::contract(g_uu_L, F_d_L, &F_u_L);
          tensor::contract(g_uu_L, P_dd_L, &P_ud_L);
          assemble_fnu(u_u_L, rJ[i4D_L], H_u_L, &fnu_u_L, rad_E_floor);

          const CCTK_REAL Gamma_L =
              compute_Gamma(fidu_w_lorentz[idxL], v_u_L, rJ[i4D_L], rE[i4D_L],
                            F_d_L, rad_E_floor, rad_eps);
          const CCTK_REAL nnu_L = rN[i4D_L] / fmax(Gamma_L, 1.0);

          CCTK_REAL flux_L[5];
          flux_L[0] = alphaF * nnu_L * fnu_u_L(dir + 1);
          flux_L[1] = calc_F_flux(alphaF, beta_uF, F_d_L, P_ud_L, dir + 1, 1);
          flux_L[2] = calc_F_flux(alphaF, beta_uF, F_d_L, P_ud_L, dir + 1, 2);
          flux_L[3] = calc_F_flux(alphaF, beta_uF, F_d_L, P_ud_L, dir + 1, 3);
          flux_L[4] = calc_E_flux(alphaF, beta_uF, rE[i4D_L], F_u_L, dir + 1);

          // Assemble RIGHT tensors and physical fluxes
          const int i4D_R = layout_cc.linear(iR, jR, kR, ig);

          tensor::generic<CCTK_REAL, 4, 1> H_d_R, H_u_R, F_u_R, F_d_R, fnu_u_R;
          tensor::symmetric2<CCTK_REAL, 4, 2> P_dd_R;
          tensor::generic<CCTK_REAL, 4, 2> P_ud_R;

          pack_F_d(betax[idxR], betay[idxR], betaz[idxR], rFx[i4D_R],
                   rFy[i4D_R], rFz[i4D_R], &F_d_R);
          pack_H_d(rHt[i4D_R], rHx[i4D_R], rHy[i4D_R], rHz[i4D_R], &H_d_R);
          pack_P_dd(betax[idxR], betay[idxR], betaz[idxR], rPxx[i4D_R],
                    rPxy[i4D_R], rPxz[i4D_R], rPyy[i4D_R], rPyz[i4D_R],
                    rPzz[i4D_R], &P_dd_R);

          tensor::contract(g_uu_R, H_d_R, &H_u_R);
          tensor::contract(g_uu_R, F_d_R, &F_u_R);
          tensor::contract(g_uu_R, P_dd_R, &P_ud_R);
          assemble_fnu(u_u_R, rJ[i4D_R], H_u_R, &fnu_u_R, rad_E_floor);

          const CCTK_REAL Gamma_R =
              compute_Gamma(fidu_w_lorentz[idxR], v_u_R, rJ[i4D_R], rE[i4D_R],
                            F_d_R, rad_E_floor, rad_eps);
          const CCTK_REAL nnu_R = rN[i4D_R] / fmax(Gamma_R, 1.0);

          CCTK_REAL flux_R[5];
          flux_R[0] = alphaF * nnu_R * fnu_u_R(dir + 1);
          flux_R[1] = calc_F_flux(alphaF, beta_uF, F_d_R, P_ud_R, dir + 1, 1);
          flux_R[2] = calc_F_flux(alphaF, beta_uF, F_d_R, P_ud_R, dir + 1, 2);
          flux_R[3] = calc_F_flux(alphaF, beta_uF, F_d_R, P_ud_R, dir + 1, 3);
          flux_R[4] = calc_E_flux(alphaF, beta_uF, rE[i4D_R], F_u_R, dir + 1);

          // Limiter φ per component using stencil (LL, L | R, RR)
          CCTK_REAL phi[5] = {0, 0, 0, 0, 0};
          CCTK_REAL sawtooth_mask[5] = {0, 0, 0, 0, 0};

          for (int iv = 0; iv < 5; ++iv) {
            const CCTK_REAL uLL =
                read_cons_by_iv(iv, rN, rFx, rFy, rFz, rE, idxLL);
            const CCTK_REAL uL =
                read_cons_by_iv(iv, rN, rFx, rFy, rFz, rE, idxL);
            const CCTK_REAL uR =
                read_cons_by_iv(iv, rN, rFx, rFy, rFz, rE, idxR);
            const CCTK_REAL uRR =
                read_cons_by_iv(iv, rN, rFx, rFy, rFz, rE, idxRR);

            const CCTK_REAL dup = uRR - uR;
            const CCTK_REAL duc = uR - uL;
            const CCTK_REAL dum = uL - uLL;

            CCTK_REAL ph = 0.0;
            bool saw = false;

            if (dup * duc > 0 && dum * duc > 0) {
              const CCTK_REAL rl = (duc != 0.0) ? (dum / duc) : 0.0;
              const CCTK_REAL rp = (duc != 0.0) ? (dup / duc) : 0.0;
              ph = minmod2(rl, rp, minmod_theta);
            } else if (dup * duc < 0 && dum * duc < 0) {
              saw = true;
            }
            phi[iv] = ph;
            sawtooth_mask[iv] = saw ? 1.0 : 0.0;
          }

          // LF and high-order fluxes; blend with φ and A (paper §3.1)
          for (int iv = 0; iv < 5; ++iv) {
            const CCTK_REAL U_L =
                read_cons_by_iv(iv, rN, rFx, rFy, rFz, rE, idxL);
            const CCTK_REAL U_R =
                read_cons_by_iv(iv, rN, rFx, rFy, rFz, rE, idxR);

            const CCTK_REAL fL = (iv == 0   ? flux_L[0]
                                  : iv == 1 ? flux_L[1]
                                  : iv == 2 ? flux_L[2]
                                  : iv == 3 ? flux_L[3]
                                            : flux_L[4]);
            const CCTK_REAL fR = (iv == 0   ? flux_R[0]
                                  : iv == 1 ? flux_R[1]
                                  : iv == 2 ? flux_R[2]
                                  : iv == 3 ? flux_R[3]
                                            : flux_R[4]);

            const CCTK_REAL flow = 0.5 * (fL + fR) - 0.5 * c_face * (U_R - U_L);
            const CCTK_REAL fhigh = 0.5 * (fL + fR);

            const CCTK_REAL limiter = 1.0 - phi[iv];
            const CCTK_REAL Aeff = (sawtooth_mask[iv] > 0.5) ? 1.0 : A;
            const CCTK_REAL fnum = fhigh - Aeff * limiter * (fhigh - flow);

            const int comp = PINDEX1D(ig, iv);
            nu_flux_dir[comp * STRIDE + ijk_fc] = fnum;

#ifndef NDEBUG
            assert(isfinite(nu_flux_dir[comp * STRIDE + ijk_fc]));
#endif
          }
        } // ig
      } // lambda
  ); // loop_int_device (faces)
}

//---------------------------------------------------------------------
// (B) RHS update kernel (cell-centred divergence of face fluxes)
//---------------------------------------------------------------------
template <int dir> void M1_UpdateRHSFromFluxes(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_M1_UpdateRHSFromFluxes;
  DECLARE_CCTK_PARAMETERS;

  const GridDescBaseDevice grid(cctkGH);

  // layouts and stride
  const GF3D2layout layout_cc(cctkGH, {1, 1, 1});
  const GF3D2layout layout_fc(cctkGH, face_centering(dir));
  const int STRIDE = face_stride(dir, cctk_lsh);

  // read-only view of packed face buffer for this direction
  const CCTK_REAL *nu_flux_dir = (dir == 0   ? nu_flux_x
                                  : dir == 1 ? nu_flux_y
                                             : nu_flux_z);

  // RHS component pointers
  CCTK_REAL *r_rhs[5] = {rN_rhs, rFx_rhs, rFy_rhs, rFz_rhs, rE_rhs};

  // grid spacings
  const CCTK_REAL delta[3] = {CCTK_DELTA_SPACE(0), CCTK_DELTA_SPACE(1),
                              CCTK_DELTA_SPACE(2)};
  const CCTK_REAL idelta[3] = {1.0 / delta[0], 1.0 / delta[1], 1.0 / delta[2]};

  // Loop over cell centres
  grid.loop_int_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        const int ijk_cc = layout_cc.linear(p.i, p.j, p.k);

        // indices of faces bounding this cell in 'dir'
        const int fL = layout_fc.linear(p.i, p.j, p.k);
        const int fR = (dir == 0)   ? layout_fc.linear(p.i + 1, p.j, p.k)
                       : (dir == 1) ? layout_fc.linear(p.i, p.j + 1, p.k)
                                    : layout_fc.linear(p.i, p.j, p.k + 1);

        const int groupspec = ngroups * nspecies;

        for (int ig = 0; ig < groupspec; ++ig) {
          const int i4D = layout_cc.linear(p.i, p.j, p.k, ig);

          for (int iv = 0; iv < 5; ++iv) {
            const int comp = PINDEX1D(ig, iv);
            const CCTK_REAL flux_L = nu_flux_dir[comp * STRIDE + fL];
            const CCTK_REAL flux_R = nu_flux_dir[comp * STRIDE + fR];

            if (!nuX_m1_mask[ijk_cc]) {
              r_rhs[iv][i4D] += idelta[dir] * (flux_L - flux_R);
#ifndef NDEBUG
              assert(isfinite(r_rhs[iv][i4D]));
#endif
            }
          } // iv
        } // ig
      } // lambda
  ); // loop_int_device (cells)
}

//---------------------------------------------------------------------
// Wrappers
//---------------------------------------------------------------------
extern "C" void nuX_M1_CalcFluxes(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_M1_CalcFluxes;
  DECLARE_CCTK_PARAMETERS;
  if (verbose) {
    CCTK_INFO(
        "nuX_M1_CalcFluxes (face-centered, LF+HO blend with limiter and A)");
  }
  M1_CalcFlux<0>(cctkGH);
  M1_CalcFlux<1>(cctkGH);
  M1_CalcFlux<2>(cctkGH);
}

extern "C" void nuX_M1_UpdateRHSFromFluxes(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_M1_UpdateRHSFromFluxes;
  DECLARE_CCTK_PARAMETERS;
  if (verbose) {
    CCTK_INFO(
        "nuX_M1_UpdateRHSFromFluxes (divergence from face-centered fluxes)");
  }
  M1_UpdateRHSFromFluxes<0>(cctkGH);
  M1_UpdateRHSFromFluxes<1>(cctkGH);
  M1_UpdateRHSFromFluxes<2>(cctkGH);
}

} // namespace nuX_M1
