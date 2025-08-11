#include <cassert>
#include <cstring>
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

#define NDIM 3

#define GFINDEX1D(__k, ig, iv)                                                 \
  ((iv) + (ig) * 5 + (__k) * (5 * ngroups * nspecies))

#define PINDEX1D(ig, iv) ((iv) + (ig) * 5)

// TODO: User–defined limits; adjust as needed
#ifndef nuX_M1_NGHOST
#define nuX_M1_NGHOST 2
#endif

#ifndef nuX_M1_MAX_K
#define nuX_M1_MAX_K 512
#endif

#ifndef MAX_GROUPSPECIES
#define MAX_GROUPSPECIES 64
#endif

namespace {

CCTK_HOST CCTK_DEVICE inline CCTK_REAL minmod2(CCTK_REAL rl, CCTK_REAL rp, CCTK_REAL th) {
  return min(1.0, min(th * rl, th * rp));
}

} // namespace

namespace nuX_M1 {

using namespace Loop;
using namespace Arith;
using namespace AsterUtils;
using namespace nuX_Utils;
using namespace std;

// Compute the numerical fluxes using a simple 2nd order flux-limited method
// with high-Peclet limit fix. The fluxes are then added to the RHSs.
template <int dir> void M1_CalcFlux(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_M1_CalcFluxes;
  DECLARE_CCTK_PARAMETERS;

  /* grid functions for fluxes */

  // const vec<GF3D2<CCTK_REAL>, 5> gf_rhs = {rN_rhs, rFx_rhs, rFy_rhs, rFz_rhs,
  // rE_rhs};

  /* grid functions */
  const GF3D2layout layout(cctkGH, {0, 0, 0});

  const vec<GF3D2<const CCTK_REAL8>, dim> gf_beta{
      GF3D2<const CCTK_REAL8>(layout, betax),
      GF3D2<const CCTK_REAL8>(layout, betay),
      GF3D2<const CCTK_REAL8>(layout, betaz)};

  const smat<GF3D2<const CCTK_REAL8>, dim> gf_g{
      GF3D2<const CCTK_REAL8>(layout, gxx),
      GF3D2<const CCTK_REAL8>(layout, gxy),
      GF3D2<const CCTK_REAL8>(layout, gxz),
      GF3D2<const CCTK_REAL8>(layout, gyy),
      GF3D2<const CCTK_REAL8>(layout, gyz),
      GF3D2<const CCTK_REAL8>(layout, gzz)};

  const GridDescBaseDevice grid(cctkGH);
  const GF3D2layout layout2(cctkGH, {1, 1, 1});

  tensor::slicing_geometry_const geom(alp, betax, betay, betaz, gxx, gxy, gxz,
                                      gyy, gyz, gzz, kxx, kxy, kxz, kyy, kyz,
                                      kzz, volform);
  tensor::fluid_velocity_field_const fidu(alp, betax, betay, betaz,
                                          fidu_w_lorentz, fidu_velx, fidu_vely,
                                          fidu_velz);

  const int GFNSIZE = cctk_lsh[0] * cctk_lsh[1] * cctk_lsh[2];

  auto &nu_flux_dir = (dir == 0 ? nu_flux_x : dir == 1 ? nu_flux_y : nu_flux_z);

  // Reconstruct flux quantities at interfaces
  grid.loop_int_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        int const ijk = layout2.linear(p.i, p.j, p.k);

        /* Interpolate metric components from vertices to centers */
/*
       	const CCTK_REAL alp_avg = calc_avg_v2c(alp, p, dir);
        const vec<CCTK_REAL, 3> betas_avg([&](int i) ARITH_INLINE {
          return calc_avg_v2c(gf_beta(i), p, dir);
        });
*/
        // Loop over groups/species.
        int groupspec = ngroups * nspecies;

        for (int ig = 0; ig < groupspec; ++ig) {
          int const i4D = layout2.linear(p.i, p.j, p.k, ig);

          //--- Extract geometry at the face ---
          tensor::inv_metric<3> gamma_uu;
          tensor::inv_metric<4> g_uu;
          tensor::generic<CCTK_REAL, 4, 1> beta_u;

          geom.get_inv_metric(ijk, &gamma_uu);
          geom.get_inv_metric(ijk, &g_uu);
          geom.get_shift_vec(ijk, &beta_u);

          //--- Extract fluid velocity and pack the three-velocity ---
          tensor::generic<CCTK_REAL, 4, 1> u_u;
          fidu.get(ijk, &u_u);
          tensor::generic<CCTK_REAL, 4, 1> v_u;
          pack_v_u(fidu_velx[ijk], fidu_vely[ijk], fidu_velz[ijk], &v_u);

          //--- Reconstruct closure tensors ---
          tensor::generic<CCTK_REAL, 4, 1> H_d, H_u, F_u;
          tensor::generic<CCTK_REAL, 4, 1> F_d;
          tensor::symmetric2<CCTK_REAL, 4, 2> P_dd;
          tensor::generic<CCTK_REAL, 4, 2> P_ud;
          tensor::generic<CCTK_REAL, 4, 1> fnu_u;

          // Get the fluid variables from the corresponding cell (for simplicity
          // we use the face point).
          pack_F_d(betax[ijk], betay[ijk], betaz[ijk], rFx[i4D], rFy[i4D],
                   rFz[i4D], &F_d);
          pack_H_d(rHt[i4D], rHx[i4D], rHy[i4D], rHz[i4D], &H_d);
          pack_P_dd(betax[ijk], betay[ijk], betaz[ijk], rPxx[i4D], rPxy[i4D],
                    rPxz[i4D], rPyy[i4D], rPyz[i4D], rPzz[i4D], &P_dd);
          tensor::contract(g_uu, H_d, &H_u);
          tensor::contract(g_uu, F_d, &F_u);
          tensor::contract(g_uu, P_dd, &P_ud);
          assemble_fnu(u_u, rJ[i4D], H_u, &fnu_u);

          // Compute Gamma and nnu (densitized neutrino number).
          CCTK_REAL Gamma =
              compute_Gamma(fidu_w_lorentz[ijk], v_u, rJ[i4D], rE[i4D], F_d);
          CCTK_REAL nnu = rN[i4D] / Gamma;

          //--- Compute the flux components at the face ---
          // Note: dir+1 converts from 0-indexed C++ to the 1-indexed convention
          // used in the physics routines.
          CCTK_REAL flux[5];
          flux[0] = alp[ijk] * nnu * fnu_u(dir + 1);
          flux[1] = calc_F_flux(alp[ijk], beta_u, F_d, P_ud, dir + 1, 1);
          flux[2] = calc_F_flux(alp[ijk], beta_u, F_d, P_ud, dir + 1, 2);
          flux[3] = calc_F_flux(alp[ijk], beta_u, F_d, P_ud, dir + 1, 3);
          flux[4] = calc_E_flux(alp[ijk], beta_u, rE[i4D], F_u, dir + 1);

          // Verify flux values are finite
          for (int iv = 0; iv < 5; ++iv) {
            assert(isfinite(flux[iv]));
          }

          //--- Store the computed flux into the proper face-centered grid
          // function ---
          // TODO: need to define the following grid functions
          /*
                    if (dir == 0) {
                      for (int iv = 0; iv < 5; ++iv)
                        nu_flux_x[PINDEX1D(ig, iv)][ijk] = flux[iv];
                    } else if (dir == 1) {
                      for (int iv = 0; iv < 5; ++iv)
                        nu_flux_y[PINDEX1D(ig, iv)][ijk] = flux[iv];
                    } else if (dir == 2) {
                      for (int iv = 0; iv < 5; ++iv)
                        nu_flux_z[PINDEX1D(ig, iv)][ijk] = flux[iv];
                    }
          */
          for (int iv = 0; iv < 5; ++iv) {
            int comp = PINDEX1D(ig, iv);
            nu_flux_dir[comp * GFNSIZE + ijk] = flux[iv];
          }

        } // end loop over ig
      } // end lambda
  ); // end grid.loop_int_device
} // end M1_CalcFlux

//---------------------------------------------------------------------
// (B) M1_UpdateRHSFromFluxes<dir>
// For a given direction, loop over cell centers (centering <1,1,1>) to update
// the RHS using the flux differences computed at cell faces. For each cell the
// left face (at i_d-1) and right face (at i_d) are used. A flux limiter and
// high-Peclet dissipation fix is applied.
//---------------------------------------------------------------------
template <int dir> void M1_UpdateRHSFromFluxes(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_M1_UpdateRHSFromFluxes;
  DECLARE_CCTK_PARAMETERS;

  const GridDescBaseDevice grid(cctkGH);
  const GF3D2layout layout2(cctkGH, {1, 1, 1});
  const int GFNSIZE = cctk_lsh[0] * cctk_lsh[1] * cctk_lsh[2];

  // Grid data
  CCTK_REAL const delta[3] = {
      CCTK_DELTA_SPACE(0),
      CCTK_DELTA_SPACE(1),
      CCTK_DELTA_SPACE(2),
  };
  CCTK_REAL const idelta[3] = {
      1.0 / CCTK_DELTA_SPACE(0),
      1.0 / CCTK_DELTA_SPACE(1),
      1.0 / CCTK_DELTA_SPACE(2),
  };

  auto &nu_flux_dir = (dir == 0 ? nu_flux_x : dir == 1 ? nu_flux_y : nu_flux_z);
  CCTK_REAL *r_rhs[5] = {rN_rhs, rFx_rhs, rFy_rhs, rFz_rhs, rE_rhs};

  // Loop over cell centers.
  grid.loop_int_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        int const ijk = layout2.linear(p.i, p.j, p.k);

        // Determine indices of the two faces in direction dir.
        int face_left_idx, face_right_idx;
        if (dir == 0) {
          face_left_idx = layout2.linear(p.i - 1, p.j, p.k);
          face_right_idx = ijk;
        } else if (dir == 1) {
          face_left_idx = layout2.linear(p.i, p.j - 1, p.k);
          face_right_idx = ijk;
        } else { // dir == 2
          face_left_idx = layout2.linear(p.i, p.j, p.k - 1);
          face_right_idx = ijk;
        }

        int groupspec = ngroups * nspecies;

        // Loop over groups/species and over the 5 flux components.
        for (int ig = 0; ig < groupspec; ++ig) {
          for (int iv = 0; iv < 5; ++iv) {

            int comp = PINDEX1D(ig, iv);
            CCTK_REAL flux_left = nu_flux_dir[comp * GFNSIZE + face_left_idx];
            CCTK_REAL flux_right = nu_flux_dir[comp * GFNSIZE + face_right_idx];
            /*
                        CCTK_REAL flux_left, flux_right;
                        if (dir == 0) {
                          flux_left = nu_flux_x[PINDEX1D(ig,
               iv)][face_left_idx]; flux_right = nu_flux_x[PINDEX1D(ig,
               iv)][face_right_idx]; } else if (dir == 1) { flux_left =
               nu_flux_y[PINDEX1D(ig, iv)][face_left_idx]; flux_right =
               nu_flux_y[PINDEX1D(ig, iv)][face_right_idx]; } else { // dir == 2
                          flux_left = nu_flux_z[PINDEX1D(ig,
               iv)][face_left_idx]; flux_right = nu_flux_z[PINDEX1D(ig,
               iv)][face_right_idx];
                        }
            */

            // --- Optional: Apply flux-limiter/dissipation fix ---
            // For demonstration we recompute a limited flux using neighboring
            // cell-centered conserved data. Here we read four neighboring cell
            // values along the direction: For simplicity we use rN, rFx, etc.
            // (depending on iv) as representative of the conserved variable. In
            // a full implementation, you would select the proper variable
            // depending on iv.
            int shift[3];
            shift[0] = (dir == 0) ? 1 : 0;
            shift[1] = (dir == 1) ? 1 : 0;
            shift[2] = (dir == 2) ? 1 : 0;
            int idx_m, idx, idx_p, idx_pp;
            // Left neighbor:
            idx_m =
                layout2.linear(p.i - shift[0], p.j - shift[1], p.k - shift[2]);
            // Current cell:
            idx = ijk;
            // Right neighbor:
            idx_p =
                layout2.linear(p.i + shift[0], p.j + shift[1], p.k + shift[2]);
            // Next to right:
            idx_pp = layout2.linear(p.i + 2 * shift[0], p.j + 2 * shift[1],
                                    p.k + 2 * shift[2]);

            // Read the cell-centered conserved variable.
            // For demonstration we use an array "cons_center" selected by iv:
            // Let cons_center = rN for iv==0, rFx for iv==1, etc.
            CCTK_REAL ujm =
                (iv == 0
                     ? rN[idx_m]
                     : (iv == 1
                            ? rFx[idx_m]
                            : (iv == 2 ? rFy[idx_m]
                                       : (iv == 3 ? rFz[idx_m] : rE[idx_m]))));
            CCTK_REAL uj =
                (iv == 0
                     ? rN[idx]
                     : (iv == 1 ? rFx[idx]
                                : (iv == 2 ? rFy[idx]
                                           : (iv == 3 ? rFz[idx] : rE[idx]))));
            CCTK_REAL ujp =
                (iv == 0
                     ? rN[idx_p]
                     : (iv == 1
                            ? rFx[idx_p]
                            : (iv == 2 ? rFy[idx_p]
                                       : (iv == 3 ? rFz[idx_p] : rE[idx_p]))));
            CCTK_REAL ujpp =
                (iv == 0 ? rN[idx_pp]
                         : (iv == 1 ? rFx[idx_pp]
                                    : (iv == 2 ? rFy[idx_pp]
                                               : (iv == 3 ? rFz[idx_pp]
                                                          : rE[idx_pp]))));

            // Compute dissipation factor from nearby "abs" and "scat" values.
            CCTK_REAL kapp =
                0.5 * (abs_1[idx] + abs_1[idx_p] + scat_1[idx] + scat_1[idx_p]);
            CCTK_REAL A = 1.0;
            if (kapp * delta[dir] > 1.0) {
              A = min(1.0, 1.0 / (delta[dir] * kapp));
              A = max(A, mindiss);
            }

            // Compute difference ratios for limiter.
            CCTK_REAL dup = ujpp - ujp;
            CCTK_REAL duc = ujp - uj;
            CCTK_REAL dum = uj - ujm;
            bool sawtooth = false;
            CCTK_REAL phi = 0.0;
            if (dup * duc > 0 && dum * duc > 0) {
              phi = minmod2(dum / duc, dup / duc, minmod_theta);
            } else if (dup * duc < 0 && dum * duc < 0) {
              sawtooth = true;
            }
            assert(isfinite(phi));
            // Compute limited flux: we mimic the CPU formulation.
            // Here we combine the left and right fluxes.
            CCTK_REAL cmx = max(abs(flux_left), abs(flux_right));
            CCTK_REAL flux_low =
                0.5 * (flux_left + flux_right - cmx * (ujp - uj));
            CCTK_REAL flux_high = 0.5 * (flux_left + flux_right);
            CCTK_REAL flux_num = flux_high - (sawtooth ? 1.0 : A) * (1 - phi) *
                                                 (flux_high - flux_low);

            // Update the RHS from the difference of flux between left face
            // (flux_left) and the limited flux (flux_num). Here we use the
            // standard finite volume update:
            int const i4D = layout2.linear(p.i, p.j, p.k, ig);
            // We assume local_rhs is available as in the original CPU code;
            // here we use r_rhs alias. For demonstration, assume r_rhs for
            // component iv is:
            //   r_rhs[0] corresponds to rN_rhs, r_rhs[1] to rFx_rhs, etc.
            if (!nuX_m1_mask[ijk]) {
              r_rhs[iv][i4D] += idelta[dir] * (flux_left - flux_num);
              assert(isfinite(r_rhs[iv][i4D]));
            }
          } // end loop over iv
        } // end loop over ig
      } // end lambda
  ); // end grid.loop_int_device (cell centers)
} // end M1_UpdateRHSFromFluxes_GPU

//---------------------------------------------------------------------
// Host wrapper: This is the extern "C" function called by CarpetX.
// It launches the face–flux kernel and the cell–centered RHS update kernel for
// each spatial direction.
//---------------------------------------------------------------------
extern "C" void nuX_M1_CalcFluxes(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_M1_CalcFluxes;
  DECLARE_CCTK_PARAMETERS;

  if (verbose) {
    CCTK_INFO("nuX_M1_CalcFluxes");
  }
  // For each spatial direction, first compute the fluxes

  M1_CalcFlux<0>(cctkGH);
  M1_CalcFlux<1>(cctkGH);
  M1_CalcFlux<2>(cctkGH);
}

extern "C" void nuX_M1_UpdateRHSFromFluxes(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_M1_UpdateRHSFromFluxes;
  DECLARE_CCTK_PARAMETERS;

  if (verbose) {
    CCTK_INFO("nuX_M1_UpdateRHSFromFluxes");
  }
  // next, update the cell-centered RHS using those fluxes.

  M1_UpdateRHSFromFluxes<0>(cctkGH);
  M1_UpdateRHSFromFluxes<1>(cctkGH);
  M1_UpdateRHSFromFluxes<2>(cctkGH);
}

} // end namespace nuX_M1
