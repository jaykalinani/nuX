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

#define PINDEX1D(ig, iv) ((iv) + (ig) * 5)

namespace nuX_M1 {

using namespace Loop;
using namespace Arith;
using namespace AsterUtils;
using namespace nuX_Utils;
using namespace std;

//-------------------------------------------------------------------------
// Limiter helper
//-------------------------------------------------------------------------
CCTK_HOST CCTK_DEVICE inline CCTK_REAL minmod2(CCTK_REAL rl,
                                               CCTK_REAL rp,
                                               CCTK_REAL th) {
  return fmin(1.0, fmin(th * rl, th * rp));
}

//-------------------------------------------------------------------------
// Read a “conserved” scalar U^{(iv)} from the 5-tuple.
//
// iv: 0 -> rN, 1 -> rFx, 2 -> rFy, 3 -> rFz, 4 -> rE
//-------------------------------------------------------------------------
CCTK_HOST CCTK_DEVICE inline CCTK_REAL read_cons_by_iv(
    int iv,
    const CCTK_REAL *__restrict__ rN,
    const CCTK_REAL *__restrict__ rFx,
    const CCTK_REAL *__restrict__ rFy,
    const CCTK_REAL *__restrict__ rFz,
    const CCTK_REAL *__restrict__ rE,
    int idx) {
  return (iv == 0)   ? rN[idx]
         : (iv == 1) ? rFx[idx]
         : (iv == 2) ? rFy[idx]
         : (iv == 3) ? rFz[idx]
                     : rE[idx];
}

//-------------------------------------------------------------------------
// Characteristic speed in a given spatial direction:
//
//   lambda = -beta_dir ± alpha sqrt(gamma^{dir,dir})
//
// dir = 0,1,2 corresponds to x,y,z direction.
//-------------------------------------------------------------------------
CCTK_HOST CCTK_DEVICE inline CCTK_REAL face_speed(
    int dir,
    const tensor::inv_metric<3> &gamma,
    const tensor::generic<CCTK_REAL, 4, 1> &beta_u,
    CCTK_REAL alpha_cell) {

  // NOTE: gamma is a 3-metric (indices 0..2) like in THC.
  const CCTK_REAL gdd  = gamma(dir, dir);
  const CCTK_REAL root = alpha_cell * sqrt(fabs(gdd));

  const CCTK_REAL lamP = -beta_u(dir + 1) + root;
  const CCTK_REAL lamM = -beta_u(dir + 1) - root;

  return fmax(fabs(lamP), fabs(lamM));
}

//---------------------------------------------------------------------
// (A) Cell-centred physical flux kernel.
//
// For each direction `dir`, compute the *physical* radiative fluxes at
// cell centres and store them into `nu_flux_{x,y,z}`. No limiting or
// numerical flux construction happens here — that is done in the RHS
// update kernel using a 1D stencil along each direction.
//---------------------------------------------------------------------
template <int dir>
CCTK_DEVICE CCTK_HOST void M1_ComputePhysicalFluxes(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_M1_CalcFluxes;
  DECLARE_CCTK_PARAMETERS;

  const GridDescBaseDevice grid(cctkGH);
  const GF3D2layout layout_cc(cctkGH, {1, 1, 1});

  // These are treated as *cell-centred* storage for physical fluxes
  CCTK_REAL *restrict nu_flux_dir =
      (dir == 0 ? nu_flux_x : (dir == 1 ? nu_flux_y : nu_flux_z));

  tensor::slicing_geometry_const geom(alp, betax, betay, betaz,
                                      gxx, gxy, gxz, gyy, gyz, gzz,
                                      kxx, kxy, kxz, kyy, kyz, kzz,
                                      volform);

  tensor::fluid_velocity_field_const fidu(alp, betax, betay, betaz,
                                          fidu_w_lorentz,
                                          fidu_velx, fidu_vely, fidu_velz);

  const int groupspec = ngroups * nspecies;

  //--------------------------------------------------------------------
  // GPU loop over all cells (INCLUDING ghost zones)
  //
  // This matches THC's "1st pass compute the fluxes", which computes
  // fluxes over the full 1D line including ghosts so that the subsequent
  // RHS update can safely reference fluxes near the interior boundary.
  //--------------------------------------------------------------------
  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) {
        const int i = p.i;
        const int j = p.j;
        const int k = p.k;

        const int idx = layout_cc.linear(i, j, k);

        // Local geometry at this cell
        tensor::inv_metric<3> gamma_uu;
        tensor::inv_metric<4> g_uu;
        tensor::generic<CCTK_REAL, 4, 1> beta_u;

        geom.get_inv_metric(idx, &gamma_uu);
        geom.get_inv_metric(idx, &g_uu);
        geom.get_shift_vec(idx, &beta_u);

        const CCTK_REAL alpha = alp[idx];

        // Fiducial 4-velocity and 3-velocity at this cell
        tensor::generic<CCTK_REAL, 4, 1> u_u, v_u;
        fidu.get(idx, &u_u);
        pack_v_u(fidu_velx[idx], fidu_vely[idx], fidu_velz[idx], &v_u);

        // Temporary tensors used for flux construction
        tensor::generic<CCTK_REAL, 4, 1> H_d, H_u, F_d, F_u, fnu_u;
        tensor::symmetric2<CCTK_REAL, 4, 2> P_dd;
        tensor::generic<CCTK_REAL, 4, 2> P_ud;

        for (int ig = 0; ig < groupspec; ++ig) {
          const int i4D = layout_cc.linear(i, j, k, ig);

          // Assemble radiation tensors at cell centre
          pack_F_d(betax[idx], betay[idx], betaz[idx],
                   rFx[i4D], rFy[i4D], rFz[i4D], &F_d);
          pack_H_d(rHt[i4D], rHx[i4D], rHy[i4D], rHz[i4D], &H_d);
          pack_P_dd(betax[idx], betay[idx], betaz[idx],
                    rPxx[i4D], rPxy[i4D], rPxz[i4D],
                    rPyy[i4D], rPyz[i4D], rPzz[i4D], &P_dd);

          tensor::contract(g_uu, H_d, &H_u);
          tensor::contract(g_uu, F_d, &F_u);
          tensor::contract(g_uu, P_dd, &P_ud);

          assemble_fnu(u_u, rJ[i4D], H_u, &fnu_u, rad_E_floor);

          const CCTK_REAL Gamma =
              compute_Gamma(fidu_w_lorentz[idx], v_u,
                            rJ[i4D], rE[i4D], F_d,
                            rad_E_floor, rad_eps);

          const CCTK_REAL nnu = rN[i4D] / Gamma; //fmax(Gamma, 1.0);

          // Store physical flux components at this cell
          for (int iv = 0; iv < 5; ++iv) {
            const int comp      = PINDEX1D(ig, iv);
            const int idx_flux  = layout_cc.linear(i, j, k, comp);

            CCTK_REAL fval = 0.0;

            switch (iv) {
            case 0: // number density flux
              fval = alpha * nnu * fnu_u(dir + 1);
              break;
            case 1: // momentum flux (x-like component)
              fval = calc_F_flux(alpha, beta_u, F_d, P_ud, dir + 1, 1);
              break;
            case 2: // momentum flux (y-like component)
              fval = calc_F_flux(alpha, beta_u, F_d, P_ud, dir + 1, 2);
              break;
            case 3: // momentum flux (z-like component)
              fval = calc_F_flux(alpha, beta_u, F_d, P_ud, dir + 1, 3);
              break;
            case 4: // energy flux
              fval = calc_E_flux(alpha, beta_u, rE[i4D], F_u, dir + 1);
              break;
            }

            nu_flux_dir[idx_flux] = fval;
          }
        } // ig
      });
}

//---------------------------------------------------------------------
// (B) RHS update kernel (cell-centred numerical flux & divergence)
//
// Uses a 4-point stencil along direction `dir`:
//
//   LL, L | R, RR  = u_{j-1}, u_j, u_{j+1}, u_{j+2}
//
// together with physical fluxes at cell centres to build a limited,
// Lax–Friedrichs–type numerical flux at the left and right interfaces
// of each cell. The divergence (F_L − F_R)/dx is then added to the RHS.
//
// This mimics the original 1D sweep algorithm, but recast into a
// per-cell GPU kernel using local stencils (no dynamic allocations).
//---------------------------------------------------------------------
template <int dir>
void M1_UpdateRHSFromFluxes(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_M1_UpdateRHSFromFluxes;
  DECLARE_CCTK_PARAMETERS;

  const GridDescBaseDevice grid(cctkGH);
  const GF3D2layout layout_cc(cctkGH, {1, 1, 1});

  const CCTK_REAL dx  = CCTK_DELTA_SPACE(dir);
  const CCTK_REAL idx = 1.0 / dx;

  const int groupspec = ngroups * nspecies;

  // RHS components (cell-centred)
  CCTK_REAL *r_rhs[5] = {rN_rhs, rFx_rhs, rFy_rhs, rFz_rhs, rE_rhs};

  // Physical fluxes at cell centres in this direction
  const CCTK_REAL *nu_flux_dir =
      (dir == 0 ? nu_flux_x : (dir == 1 ? nu_flux_y : nu_flux_z));

  tensor::slicing_geometry_const geom(alp, betax, betay, betaz,
                                      gxx, gxy, gxz, gyy, gyz, gzz,
                                      kxx, kxy, kxz, kyy, kyz, kzz,
                                      volform);

  // THC needs at least 2 ghost zones along each direction
  assert(cctk_nghostzones[dir] >= 2);

  grid.loop_int_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        const int i = p.i;
        const int j = p.j;
        const int k = p.k;

        const int idxC = layout_cc.linear(i, j, k);

        // Optional mask: skip cells where transport should not act
        if (nuX_m1_mask[idxC]) {
          return;
        }

        // Directional unit offset (di,dj,dk)
        const int di = (dir == 0);
        const int dj = (dir == 1);
        const int dk = (dir == 2);

        auto idx_cell = [&](int ii, int jj, int kk) {
          return layout_cc.linear(ii, jj, kk);
        };

        // 1D indices around this cell along direction `dir`
        // Treat this cell as "j".
        const int i_jm2 = i - 2 * di;
        const int j_jm2 = j - 2 * dj;
        const int k_jm2 = k - 2 * dk;

        const int i_jm1 = i - 1 * di;
        const int j_jm1 = j - 1 * dj;
        const int k_jm1 = k - 1 * dk;

        const int i_jp1 = i + 1 * di;
        const int j_jp1 = j + 1 * dj;
        const int k_jp1 = k + 1 * dk;

        const int i_jp2 = i + 2 * di;
        const int j_jp2 = j + 2 * dj;
        const int k_jp2 = k + 2 * dk;

        const int idx_jm2 = idx_cell(i_jm2, j_jm2, k_jm2);
        const int idx_jm1 = idx_cell(i_jm1, j_jm1, k_jm1);
        const int idx_j   = idx_cell(i,      j,      k     );
        const int idx_jp1 = idx_cell(i_jp1, j_jp1, k_jp1);
        const int idx_jp2 = idx_cell(i_jp2, j_jp2, k_jp2);

        // Precompute characteristic speeds at j-1, j, j+1
        tensor::inv_metric<3> gamma_jm1, gamma_j, gamma_jp1;
        tensor::generic<CCTK_REAL, 4, 1> beta_jm1, beta_j, beta_jp1;

        geom.get_inv_metric(idx_jm1, &gamma_jm1);
        geom.get_shift_vec(idx_jm1, &beta_jm1);

        geom.get_inv_metric(idx_j, &gamma_j);
        geom.get_shift_vec(idx_j, &beta_j);

        geom.get_inv_metric(idx_jp1, &gamma_jp1);
        geom.get_shift_vec(idx_jp1, &beta_jp1);

        const CCTK_REAL alpha_jm1 = alp[idx_jm1];
        const CCTK_REAL alpha_j   = alp[idx_j];
        const CCTK_REAL alpha_jp1 = alp[idx_jp1];

        const CCTK_REAL clight_jm1 =
            face_speed(dir, gamma_jm1, beta_jm1, alpha_jm1);
        const CCTK_REAL clight_j   =
            face_speed(dir, gamma_j,   beta_j,   alpha_j);
        const CCTK_REAL clight_jp1 =
            face_speed(dir, gamma_jp1, beta_jp1, alpha_jp1);

        // Loop over groups/species and conserved components
        for (int ig = 0; ig < groupspec; ++ig) {

          //------------------------------------------------------------------
          // Opacity factors for right (j | j+1) and left (j-1 | j) interfaces
          //------------------------------------------------------------------

          // Right interface j+1/2: use cells j and j+1
          const int idx4_j   = layout_cc.linear(i,      j,      k,      ig);
          const int idx4_jp1 = layout_cc.linear(i_jp1, j_jp1, k_jp1, ig);

          const CCTK_REAL kappa_j   = abs_1[idx4_j]   + scat_1[idx4_j];
          const CCTK_REAL kappa_jp1 = abs_1[idx4_jp1] + scat_1[idx4_jp1];

          const CCTK_REAL kapa_R = 0.5 * (kappa_j + kappa_jp1);

          CCTK_REAL A_R = 1.0;
          if (kapa_R * dx > 1.0) {
            A_R = 1.0 / (dx * kapa_R);
            if (A_R < mindiss) A_R = mindiss;
            if (A_R > 1.0)     A_R = 1.0;
          }

          // Left interface j-1/2: use cells j-1 and j
          const int idx4_jm1 = layout_cc.linear(i_jm1, j_jm1, k_jm1, ig);

          const CCTK_REAL kappa_jm1 = abs_1[idx4_jm1] + scat_1[idx4_jm1];
          const CCTK_REAL kapa_L    = 0.5 * (kappa_jm1 + kappa_j);

          CCTK_REAL A_L = 1.0;
          if (kapa_L * dx > 1.0) {
            A_L = 1.0 / (dx * kapa_L);
            if (A_L < mindiss) A_L = mindiss;
            if (A_L > 1.0)     A_L = 1.0;
          }

          for (int iv = 0; iv < 5; ++iv) {
            const int comp = PINDEX1D(ig, iv);

            //----------------------------------------------------------------
            // Conserved variables on stencil
            //
            // For right interface j+1/2:
            //   {u_{j-1}, u_j, u_{j+1}, u_{j+2}}
            //
            // NOTE: These must be indexed with (i,j,k,ig), matching THC.
            //----------------------------------------------------------------
            const int idx4_jm2 = layout_cc.linear(i_jm2, j_jm2, k_jm2, ig);
            const int idx4_j   = layout_cc.linear(i,      j,      k,      ig);
            const int idx4_jp2 = layout_cc.linear(i_jp2, j_jp2, k_jp2, ig);

            const CCTK_REAL u_jm1 =
                read_cons_by_iv(iv, rN, rFx, rFy, rFz, rE, idx4_jm1);
            const CCTK_REAL u_j =
                read_cons_by_iv(iv, rN, rFx, rFy, rFz, rE, idx4_j);
            const CCTK_REAL u_jp1 =
                read_cons_by_iv(iv, rN, rFx, rFy, rFz, rE, idx4_jp1);
            const CCTK_REAL u_jp2 =
                read_cons_by_iv(iv, rN, rFx, rFy, rFz, rE, idx4_jp2);

            //----------------------------------------------------------------
            // 1) Right interface F_R at j+1/2 (between j and j+1)
            //    This exactly mirrors the original scheme with:
            //
            //    ujm  = u_{j-1}, uj = u_j, ujp = u_{j+1}, ujpp = u_{j+2}
            //----------------------------------------------------------------
            const CCTK_REAL dup_R = u_jp2 - u_jp1;  // ujpp - ujp
            const CCTK_REAL duc_R = u_jp1 - u_j;    // ujp  - uj
            const CCTK_REAL dum_R = u_j   - u_jm1;  // uj   - ujm

            CCTK_REAL phi_R = 0.0;
            bool saw_R      = false;

            if (dup_R * duc_R > 0.0 && dum_R * duc_R > 0.0 && duc_R != 0.0) {
              phi_R = minmod2(dum_R / duc_R, dup_R / duc_R, minmod_theta);
            } else if (dup_R * duc_R < 0.0 && dum_R * duc_R < 0.0) {
              saw_R = true;
            }

            // Physical fluxes at j and j+1
            const CCTK_REAL f_j =
                nu_flux_dir[layout_cc.linear(i,      j,      k,      comp)];
            const CCTK_REAL f_jp1 =
                nu_flux_dir[layout_cc.linear(i_jp1, j_jp1, k_jp1, comp)];

            // Characteristic speed at the interface: max(c_j, c_{j+1})
            const CCTK_REAL cmx_R =
                fmax(clight_j, clight_jp1);

            const CCTK_REAL flux_low_R =
                0.5 * (f_j + f_jp1 - cmx_R * (u_jp1 - u_j));
            const CCTK_REAL flux_high_R =
                0.5 * (f_j + f_jp1);

            const CCTK_REAL Aeff_R = (saw_R ? 1.0 : A_R);

            const CCTK_REAL F_R =
                flux_high_R - Aeff_R * (1.0 - phi_R) *
                              (flux_high_R - flux_low_R);

            //----------------------------------------------------------------
            // 2) Left interface F_L at j-1/2 (between j-1 and j)
            //
            //    In the original algorithm this is computed one step back
            //    using the stencil {u_{j-2}, u_{j-1}, u_j, u_{j+1}}.
            //----------------------------------------------------------------
            const CCTK_REAL u_jm2 =
                read_cons_by_iv(iv, rN, rFx, rFy, rFz, rE, idx4_jm2);

            const CCTK_REAL dup_L = u_jp1 - u_j;     // ujpp - ujp (j -> j-1 shift)
            const CCTK_REAL duc_L = u_j   - u_jm1;   // ujp  - uj
            const CCTK_REAL dum_L = u_jm1 - u_jm2;   // uj   - ujm

            CCTK_REAL phi_L = 0.0;
            bool saw_L      = false;

            if (dup_L * duc_L > 0.0 && dum_L * duc_L > 0.0 && duc_L != 0.0) {
              phi_L = minmod2(dum_L / duc_L, dup_L / duc_L, minmod_theta);
            } else if (dup_L * duc_L < 0.0 && dum_L * duc_L < 0.0) {
              saw_L = true;
            }

            // Physical fluxes at j-1 and j
            const CCTK_REAL f_jm1 =
                nu_flux_dir[layout_cc.linear(i_jm1, j_jm1, k_jm1, comp)];
            const CCTK_REAL f_jL =
                nu_flux_dir[layout_cc.linear(i,      j,      k,      comp)];

            // Characteristic speed at left interface: max(c_{j-1}, c_j)
            const CCTK_REAL cmx_L =
                fmax(clight_jm1, clight_j);

            const CCTK_REAL flux_low_L =
                0.5 * (f_jm1 + f_jL - cmx_L * (u_j - u_jm1));
            const CCTK_REAL flux_high_L =
                0.5 * (f_jm1 + f_jL);

            const CCTK_REAL Aeff_L = (saw_L ? 1.0 : A_L);

            const CCTK_REAL F_L =
                flux_high_L - Aeff_L * (1.0 - phi_L) *
                              (flux_high_L - flux_low_L);

            //----------------------------------------------------------------
            // Discrete divergence contribution for this cell/component
            //----------------------------------------------------------------
            const int idx_rhs = layout_cc.linear(i, j, k, ig);
            r_rhs[iv][idx_rhs] += idx * (F_L - F_R);
          } // iv
        }   // ig
      });
}

//---------------------------------------------------------------------
// Wrappers
//---------------------------------------------------------------------
extern "C" void nuX_M1_CalcFluxes(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_M1_CalcFluxes;
  DECLARE_CCTK_PARAMETERS;

  if (verbose) {
    CCTK_INFO("nuX_M1_CalcFluxes");
  }

  // Compute cell-centred physical fluxes in each direction
  M1_ComputePhysicalFluxes<0>(cctkGH);
  M1_ComputePhysicalFluxes<1>(cctkGH);
  M1_ComputePhysicalFluxes<2>(cctkGH);
}

extern "C" void nuX_M1_UpdateRHSFromFluxes(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_M1_UpdateRHSFromFluxes;
  DECLARE_CCTK_PARAMETERS;

  if (verbose) {
    CCTK_INFO("nuX_M1_UpdateRHSFromFluxes");
  }

  // Build limited numerical fluxes and apply divergence to the RHS
  M1_UpdateRHSFromFluxes<0>(cctkGH);
  M1_UpdateRHSFromFluxes<1>(cctkGH);
  M1_UpdateRHSFromFluxes<2>(cctkGH);
}

} // namespace nuX_M1

