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
CCTK_HOST CCTK_DEVICE inline CCTK_REAL minmod2(CCTK_REAL rl, CCTK_REAL rp,
                                               CCTK_REAL th) {
  return fmin(1.0, fmin(th * rl, th * rp));
}

//-------------------------------------------------------------------------
// Characteristic speed in a given spatial direction:
//
//   lambda = -beta_dir ± alpha sqrt(gamma^{dir,dir})
//
// dir = 0,1,2 corresponds to x,y,z direction.
//-------------------------------------------------------------------------
CCTK_HOST CCTK_DEVICE inline CCTK_REAL
face_speed(int dir, const tensor::inv_metric<3> &gamma,
           const tensor::generic<CCTK_REAL, 4, 1> &beta_u,
           CCTK_REAL alpha_cell) {

  // NOTE: gamma is a 3-metric (indices 0..2) like in THC.
  const CCTK_REAL root = alpha_cell * sqrt(gamma(dir, dir));

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
template <int dir> void M1_ComputePhysicalFluxes(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_M1_CalcFluxes;
  DECLARE_CCTK_PARAMETERS;

  const GridDescBaseDevice grid(cctkGH);
  const GF3D2layout layout_cc(cctkGH, {1, 1, 1});
  const GF3D2layout layout_vc(cctkGH, {0, 0, 0});

  // These are treated as *cell-centred* storage for physical fluxes
  CCTK_REAL *restrict nu_flux_dir =
      (dir == 0 ? nu_flux_x : (dir == 1 ? nu_flux_y : nu_flux_z));
  CCTK_REAL *restrict nu_cmax_dir =
      (dir == 0 ? nu_cmax_x : (dir == 1 ? nu_cmax_y : nu_cmax_z));

  tensor::slicing_geometry_const geom(layout_vc, layout_cc, alp, betax, betay,
                                      betaz, gxx, gxy, gxz, gyy, gyz, gzz, kxx,
                                      kxy, kxz, kyy, kyz, kzz);

  tensor::fluid_velocity_field_const fidu(layout_vc, layout_cc, alp, betax,
                                          betay, betaz, fidu_w_lorentz,
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
      grid.nghostzones, [=] CCTK_DEVICE(const PointDesc &p) {
        const int i = p.i;
        const int j = p.j;
        const int k = p.k;

        const int idx = layout_cc.linear(i, j, k);

        // Local geometry at this cell
        tensor::inv_metric<3> gamma_uu;
        tensor::inv_metric<4> g_uu;
        tensor::generic<CCTK_REAL, 4, 1> beta_u;

        geom.get_inv_metric(p, &gamma_uu);
        geom.get_inv_metric(p, &g_uu);
        geom.get_shift_vec(p, &beta_u);

        const CCTK_REAL alpha = geom.get_lapse(p);
        const CCTK_REAL beta_x = beta_u(1);
        const CCTK_REAL beta_y = beta_u(2);
        const CCTK_REAL beta_z = beta_u(3);
        const CCTK_REAL W = fidu_w_lorentz[idx];
        const CCTK_REAL vel_x = fidu_velx[idx];
        const CCTK_REAL vel_y = fidu_vely[idx];
        const CCTK_REAL vel_z = fidu_velz[idx];

        // Fiducial 4-velocity and 3-velocity at this cell
        tensor::generic<CCTK_REAL, 4, 1> u_u, v_u;
        fidu.get(p, &u_u);
        pack_v_u(vel_x, vel_y, vel_z, &v_u);

        // Temporary tensors used for flux construction
        tensor::generic<CCTK_REAL, 4, 1> H_d, H_u, F_d, F_u, fnu_u;
        tensor::symmetric2<CCTK_REAL, 4, 2> P_dd;
        tensor::generic<CCTK_REAL, 4, 2> P_ud;

        for (int ig = 0; ig < groupspec; ++ig) {
          const int i4D = layout_cc.linear(i, j, k, ig);

          // Assemble radiation tensors at cell centre
          pack_F_d(beta_x, beta_y, beta_z, rFx_transport_input[i4D],
                   rFy_transport_input[i4D], rFz_transport_input[i4D], &F_d);
          pack_H_d(rHt[i4D], rHx[i4D], rHy[i4D], rHz[i4D], &H_d);
          pack_P_dd(beta_x, beta_y, beta_z, rPxx[i4D], rPxy[i4D], rPxz[i4D],
                    rPyy[i4D], rPyz[i4D], rPzz[i4D], &P_dd);

          tensor::contract(g_uu, H_d, &H_u);
          tensor::contract(g_uu, F_d, &F_u);
          tensor::contract(g_uu, P_dd, &P_ud);

          assemble_fnu(u_u, rJ[i4D], H_u, &fnu_u, rad_E_floor);

          const CCTK_REAL Gamma =
              compute_Gamma(W, v_u, rJ[i4D], rE_transport_input[i4D], F_d,
                            rad_E_floor, rad_eps);

          const CCTK_REAL nnu = rN_transport_input[i4D] / Gamma; // fmax(Gamma, 1.0);

          const int cmax_idx = layout_cc.linear(i, j, k, ig);
          const int comp_n = PINDEX1D(ig, 0);
          const int comp_fx = PINDEX1D(ig, 1);
          const int comp_fy = PINDEX1D(ig, 2);
          const int comp_fz = PINDEX1D(ig, 3);
          const int comp_e = PINDEX1D(ig, 4);

          nu_flux_dir[layout_cc.linear(i, j, k, comp_n)] =
              alpha * nnu * fnu_u(dir + 1);
          nu_flux_dir[layout_cc.linear(i, j, k, comp_fx)] =
              calc_F_flux(alpha, beta_u, F_d, P_ud, dir + 1, 1);
          nu_flux_dir[layout_cc.linear(i, j, k, comp_fy)] =
              calc_F_flux(alpha, beta_u, F_d, P_ud, dir + 1, 2);
          nu_flux_dir[layout_cc.linear(i, j, k, comp_fz)] =
              calc_F_flux(alpha, beta_u, F_d, P_ud, dir + 1, 3);
          nu_flux_dir[layout_cc.linear(i, j, k, comp_e)] =
              calc_E_flux(alpha, beta_u, rE_transport_input[i4D], F_u, dir + 1);

          nu_cmax_dir[cmax_idx] = face_speed(dir, gamma_uu, beta_u, alpha);
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
template <int dir> void M1_UpdateRHSFromFluxes(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_M1_UpdateRHSFromFluxes;
  DECLARE_CCTK_PARAMETERS;

  const GridDescBaseDevice grid(cctkGH);
  const GF3D2layout layout_cc(cctkGH, {1, 1, 1});

  const CCTK_REAL dx = CCTK_DELTA_SPACE(dir);
  const CCTK_REAL idx = 1.0 / dx;

  const int groupspec = ngroups * nspecies;

  // RHS components (cell-centred)
  CCTK_REAL *r_rhs[5] = {rN_rhs, rFx_rhs, rFy_rhs, rFz_rhs, rE_rhs};
  CCTK_REAL *diag_cmax_L[3] = {transport_cmax_Lx, transport_cmax_Ly,
                               transport_cmax_Lz};
  CCTK_REAL *diag_cmax_R[3] = {transport_cmax_Rx, transport_cmax_Ry,
                               transport_cmax_Rz};
  CCTK_REAL *diag_A_L[3] = {transport_A_Lx, transport_A_Ly, transport_A_Lz};
  CCTK_REAL *diag_A_R[3] = {transport_A_Rx, transport_A_Ry, transport_A_Rz};
  CCTK_REAL *diag_FN_L[3] = {transport_FN_Lx, transport_FN_Ly,
                             transport_FN_Lz};
  CCTK_REAL *diag_FN_R[3] = {transport_FN_Rx, transport_FN_Ry,
                             transport_FN_Rz};
  CCTK_REAL *diag_FE_L[3] = {transport_FE_Lx, transport_FE_Ly,
                             transport_FE_Lz};
  CCTK_REAL *diag_FE_R[3] = {transport_FE_Rx, transport_FE_Ry,
                             transport_FE_Rz};
  CCTK_REAL *diag_favgN_L[3] = {transport_favgN_Lx, transport_favgN_Ly,
                                transport_favgN_Lz};
  CCTK_REAL *diag_favgN_R[3] = {transport_favgN_Rx, transport_favgN_Ry,
                                transport_favgN_Rz};
  CCTK_REAL *diag_favgE_L[3] = {transport_favgE_Lx, transport_favgE_Ly,
                                transport_favgE_Lz};
  CCTK_REAL *diag_favgE_R[3] = {transport_favgE_Rx, transport_favgE_Ry,
                                transport_favgE_Rz};
  CCTK_REAL *diag_fdissN_L[3] = {transport_fdissN_Lx, transport_fdissN_Ly,
                                 transport_fdissN_Lz};
  CCTK_REAL *diag_fdissN_R[3] = {transport_fdissN_Rx, transport_fdissN_Ry,
                                 transport_fdissN_Rz};
  CCTK_REAL *diag_fdissE_L[3] = {transport_fdissE_Lx, transport_fdissE_Ly,
                                 transport_fdissE_Lz};
  CCTK_REAL *diag_fdissE_R[3] = {transport_fdissE_Rx, transport_fdissE_Ry,
                                 transport_fdissE_Rz};
  CCTK_REAL *diag_duN_L[3] = {transport_duN_Lx, transport_duN_Ly,
                              transport_duN_Lz};
  CCTK_REAL *diag_duN_R[3] = {transport_duN_Rx, transport_duN_Ry,
                              transport_duN_Rz};
  CCTK_REAL *diag_duE_L[3] = {transport_duE_Lx, transport_duE_Ly,
                              transport_duE_Lz};
  CCTK_REAL *diag_duE_R[3] = {transport_duE_Rx, transport_duE_Ry,
                              transport_duE_Rz};
  CCTK_REAL *diag_fphysN[3][3] = {
      {transport_fphysN_Lx, transport_fphysN_Cx, transport_fphysN_Rx},
      {transport_fphysN_Ly, transport_fphysN_Cy, transport_fphysN_Ry},
      {transport_fphysN_Lz, transport_fphysN_Cz, transport_fphysN_Rz},
  };
  CCTK_REAL *diag_fphysE[3][3] = {
      {transport_fphysE_Lx, transport_fphysE_Cx, transport_fphysE_Rx},
      {transport_fphysE_Ly, transport_fphysE_Cy, transport_fphysE_Ry},
      {transport_fphysE_Lz, transport_fphysE_Cz, transport_fphysE_Rz},
  };
  CCTK_REAL *diag_flowN_L[3] = {transport_flowN_Lx, transport_flowN_Ly,
                                transport_flowN_Lz};
  CCTK_REAL *diag_flowN_R[3] = {transport_flowN_Rx, transport_flowN_Ry,
                                transport_flowN_Rz};
  CCTK_REAL *diag_flowE_L[3] = {transport_flowE_Lx, transport_flowE_Ly,
                                transport_flowE_Lz};
  CCTK_REAL *diag_flowE_R[3] = {transport_flowE_Rx, transport_flowE_Ry,
                                transport_flowE_Rz};
  CCTK_REAL *diag_fhighN_L[3] = {transport_fhighN_Lx, transport_fhighN_Ly,
                                 transport_fhighN_Lz};
  CCTK_REAL *diag_fhighN_R[3] = {transport_fhighN_Rx, transport_fhighN_Ry,
                                 transport_fhighN_Rz};
  CCTK_REAL *diag_fhighE_L[3] = {transport_fhighE_Lx, transport_fhighE_Ly,
                                 transport_fhighE_Lz};
  CCTK_REAL *diag_fhighE_R[3] = {transport_fhighE_Rx, transport_fhighE_Ry,
                                 transport_fhighE_Rz};
  CCTK_REAL *diag_phiN_L[3] = {transport_phiN_Lx, transport_phiN_Ly,
                               transport_phiN_Lz};
  CCTK_REAL *diag_phiN_R[3] = {transport_phiN_Rx, transport_phiN_Ry,
                               transport_phiN_Rz};
  CCTK_REAL *diag_phiE_L[3] = {transport_phiE_Lx, transport_phiE_Ly,
                               transport_phiE_Lz};
  CCTK_REAL *diag_phiE_R[3] = {transport_phiE_Rx, transport_phiE_Ry,
                               transport_phiE_Rz};
  CCTK_REAL *diag_sawN_L[3] = {transport_sawN_Lx, transport_sawN_Ly,
                               transport_sawN_Lz};
  CCTK_REAL *diag_sawN_R[3] = {transport_sawN_Rx, transport_sawN_Ry,
                               transport_sawN_Rz};
  CCTK_REAL *diag_sawE_L[3] = {transport_sawE_Lx, transport_sawE_Ly,
                               transport_sawE_Lz};
  CCTK_REAL *diag_sawE_R[3] = {transport_sawE_Rx, transport_sawE_Ry,
                               transport_sawE_Rz};
  CCTK_REAL *diag_divN[3] = {transport_divN_x, transport_divN_y,
                             transport_divN_z};
  CCTK_REAL *diag_divE[3] = {transport_divE_x, transport_divE_y,
                             transport_divE_z};

  // Physical fluxes at cell centres in this direction
  const CCTK_REAL *nu_flux_dir =
      (dir == 0 ? nu_flux_x : (dir == 1 ? nu_flux_y : nu_flux_z));
  const CCTK_REAL *nu_cmax_dir =
      (dir == 0 ? nu_cmax_x : (dir == 1 ? nu_cmax_y : nu_cmax_z));

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

        // Loop over groups/species and conserved components
        for (int ig = 0; ig < groupspec; ++ig) {

          //------------------------------------------------------------------
          // Opacity factors for right (j | j+1) and left (j-1 | j) interfaces
          //------------------------------------------------------------------

          // Right interface j+1/2: use cells j and j+1
          const int idx4_j = layout_cc.linear(i, j, k, ig);
          const int idx_rhs = idx4_j;
          const int idx4_jm2 = layout_cc.linear(i_jm2, j_jm2, k_jm2, ig);
          const int idx4_jm1 = layout_cc.linear(i_jm1, j_jm1, k_jm1, ig);
          const int idx4_jp1 = layout_cc.linear(i_jp1, j_jp1, k_jp1, ig);
          const int idx4_jp2 = layout_cc.linear(i_jp2, j_jp2, k_jp2, ig);

          const CCTK_REAL kappa_j = abs_1[idx4_j] + scat_1[idx4_j];
          const CCTK_REAL kappa_jp1 = abs_1[idx4_jp1] + scat_1[idx4_jp1];

          const CCTK_REAL kapa_R = 0.5 * (kappa_j + kappa_jp1);

          CCTK_REAL A_R = 1.0;
          if (kapa_R * dx > 1.0) {
            A_R = 1.0 / (dx * kapa_R);
            if (A_R < mindiss)
              A_R = mindiss;
            if (A_R > 1.0)
              A_R = 1.0;
          }

          // Left interface j-1/2: use cells j-1 and j
          const CCTK_REAL kappa_jm1 = abs_1[idx4_jm1] + scat_1[idx4_jm1];
          const CCTK_REAL kapa_L = 0.5 * (kappa_jm1 + kappa_j);

          CCTK_REAL A_L = 1.0;
          if (kapa_L * dx > 1.0) {
            A_L = 1.0 / (dx * kapa_L);
            if (A_L < mindiss)
              A_L = mindiss;
            if (A_L > 1.0)
              A_L = 1.0;
          }

          const CCTK_REAL c_jm1 =
              nu_cmax_dir[layout_cc.linear(i_jm1, j_jm1, k_jm1, ig)];
          const CCTK_REAL c_j = nu_cmax_dir[idx4_j];
          const CCTK_REAL c_jp1 =
              nu_cmax_dir[layout_cc.linear(i_jp1, j_jp1, k_jp1, ig)];

          const CCTK_REAL cmx_R = fmax(c_j, c_jp1);
          const CCTK_REAL cmx_L = fmax(c_jm1, c_j);

          if (transport_flux_diagnostics) {
            diag_cmax_L[dir][idx4_j] = cmx_L;
            diag_cmax_R[dir][idx4_j] = cmx_R;
            diag_A_L[dir][idx4_j] = A_L;
            diag_A_R[dir][idx4_j] = A_R;
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
            const CCTK_REAL *const u_cons =
                (iv == 0)   ? rN_transport_input
                : (iv == 1) ? rFx_transport_input
                : (iv == 2) ? rFy_transport_input
                : (iv == 3) ? rFz_transport_input
                            : rE_transport_input;

            const CCTK_REAL u_jm1 = u_cons[idx4_jm1];
            const CCTK_REAL u_j = u_cons[idx4_j];
            const CCTK_REAL u_jp1 = u_cons[idx4_jp1];
            const CCTK_REAL u_jp2 = u_cons[idx4_jp2];
            const CCTK_REAL u_jm2 = u_cons[idx4_jm2];
            const CCTK_REAL f_jm1 =
                nu_flux_dir[layout_cc.linear(i_jm1, j_jm1, k_jm1, comp)];
            const CCTK_REAL f_j =
                nu_flux_dir[layout_cc.linear(i, j, k, comp)];
            const CCTK_REAL f_jp1 =
                nu_flux_dir[layout_cc.linear(i_jp1, j_jp1, k_jp1, comp)];

            CCTK_REAL F_R;
            CCTK_REAL phi_R = 0.0;
            bool saw_R = false;
            CCTK_REAL flux_low_R = 0.0;
            CCTK_REAL flux_high_R = 0.0;
            CCTK_REAL flux_avg_R = 0.0;
            CCTK_REAL flux_diss_R = 0.0;
            CCTK_REAL du_state_R = 0.0;
            {
              //----------------------------------------------------------------
              // 1) Right interface F_R at j+1/2 (between j and j+1)
              //    This exactly mirrors the original scheme with:
              //
              //    ujm  = u_{j-1}, uj = u_j, ujp = u_{j+1}, ujpp = u_{j+2}
              //----------------------------------------------------------------
              const CCTK_REAL dup_R = u_jp2 - u_jp1; // ujpp - ujp
              const CCTK_REAL duc_R = u_jp1 - u_j;   // ujp  - uj
              const CCTK_REAL dum_R = u_j - u_jm1;   // uj   - ujm
              du_state_R = duc_R;

              if (dup_R * duc_R > 0.0 && dum_R * duc_R > 0.0 && duc_R != 0.0) {
                phi_R = minmod2(dum_R / duc_R, dup_R / duc_R, minmod_theta);
              } else if (dup_R * duc_R < 0.0 && dum_R * duc_R < 0.0) {
                saw_R = true;
              }

              flux_avg_R = 0.5 * (f_j + f_jp1);
              flux_diss_R = 0.5 * cmx_R * du_state_R;
              flux_low_R = flux_avg_R - flux_diss_R;
              flux_high_R = 0.5 * (f_j + f_jp1);

              const CCTK_REAL Aeff_R = (saw_R ? 1.0 : A_R);

              F_R = flux_high_R -
                    Aeff_R * (1.0 - phi_R) * (flux_high_R - flux_low_R);
            }

            //----------------------------------------------------------------
            // 2) Left interface F_L at j-1/2 (between j-1 and j)
            //
            //    In the original algorithm this is computed one step back
            //    using the stencil {u_{j-2}, u_{j-1}, u_j, u_{j+1}}.
            //----------------------------------------------------------------
            CCTK_REAL F_L;
            CCTK_REAL phi_L = 0.0;
            bool saw_L = false;
            CCTK_REAL flux_low_L = 0.0;
            CCTK_REAL flux_high_L = 0.0;
            CCTK_REAL flux_avg_L = 0.0;
            CCTK_REAL flux_diss_L = 0.0;
            CCTK_REAL du_state_L = 0.0;
            {
              const CCTK_REAL dup_L =
                  u_jp1 - u_j; // ujpp - ujp (j -> j-1 shift)
              const CCTK_REAL duc_L = u_j - u_jm1;   // ujp  - uj
              const CCTK_REAL dum_L = u_jm1 - u_jm2; // uj   - ujm
              du_state_L = duc_L;

              if (dup_L * duc_L > 0.0 && dum_L * duc_L > 0.0 && duc_L != 0.0) {
                phi_L = minmod2(dum_L / duc_L, dup_L / duc_L, minmod_theta);
              } else if (dup_L * duc_L < 0.0 && dum_L * duc_L < 0.0) {
                saw_L = true;
              }

              flux_avg_L = 0.5 * (f_jm1 + f_j);
              flux_diss_L = 0.5 * cmx_L * du_state_L;
              flux_low_L = flux_avg_L - flux_diss_L;
              flux_high_L = 0.5 * (f_jm1 + f_j);

              const CCTK_REAL Aeff_L = (saw_L ? 1.0 : A_L);

              F_L = flux_high_L -
                    Aeff_L * (1.0 - phi_L) * (flux_high_L - flux_low_L);
            }

            //----------------------------------------------------------------
            // Discrete divergence contribution for this cell/component
            //----------------------------------------------------------------
            const CCTK_REAL div_contrib = idx * (F_L - F_R);
            r_rhs[iv][idx_rhs] += div_contrib;

            if (transport_flux_diagnostics) {
              if (iv == 0) {
                diag_fphysN[dir][0][idx4_j] = f_jm1;
                diag_fphysN[dir][1][idx4_j] = f_j;
                diag_fphysN[dir][2][idx4_j] = f_jp1;
                diag_favgN_L[dir][idx4_j] = flux_avg_L;
                diag_favgN_R[dir][idx4_j] = flux_avg_R;
                diag_fdissN_L[dir][idx4_j] = flux_diss_L;
                diag_fdissN_R[dir][idx4_j] = flux_diss_R;
                diag_duN_L[dir][idx4_j] = du_state_L;
                diag_duN_R[dir][idx4_j] = du_state_R;
                diag_flowN_L[dir][idx4_j] = flux_low_L;
                diag_flowN_R[dir][idx4_j] = flux_low_R;
                diag_fhighN_L[dir][idx4_j] = flux_high_L;
                diag_fhighN_R[dir][idx4_j] = flux_high_R;
                diag_phiN_L[dir][idx4_j] = phi_L;
                diag_phiN_R[dir][idx4_j] = phi_R;
                diag_sawN_L[dir][idx4_j] = saw_L ? 1.0 : 0.0;
                diag_sawN_R[dir][idx4_j] = saw_R ? 1.0 : 0.0;
                diag_FN_L[dir][idx4_j] = F_L;
                diag_FN_R[dir][idx4_j] = F_R;
                diag_divN[dir][idx4_j] = div_contrib;
              } else if (iv == 4) {
                diag_fphysE[dir][0][idx4_j] = f_jm1;
                diag_fphysE[dir][1][idx4_j] = f_j;
                diag_fphysE[dir][2][idx4_j] = f_jp1;
                diag_favgE_L[dir][idx4_j] = flux_avg_L;
                diag_favgE_R[dir][idx4_j] = flux_avg_R;
                diag_fdissE_L[dir][idx4_j] = flux_diss_L;
                diag_fdissE_R[dir][idx4_j] = flux_diss_R;
                diag_duE_L[dir][idx4_j] = du_state_L;
                diag_duE_R[dir][idx4_j] = du_state_R;
                diag_flowE_L[dir][idx4_j] = flux_low_L;
                diag_flowE_R[dir][idx4_j] = flux_low_R;
                diag_fhighE_L[dir][idx4_j] = flux_high_L;
                diag_fhighE_R[dir][idx4_j] = flux_high_R;
                diag_phiE_L[dir][idx4_j] = phi_L;
                diag_phiE_R[dir][idx4_j] = phi_R;
                diag_sawE_L[dir][idx4_j] = saw_L ? 1.0 : 0.0;
                diag_sawE_R[dir][idx4_j] = saw_R ? 1.0 : 0.0;
                diag_FE_L[dir][idx4_j] = F_L;
                diag_FE_R[dir][idx4_j] = F_R;
                diag_divE[dir][idx4_j] = div_contrib;

                if (dir == 2 && ig == 0) {
                  const int idx_trace = layout_cc.linear(p.i, p.j, p.k, 0);
                  transport_uE_jm2_z[idx_trace] = u_jm2;
                  transport_uE_jm1_z[idx_trace] = u_jm1;
                  transport_uE_j_z[idx_trace] = u_j;
                  transport_uE_jp1_z[idx_trace] = u_jp1;
                  transport_uE_jp2_z[idx_trace] = u_jp2;
                  transport_idxE_jm2_z[idx_trace] = static_cast<CCTK_REAL>(idx4_jm2);
                  transport_idxE_jm1_z[idx_trace] = static_cast<CCTK_REAL>(idx4_jm1);
                  transport_idxE_j_z[idx_trace] = static_cast<CCTK_REAL>(idx4_j);
                  transport_idxE_jp1_z[idx_trace] = static_cast<CCTK_REAL>(idx4_jp1);
                  transport_idxE_jp2_z[idx_trace] = static_cast<CCTK_REAL>(idx4_jp2);
                }
              }
            }
          } // iv
        } // ig
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
  if (cctk_nghostzones[0] < 2 || cctk_nghostzones[1] < 2 ||
      cctk_nghostzones[2] < 2) {
    CCTK_VERROR("nuX_M1_UpdateRHSFromFluxes requires at least 2 ghost zones "
                "in every direction, got (%d, %d, %d)",
                int(cctk_nghostzones[0]), int(cctk_nghostzones[1]),
                int(cctk_nghostzones[2]));
  }

  // Build limited numerical fluxes and apply divergence to the RHS
  M1_UpdateRHSFromFluxes<0>(cctkGH);
  M1_UpdateRHSFromFluxes<1>(cctkGH);
  M1_UpdateRHSFromFluxes<2>(cctkGH);
}

} // namespace nuX_M1
