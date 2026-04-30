#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include "loop_device.hxx"
#include "vect.hxx"

namespace nuX_M1 {

using namespace Loop;

extern "C" void nuX_M1_InitFluxesRHS(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  const int groupspec = ngroups * nspecies;
  const int flux_comps = 5 * groupspec;
  CCTK_REAL *diag6[18][6] = {
      {transport_cmax_Lx, transport_cmax_Ly, transport_cmax_Lz,
       transport_cmax_Rx, transport_cmax_Ry, transport_cmax_Rz},
      {transport_A_Lx, transport_A_Ly, transport_A_Lz, transport_A_Rx,
       transport_A_Ry, transport_A_Rz},
      {transport_FN_Lx, transport_FN_Ly, transport_FN_Lz, transport_FN_Rx,
       transport_FN_Ry, transport_FN_Rz},
      {transport_FE_Lx, transport_FE_Ly, transport_FE_Lz, transport_FE_Rx,
       transport_FE_Ry, transport_FE_Rz},
      {transport_favgN_Lx, transport_favgN_Ly, transport_favgN_Lz,
       transport_favgN_Rx, transport_favgN_Ry, transport_favgN_Rz},
      {transport_favgE_Lx, transport_favgE_Ly, transport_favgE_Lz,
       transport_favgE_Rx, transport_favgE_Ry, transport_favgE_Rz},
      {transport_fdissN_Lx, transport_fdissN_Ly, transport_fdissN_Lz,
       transport_fdissN_Rx, transport_fdissN_Ry, transport_fdissN_Rz},
      {transport_fdissE_Lx, transport_fdissE_Ly, transport_fdissE_Lz,
       transport_fdissE_Rx, transport_fdissE_Ry, transport_fdissE_Rz},
      {transport_duN_Lx, transport_duN_Ly, transport_duN_Lz, transport_duN_Rx,
       transport_duN_Ry, transport_duN_Rz},
      {transport_duE_Lx, transport_duE_Ly, transport_duE_Lz, transport_duE_Rx,
       transport_duE_Ry, transport_duE_Rz},
      {transport_flowN_Lx, transport_flowN_Ly, transport_flowN_Lz,
       transport_flowN_Rx, transport_flowN_Ry, transport_flowN_Rz},
      {transport_flowE_Lx, transport_flowE_Ly, transport_flowE_Lz,
       transport_flowE_Rx, transport_flowE_Ry, transport_flowE_Rz},
      {transport_fhighN_Lx, transport_fhighN_Ly, transport_fhighN_Lz,
       transport_fhighN_Rx, transport_fhighN_Ry, transport_fhighN_Rz},
      {transport_fhighE_Lx, transport_fhighE_Ly, transport_fhighE_Lz,
       transport_fhighE_Rx, transport_fhighE_Ry, transport_fhighE_Rz},
      {transport_phiN_Lx, transport_phiN_Ly, transport_phiN_Lz,
       transport_phiN_Rx, transport_phiN_Ry, transport_phiN_Rz},
      {transport_phiE_Lx, transport_phiE_Ly, transport_phiE_Lz,
       transport_phiE_Rx, transport_phiE_Ry, transport_phiE_Rz},
      {transport_sawN_Lx, transport_sawN_Ly, transport_sawN_Lz,
       transport_sawN_Rx, transport_sawN_Ry, transport_sawN_Rz},
      {transport_sawE_Lx, transport_sawE_Ly, transport_sawE_Lz,
       transport_sawE_Rx, transport_sawE_Ry, transport_sawE_Rz},
  };
  CCTK_REAL *diag_physN[3][3] = {
      {transport_fphysN_Lx, transport_fphysN_Cx, transport_fphysN_Rx},
      {transport_fphysN_Ly, transport_fphysN_Cy, transport_fphysN_Ry},
      {transport_fphysN_Lz, transport_fphysN_Cz, transport_fphysN_Rz},
  };
  CCTK_REAL *diag_physE[3][3] = {
      {transport_fphysE_Lx, transport_fphysE_Cx, transport_fphysE_Rx},
      {transport_fphysE_Ly, transport_fphysE_Cy, transport_fphysE_Ry},
      {transport_fphysE_Lz, transport_fphysE_Cz, transport_fphysE_Rz},
  };
  CCTK_REAL *diag3[2][3] = {
      {transport_divN_x, transport_divN_y, transport_divN_z},
      {transport_divE_x, transport_divE_y, transport_divE_z},
  };
  CCTK_REAL *traceEz[10] = {
      transport_uE_jm2_z, transport_uE_jm1_z, transport_uE_j_z,
      transport_uE_jp1_z, transport_uE_jp2_z, transport_idxE_jm2_z,
      transport_idxE_jm1_z, transport_idxE_j_z, transport_idxE_jp1_z,
      transport_idxE_jp2_z,
  };

  const auto &grid = GridDescBaseDevice(cctkGH);
  GF3D2layout layout(cctkGH, {1, 1, 1});

  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones, [=] CCTK_DEVICE(const PointDesc &p) {
        const int idx_trace = layout.linear(p.i, p.j, p.k, 0);
        for (int ig = 0; ig < groupspec; ++ig) {
          const int i4D = layout.linear(p.i, p.j, p.k, ig);
          rN_rhs[i4D] = 0.0;
          rE_rhs[i4D] = 0.0;
          rFx_rhs[i4D] = 0.0;
          rFy_rhs[i4D] = 0.0;
          rFz_rhs[i4D] = 0.0;
          for (int group = 0; group < 18; ++group) {
            for (int comp = 0; comp < 6; ++comp) {
              diag6[group][comp][i4D] = 0.0;
            }
          }
          for (int dir = 0; dir < 3; ++dir) {
            for (int comp = 0; comp < 3; ++comp) {
              diag_physN[dir][comp][i4D] = 0.0;
              diag_physE[dir][comp][i4D] = 0.0;
            }
          }
          for (int group = 0; group < 2; ++group) {
            for (int comp = 0; comp < 3; ++comp) {
              diag3[group][comp][i4D] = 0.0;
            }
          }
        }
        for (int comp = 0; comp < 10; ++comp) {
          traceEz[comp][idx_trace] = 0.0;
        }
      });

  for (int dir = 0; dir < 3; ++dir) {
    CCTK_REAL *nu_flux_dir =
        (dir == 0 ? nu_flux_x : (dir == 1 ? nu_flux_y : nu_flux_z));
    CCTK_REAL *nu_cmax_dir =
        (dir == 0 ? nu_cmax_x : (dir == 1 ? nu_cmax_y : nu_cmax_z));

    grid.loop_all_device<1, 1, 1>(
        grid.nghostzones, [=] CCTK_DEVICE(const PointDesc &p) {
          for (int ig = 0; ig < groupspec; ++ig) {
            const int idx_cmax = layout.linear(p.i, p.j, p.k, ig);
            nu_cmax_dir[idx_cmax] = 0.0;
          }
          for (int comp = 0; comp < flux_comps; ++comp) {
            const int idx_flux = layout.linear(p.i, p.j, p.k, comp);
            nu_flux_dir[idx_flux] = 0.0;
          }
        });
  }
}

} // namespace nuX_M1
