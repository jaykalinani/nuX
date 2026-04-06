#include <cassert>
#include <cmath>
#include <memory>

#include <AMReX_Gpu.H>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "nuX_M0_kernel.hxx"

namespace nuX_M0 {

std::unique_ptr<SphericalGrid> M0Grid;

using namespace nuX_Utils::sph_grid;

extern "C" void nuX_M0_SetupGrid(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_M0_SetupGrid;
  DECLARE_CCTK_PARAMETERS;

  int const ntheta = static_cast<int>(std::round(std::sqrt(nray / 2.0)));
  int const nphi = 2 * ntheta;
  assert(ntheta * nphi == nray);

  CCTK_REAL const zero[3] = {0.0, 0.0, 0.0};
  M0Grid = std::make_unique<SphericalGrid>(
      make_grid(zero, rmax, nrad, ntheta, nphi, true));

  int const group_id = CCTK_GroupIndex("nuX_M0::nuX_M0_grid_vars");
  cGroupDynamicData group_data;
  int const ierr = CCTK_GroupDynamicData(cctkGH, group_id, &group_data);
  assert(!ierr);

  SphericalGrid const grid = *M0Grid;
  int const npts = group_data.ash[0] * group_data.ash[1];
  amrex::ParallelFor(npts, [=] CCTK_DEVICE(int const idx) {
    int const irad = idx % group_data.ash[0];
    int const iray = idx / group_data.ash[0];
    nuX_M0_x[idx] = 0.0;
    nuX_M0_y[idx] = 0.0;
    nuX_M0_z[idx] = 0.0;
    nuX_M0_alp[idx] = 0.0;
    nuX_M0_betax[idx] = 0.0;
    nuX_M0_betay[idx] = 0.0;
    nuX_M0_betaz[idx] = 0.0;
    nuX_M0_gxx[idx] = 0.0;
    nuX_M0_gxy[idx] = 0.0;
    nuX_M0_gxz[idx] = 0.0;
    nuX_M0_gyy[idx] = 0.0;
    nuX_M0_gyz[idx] = 0.0;
    nuX_M0_gzz[idx] = 0.0;
    nuX_M0_rho[idx] = 0.0;
    nuX_M0_zvecx[idx] = 0.0;
    nuX_M0_zvecy[idx] = 0.0;
    nuX_M0_zvecz[idx] = 0.0;
    nuX_M0_temp[idx] = 0.0;
    nuX_M0_Ye[idx] = 0.0;
    nuX_M0_kappa_0_nue[idx] = 0.0;
    nuX_M0_kappa_0_nua[idx] = 0.0;
    nuX_M0_kappa_0_nux[idx] = 0.0;
    nuX_M0_optd_0_nue[idx] = 0.0;
    nuX_M0_optd_0_nua[idx] = 0.0;
    nuX_M0_optd_0_nux[idx] = 0.0;
    nuX_M0_optd_1_nue[idx] = 0.0;
    nuX_M0_optd_1_nua[idx] = 0.0;
    nuX_M0_optd_1_nux[idx] = 0.0;
    nuX_M0_R_nue[idx] = 0.0;
    nuX_M0_R_nua[idx] = 0.0;
    nuX_M0_R_nux[idx] = 0.0;
    nuX_M0_Q_nue[idx] = 0.0;
    nuX_M0_Q_nua[idx] = 0.0;
    nuX_M0_Q_nux[idx] = 0.0;
    get_x_y_z(grid, irad, iray, nuX_M0_x[idx], nuX_M0_y[idx], nuX_M0_z[idx]);
  });
}

extern "C" void nuX_M0_RecoverGrid(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_M0_RecoverGrid;
  DECLARE_CCTK_PARAMETERS;

  int const ntheta = static_cast<int>(std::round(std::sqrt(nray / 2.0)));
  int const nphi = 2 * ntheta;
  CCTK_REAL const zero[3] = {0.0, 0.0, 0.0};
  M0Grid = std::make_unique<SphericalGrid>(
      make_grid(zero, rmax, nrad, ntheta, nphi, true));
}

extern "C" void nuX_M0_FreeGrid(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_M0_FreeGrid;
  DECLARE_CCTK_PARAMETERS;

  M0Grid.reset();
}

} // namespace nuX_M0
