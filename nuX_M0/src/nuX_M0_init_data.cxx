#include <cassert>

#include <AMReX_Gpu.H>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "nuX_M0_kernel.hxx"

namespace nuX_M0 {

extern "C" void nuX_M0_InitData(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_M0_InitData;
  DECLARE_CCTK_PARAMETERS;

  if (verbose && CCTK_MyProc(cctkGH) == 0) {
    CCTK_INFO("nuX_M0_InitData");
  }

  int const group_id = CCTK_GroupIndex("nuX_M0::nuX_M0_grid_vars");
  cGroupDynamicData group_data;
  int const ierr = CCTK_GroupDynamicData(cctkGH, group_id, &group_data);
  assert(!ierr);

  int const lsiz = group_data.lsh[0] * group_data.lsh[1];

  *nuX_M0_is_on = 0;
  *nuX_M0_nue_num_flux = 0.0;
  *nuX_M0_nua_num_flux = 0.0;
  *nuX_M0_nux_num_flux = 0.0;
  *nuX_M0_nue_ene_flux = 0.0;
  *nuX_M0_nua_ene_flux = 0.0;
  *nuX_M0_nux_ene_flux = 0.0;

  assert(lsiz == group_data.ash[0] * group_data.ash[1]);

  amrex::ParallelFor(lsiz, [=] CCTK_DEVICE(int const ij) {
    nuX_M0_alp[ij] = 0.0;
    nuX_M0_betax[ij] = 0.0;
    nuX_M0_betay[ij] = 0.0;
    nuX_M0_betaz[ij] = 0.0;
    nuX_M0_gxx[ij] = 0.0;
    nuX_M0_gxy[ij] = 0.0;
    nuX_M0_gxz[ij] = 0.0;
    nuX_M0_gyy[ij] = 0.0;
    nuX_M0_gyz[ij] = 0.0;
    nuX_M0_gzz[ij] = 0.0;
    nuX_M0_rho[ij] = 0.0;
    nuX_M0_zvecx[ij] = 0.0;
    nuX_M0_zvecy[ij] = 0.0;
    nuX_M0_zvecz[ij] = 0.0;
    nuX_M0_temp[ij] = 0.0;
    nuX_M0_Ye[ij] = 0.0;
    nuX_M0_kappa_0_nue[ij] = 0.0;
    nuX_M0_kappa_0_nua[ij] = 0.0;
    nuX_M0_kappa_0_nux[ij] = 0.0;
    nuX_M0_optd_0_nue[ij] = 0.0;
    nuX_M0_optd_0_nua[ij] = 0.0;
    nuX_M0_optd_0_nux[ij] = 0.0;
    nuX_M0_optd_1_nue[ij] = 0.0;
    nuX_M0_optd_1_nua[ij] = 0.0;
    nuX_M0_optd_1_nux[ij] = 0.0;
    nuX_M0_R_nue[ij] = 0.0;
    nuX_M0_R_nua[ij] = 0.0;
    nuX_M0_R_nux[ij] = 0.0;
    nuX_M0_Q_nue[ij] = 0.0;
    nuX_M0_Q_nua[ij] = 0.0;
    nuX_M0_Q_nux[ij] = 0.0;

    nuX_M0_kt[ij] = 0.0;
    nuX_M0_chi[ij] = 0.0;
    nuX_M0_sdetg[ij] = 0.0;
    nuX_M0_abs_nue[ij] = 0.0;
    nuX_M0_abs_nua[ij] = 0.0;
    nuX_M0_ndens_nue[ij] = 0.0;
    nuX_M0_ndens_nua[ij] = 0.0;
    nuX_M0_ndens_nux[ij] = 0.0;
    nuX_M0_eave_nue[ij] = 0.0;
    nuX_M0_eave_nua[ij] = 0.0;
    nuX_M0_eave_nux[ij] = 0.0;
    nuX_M0_abs_number[ij] = 0.0;
    nuX_M0_abs_energy[ij] = 0.0;

    nuX_M0_flux_fac[ij] = 0.0;

    nuX_M0_N_nue[ij] = 0.0;
    nuX_M0_N_nua[ij] = 0.0;
    nuX_M0_N_nux[ij] = 0.0;
    nuX_M0_E_nue[ij] = 0.0;
    nuX_M0_E_nua[ij] = 0.0;
    nuX_M0_E_nux[ij] = 0.0;

    nuX_M0_N_nue_old[ij] = 0.0;
    nuX_M0_N_nua_old[ij] = 0.0;
    nuX_M0_N_nux_old[ij] = 0.0;
    nuX_M0_E_nue_old[ij] = 0.0;
    nuX_M0_E_nua_old[ij] = 0.0;
    nuX_M0_E_nux_old[ij] = 0.0;

    nuX_M0_mask[ij] = 0;
  });

  *nuX_M0_time = cctk_time;
}

} // namespace nuX_M0
