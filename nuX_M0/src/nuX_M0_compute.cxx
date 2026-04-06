#include <AMReX_Gpu.H>

#include <array>
#include <cassert>

#include <mpi.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "nuX_M0_kernel.hxx"
#include "avg_baryon_mass.hpp"

namespace nuX_M0 {

extern "C" void nuX_M0_Compute(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_M0_Compute;
  DECLARE_CCTK_PARAMETERS;

  if ((cctk_iteration - 1) % compute_every != 0 || !*nuX_M0_is_on) {
    return;
  }

  if (verbose && CCTK_MyProc(cctkGH) == 0) {
    CCTK_INFO("nuX_M0_Compute");
  }

  assert(sizeof(CCTK_REAL) == sizeof(double));

  CCTK_REAL const dt = cctk_time - *nuX_M0_time;
  assert(dt >= 0.0);

  int const group_id = CCTK_GroupIndex("nuX_M0::nuX_M0_grid_vars");
  cGroupDynamicData group_data;
  int const ierr = CCTK_GroupDynamicData(cctkGH, group_id, &group_data);
  assert(!ierr);

  int const nrad_local = group_data.ash[0];
  int const nray_local = group_data.ash[1];
  int const npts = nrad_local * nray_local;

  CCTK_REAL const mb = AverageBaryonMass(particle_mass);
  SphericalGrid const grid = *M0Grid;

  amrex::ParallelFor(npts, [=] CCTK_DEVICE(int const idx) {
    nuX_M0_N_nue_old[idx] = nuX_M0_N_nue[idx];
    nuX_M0_N_nua_old[idx] = nuX_M0_N_nua[idx];
    nuX_M0_N_nux_old[idx] = nuX_M0_N_nux[idx];
    nuX_M0_E_nue_old[idx] = nuX_M0_E_nue[idx];
    nuX_M0_E_nua_old[idx] = nuX_M0_E_nua[idx];
    nuX_M0_E_nux_old[idx] = nuX_M0_E_nux[idx];
  });

  CCTK_REAL *ray_nue_num = static_cast<CCTK_REAL *>(
      amrex::The_Arena()->alloc(nray_local * sizeof(CCTK_REAL)));
  CCTK_REAL *ray_nua_num = static_cast<CCTK_REAL *>(
      amrex::The_Arena()->alloc(nray_local * sizeof(CCTK_REAL)));
  CCTK_REAL *ray_nux_num = static_cast<CCTK_REAL *>(
      amrex::The_Arena()->alloc(nray_local * sizeof(CCTK_REAL)));
  CCTK_REAL *ray_nue_ene = static_cast<CCTK_REAL *>(
      amrex::The_Arena()->alloc(nray_local * sizeof(CCTK_REAL)));
  CCTK_REAL *ray_nua_ene = static_cast<CCTK_REAL *>(
      amrex::The_Arena()->alloc(nray_local * sizeof(CCTK_REAL)));
  CCTK_REAL *ray_nux_ene = static_cast<CCTK_REAL *>(
      amrex::The_Arena()->alloc(nray_local * sizeof(CCTK_REAL)));

  amrex::ParallelFor(nray_local, [=] CCTK_DEVICE(int const iray) {
    int const offset = index(nrad_local, 0, iray);

    for (int irad = 0; irad < nrad_local; ++irad) {
      int const ij = offset + irad;
      nuX_M0_mask[ij] = 0;
      if (excision && nuX_M0_alp[ij] <= 0.0) {
        nuX_M0_mask[ij] = 1;
      }

      CCTK_REAL kr = 0.0;
      rad_null(grid, irad, iray, nuX_M0_alp[ij], nuX_M0_betax[ij],
               nuX_M0_betay[ij], nuX_M0_betaz[ij], nuX_M0_gxx[ij],
               nuX_M0_gxy[ij], nuX_M0_gxz[ij], nuX_M0_gyy[ij], nuX_M0_gyz[ij],
               nuX_M0_gzz[ij], nuX_M0_zvecx[ij], nuX_M0_zvecy[ij],
               nuX_M0_zvecz[ij], nuX_M0_mask[ij], nuX_M0_kt[ij], kr,
               nuX_M0_chi[ij], nuX_M0_sdetg[ij]);
      nuX_M0_flux_fac[ij] = nuX_M0_kt[ij] > 0.0 ? kr / nuX_M0_kt[ij] : 0.0;
    }

    for (int irad = 0; irad < nrad_local; ++irad) {
      int const ij = offset + irad;
      if (nuX_M0_mask[ij]) {
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
        nuX_M0_N_nue_old[ij] = 0.0;
        nuX_M0_N_nua_old[ij] = 0.0;
        nuX_M0_N_nux_old[ij] = 0.0;
        nuX_M0_E_nue_old[ij] = 0.0;
        nuX_M0_E_nua_old[ij] = 0.0;
        nuX_M0_E_nux_old[ij] = 0.0;
        nuX_M0_N_nue[ij] = 0.0;
        nuX_M0_N_nua[ij] = 0.0;
        nuX_M0_N_nux[ij] = 0.0;
        nuX_M0_E_nue[ij] = 0.0;
        nuX_M0_E_nua[ij] = 0.0;
        nuX_M0_E_nux[ij] = 0.0;
        continue;
      }

      CCTK_REAL abs_nue = nuX_M0_kappa_0_nue[ij];
      CCTK_REAL abs_nua = nuX_M0_kappa_0_nua[ij];
      if (use_reduced_opacity) {
        abs_nue *= ::exp(-nuX_M0_optd_0_nue[ij]);
        abs_nua *= ::exp(-nuX_M0_optd_0_nua[ij]);
      }
      nuX_M0_abs_nue[ij] = abs_nue;
      nuX_M0_abs_nua[ij] = abs_nua;
    }

    evol_density_slice(grid, dt, &nuX_M0_mask[offset], &nuX_M0_flux_fac[offset],
                       &nuX_M0_sdetg[offset], &nuX_M0_kt[offset],
                       &nuX_M0_R_nue[offset], &nuX_M0_abs_nue[offset],
                       &nuX_M0_N_nue_old[offset], &nuX_M0_N_nue[offset],
                       &nuX_M0_ndens_nue[offset]);
    evol_density_slice(grid, dt, &nuX_M0_mask[offset], &nuX_M0_flux_fac[offset],
                       &nuX_M0_sdetg[offset], &nuX_M0_kt[offset],
                       &nuX_M0_R_nua[offset], &nuX_M0_abs_nua[offset],
                       &nuX_M0_N_nua_old[offset], &nuX_M0_N_nua[offset],
                       &nuX_M0_ndens_nua[offset]);
    evol_density_slice(grid, dt, &nuX_M0_mask[offset], &nuX_M0_flux_fac[offset],
                       &nuX_M0_sdetg[offset], &nuX_M0_kt[offset],
                       &nuX_M0_R_nux[offset], nullptr,
                       &nuX_M0_N_nux_old[offset], &nuX_M0_N_nux[offset],
                       &nuX_M0_ndens_nux[offset]);

    evol_energy_slice(grid, dt, &nuX_M0_mask[offset], &nuX_M0_flux_fac[offset],
                      &nuX_M0_kt[offset], &nuX_M0_chi[offset],
                      &nuX_M0_R_nue[offset], &nuX_M0_Q_nue[offset],
                      &nuX_M0_ndens_nue[offset], &nuX_M0_E_nue_old[offset],
                      &nuX_M0_E_nue[offset], &nuX_M0_eave_nue[offset]);
    evol_energy_slice(grid, dt, &nuX_M0_mask[offset], &nuX_M0_flux_fac[offset],
                      &nuX_M0_kt[offset], &nuX_M0_chi[offset],
                      &nuX_M0_R_nua[offset], &nuX_M0_Q_nua[offset],
                      &nuX_M0_ndens_nua[offset], &nuX_M0_E_nua_old[offset],
                      &nuX_M0_E_nua[offset], &nuX_M0_eave_nua[offset]);
    evol_energy_slice(grid, dt, &nuX_M0_mask[offset], &nuX_M0_flux_fac[offset],
                      &nuX_M0_kt[offset], &nuX_M0_chi[offset],
                      &nuX_M0_R_nux[offset], &nuX_M0_Q_nux[offset],
                      &nuX_M0_ndens_nux[offset], &nuX_M0_E_nux_old[offset],
                      &nuX_M0_E_nux[offset], &nuX_M0_eave_nux[offset]);

    for (int irad = 0; irad < nrad_local; ++irad) {
      int const ij = offset + irad;
      nuX_M0_abs_number[ij] =
          mb * (nuX_M0_abs_nue[ij] * nuX_M0_ndens_nue[ij] -
                nuX_M0_abs_nua[ij] * nuX_M0_ndens_nua[ij]);
      nuX_M0_abs_energy[ij] =
          nuX_M0_abs_nue[ij] * nuX_M0_ndens_nue[ij] * nuX_M0_eave_nue[ij] +
          nuX_M0_abs_nua[ij] * nuX_M0_ndens_nua[ij] * nuX_M0_eave_nua[ij];
    }

    int const itheta = nuX_Utils::sph_grid::get_itheta(grid, iray);
    int const iphi = nuX_Utils::sph_grid::get_iphi(grid, iray);
    CCTK_REAL dS = grid.dtheta * grid.dphi;
    if (iphi == 0 || iphi == grid.nphi - 1)
      dS *= 0.5;
    if (itheta == 0 || itheta == grid.ntheta - 1)
      dS /= grid.nphi;

    int const iout = offset + (nrad_local - 1);
    CCTK_REAL const wt = nuX_M0_flux_fac[iout] * dS;
    ray_nue_num[iray] = nuX_M0_N_nue[iout] * wt;
    ray_nua_num[iray] = nuX_M0_N_nua[iout] * wt;
    ray_nux_num[iray] = nuX_M0_N_nux[iout] * wt;
    ray_nue_ene[iray] = nuX_M0_N_nue[iout] * nuX_M0_eave_nue[iout] * wt;
    ray_nua_ene[iray] = nuX_M0_N_nua[iout] * nuX_M0_eave_nua[iout] * wt;
    ray_nux_ene[iray] = nuX_M0_N_nux[iout] * nuX_M0_eave_nux[iout] * wt;
  });

  amrex::Gpu::Device::streamSynchronize();

  CCTK_REAL my_nue_num = 0.0;
  CCTK_REAL my_nua_num = 0.0;
  CCTK_REAL my_nux_num = 0.0;
  CCTK_REAL my_nue_ene = 0.0;
  CCTK_REAL my_nua_ene = 0.0;
  CCTK_REAL my_nux_ene = 0.0;
  for (int iray = 0; iray < nray_local; ++iray) {
    my_nue_num += ray_nue_num[iray];
    my_nua_num += ray_nua_num[iray];
    my_nux_num += ray_nux_num[iray];
    my_nue_ene += ray_nue_ene[iray];
    my_nua_ene += ray_nua_ene[iray];
    my_nux_ene += ray_nux_ene[iray];
  }

  MPI_Allreduce(&my_nue_num, nuX_M0_nue_num_flux, 1, MPI_DOUBLE, MPI_SUM,
                MPI_COMM_WORLD);
  MPI_Allreduce(&my_nua_num, nuX_M0_nua_num_flux, 1, MPI_DOUBLE, MPI_SUM,
                MPI_COMM_WORLD);
  MPI_Allreduce(&my_nux_num, nuX_M0_nux_num_flux, 1, MPI_DOUBLE, MPI_SUM,
                MPI_COMM_WORLD);
  MPI_Allreduce(&my_nue_ene, nuX_M0_nue_ene_flux, 1, MPI_DOUBLE, MPI_SUM,
                MPI_COMM_WORLD);
  MPI_Allreduce(&my_nua_ene, nuX_M0_nua_ene_flux, 1, MPI_DOUBLE, MPI_SUM,
                MPI_COMM_WORLD);
  MPI_Allreduce(&my_nux_ene, nuX_M0_nux_ene_flux, 1, MPI_DOUBLE, MPI_SUM,
                MPI_COMM_WORLD);

  amrex::The_Arena()->free(ray_nue_num);
  amrex::The_Arena()->free(ray_nua_num);
  amrex::The_Arena()->free(ray_nux_num);
  amrex::The_Arena()->free(ray_nue_ene);
  amrex::The_Arena()->free(ray_nua_ene);
  amrex::The_Arena()->free(ray_nux_ene);

  *nuX_M0_time = cctk_time;
}

} // namespace nuX_M0
