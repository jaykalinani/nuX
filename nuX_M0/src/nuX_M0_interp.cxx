#include <AMReX_Gpu.H>

#include <array>
#include <cassert>

#include <mpi.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "nuX_M0_kernel.hxx"

namespace nuX_M0 {

using namespace Loop;

extern "C" void nuX_M0_InterpToSph(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_M0_InterpToSph;
  DECLARE_CCTK_PARAMETERS;

  if ((cctk_iteration - 1) % compute_every != 0 || !*nuX_M0_is_on) {
    return;
  }

  if (verbose && CCTK_MyProc(cctkGH) == 0) {
    CCTK_INFO("nuX_M0_InterpToSph");
  }

  assert(sizeof(CCTK_REAL) == sizeof(double));

  const GridDescBaseDevice cart_grid(cctkGH);
  const GF3D2layout layout_cc(cctkGH, {1, 1, 1});
  const GF3D2layout layout_vc(cctkGH, {0, 0, 0});

  const GF3D2<const CCTK_REAL> gf_alp(layout_cc, alp);
  const GF3D2<const CCTK_REAL> gf_betax(layout_cc, betax);
  const GF3D2<const CCTK_REAL> gf_betay(layout_cc, betay);
  const GF3D2<const CCTK_REAL> gf_betaz(layout_cc, betaz);
  const GF3D2<const CCTK_REAL> gf_gxx(layout_vc, gxx);
  const GF3D2<const CCTK_REAL> gf_gxy(layout_vc, gxy);
  const GF3D2<const CCTK_REAL> gf_gxz(layout_vc, gxz);
  const GF3D2<const CCTK_REAL> gf_gyy(layout_vc, gyy);
  const GF3D2<const CCTK_REAL> gf_gyz(layout_vc, gyz);
  const GF3D2<const CCTK_REAL> gf_gzz(layout_vc, gzz);
  const GF3D2<const CCTK_REAL> gf_rho(layout_cc, rho);
  const GF3D2<const CCTK_REAL> gf_temp(layout_cc, temperature);
  const GF3D2<const CCTK_REAL> gf_ye(layout_cc, Ye);
  const GF3D2<const CCTK_REAL> gf_zvecx(layout_cc, zvec_x);
  const GF3D2<const CCTK_REAL> gf_zvecy(layout_cc, zvec_y);
  const GF3D2<const CCTK_REAL> gf_zvecz(layout_cc, zvec_z);
  const GF3D2<const CCTK_REAL> gf_kappa_0_nue(layout_cc, kappa_0_nue);
  const GF3D2<const CCTK_REAL> gf_kappa_0_nua(layout_cc, kappa_0_nua);
  const GF3D2<const CCTK_REAL> gf_kappa_0_nux(layout_cc, kappa_0_nux);
  const GF3D2<const CCTK_REAL> gf_optd_0_nue(layout_cc, optd_0_nue);
  const GF3D2<const CCTK_REAL> gf_optd_0_nua(layout_cc, optd_0_nua);
  const GF3D2<const CCTK_REAL> gf_optd_0_nux(layout_cc, optd_0_nux);
  const GF3D2<const CCTK_REAL> gf_optd_1_nue(layout_cc, optd_1_nue);
  const GF3D2<const CCTK_REAL> gf_optd_1_nua(layout_cc, optd_1_nua);
  const GF3D2<const CCTK_REAL> gf_optd_1_nux(layout_cc, optd_1_nux);
  const GF3D2<const CCTK_REAL> gf_R_nue(layout_cc, R_eff_nue);
  const GF3D2<const CCTK_REAL> gf_R_nua(layout_cc, R_eff_nua);
  const GF3D2<const CCTK_REAL> gf_R_nux(layout_cc, R_eff_nux);
  const GF3D2<const CCTK_REAL> gf_Q_nue(layout_cc, Q_eff_nue);
  const GF3D2<const CCTK_REAL> gf_Q_nua(layout_cc, Q_eff_nua);
  const GF3D2<const CCTK_REAL> gf_Q_nux(layout_cc, Q_eff_nux);

  OwnedBox owned_box{};
  for (int d = 0; d < 3; ++d) {
    int const imin = cart_grid.nghostzones[d];
    int const imax = cart_grid.lsh[d] - cart_grid.nghostzones[d];
    if (imax <= imin) {
      owned_box.xmin[d] = 1.0;
      owned_box.xmax[d] = 0.0;
      owned_box.upper_inclusive[d] = 0;
      continue;
    }
    owned_box.xmin[d] =
        cart_grid.x0[d] + (cart_grid.lbnd[d] + imin) * cart_grid.dx[d] -
        0.5 * cart_grid.dx[d];
    owned_box.xmax[d] =
        cart_grid.x0[d] + (cart_grid.lbnd[d] + (imax - 1)) * cart_grid.dx[d] +
        0.5 * cart_grid.dx[d];
    owned_box.upper_inclusive[d] = cart_grid.bbox[1][d] ? 1 : 0;
  }

  int const group_id = CCTK_GroupIndex("nuX_M0::nuX_M0_grid_vars");
  cGroupDynamicData group_data;
  int const ierr = CCTK_GroupDynamicData(cctkGH, group_id, &group_data);
  assert(!ierr);
  int const npts = group_data.ash[0] * group_data.ash[1];

  CCTK_REAL *hits =
      static_cast<CCTK_REAL *>(amrex::The_Arena()->alloc(npts * sizeof(*hits)));

  amrex::ParallelFor(npts, [=] CCTK_DEVICE(int const idx) {
    CCTK_REAL const xp = nuX_M0_x[idx];
    CCTK_REAL const yp = nuX_M0_y[idx];
    CCTK_REAL const zp = nuX_M0_z[idx];

    if (!inside_owned_box(owned_box, xp, yp, zp)) {
      hits[idx] = 0.0;
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
      return;
    }

    hits[idx] = 1.0;
    nuX_M0_alp[idx] = sample_cart_linear(gf_alp, cart_grid, xp, yp, zp, true);
    nuX_M0_betax[idx] =
        sample_cart_linear(gf_betax, cart_grid, xp, yp, zp, true);
    nuX_M0_betay[idx] =
        sample_cart_linear(gf_betay, cart_grid, xp, yp, zp, true);
    nuX_M0_betaz[idx] =
        sample_cart_linear(gf_betaz, cart_grid, xp, yp, zp, true);
    nuX_M0_gxx[idx] = sample_cart_linear(gf_gxx, cart_grid, xp, yp, zp, false);
    nuX_M0_gxy[idx] = sample_cart_linear(gf_gxy, cart_grid, xp, yp, zp, false);
    nuX_M0_gxz[idx] = sample_cart_linear(gf_gxz, cart_grid, xp, yp, zp, false);
    nuX_M0_gyy[idx] = sample_cart_linear(gf_gyy, cart_grid, xp, yp, zp, false);
    nuX_M0_gyz[idx] = sample_cart_linear(gf_gyz, cart_grid, xp, yp, zp, false);
    nuX_M0_gzz[idx] = sample_cart_linear(gf_gzz, cart_grid, xp, yp, zp, false);
    nuX_M0_rho[idx] = sample_cart_linear(gf_rho, cart_grid, xp, yp, zp, true);
    nuX_M0_zvecx[idx] =
        sample_cart_linear(gf_zvecx, cart_grid, xp, yp, zp, true);
    nuX_M0_zvecy[idx] =
        sample_cart_linear(gf_zvecy, cart_grid, xp, yp, zp, true);
    nuX_M0_zvecz[idx] =
        sample_cart_linear(gf_zvecz, cart_grid, xp, yp, zp, true);
    nuX_M0_temp[idx] =
        sample_cart_linear(gf_temp, cart_grid, xp, yp, zp, true);
    nuX_M0_Ye[idx] = sample_cart_linear(gf_ye, cart_grid, xp, yp, zp, true);
    nuX_M0_kappa_0_nue[idx] =
        sample_cart_linear(gf_kappa_0_nue, cart_grid, xp, yp, zp, true);
    nuX_M0_kappa_0_nua[idx] =
        sample_cart_linear(gf_kappa_0_nua, cart_grid, xp, yp, zp, true);
    nuX_M0_kappa_0_nux[idx] =
        sample_cart_linear(gf_kappa_0_nux, cart_grid, xp, yp, zp, true);
    nuX_M0_optd_0_nue[idx] =
        sample_cart_linear(gf_optd_0_nue, cart_grid, xp, yp, zp, true);
    nuX_M0_optd_0_nua[idx] =
        sample_cart_linear(gf_optd_0_nua, cart_grid, xp, yp, zp, true);
    nuX_M0_optd_0_nux[idx] =
        sample_cart_linear(gf_optd_0_nux, cart_grid, xp, yp, zp, true);
    nuX_M0_optd_1_nue[idx] =
        sample_cart_linear(gf_optd_1_nue, cart_grid, xp, yp, zp, true);
    nuX_M0_optd_1_nua[idx] =
        sample_cart_linear(gf_optd_1_nua, cart_grid, xp, yp, zp, true);
    nuX_M0_optd_1_nux[idx] =
        sample_cart_linear(gf_optd_1_nux, cart_grid, xp, yp, zp, true);
    nuX_M0_R_nue[idx] =
        sample_cart_linear(gf_R_nue, cart_grid, xp, yp, zp, true);
    nuX_M0_R_nua[idx] =
        sample_cart_linear(gf_R_nua, cart_grid, xp, yp, zp, true);
    nuX_M0_R_nux[idx] =
        sample_cart_linear(gf_R_nux, cart_grid, xp, yp, zp, true);
    nuX_M0_Q_nue[idx] =
        sample_cart_linear(gf_Q_nue, cart_grid, xp, yp, zp, true);
    nuX_M0_Q_nua[idx] =
        sample_cart_linear(gf_Q_nua, cart_grid, xp, yp, zp, true);
    nuX_M0_Q_nux[idx] =
        sample_cart_linear(gf_Q_nux, cart_grid, xp, yp, zp, true);
  });

  amrex::Gpu::Device::streamSynchronize();

  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  MPI_Allreduce(MPI_IN_PLACE, hits, npts, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  std::array<CCTK_REAL *, 25> reduce_arrays = {
      nuX_M0_alp,         nuX_M0_betax,      nuX_M0_betay,
      nuX_M0_betaz,       nuX_M0_gxx,        nuX_M0_gxy,
      nuX_M0_gxz,         nuX_M0_gyy,        nuX_M0_gyz,
      nuX_M0_gzz,         nuX_M0_rho,        nuX_M0_zvecx,
      nuX_M0_zvecy,       nuX_M0_zvecz,      nuX_M0_temp,
      nuX_M0_Ye,          nuX_M0_kappa_0_nue, nuX_M0_kappa_0_nua,
      nuX_M0_kappa_0_nux, nuX_M0_optd_0_nue, nuX_M0_optd_0_nua,
      nuX_M0_optd_0_nux,  nuX_M0_optd_1_nue, nuX_M0_optd_1_nua,
      nuX_M0_optd_1_nux};
  for (CCTK_REAL *ptr : reduce_arrays)
    MPI_Allreduce(MPI_IN_PLACE, ptr, npts, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  std::array<CCTK_REAL *, 6> reduce_rates = {nuX_M0_R_nue, nuX_M0_R_nua,
                                             nuX_M0_R_nux, nuX_M0_Q_nue,
                                             nuX_M0_Q_nua, nuX_M0_Q_nux};
  for (CCTK_REAL *ptr : reduce_rates)
    MPI_Allreduce(MPI_IN_PLACE, ptr, npts, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  if (rank == 0) {
    int miss_count = 0;
    int overlap_count = 0;
    for (int idx = 0; idx < npts; ++idx) {
      if (hits[idx] < 0.5)
        ++miss_count;
      if (hits[idx] > 1.5)
        ++overlap_count;
    }
    if (miss_count > 0) {
      CCTK_VWARN(CCTK_WARN_ALERT,
                 "nuX_M0_InterpToSph missed %d ray points on the Cartesian grid",
                 miss_count);
    }
    if (overlap_count > 0) {
      CCTK_VWARN(CCTK_WARN_ALERT,
                 "nuX_M0_InterpToSph found %d multiply-owned ray points",
                 overlap_count);
    }
  }

  amrex::The_Arena()->free(hits);
}

extern "C" void nuX_M0_InterpToCart(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_nuX_M0_InterpToCart;
  DECLARE_CCTK_PARAMETERS;

  if ((cctk_iteration - 1) % compute_every != 0) {
    return;
  }

  if (verbose && CCTK_MyProc(cctkGH) == 0) {
    CCTK_INFO("nuX_M0_InterpToCart");
  }

  if (!*nuX_M0_is_on) {
    grid.loop_all_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          abs_number(p.I) = 0.0;
          abs_energy(p.I) = 0.0;
        });
    return;
  }

  SphericalGrid const sph_grid = *M0Grid;
  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        CCTK_REAL r, theta, phi;
        nuX_Utils::sph_grid::coord_cart_to_sph(p.x, p.y, p.z, r, theta, phi);
        if (r > sph_grid.rmax)
          r = sph_grid.rmax;
        abs_number(p.I) =
            sample_spherical_linear(sph_grid, nuX_M0_abs_number, r, theta, phi);
        abs_energy(p.I) =
            sample_spherical_linear(sph_grid, nuX_M0_abs_energy, r, theta, phi);
      });
}

} // namespace nuX_M0
