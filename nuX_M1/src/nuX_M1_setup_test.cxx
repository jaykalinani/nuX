#include <cmath>
#include <cassert>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "nuX_M1_macro.hxx"

namespace nuX_M1 {

using namespace Loop;

CCTK_DEVICE CCTK_REAL volume(CCTK_REAL R, CCTK_REAL xp, CCTK_REAL yp,
                             CCTK_REAL zp, CCTK_REAL dx, CCTK_REAL dy,
                             CCTK_REAL dz) {
  const int NPOINTS = 10;
  int inside = 0, count = 0;
  for (int ii = 0; ii < NPOINTS; ++ii) {
    CCTK_REAL const myx = (xp - dx / 2.) + (ii + 0.5) * (dx / NPOINTS);
    for (int jj = 0; jj < NPOINTS; ++jj) {
      CCTK_REAL const myy = (yp - dy / 2.) + (jj + 0.5) * (dy / NPOINTS);
      for (int kk = 0; kk < NPOINTS; ++kk) {
        CCTK_REAL const myz = (zp - dz / 2.) + (kk + 0.5) * (dz / NPOINTS);
        ++count;
        if (myx * myx + myy * myy + myz * myz < R * R)
          ++inside;
      }
    }
  }
  return static_cast<CCTK_REAL>(inside) / static_cast<CCTK_REAL>(count);
}

// -----------------------------------------------------------------------------
// Main setup routine
// -----------------------------------------------------------------------------
extern "C" void nuX_M1_SetupTest(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_M1_SetupTest;
  DECLARE_CCTK_PARAMETERS;

  if (verbose)
    CCTK_INFO("nuX_M1_SetupTest");

  const GridDescBaseDevice grid(cctkGH);
  const GF3D2layout layout2(cctkGH, {1, 1, 1});

  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones, [=] CCTK_DEVICE(const PointDesc &p) {
        const int ijk = layout2.linear(p.i, p.j, p.k);

        if (CCTK_Equals(thc_m1_test, "beam")) {
          CCTK_REAL nx = beam_test_dir[0];
          CCTK_REAL ny = beam_test_dir[1];
          CCTK_REAL nz = beam_test_dir[2];
          CCTK_REAL n2 = nx * nx + ny * ny + nz * nz;
          if (n2 > 0) {
            CCTK_REAL nn = sqrt(n2);
            nx /= nn;
            ny /= nn;
            nz /= nn;
          } else {
            nx = 1.0;
            ny = 0.0;
            nz = 0.0;
          }
          CCTK_REAL proj = nx * x[ijk] + ny * y[ijk] + nz * z[ijk];
          CCTK_REAL offset2 = (x[ijk] - nx * x[ijk]) * (x[ijk] - nx * x[ijk]) +
                              (y[ijk] - ny * y[ijk]) * (y[ijk] - ny * y[ijk]) +
                              (z[ijk] - nz * z[ijk]) * (z[ijk] - nz * z[ijk]);
          for (int ig = 0; ig < ngroups * nspecies; ++ig) {
            int const i4D = layout2.linear(p.i, p.j, p.k, ig);
            if (proj < beam_position && offset2 < beam_width * beam_width) {
              rE[i4D] = 1.0;
              rN[i4D] = 1.0;
              rFx[i4D] = nx;
              rFy[i4D] = ny;
              rFz[i4D] = nz;
            } else {
              rE[i4D] = 0.0;
              rN[i4D] = 0.0;
              rFx[i4D] = 0.0;
              rFy[i4D] = 0.0;
              rFz[i4D] = 0.0;
            }
          }
        } else if (CCTK_Equals(thc_m1_test, "diff")) {
          for (int ig = 0; ig < ngroups * nspecies; ++ig) {
            int const i4D = layout2.linear(p.i, p.j, p.k, ig);
            if (CCTK_Equals(diff_profile, "step")) {
              rE[i4D] = (x[ijk] > -0.5 && x[ijk] < 0.5) ? 1.0 : 0.0;
            } else if (CCTK_Equals(diff_profile, "gaussian")) {
              rE[i4D] = exp(-(3 * x[ijk]) * (3 * x[ijk]));
            }
            rN[i4D] = rE[i4D];
            CCTK_REAL const W = fidu_w_lorentz[ijk];
            CCTK_REAL const Jo3 = rE[i4D] / (4 * W * W - 1);
            rFx[i4D] = 4 * W * W * fidu_velx[ijk] * Jo3;
            rFy[i4D] = 4 * W * W * fidu_vely[ijk] * Jo3;
            rFz[i4D] = 4 * W * W * fidu_velz[ijk] * Jo3;
          }
        } else if (CCTK_Equals(thc_m1_test, "equil")) {
          assert(ngroups == 1 && nspecies == 3);
          CCTK_REAL const W = fidu_w_lorentz[ijk];
          for (int is = 0; is < nspecies; ++is) {
            int const i4D = layout2.linear(p.i, p.j, p.k, is);
            CCTK_REAL const Jnu = equil_nudens_1[is];
            rE[i4D] = (4. * W * W - 1.) / 3. * Jnu;
            rFx[i4D] = 4. / 3. * W * W * fidu_velx[ijk] * Jnu;
            rFy[i4D] = 4. / 3. * W * W * fidu_vely[ijk] * Jnu;
            rFz[i4D] = 4. / 3. * W * W * fidu_velz[ijk] * Jnu;
            rN[i4D] = equil_nudens_0[is] * W;
          }
        } else if (CCTK_Equals(thc_m1_test, "kerrschild") ||
                   CCTK_Equals(thc_m1_test, "shadow") ||
                   CCTK_Equals(thc_m1_test, "sphere")) {
          for (int ig = 0; ig < ngroups * nspecies; ++ig) {
            int const i4D = layout2.linear(p.i, p.j, p.k, ig);
            rE[i4D] = rN[i4D] = rFx[i4D] = rFy[i4D] = rFz[i4D] = 0.0;
          }
        } else {
          CCTK_ERROR("Unknown test problem");
        }
      });
}

// -----------------------------------------------------------------------------
// KerrSchild Mask
// -----------------------------------------------------------------------------
extern "C" void nuX_M1_KerrSchild_Mask(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_M1_KerrSchild_Mask;
  DECLARE_CCTK_PARAMETERS;

  if (verbose)
    CCTK_INFO("nuX_M1_KerrSchild_Mask");

  const GridDescBaseDevice grid(cctkGH);
  const GF3D2layout layout2(cctkGH, {1, 1, 1});

  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones, [=] CCTK_DEVICE(const PointDesc &p) {
        const int ijk = layout2.linear(p.i, p.j, p.k);
        if (x[ijk] * x[ijk] + y[ijk] * y[ijk] + z[ijk] * z[ijk] <
            kerr_mask_radius * kerr_mask_radius) {
          nuX_m1_mask[ijk] = 1;
          for (int ig = 0; ig < nspecies * ngroups; ++ig) {
            int const i4D = layout2.linear(p.i, p.j, p.k, ig);
            rN[i4D] = rE[i4D] = rFx[i4D] = rFy[i4D] = rFz[i4D] = 0.0;
          }
        } else {
          nuX_m1_mask[ijk] = 0;
        }
      });
}

// -----------------------------------------------------------------------------
// Hydro setup for shadow/sphere
// -----------------------------------------------------------------------------
extern "C" void nuX_M1_SetupTest_Hydro(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_M1_SetupTest_Hydro;
  DECLARE_CCTK_PARAMETERS;

  if (verbose)
    CCTK_INFO("nuX_M1_SetupTest_Hydro");

  CCTK_REAL dx = CCTK_DELTA_SPACE(0);
  CCTK_REAL dy = CCTK_DELTA_SPACE(1);
  CCTK_REAL dz = CCTK_DELTA_SPACE(2);

  const GridDescBaseDevice grid(cctkGH);
  const GF3D2layout layout2(cctkGH, {1, 1, 1});

  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones, [=] CCTK_DEVICE(const PointDesc &p) {
        const int ijk = layout2.linear(p.i, p.j, p.k);
        if (CCTK_Equals(thc_m1_test, "shadow") ||
            CCTK_Equals(thc_m1_test, "sphere")) {
          rho[ijk] = volume(1.0, x[ijk], y[ijk], z[ijk], dx, dy, dz);
          velx[ijk] = 0.0;
          vely[ijk] = 0.0;
          velz[ijk] = 0.0;
        } else {
          CCTK_ERROR("Unknown test problem in SetupTest_Hydro");
        }
      });
}

} // namespace nuX_M1
