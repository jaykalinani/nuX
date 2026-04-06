#include <cassert>
#include <cmath>
#include <loop_device.hxx>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "nuX_utils.hxx"
#include "nuX_M1_closure.hxx"

namespace nuX_M1 {

using namespace nuX_Utils;
using namespace Loop;
using namespace std;

extern "C" void nuX_M1_CalcRadialFluxes(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_M1_CalcRadialFluxes;
  DECLARE_CCTK_PARAMETERS

  if (verbose && CCTK_MyProc(cctkGH) == 0) {
    CCTK_INFO("nuX_M1_CalcRadialFluxes");
  }

  const GridDescBaseDevice grid(cctkGH);
  const GF3D2layout layout_cc(cctkGH, {1, 1, 1});
  const GF3D2layout layout_vc(cctkGH, {0, 0, 0});

  tensor::slicing_geometry_const geom(layout_vc, layout_cc, alp, betax, betay,
                                      betaz, gxx, gxy, gxz, gyy, gyz, gzz, kxx,
                                      kxy, kxz, kyy, kyz, kzz);
  tensor::fluid_velocity_field_const fidu(layout_vc, layout_cc, alp, betax,
                                          betay, betaz, fidu_w_lorentz,
                                          fidu_velx, fidu_vely, fidu_velz);

  // A scale-aware radius floor: ~1e-12 * dx (dx from the grid)
  const CCTK_REAL dx = cctk_delta_space[0];
  const CCTK_REAL dy = cctk_delta_space[1];
  const CCTK_REAL dz = cctk_delta_space[2];
  const CCTK_REAL dmin = fmin(dx, fmin(dy, dz));
  const CCTK_REAL r_floor = 1.0e-12 * dmin; // very small, but nonzero
  const CCTK_REAL r_floor2 = r_floor * r_floor;

  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        const int ijk = layout_cc.linear(p.i, p.j, p.k);

        if (nuX_m1_mask[ijk]) {
          for (int ig = 0; ig < nspecies * ngroups; ++ig) {
            int const i4D = layout_cc.linear(p.i, p.j, p.k, ig);
            radial_flux_0[i4D] = 0.0;
            radial_flux_1[i4D] = 0.0;
          }
          return;
        }

        tensor::metric<4> g_dd;
        geom.get_metric(p, &g_dd);

        tensor::inv_metric<4> g_uu;
        geom.get_inv_metric(p, &g_uu);

        tensor::generic<CCTK_REAL, 4, 1> beta_u;
        geom.get_shift_vec(p, &beta_u);
        const CCTK_REAL alp_cc = geom.get_lapse(p);

        // Coordinate position vector (spatial components only)
        tensor::generic<CCTK_REAL, 4, 1> r_d;
        r_d(0) = 0.0;
        r_d(1) = p.x;
        r_d(2) = p.y;
        r_d(3) = p.z;

        // Robust, positive-definite coordinate radius
        const CCTK_REAL rr2 = p.x * p.x + p.y * p.y + p.z * p.z;

        // Guard against NaNs / zeros so we never divide by NaN/0
        // Note: (rr2 > r_floor2) is false if rr2 is NaN.
        if (!(rr2 > r_floor2)) {
          for (int ig = 0; ig < nspecies * ngroups; ++ig) {
            int const i4D = layout_cc.linear(p.i, p.j, p.k, ig);
            radial_flux_0[i4D] = 0.0;
            radial_flux_1[i4D] = 0.0;
          }
          return;
        }

        const CCTK_REAL rr = sqrt(rr2);

        // Extra guard (covers weirdness if sqrt ever returns NaN)
        if (!(rr > 0.0) || !isfinite(rr)) {
          for (int ig = 0; ig < nspecies * ngroups; ++ig) {
            int const i4D = layout_cc.linear(p.i, p.j, p.k, ig);
            radial_flux_0[i4D] = 0.0;
            radial_flux_1[i4D] = 0.0;
          }
          return;
        }

        const CCTK_REAL irr = 1.0 / rr;

        tensor::generic<CCTK_REAL, 4, 1> u_u;
        fidu.get(p, &u_u);

        tensor::generic<CCTK_REAL, 4, 1> F_d;
        tensor::generic<CCTK_REAL, 4, 1> F_u;
        tensor::generic<CCTK_REAL, 4, 1> H_d;
        tensor::generic<CCTK_REAL, 4, 1> H_u;
        tensor::generic<CCTK_REAL, 4, 1> fnu_u;

        for (int ig = 0; ig < ngroups * nspecies; ++ig) {
          int const i4D = layout_cc.linear(p.i, p.j, p.k, ig);
          pack_F_d(beta_u(1), beta_u(2), beta_u(3), rFx[i4D], rFy[i4D],
                   rFz[i4D], &F_d);
          pack_H_d(rHt[i4D], rHx[i4D], rHy[i4D], rHz[i4D], &H_d);

          tensor::contract(g_uu, F_d, &F_u);
          tensor::contract(g_uu, H_d, &H_u);

          radial_flux_1[i4D] = 0.0;
          for (int a = 1; a < 4; ++a) {
            radial_flux_1[i4D] +=
                r_d(a) * irr * calc_E_flux(alp_cc, beta_u, rE[i4D], F_u, a);
          }

          assemble_fnu(u_u, rJ[i4D], H_u, &fnu_u, rad_E_floor);
          const CCTK_REAL Gamma = alp_cc * fnu_u(0);
          const CCTK_REAL nnu = rN[i4D] / Gamma;

          radial_flux_0[i4D] = alp_cc * irr * nnu * tensor::dot(fnu_u, r_d);
        }
      });
}

} // namespace nuX_M1
