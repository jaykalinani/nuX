#include <algorithm>
#include <cassert>
#include <cmath>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "nuX_M1_closure.hxx"
#include "nuX_M1_macro.hxx"
#include "nuX_utils.hxx"

namespace nuX_M1 {

using namespace std;
using namespace Loop;
using namespace nuX_Utils;

extern "C" void nuX_M1_CalcSourceRHS(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_M1_CalcSourceRHS;
  DECLARE_CCTK_PARAMETERS;

  if (verbose)
    CCTK_INFO("nuX_M1_CalcSourceRHS");

  const GridDescBaseDevice grid(cctkGH);
  const GF3D2layout layout_cc(cctkGH, {1, 1, 1});
  const GF3D2layout layout_vc(cctkGH, {0, 0, 0});

  tensor::slicing_geometry_const geom(layout_vc, layout_cc, alp, betax, betay,
                                      betaz, gxx, gxy, gxz, gyy, gyz, gzz, kxx,
                                      kxy, kxz, kyy, kyz, kzz);
  tensor::fluid_velocity_field_const fidu(layout_vc, layout_cc, alp, betax,
                                          betay, betaz, fidu_w_lorentz,
                                          fidu_velx, fidu_vely, fidu_velz);

  grid.loop_int_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        const int ijk = layout_cc.linear(p.i, p.j, p.k);
        const int groupspec = ngroups * nspecies;

        if (nuX_m1_mask[ijk]) {
          for (int ig = 0; ig < groupspec; ++ig) {
            const int i4D = layout_cc.linear(p.i, p.j, p.k, ig);
            rN_rhs[i4D] = 0.0;
            rE_rhs[i4D] = 0.0;
            rFx_rhs[i4D] = 0.0;
            rFy_rhs[i4D] = 0.0;
            rFz_rhs[i4D] = 0.0;
          }
          return;
        }

        tensor::generic<CCTK_REAL, 4, 1> beta_u;
        geom.get_shift_vec(p, &beta_u);
        const CCTK_REAL alp_ijk = geom.get_lapse(p);
        const CCTK_REAL betax_ijk = beta_u(1);
        const CCTK_REAL betay_ijk = beta_u(2);
        const CCTK_REAL betaz_ijk = beta_u(3);
        const CCTK_REAL W_ijk = fidu_w_lorentz[ijk];

        tensor::metric<4> g_dd;
        tensor::generic<CCTK_REAL, 4, 1> n_u;
        tensor::generic<CCTK_REAL, 4, 2> gamma_ud;
        geom.get_metric(p, &g_dd);
        geom.get_normal(p, &n_u);
        geom.get_space_proj(p, &gamma_ud);
        const CCTK_REAL volform_ijk = sqrt(
            nuX_Utils::metric::spatial_det(g_dd(1, 1), g_dd(1, 2), g_dd(1, 3),
                                           g_dd(2, 2), g_dd(2, 3), g_dd(3, 3)));

        tensor::generic<CCTK_REAL, 4, 1> u_u, u_d;
        fidu.get(p, &u_u);
        tensor::contract(g_dd, u_u, &u_d);

        tensor::generic<CCTK_REAL, 4, 1> v_u;
        pack_v_u(fidu_velx[ijk], fidu_vely[ijk], fidu_velz[ijk], &v_u);

        for (int ig = 0; ig < groupspec; ++ig) {
          const int i4D = layout_cc.linear(p.i, p.j, p.k, ig);
          assert(isfinite(rN[i4D]));
          assert(isfinite(rE[i4D]));
          assert(isfinite(rFx[i4D]));
          assert(isfinite(rFy[i4D]));
          assert(isfinite(rFz[i4D]));

          tensor::generic<CCTK_REAL, 4, 1> F_d;
          pack_F_d(betax_ijk, betay_ijk, betaz_ijk, rFx[i4D], rFy[i4D],
                   rFz[i4D], &F_d);

          const CCTK_REAL Gamma = compute_Gamma(
              W_ijk, v_u, rJ[i4D], rE[i4D], F_d, rad_E_floor, rad_eps);

          tensor::generic<CCTK_REAL, 4, 1> H_d;
          pack_H_d(rHt[i4D], rHx[i4D], rHy[i4D], rHz[i4D], &H_d);

          tensor::generic<CCTK_REAL, 4, 1> S_d, tS_d;
          calc_rad_sources(eta_1[i4D] * volform_ijk, abs_1[i4D], scat_1[i4D],
                           u_d, rJ[i4D], H_d, &S_d);

          rE_rhs[i4D] = calc_rE_source(alp_ijk, n_u, S_d);
          calc_rF_source(alp_ijk, gamma_ud, S_d, &tS_d);
          rFx_rhs[i4D] = tS_d(1);
          rFy_rhs[i4D] = tS_d(2);
          rFz_rhs[i4D] = tS_d(3);
          rN_rhs[i4D] =
              alp_ijk *
              (volform_ijk * eta_0[i4D] - abs_0[i4D] * rN[i4D] / Gamma);

          assert(isfinite(rN_rhs[i4D]));
          assert(isfinite(rE_rhs[i4D]));
          assert(isfinite(rFx_rhs[i4D]));
          assert(isfinite(rFy_rhs[i4D]));
          assert(isfinite(rFz_rhs[i4D]));
        }
      });
}

} // namespace nuX_M1
