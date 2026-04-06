#include <loop_device.hxx>
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "nuX_M1_closure.hxx"
#include "nuX_M1_macro.hxx"
#include "aster_utils.hxx"
#include "nuX_utils.hxx"

#define FDORDER 2

namespace nuX_M1 {

using namespace AsterUtils; // calc_fd_v2c<>
using namespace nuX_Utils;  // tensor helpers
using namespace Loop;

// NOTE: unlike the rest of M1, spatial indices run over 0,1,2 here
extern "C" void nuX_M1_CalcGRSources(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_M1_CalcGRSources;
  DECLARE_CCTK_PARAMETERS;

  if (verbose && CCTK_MyProc(cctkGH) == 0) {
    CCTK_INFO("nuX_M1_CalcGRSources");
  }

  const GridDescBaseDevice grid(cctkGH);
  const GF3D2layout layout_cc(cctkGH, {1, 1, 1});
  const GF3D2layout layout_vc(cctkGH, {0, 0, 0});
  tensor::slicing_geometry_const geom(layout_vc, layout_cc, alp, betax, betay,
                                      betaz, gxx, gxy, gxz, gyy, gyz, gzz, kxx,
                                      kxy, kxz, kyy, kyz, kzz);

  const GF3D2<const CCTK_REAL> gf_alp(layout_vc, alp);
  const GF3D2<const CCTK_REAL> gf_betax(layout_vc, betax);
  const GF3D2<const CCTK_REAL> gf_betay(layout_vc, betay);
  const GF3D2<const CCTK_REAL> gf_betaz(layout_vc, betaz);
  const GF3D2<const CCTK_REAL> gf_gxx(layout_vc, gxx);
  const GF3D2<const CCTK_REAL> gf_gxy(layout_vc, gxy);
  const GF3D2<const CCTK_REAL> gf_gxz(layout_vc, gxz);
  const GF3D2<const CCTK_REAL> gf_gyy(layout_vc, gyy);
  const GF3D2<const CCTK_REAL> gf_gyz(layout_vc, gyz);
  const GF3D2<const CCTK_REAL> gf_gzz(layout_vc, gzz);

  grid.loop_int_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        const int ijk = layout_cc.linear(p.i, p.j, p.k);
        if (nuX_m1_mask[ijk])
          return; // skip masked zones

        tensor::inv_metric<3> g_uu;
        geom.get_inv_metric(p, &g_uu);

        tensor::symmetric2<CCTK_REAL, 3, 2> K_dd;
        geom.get_extr_curv(p, &K_dd);
        const CCTK_REAL alp_ijk = geom.get_lapse(p);

        tensor::generic<CCTK_REAL, 3, 1> dalp_d;
        for (int a = 0; a < 3; ++a) {
          dalp_d(a) = calc_fd_v2c<FDORDER>(gf_alp, p, a);
        }

        tensor::generic<CCTK_REAL, 3, 2> dbeta_du;
        for (int a = 0; a < 3; ++a) {
          dbeta_du(a, 0) = calc_fd_v2c<FDORDER>(gf_betax, p, a);
          dbeta_du(a, 1) = calc_fd_v2c<FDORDER>(gf_betay, p, a);
          dbeta_du(a, 2) = calc_fd_v2c<FDORDER>(gf_betaz, p, a);
        }

        tensor::symmetric2<CCTK_REAL, 3, 3> dg_ddd;
        for (int a = 0; a < 3; ++a) {
          dg_ddd(a, 0, 0) = calc_fd_v2c<FDORDER>(gf_gxx, p, a);
          dg_ddd(a, 0, 1) = calc_fd_v2c<FDORDER>(gf_gxy, p, a);
          dg_ddd(a, 0, 2) = calc_fd_v2c<FDORDER>(gf_gxz, p, a);
          dg_ddd(a, 1, 1) = calc_fd_v2c<FDORDER>(gf_gyy, p, a);
          dg_ddd(a, 1, 2) = calc_fd_v2c<FDORDER>(gf_gyz, p, a);
          dg_ddd(a, 2, 2) = calc_fd_v2c<FDORDER>(gf_gzz, p, a);
        }

        // ------------------------------------------------------------------
        //   Radiation-moment source terms
        // ------------------------------------------------------------------
        for (int ig = 0; ig < ngroups * nspecies; ++ig) {

          int const i4D = layout_cc.linear(p.i, p.j, p.k, ig);

          // Pack Fᵢ and P_{ij}
          tensor::generic<CCTK_REAL, 3, 1> F_d;
          pack_F_d(rFx[i4D], rFy[i4D], rFz[i4D], &F_d);

          tensor::symmetric2<CCTK_REAL, 3, 2> P_dd;
          pack_P_dd(rPxx[i4D], rPxy[i4D], rPxz[i4D], rPyy[i4D], rPyz[i4D],
                    rPzz[i4D], &P_dd);

          // P^{ij} = γ^{ik} γ^{jl} P_{kl}
          tensor::symmetric2<CCTK_REAL, 3, 2> P_uu;
          tensor::contract2(g_uu, P_dd, &P_uu);

          CCTK_REAL rF_rhs_local[3] = {rFx_rhs[i4D], rFy_rhs[i4D],
                                       rFz_rhs[i4D]};

          // ------------------------------------------------------------
          //   Energy source   (densitised already)
          // ------------------------------------------------------------
          rE_rhs[i4D] += alp_ijk * tensor::dot(P_uu, K_dd) -
                         tensor::dot(g_uu, F_d, dalp_d);

          // ------------------------------------------------------------
          //   Flux source
          // ------------------------------------------------------------
          for (int a = 0; a < 3; ++a) {

            rF_rhs_local[a] -= rE[i4D] * dalp_d(a);

            for (int b = 0; b < 3; ++b)
              rF_rhs_local[a] += F_d(b) * dbeta_du(a, b);

            for (int b = 0; b < 3; ++b)
              for (int c = 0; c < 3; ++c)
                rF_rhs_local[a] += 0.5 * alp_ijk * P_uu(b, c) * dg_ddd(a, b, c);
          }
          rFx_rhs[i4D] = rF_rhs_local[0];
          rFy_rhs[i4D] = rF_rhs_local[1];
          rFz_rhs[i4D] = rF_rhs_local[2];
        }
      }); // end loop_all_device
} // nuX_M1_CalcGRSources

} // namespace nuX_M1
