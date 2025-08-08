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

  if (verbose) {
    CCTK_INFO("nuX_M1_CalcGRSources");
  }

  tensor::slicing_geometry_const geom(alp, betax, betay, betaz, gxx, gxy, gxz,
                                      gyy, gyz, gzz, kxx, kxy, kxz, kyy, kyz,
                                      kzz, volform);

  const GridDescBaseDevice grid(cctkGH);
  const GF3D2layout layout2(cctkGH, {1, 1, 1});

  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        const int ijk = layout2.linear(p.i, p.j, p.k);
        if (nuX_m1_mask[ijk])
          return; // skip masked zones

        tensor::inv_metric<3> g_uu;
        geom.get_inv_metric(ijk, &g_uu);

        tensor::symmetric2<CCTK_REAL, 3, 2> K_dd;
        geom.get_extr_curv(ijk, &K_dd);

        GF3D2<const CCTK_REAL> gf_alp{alp};
        tensor::generic<CCTK_REAL, 3, 1> dalp_d;
        for (int a = 0; a < 3; ++a) {
          dalp_d(a) = calc_fd_v2c<FDORDER>(gf_alp, p, a);
        }

        tensor::generic<CCTK_REAL, 3, 2> dbeta_du;
        for (int b = 0; b < 3; ++b) {
          GF3D2<const CCTK_REAL> gf_beta_b{geom.get_shift_comp(b)};
          for (int a = 0; a < 3; ++a) {
            dbeta_du(a, b) = calc_fd_v2c<FDORDER>(gf_beta_b, p, a);
          }
        }

        tensor::symmetric2<CCTK_REAL, 3, 3> dg_ddd;
        for (int b = 0; b < 3; ++b)
          for (int c = b; c < 3; ++c) {
            GF3D2<const CCTK_REAL> gf_g_bc{geom.get_space_metric_comp(b, c)};
            for (int a = 0; a < 3; ++a) {
              dg_ddd(a, b, c) = calc_fd_v2c<FDORDER>(gf_g_bc, p, a);
            }
          }

        // ------------------------------------------------------------------
        //   Radiation-moment source terms
        // ------------------------------------------------------------------
        for (int ig = 0; ig < ngroups * nspecies; ++ig) {

          int const i4D = layout2.linear(p.i, p.j, p.k, ig);

          // Pack Fᵢ and P_{ij}
          tensor::generic<CCTK_REAL, 3, 1> F_d;
          pack_F_d(rFx[i4D], rFy[i4D], rFz[i4D], &F_d);

          tensor::symmetric2<CCTK_REAL, 3, 2> P_dd;
          pack_P_dd(rPxx[i4D], rPxy[i4D], rPxz[i4D], rPyy[i4D], rPyz[i4D],
                    rPzz[i4D], &P_dd);

          // P^{ij} = γ^{ik} γ^{jl} P_{kl}
          tensor::symmetric2<CCTK_REAL, 3, 2> P_uu;
          tensor::contract2(g_uu, P_dd, &P_uu);

          // Pointers to flux RHS components
          tensor::generic<CCTK_REAL *, 3, 1> rF_rhs;
          rF_rhs(0) = &rFx_rhs[i4D];
          rF_rhs(1) = &rFy_rhs[i4D];
          rF_rhs(2) = &rFz_rhs[i4D];

          // ------------------------------------------------------------
          //   Energy source   (densitised already)
          // ------------------------------------------------------------
          rE_rhs[i4D] += alp[ijk] * tensor::dot(P_uu, K_dd) -
                         tensor::dot(g_uu, F_d, dalp_d);

          // ------------------------------------------------------------
          //   Flux source
          // ------------------------------------------------------------
          for (int a = 0; a < 3; ++a) {

            *rF_rhs(a) -= rE[i4D] * dalp_d(a);

            for (int b = 0; b < 3; ++b)
              *rF_rhs(a) += F_d(b) * dbeta_du(a, b);

            for (int b = 0; b < 3; ++b)
              for (int c = 0; c < 3; ++c)
                *rF_rhs(a) += 0.5 * alp[ijk] * P_uu(b, c) * dg_ddd(a, b, c);
          }
        }
      }); // end loop_all_device
} // nuX_M1_CalcGRSources

} // namespace nuX_M1
