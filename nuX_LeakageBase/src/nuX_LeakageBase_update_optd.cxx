#include <array>
#include <cfloat>
#include <cmath>
#include <loop_device.hxx>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "nuX_utils.hxx"

namespace nuX_LeakageBase {

using namespace Loop;
using namespace nuX_Utils;

extern "C" void nuX_LeakageBase_UpdateOpticalDepth(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_LeakageBase_UpdateOpticalDepth;
  DECLARE_CCTK_PARAMETERS;
  using std::isfinite;

  if (verbose && CCTK_MyProc(cctkGH) == 0) {
    CCTK_INFO("nuX_LeakageBase_UpdateOpticalDepth");
  }

  const GridDescBaseDevice grid(cctkGH);
  const GF3D2layout layout_cc(cctkGH, {1, 1, 1});
  const GF3D2layout layout_vc(cctkGH, {0, 0, 0});
  const GF3D2<const CCTK_REAL> gf_gxx(layout_vc, gxx);
  const GF3D2<const CCTK_REAL> gf_gxy(layout_vc, gxy);
  const GF3D2<const CCTK_REAL> gf_gxz(layout_vc, gxz);
  const GF3D2<const CCTK_REAL> gf_gyy(layout_vc, gyy);
  const GF3D2<const CCTK_REAL> gf_gyz(layout_vc, gyz);
  const GF3D2<const CCTK_REAL> gf_gzz(layout_vc, gzz);
  const std::array<GF3D2<CCTK_REAL>, 6> tau_vec = {
      GF3D2<CCTK_REAL>(layout_cc, optd_0_nue),
      GF3D2<CCTK_REAL>(layout_cc, optd_0_nua),
      GF3D2<CCTK_REAL>(layout_cc, optd_0_nux),
      GF3D2<CCTK_REAL>(layout_cc, optd_1_nue),
      GF3D2<CCTK_REAL>(layout_cc, optd_1_nua),
      GF3D2<CCTK_REAL>(layout_cc, optd_1_nux)};
  const std::array<GF3D2<const CCTK_REAL>, 6> tau_vec_p = {
      GF3D2<const CCTK_REAL>(layout_cc, optd_0_nue_p),
      GF3D2<const CCTK_REAL>(layout_cc, optd_0_nua_p),
      GF3D2<const CCTK_REAL>(layout_cc, optd_0_nux_p),
      GF3D2<const CCTK_REAL>(layout_cc, optd_1_nue_p),
      GF3D2<const CCTK_REAL>(layout_cc, optd_1_nua_p),
      GF3D2<const CCTK_REAL>(layout_cc, optd_1_nux_p)};
  const std::array<GF3D2<const CCTK_REAL>, 6> kappa_vec = {
      GF3D2<const CCTK_REAL>(layout_cc, kappa_0_nue),
      GF3D2<const CCTK_REAL>(layout_cc, kappa_0_nua),
      GF3D2<const CCTK_REAL>(layout_cc, kappa_0_nux),
      GF3D2<const CCTK_REAL>(layout_cc, kappa_1_nue),
      GF3D2<const CCTK_REAL>(layout_cc, kappa_1_nua),
      GF3D2<const CCTK_REAL>(layout_cc, kappa_1_nux)};

  CCTK_REAL const dx = CCTK_DELTA_SPACE(0);
  CCTK_REAL const dy = CCTK_DELTA_SPACE(1);
  CCTK_REAL const dz = CCTK_DELTA_SPACE(2);

  grid.loop_int_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        const CCTK_REAL gxx_cc = tensor::interp_v2c(gf_gxx, p);
        const CCTK_REAL gxy_cc = tensor::interp_v2c(gf_gxy, p);
        const CCTK_REAL gxz_cc = tensor::interp_v2c(gf_gxz, p);
        const CCTK_REAL gyy_cc = tensor::interp_v2c(gf_gyy, p);
        const CCTK_REAL gyz_cc = tensor::interp_v2c(gf_gyz, p);
        const CCTK_REAL gzz_cc = tensor::interp_v2c(gf_gzz, p);

        for (int v = 0; v < 6; ++v) {
          tau_vec[v](p.I) = DBL_MAX;
          for (int nk = -1; nk <= 1; ++nk)
            for (int nj = -1; nj <= 1; ++nj)
              for (int ni = -1; ni <= 1; ++ni) {
                if (ni == 0 && nj == 0 && nk == 0) {
                  continue;
                }

                const auto nidx =
                    p.I + p.DI[0] * ni + p.DI[1] * nj + p.DI[2] * nk;
                const CCTK_REAL gxx_n =
                    tensor::interp_v2c(gf_gxx, p, ni, nj, nk);
                const CCTK_REAL gxy_n =
                    tensor::interp_v2c(gf_gxy, p, ni, nj, nk);
                const CCTK_REAL gxz_n =
                    tensor::interp_v2c(gf_gxz, p, ni, nj, nk);
                const CCTK_REAL gyy_n =
                    tensor::interp_v2c(gf_gyy, p, ni, nj, nk);
                const CCTK_REAL gyz_n =
                    tensor::interp_v2c(gf_gyz, p, ni, nj, nk);
                const CCTK_REAL gzz_n =
                    tensor::interp_v2c(gf_gzz, p, ni, nj, nk);
                CCTK_REAL const delta[3] = {ni * dx, nj * dy, nk * dz};
                CCTK_REAL const metric[9] = {
                    0.5 * (gxx_cc + gxx_n), 0.5 * (gxy_cc + gxy_n),
                    0.5 * (gxz_cc + gxz_n), 0.5 * (gxy_cc + gxy_n),
                    0.5 * (gyy_cc + gyy_n), 0.5 * (gyz_cc + gyz_n),
                    0.5 * (gxz_cc + gxz_n), 0.5 * (gyz_cc + gyz_n),
                    0.5 * (gzz_cc + gzz_n)};
                CCTK_REAL dl2 = 0.0;
                for (int a = 0; a < 3; ++a)
                  for (int b = 0; b < 3; ++b) {
                    dl2 += metric[3 * a + b] * delta[a] * delta[b];
                  }
                if (!(dl2 > 0.0) || !isfinite(dl2)) {
                  continue;
                }

                CCTK_REAL const dl = ::sqrt(dl2);
                CCTK_REAL const kappa_avg =
                    0.5 * (kappa_vec[v](p.I) + kappa_vec[v](nidx));
                CCTK_REAL const tau_n = tau_vec_p[v](nidx);
                if (!isfinite(dl) || !isfinite(kappa_avg) ||
                    !isfinite(tau_n)) {
                  continue;
                }

                CCTK_REAL const path_optd = kappa_avg * dl + tau_n;
                if (!isfinite(path_optd)) {
                  continue;
                }

                if (tau_vec[v](p.I) > path_optd) {
                  tau_vec[v](p.I) = path_optd;
                }
              }
          if (tau_vec[v](p.I) == DBL_MAX) {
            tau_vec[v](p.I) = tau_vec_p[v](p.I);
          }
        }
      });
}

} // namespace nuX_LeakageBase
