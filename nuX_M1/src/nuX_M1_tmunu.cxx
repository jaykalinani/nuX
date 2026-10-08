#include <cmath>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

#define NUX_M1_CLOSURE_IMPLEMENTATION
#include "nuX_M1_closure.hxx"
#include "nuX_utils.hxx"

namespace nuX_M1 {

using namespace Loop;
using namespace nuX_Utils;

template <int interp_order>
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline bool
mask_active_c2v(const GF3D2layout &layout, const CCTK_REAL *mask,
                const PointDesc &p) {
  if constexpr (interp_order == 2) {
    for (int dk = 0; dk < 2; ++dk)
      for (int dj = 0; dj < 2; ++dj)
        for (int di = 0; di < 2; ++di)
          if (mask[layout.linear(p.i - di, p.j - dj, p.k - dk)])
            return true;
  } else {
    for (int dk = 0; dk < 4; ++dk)
      for (int dj = 0; dj < 4; ++dj)
        for (int di = 0; di < 4; ++di)
          if (mask[layout.linear(p.i + di - 2, p.j + dj - 2, p.k + dk - 2)])
            return true;
  }
  return false;
}

template <int interp_order> void AddToTmunu(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_M1_AddToTmunu;
  DECLARE_CCTK_PARAMETERS;

  closure_t closure_fun = CLOSURE_EDDINGTON;
  if (CCTK_Equals(closure, "Eddington")) {
    closure_fun = CLOSURE_EDDINGTON;
  } else if (CCTK_Equals(closure, "Kershaw")) {
    closure_fun = CLOSURE_KERSHAW;
  } else if (CCTK_Equals(closure, "Minerbo")) {
    closure_fun = CLOSURE_MINERBO;
  } else if (CCTK_Equals(closure, "thin")) {
    closure_fun = CLOSURE_THIN;
  } else {
    CCTK_VINFO("Unknown closure \"%s\"", closure);
    CCTK_ERROR("Unsupported closure");
  }

  const GridDescBaseDevice grid(cctkGH);
  const GF3D2layout layout_cc(cctkGH, {1, 1, 1});
  const GF3D2layout layout_vc(cctkGH, {0, 0, 0});

  grid.loop_int_device<0, 0, 0>(grid.nghostzones, [=] CCTK_DEVICE(
                                                      const PointDesc &p) {
    if (mask_active_c2v<interp_order>(layout_cc, nuX_m1_mask, p)) {
      return;
    }

    const int ijk = layout_vc.linear(p.i, p.j, p.k);
    const CCTK_REAL alpha = alp[ijk];
    const CCTK_REAL beta_x = betax[ijk];
    const CCTK_REAL beta_y = betay[ijk];
    const CCTK_REAL beta_z = betaz[ijk];
    const CCTK_REAL spatial_det = nuX_Utils::metric::spatial_det(
        gxx[ijk], gxy[ijk], gxz[ijk], gyy[ijk], gyz[ijk], gzz[ijk]);
    if (!(isfinite(alpha) && alpha > 0.0 && isfinite(beta_x) &&
          isfinite(beta_y) && isfinite(beta_z) && isfinite(gxx[ijk]) &&
          isfinite(gxy[ijk]) && isfinite(gxz[ijk]) && isfinite(gyy[ijk]) &&
          isfinite(gyz[ijk]) && isfinite(gzz[ijk]) && isfinite(spatial_det) &&
          spatial_det > 0.0))
      device_abort();

    tensor::metric<4> g_dd;
    tensor::inv_metric<4> g_uu;
    tensor::generic<CCTK_REAL, 4, 1> n_d;
    g_dd.from_adm(alpha, beta_x, beta_y, beta_z, gxx[ijk], gxy[ijk], gxz[ijk],
                  gyy[ijk], gyz[ijk], gzz[ijk]);
    g_uu.from_adm(alpha, beta_x, beta_y, beta_z, gxx[ijk], gxy[ijk], gxz[ijk],
                  gyy[ijk], gyz[ijk], gzz[ijk]);
    n_d(0) = -alpha;
    n_d(1) = 0.0;
    n_d(2) = 0.0;
    n_d(3) = 0.0;

    CCTK_REAL vx = tensor::interp_c2v<interp_order>(layout_cc, fidu_velx, p);
    CCTK_REAL vy = tensor::interp_c2v<interp_order>(layout_cc, fidu_vely, p);
    CCTK_REAL vz = tensor::interp_c2v<interp_order>(layout_cc, fidu_velz, p);

    // Interpolating W independently from v^i does not preserve
    // u^a u_a=-1.  Interpolate the velocity and metric, repair a possible
    // high-order superluminal overshoot, and reconstruct W at the vertex.
    const CCTK_REAL W = repair_velocity_and_compute_W(g_dd, &vx, &vy, &vz);

    tensor::generic<CCTK_REAL, 4, 1> u_u;
    tensor::generic<CCTK_REAL, 4, 1> u_d;
    tensor::generic<CCTK_REAL, 4, 2> proj_ud;
    valencia::uvel(alpha, beta_x, beta_y, beta_z, W, vx, vy, vz, u_u.data());
    tensor::contract(g_dd, u_u, &u_d);
    calc_proj(u_d, u_u, &proj_ud);

    tensor::generic<CCTK_REAL, 4, 1> v_u;
    tensor::generic<CCTK_REAL, 4, 1> v_d;
    pack_v_u(vx, vy, vz, &v_u);
    tensor::contract(g_dd, v_u, &v_d);

    tensor::generic<CCTK_REAL, 4, 1> F_d;
    tensor::generic<CCTK_REAL, 4, 1> beta_u;
    tensor::symmetric2<CCTK_REAL, 4, 2> P_dd;
    tensor::symmetric2<CCTK_REAL, 4, 2> rT_dd;

    CCTK_REAL const volform = sqrt(spatial_det);
    CCTK_REAL const iV = 1.0 / volform;
    beta_u(0) = 0.0;
    beta_u(1) = beta_x;
    beta_u(2) = beta_y;
    beta_u(3) = beta_z;

    for (int ig = 0; ig < nspecies * ngroups; ++ig) {
      CCTK_REAL rE_v = tensor::interp_c2v<interp_order>(layout_cc, rE, p, ig);
      CCTK_REAL const rFx_v =
          tensor::interp_c2v<interp_order>(layout_cc, rFx, p, ig);
      CCTK_REAL const rFy_v =
          tensor::interp_c2v<interp_order>(layout_cc, rFy, p, ig);
      CCTK_REAL const rFz_v =
          tensor::interp_c2v<interp_order>(layout_cc, rFz, p, ig);

      pack_F_d(beta_u(1), beta_u(2), beta_u(3), rFx_v, rFy_v, rFz_v, &F_d);
      repair_moments(g_uu, &rE_v, &F_d, rad_E_floor, rad_eps);

      CCTK_REAL mychi = 0.5;

      calc_closure(cctkGH, p.i, p.j, p.k, ig, closure_fun, g_dd, g_uu, n_d, W,
                   u_u, v_d, proj_ud, rE_v, F_d, &mychi, &P_dd, closure_epsilon,
                   closure_maxiter, use_fallback != 0);

      assemble_rT(n_d, rE_v, F_d, P_dd, &rT_dd);

      eTtt[ijk] += rT_dd(0, 0) * iV;
      eTtx[ijk] += rT_dd(0, 1) * iV;
      eTty[ijk] += rT_dd(0, 2) * iV;
      eTtz[ijk] += rT_dd(0, 3) * iV;
      eTxx[ijk] += rT_dd(1, 1) * iV;
      eTxy[ijk] += rT_dd(1, 2) * iV;
      eTxz[ijk] += rT_dd(1, 3) * iV;
      eTyy[ijk] += rT_dd(2, 2) * iV;
      eTyz[ijk] += rT_dd(2, 3) * iV;
      eTzz[ijk] += rT_dd(3, 3) * iV;
    }
  });
}

extern "C" void nuX_M1_AddToTmunu(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_M1_AddToTmunu;
  DECLARE_CCTK_PARAMETERS;

  if (ODESolvers_InPreImplicitStage())
    return;

  if (verbose) {
    CCTK_INFO("nuX_M1_AddToTmunu");
  }

  switch (tmunu_interp_order) {
  case 2:
    AddToTmunu<2>(cctkGH);
    break;
  case 4:
    AddToTmunu<4>(cctkGH);
    break;
  default:
    CCTK_ERROR("tmunu_interp_order must be set to 2 or 4.");
  }
}

} // namespace nuX_M1
