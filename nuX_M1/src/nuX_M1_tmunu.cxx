#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include <cmath>

#include "nuX_M1_macro.hxx"
#include "nuX_M1_closure.hxx"
#include "nuX_utils.hxx"

namespace nuX_M1 {

using namespace Loop;
using namespace nuX_Utils;

extern "C" void nuX_M1_AddToTmunu(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_M1_AddToTmunu;
  DECLARE_CCTK_PARAMETERS;

  if (verbose && CCTK_MyProc(cctkGH) == 0) {
    CCTK_INFO("nuX_M1_AddToTmunu");
  }

  // Pick closure
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
    if (CCTK_MyProc(cctkGH) == 0) {
      CCTK_VINFO("Unknown closure \"%s\"", closure);
    }
    CCTK_ERROR("Unsupported closure");
  }
  // Setup grid layout
  const GridDescBaseDevice grid(cctkGH);
  const GF3D2layout layout_cc(cctkGH, {1, 1, 1});
  const GF3D2layout layout_vc(cctkGH, {0, 0, 0});
  tensor::slicing_geometry_const geom(layout_vc, layout_cc, alp, betax, betay,
                                      betaz, gxx, gxy, gxz, gyy, gyz, gzz, kxx,
                                      kxy, kxz, kyy, kyz, kzz);
  tensor::fluid_velocity_field_const fidu(layout_vc, layout_cc, alp, betax,
                                          betay, betaz, fidu_w_lorentz,
                                          fidu_velx, fidu_vely, fidu_velz);

  // Loop over all points
  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones, [=] CCTK_DEVICE(const PointDesc &p) {
        const int ijk = layout_cc.linear(p.i, p.j, p.k);

        if (nuX_m1_mask[ijk]) {
          return;
        }

        // Geometry and velocity fields
        tensor::metric<4> g_dd;
        tensor::inv_metric<4> g_uu;
        tensor::generic<CCTK_REAL, 4, 1> n_d;
        geom.get_metric(p, &g_dd);
        geom.get_inv_metric(p, &g_uu);
        geom.get_normal_form(p, &n_d);

        CCTK_REAL const W = fidu_w_lorentz[ijk];
        tensor::generic<CCTK_REAL, 4, 1> u_u;
        tensor::generic<CCTK_REAL, 4, 1> u_d;
        tensor::generic<CCTK_REAL, 4, 2> proj_ud;
        fidu.get(p, &u_u);
        tensor::contract(g_dd, u_u, &u_d);
        calc_proj(u_d, u_u, &proj_ud);

        tensor::generic<CCTK_REAL, 4, 1> v_u;
        tensor::generic<CCTK_REAL, 4, 1> v_d;
        pack_v_u(fidu_velx[ijk], fidu_vely[ijk], fidu_velz[ijk], &v_u);
        tensor::contract(g_dd, v_u, &v_d);

        tensor::generic<CCTK_REAL, 4, 1> F_d;
        tensor::generic<CCTK_REAL, 4, 1> beta_u;
        tensor::symmetric2<CCTK_REAL, 4, 2> P_dd;
        tensor::symmetric2<CCTK_REAL, 4, 2> rT_dd;

        CCTK_REAL const volform_ijk = sqrt(
            nuX_Utils::metric::spatial_det(g_dd(1, 1), g_dd(1, 2), g_dd(1, 3),
                                           g_dd(2, 2), g_dd(2, 3), g_dd(3, 3)));
        CCTK_REAL const iV = 1.0 / volform_ijk;
        geom.get_shift_vec(p, &beta_u);

        for (int ig = 0; ig < nspecies * ngroups; ++ig) {
          int const i4D = layout_cc.linear(p.i, p.j, p.k, ig);

          pack_F_d(beta_u(1), beta_u(2), beta_u(3), rFx[i4D], rFy[i4D],
                   rFz[i4D], &F_d);

          CCTK_REAL mychi = 0.5;

          calc_closure(cctkGH, p.i, p.j, p.k, ig, closure_fun, g_dd, g_uu, n_d,
                       W, u_u, v_d, proj_ud, rE[i4D], F_d, &mychi, &P_dd,
                       closure_epsilon, closure_maxiter, use_fallback != 0);

          assemble_rT(n_d, rE[i4D], F_d, P_dd, &rT_dd);

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

} // namespace nuX_M1
