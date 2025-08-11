#include <cassert>
#include <loop_device.hxx>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "nuX_utils.hxx"
#include "nuX_M1_closure.hxx"

namespace nuX_M1 {

using namespace nuX_Utils;
using namespace std;
using namespace Loop;

extern "C" void nuX_M1_CalcClosure(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_M1_CalcClosure;
  DECLARE_CCTK_PARAMETERS;

  if (verbose) {
    CCTK_INFO("nuX_M1_CalcClosure");
  }

  // Disable GSL error handler
  //gsl_error_handler_t *gsl_err = gsl_set_error_handler_off();

  closure_t closure_fun;
  if (CCTK_Equals(closure, "Eddington")) {
    closure_fun = eddington;
  } else if (CCTK_Equals(closure, "Kershaw")) {
    closure_fun = kershaw;
  } else if (CCTK_Equals(closure, "Minerbo")) {
    closure_fun = minerbo;
  } else if (CCTK_Equals(closure, "thin")) {
    closure_fun = thin;
  } else {
    char msg[BUFSIZ];
    snprintf(msg, BUFSIZ, "Unknown closure \"%s\"", closure);
    CCTK_ERROR(msg);
  }

  tensor::slicing_geometry_const geom(alp, betax, betay, betaz, gxx, gxy, gxz,
                                      gyy, gyz, gzz, kxx, kxy, kxz, kyy, kyz,
                                      kzz, volform);
  tensor::fluid_velocity_field_const fidu(alp, betax, betay, betaz,
                                          fidu_w_lorentz, fidu_velx, fidu_vely,
                                          fidu_velz);

  {
    //gsl_root_fsolver *gsl_solver =
    //    gsl_root_fsolver_alloc(gsl_root_fsolver_brent);

    const GridDescBaseDevice grid(cctkGH);
    const GF3D2layout layout2(cctkGH, {1, 1, 1});

    grid.loop_all_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          const int ijk = layout2.linear(p.i, p.j, p.k);
          if (nuX_m1_mask[ijk]) {
            for (int ig = 0; ig < nspecies * ngroups; ++ig) {
              int const i4D = layout2.linear(p.i, p.j, p.k, ig);
              rJ[i4D] = 0;
              rHt[i4D] = 0;
              rHx[i4D] = 0;
              rHy[i4D] = 0;
              rHz[i4D] = 0;
              rPxx[i4D] = 0;
              rPxy[i4D] = 0;
              rPxz[i4D] = 0;
              rPyy[i4D] = 0;
              rPyz[i4D] = 0;
              rPzz[i4D] = 0;
              rnnu[i4D] = 0;
              continue;
            }
          }

          tensor::metric<4> g_dd;
          tensor::inv_metric<4> g_uu;
          tensor::generic<CCTK_REAL, 4, 1> n_d;
          geom.get_metric(ijk, &g_dd);
          geom.get_inv_metric(ijk, &g_uu);
          geom.get_normal_form(ijk, &n_d);

          CCTK_REAL const W = fidu_w_lorentz[ijk];
          tensor::generic<CCTK_REAL, 4, 1> u_u;
          tensor::generic<CCTK_REAL, 4, 1> u_d;
          tensor::generic<CCTK_REAL, 4, 2> proj_ud;
          fidu.get(ijk, &u_u);
          tensor::contract(g_dd, u_u, &u_d);
          calc_proj(u_d, u_u, &proj_ud);

          tensor::generic<CCTK_REAL, 4, 1> v_u;
          tensor::generic<CCTK_REAL, 4, 1> v_d;
          pack_v_u(fidu_velx[ijk], fidu_vely[ijk], fidu_velz[ijk], &v_u);
          tensor::contract(g_dd, v_u, &v_d);

          tensor::generic<CCTK_REAL, 4, 1> H_d;
          tensor::generic<CCTK_REAL, 4, 1> H_u;
          tensor::generic<CCTK_REAL, 4, 1> fnu_u;
          tensor::generic<CCTK_REAL, 4, 1> F_d;
          tensor::symmetric2<CCTK_REAL, 4, 2> P_dd;
          tensor::symmetric2<CCTK_REAL, 4, 2> rT_dd;

          for (int ig = 0; ig < nspecies * ngroups; ++ig) {
            int const i4D = layout2.linear(p.i, p.j, p.k, ig);

            pack_F_d(betax[ijk], betay[ijk], betaz[ijk], rFx[i4D], rFy[i4D],
                     rFz[i4D], &F_d);

            calc_closure(cctkGH, p.i, p.j, p.k, ig, closure_fun,
                         g_dd, g_uu, n_d, W, u_u, v_d, proj_ud, rE[i4D], F_d,
                         &chi[i4D], &P_dd, closure_epsilon, closure_maxiter);

            unpack_P_dd(P_dd, &rPxx[i4D], &rPxy[i4D], &rPxz[i4D], &rPyy[i4D],
                        &rPyz[i4D], &rPzz[i4D]);
            assert(isfinite(rPxx[i4D]));
            assert(isfinite(rPxy[i4D]));
            assert(isfinite(rPxz[i4D]));
            assert(isfinite(rPyy[i4D]));
            assert(isfinite(rPyz[i4D]));
            assert(isfinite(rPzz[i4D]));

            assemble_rT(n_d, rE[i4D], F_d, P_dd, &rT_dd);

            rJ[i4D] = calc_J_from_rT(rT_dd, u_u);
            calc_H_from_rT(rT_dd, u_u, proj_ud, &H_d);
            apply_floor(g_uu, &rJ[i4D], &H_d, rad_E_floor, rad_eps);

            unpack_H_d(H_d, &rHt[i4D], &rHx[i4D], &rHy[i4D], &rHz[i4D]);
            assert(isfinite(rHt[i4D]));
            assert(isfinite(rHx[i4D]));
            assert(isfinite(rHy[i4D]));
            assert(isfinite(rHz[i4D]));

            tensor::contract(g_uu, H_d, &H_u);
            assemble_fnu(u_u, rJ[i4D], H_u, &fnu_u, rad_E_floor);
            CCTK_REAL const Gamma =
                compute_Gamma(fidu_w_lorentz[ijk], v_u, rJ[i4D], rE[i4D], F_d,
                              rad_E_floor, rad_eps);
            assert(Gamma > 0);
            rnnu[i4D] = rN[i4D] / Gamma;
          }
        });
    //gsl_root_fsolver_free(gsl_solver);
  }

  // Restore GSL error handler
  //gsl_set_error_handler(gsl_err);
}

} // namespace nuX_M1
