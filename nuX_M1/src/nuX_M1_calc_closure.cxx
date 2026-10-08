#include <cassert>
#include <loop_device.hxx>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#define NUX_M1_CLOSURE_IMPLEMENTATION
#include "nuX_M1_closure.hxx"
#include "nuX_utils.hxx"

namespace nuX_M1 {

using namespace nuX_Utils;
using namespace std;
using namespace Loop;

void CalcClosure(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_M1_CalcClosure;
  DECLARE_CCTK_PARAMETERS;

  if (verbose) {
    CCTK_INFO("nuX_M1_CalcClosure");
  }

  closure_t closure_fun;
  if (CCTK_Equals(closure, "Eddington")) {
    closure_fun = CLOSURE_EDDINGTON;
  } else if (CCTK_Equals(closure, "Kershaw")) {
    closure_fun = CLOSURE_KERSHAW;
  } else if (CCTK_Equals(closure, "Minerbo")) {
    closure_fun = CLOSURE_MINERBO;
  } else if (CCTK_Equals(closure, "thin")) {
    closure_fun = CLOSURE_THIN;
  } else {
    char msg[BUFSIZ];
    snprintf(msg, BUFSIZ, "Unknown closure \"%s\"", closure);
    CCTK_ERROR(msg);
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
  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        const int ijk = layout_cc.linear(p.i, p.j, p.k);
        if (nuX_m1_mask[ijk]) {
          for (int ig = 0; ig < nspecies * ngroups; ++ig) {
            int const i4D = layout_cc.linear(p.i, p.j, p.k, ig);
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
            chi[i4D] = 0;
          }
          return;
        }

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

        tensor::generic<CCTK_REAL, 4, 1> H_d;
        tensor::generic<CCTK_REAL, 4, 1> F_d;
        tensor::generic<CCTK_REAL, 4, 1> beta_u;
        tensor::symmetric2<CCTK_REAL, 4, 2> P_dd;
        tensor::symmetric2<CCTK_REAL, 4, 2> rT_dd;
        geom.get_shift_vec(p, &beta_u);

        for (int ig = 0; ig < nspecies * ngroups; ++ig) {
          int const i4D = layout_cc.linear(p.i, p.j, p.k, ig);

          pack_F_d(beta_u(1), beta_u(2), beta_u(3), rFx[i4D], rFy[i4D],
                   rFz[i4D], &F_d);

          CCTK_REAL E_closure = rE[i4D];
          repair_moments(g_uu, &E_closure, &F_d, rad_E_floor, rad_eps);

          assert(isfinite(E_closure));
          assert(isfinite(F_d(0)));
          assert(isfinite(F_d(1)));
          assert(isfinite(F_d(2)));
          assert(isfinite(F_d(3)));
          assert(isfinite(tensor::dot(g_uu, F_d, F_d)));
          assert(isfinite(W));

          calc_closure(cctkGH, p.i, p.j, p.k, ig, closure_fun, g_dd, g_uu, n_d,
                       W, u_u, v_d, proj_ud, E_closure, F_d, &chi[i4D], &P_dd,
                       closure_epsilon, closure_maxiter, use_fallback != 0);
          unpack_P_dd(P_dd, &rPxx[i4D], &rPxy[i4D], &rPxz[i4D], &rPyy[i4D],
                      &rPyz[i4D], &rPzz[i4D]);
          assert(isfinite(rPxx[i4D]));
          assert(isfinite(rPxy[i4D]));
          assert(isfinite(rPxz[i4D]));
          assert(isfinite(rPyy[i4D]));
          assert(isfinite(rPyz[i4D]));
          assert(isfinite(rPzz[i4D]));

          assemble_rT(n_d, E_closure, F_d, P_dd, &rT_dd);

          rJ[i4D] = calc_J_from_rT(rT_dd, u_u);
          calc_H_from_rT(rT_dd, u_u, proj_ud, &H_d);

          // J and H_a are projections of the same stress tensor used to store
          // P_ab. Repairing them independently here, as THC does, would make
          // the stored comoving moments inconsistent with that tensor. The
          // variable M1 closures must preserve the comoving moment cone. The
          // fixed Eddington closure is only a diffusion approximation and does
          // not preserve that cone for arbitrary flux-dominated lab states;
          // those can occur at the radiation floor in diffusion tests.
          if (!comoving_state_is_acceptable(closure_fun, g_uu, rJ[i4D],
                                             H_d))
            device_abort();

          unpack_H_d(H_d, &rHt[i4D], &rHx[i4D], &rHy[i4D], &rHz[i4D]);
          assert(isfinite(rHt[i4D]));
          assert(isfinite(rHx[i4D]));
          assert(isfinite(rHy[i4D]));
          assert(isfinite(rHz[i4D]));

          CCTK_REAL const Gamma = compute_Gamma(W, v_u, rJ[i4D], E_closure, F_d,
                                                rad_E_floor, rad_eps);
          assert(Gamma > 0);
          rnnu[i4D] = max(rN[i4D], rad_N_floor) / Gamma;
        }
      });
}

extern "C" void nuX_M1_CalcClosure(CCTK_ARGUMENTS) {
  CalcClosure(cctkGH);
}

} // namespace nuX_M1
