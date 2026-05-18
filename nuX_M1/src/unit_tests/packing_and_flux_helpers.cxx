#include <cctk.h>

#include "../test.hxx"
#include "nuX_M1_closure.hxx"

#include <random>

namespace {

using nuX_M1_Tests::require_close;

} // namespace

void nuX_M1_Tests::test_packing_and_flux_helpers(std::mt19937_64 &engine,
                                                 int repetitions) {
  using namespace nuX_M1;
  using namespace nuX_Utils;
  using std::uniform_real_distribution;

  uniform_real_distribution<CCTK_REAL> values{-2.0, 2.0};

  for (int n = 0; n < repetitions; ++n) {
    const CCTK_REAL betax = values(engine);
    const CCTK_REAL betay = values(engine);
    const CCTK_REAL betaz = values(engine);
    const CCTK_REAL Fx = values(engine);
    const CCTK_REAL Fy = values(engine);
    const CCTK_REAL Fz = values(engine);

    tensor::generic<CCTK_REAL, 4, 1> F_d;
    pack_F_d(betax, betay, betaz, Fx, Fy, Fz, &F_d);

    require_close(F_d(0), betax * Fx + betay * Fy + betaz * Fz, 0.0,
                  "pack_F_d stores F_0 = beta_i F^i");
    CCTK_REAL Fx_out = 0.0;
    CCTK_REAL Fy_out = 0.0;
    CCTK_REAL Fz_out = 0.0;
    unpack_F_d(F_d, &Fx_out, &Fy_out, &Fz_out);
    require_close(Fx_out, Fx, 0.0, "unpack_F_d returns Fx");
    require_close(Fy_out, Fy, 0.0, "unpack_F_d returns Fy");
    require_close(Fz_out, Fz, 0.0, "unpack_F_d returns Fz");
  }

  tensor::symmetric2<CCTK_REAL, 4, 2> P_dd;
  const CCTK_REAL betax = 0.1;
  const CCTK_REAL betay = -0.2;
  const CCTK_REAL betaz = 0.3;
  const CCTK_REAL Pxx = 1.0;
  const CCTK_REAL Pxy = 2.0;
  const CCTK_REAL Pxz = 3.0;
  const CCTK_REAL Pyy = 4.0;
  const CCTK_REAL Pyz = 5.0;
  const CCTK_REAL Pzz = 6.0;
  pack_P_dd(betax, betay, betaz, Pxx, Pxy, Pxz, Pyy, Pyz, Pzz, &P_dd);

  const CCTK_REAL Pbetax = Pxx * betax + Pxy * betay + Pxz * betaz;
  const CCTK_REAL Pbetay = Pxy * betax + Pyy * betay + Pyz * betaz;
  const CCTK_REAL Pbetaz = Pxz * betax + Pyz * betay + Pzz * betaz;
  require_close(P_dd(0, 0), Pbetax * betax + Pbetay * betay + Pbetaz * betaz,
                0.0, "pack_P_dd stores shifted P_00");
  require_close(P_dd(0, 1), Pbetax, 0.0, "pack_P_dd stores P_01");
  require_close(P_dd(0, 2), Pbetay, 0.0, "pack_P_dd stores P_02");
  require_close(P_dd(0, 3), Pbetaz, 0.0, "pack_P_dd stores P_03");

  CCTK_REAL Pxx_out = 0.0;
  CCTK_REAL Pxy_out = 0.0;
  CCTK_REAL Pxz_out = 0.0;
  CCTK_REAL Pyy_out = 0.0;
  CCTK_REAL Pyz_out = 0.0;
  CCTK_REAL Pzz_out = 0.0;
  unpack_P_dd(P_dd, &Pxx_out, &Pxy_out, &Pxz_out, &Pyy_out, &Pyz_out,
              &Pzz_out);
  require_close(Pxx_out, Pxx, 0.0, "unpack_P_dd returns Pxx");
  require_close(Pxy_out, Pxy, 0.0, "unpack_P_dd returns Pxy");
  require_close(Pxz_out, Pxz, 0.0, "unpack_P_dd returns Pxz");
  require_close(Pyy_out, Pyy, 0.0, "unpack_P_dd returns Pyy");
  require_close(Pyz_out, Pyz, 0.0, "unpack_P_dd returns Pyz");
  require_close(Pzz_out, Pzz, 0.0, "unpack_P_dd returns Pzz");

  tensor::generic<CCTK_REAL, 4, 1> u_u;
  tensor::generic<CCTK_REAL, 4, 1> u_d;
  u_u(0) = 1.0;
  u_u(1) = 0.0;
  u_u(2) = 0.0;
  u_u(3) = 0.0;
  u_d(0) = -1.0;
  u_d(1) = 0.0;
  u_d(2) = 0.0;
  u_d(3) = 0.0;
  tensor::generic<CCTK_REAL, 4, 2> proj_ud;
  calc_proj(u_d, u_u, &proj_ud);
  require_close(proj_ud(0, 0), 0.0, 0.0, "calc_proj removes time component");
  require_close(proj_ud(1, 1), 1.0, 0.0, "calc_proj keeps x component");
  require_close(proj_ud(2, 2), 1.0, 0.0, "calc_proj keeps y component");
  require_close(proj_ud(3, 3), 1.0, 0.0, "calc_proj keeps z component");

  tensor::generic<CCTK_REAL, 4, 1> beta_u;
  beta_u(0) = 0.0;
  beta_u(1) = 0.1;
  beta_u(2) = -0.2;
  beta_u(3) = 0.3;
  tensor::generic<CCTK_REAL, 4, 1> F_u;
  F_u(0) = 0.0;
  F_u(1) = 1.0;
  F_u(2) = 2.0;
  F_u(3) = 3.0;
  require_close(calc_E_flux(0.9, beta_u, 4.0, F_u, 2), 2.6, 0.0,
                "calc_E_flux computes alpha F^i - beta^i E");

  tensor::generic<CCTK_REAL, 4, 2> P_ud;
  for (int a = 0; a < 4; ++a) {
    for (int b = 0; b < 4; ++b) {
      P_ud(a, b) = 10.0 * a + b;
    }
  }
  tensor::generic<CCTK_REAL, 4, 1> F_flux_d;
  pack_F_d(beta_u(1), beta_u(2), beta_u(3), 1.0, 2.0, 3.0, &F_flux_d);
  require_close(calc_F_flux(0.9, beta_u, F_flux_d, P_ud, 3, 1),
                0.9 * P_ud(3, 1) - beta_u(3) * F_flux_d(1), 0.0,
                "calc_F_flux computes alpha P^i_j - beta^i F_j");

  tensor::generic<CCTK_REAL, 4, 1> v_u;
  pack_v_u(0.2, 0.0, 0.0, &v_u);
  require_close(compute_Gamma(1.1, v_u, 2.0, 5.0, F_flux_d, 1.0e-15,
                              1.0e-5),
                1.1 * (5.0 / 2.0) *
                    (1.0 - (tensor::dot(F_flux_d, v_u) / 5.0)),
                1.0e-14, "compute_Gamma handles finite radiation state");
  require_close(compute_Gamma(1.1, v_u, 0.0, 5.0, F_flux_d, 1.0e-15,
                              1.0e-5),
                1.0, 0.0, "compute_Gamma falls back below floor");
}
