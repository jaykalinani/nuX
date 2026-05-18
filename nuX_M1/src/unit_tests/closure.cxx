#include <cctk.h>

#include "../test.hxx"
#include "nuX_M1_closure.hxx"

#include <random>

namespace {

using namespace nuX_Utils;

struct FlatState {
  tensor::metric<4> g_dd;
  tensor::inv_metric<4> g_uu;
  tensor::generic<CCTK_REAL, 4, 1> n_d;
  tensor::generic<CCTK_REAL, 4, 1> u_u;
  tensor::generic<CCTK_REAL, 4, 1> v_d;
  tensor::generic<CCTK_REAL, 4, 2> proj_ud;

  FlatState() {
    for (int a = 0; a < 4; ++a) {
      n_d(a) = 0.0;
      u_u(a) = 0.0;
      v_d(a) = 0.0;
      for (int b = 0; b < 4; ++b) {
        g_dd(a, b) = 0.0;
        g_uu(a, b) = 0.0;
      }
    }

    g_dd(0, 0) = -1.0;
    g_dd(1, 1) = 1.0;
    g_dd(2, 2) = 1.0;
    g_dd(3, 3) = 1.0;
    g_uu(0, 0) = -1.0;
    g_uu(1, 1) = 1.0;
    g_uu(2, 2) = 1.0;
    g_uu(3, 3) = 1.0;

    n_d(0) = -1.0;
    u_u(0) = 1.0;
    tensor::generic<CCTK_REAL, 4, 1> u_d;
    tensor::contract(g_dd, u_u, &u_d);
    nuX_M1::calc_proj(u_d, u_u, &proj_ud);
  }
};

} // namespace

void nuX_M1_Tests::test_closure(std::mt19937_64 &engine, int repetitions) {
  using namespace nuX_M1;
  using std::uniform_real_distribution;

  uniform_real_distribution<CCTK_REAL> xi_values{0.0, 1.0};
  for (int n = 0; n < repetitions; ++n) {
    const CCTK_REAL xi = xi_values(engine);
    require_close(eval_closure(CLOSURE_EDDINGTON, xi), 1.0 / 3.0, 0.0,
                  "Eddington closure is isotropic");
    require_close(eval_closure(CLOSURE_THIN, xi), 1.0, 0.0,
                  "thin closure is free-streaming");
    require_close(eval_closure(CLOSURE_KERSHAW, xi),
                  1.0 / 3.0 + 2.0 / 3.0 * xi * xi, 0.0,
                  "Kershaw closure formula");
    require_close(eval_closure(CLOSURE_MINERBO, xi),
                  1.0 / 3.0 +
                      xi * xi * (6.0 - 2.0 * xi + 6.0 * xi * xi) / 15.0,
                  1.0e-15, "Minerbo closure formula");
  }

  FlatState flat;
  tensor::generic<CCTK_REAL, 4, 1> F_d;
  F_d(0) = 0.0;
  F_d(1) = 0.0;
  F_d(2) = 0.0;
  F_d(3) = 0.0;
  tensor::symmetric2<CCTK_REAL, 4, 2> P_dd;
  apply_closure(flat.g_dd, flat.g_uu, flat.n_d, 1.0, flat.u_u, flat.v_d,
                flat.proj_ud, 3.0, F_d, 1.0 / 3.0, &P_dd);

  require_close(P_dd(0, 0), 0.0, 0.0,
                "isotropic rest-frame closure has no P_00");
  require_close(P_dd(1, 1), 1.0, 1.0e-15,
                "isotropic rest-frame closure gives Pxx=E/3");
  require_close(P_dd(2, 2), 1.0, 1.0e-15,
                "isotropic rest-frame closure gives Pyy=E/3");
  require_close(P_dd(3, 3), 1.0, 1.0e-15,
                "isotropic rest-frame closure gives Pzz=E/3");
  require_close(P_dd(1, 2), 0.0, 0.0,
                "isotropic rest-frame closure has no off-diagonal pressure");

  F_d(1) = 0.5;
  CCTK_REAL chi = 0.0;
  CCTK_REAL status = -1.0;
  calc_closure(nullptr, 0, 0, 0, 0, CLOSURE_MINERBO, flat.g_dd, flat.g_uu,
               flat.n_d, 1.0, flat.u_u, flat.v_d, flat.proj_ud, 2.0, F_d,
               &chi, &P_dd, 1.0e-12, 128, true, &status);
  require_close(status, NUX_M1_CLOSURE_OK, 0.0,
                "calc_closure reports success for a bracketed root");
  require_close(chi, minerbo(0.25), 1.0e-12,
                "calc_closure recovers rest-frame Minerbo chi");

  F_d(1) = 0.0;
  chi = 0.0;
  status = -1.0;
  calc_closure(nullptr, 0, 0, 0, 0, CLOSURE_MINERBO, flat.g_dd, flat.g_uu,
               flat.n_d, 1.0, flat.u_u, flat.v_d, flat.proj_ud, 0.0, F_d,
               &chi, &P_dd, 1.0e-12, 128, true, &status);
  require_close(status, NUX_M1_CLOSURE_OK, 0.0,
                "zero-state closure endpoint solve reports success");
  require_close(chi, 1.0, 0.0,
                "zero-state closure follows THC/GSL endpoint handling");
}
