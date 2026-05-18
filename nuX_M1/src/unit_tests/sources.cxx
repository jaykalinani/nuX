#include <cctk.h>

#include "../test.hxx"
#include "nuX_M1_sources.hxx"

#include <random>

namespace {

using namespace nuX_Utils;

struct FlatSourceState {
  tensor::metric<4> g_dd;
  tensor::inv_metric<4> g_uu;
  tensor::generic<CCTK_REAL, 4, 1> n_d;
  tensor::generic<CCTK_REAL, 4, 1> n_u;
  tensor::generic<CCTK_REAL, 4, 2> gamma_ud;
  tensor::generic<CCTK_REAL, 4, 1> u_d;
  tensor::generic<CCTK_REAL, 4, 1> u_u;
  tensor::generic<CCTK_REAL, 4, 1> v_d;
  tensor::generic<CCTK_REAL, 4, 1> v_u;
  tensor::generic<CCTK_REAL, 4, 2> proj_ud;

  FlatSourceState() {
    for (int a = 0; a < 4; ++a) {
      n_d(a) = 0.0;
      n_u(a) = 0.0;
      u_d(a) = 0.0;
      u_u(a) = 0.0;
      v_d(a) = 0.0;
      v_u(a) = 0.0;
      for (int b = 0; b < 4; ++b) {
        g_dd(a, b) = 0.0;
        g_uu(a, b) = 0.0;
        gamma_ud(a, b) = 0.0;
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
    n_u(0) = 1.0;
    u_u(0) = 1.0;
    tensor::contract(g_dd, u_u, &u_d);
    nuX_M1::calc_proj(u_d, u_u, &proj_ud);

    for (int a = 0; a < 4; ++a) {
      for (int b = 0; b < 4; ++b) {
        gamma_ud(a, b) = tensor::delta(a, b) + n_u(a) * n_d(b);
      }
    }
  }
};

} // namespace

void nuX_M1_Tests::test_sources(std::mt19937_64 &engine, int repetitions) {
  using namespace nuX_M1;

  (void)engine;
  (void)repetitions;

  require(source_rates_are_zero(0.0, 0.0, 0.0, 0.0, 0.0),
          "source_rates_are_zero detects exactly zero rates");
  require(!source_rates_are_zero(0.0, 0.0, 0.0, 0.0, 1.0),
          "source_rates_are_zero detects scattering");
  require(source_is_nonstiff(0.5, 1.0, 1.0),
          "source_is_nonstiff accepts tau < 1");
  require(!source_is_nonstiff(2.0, 1.0, 0.0),
          "source_is_nonstiff rejects tau >= 1");
  require(source_uses_thick_limit(1.0, 10.0, 0.0, 5.0),
          "source_uses_thick_limit detects optically thick cells");
  require(!source_uses_thick_limit(1.0, 10.0, 0.0, -1.0),
          "negative thick-limit parameter disables thick limit");
  require(source_uses_scat_limit(1.0, 10.0, 5.0),
          "source_uses_scat_limit detects large scattering depth");
  require(!source_uses_scat_limit(1.0, 10.0, -1.0),
          "negative scattering-limit parameter disables scattering limit");

  tensor::generic<CCTK_REAL, 4, 1> Fstar_d;
  Fstar_d(0) = 0.0;
  Fstar_d(1) = 0.1;
  Fstar_d(2) = -0.2;
  Fstar_d(3) = 0.3;
  require(source_can_use_no_source_fallback(1.0e-4, 1.0, 1.0, 2.0, Fstar_d),
          "no-source fallback is allowed for tiny optical depth");
  require(!source_can_use_no_source_fallback(1.0e-2, 1.0, 1.0, 2.0, Fstar_d),
          "no-source fallback is rejected for finite optical depth");

  CCTK_REAL chi = 0.0;
  CCTK_REAL Enew = 0.0;
  tensor::generic<CCTK_REAL, 4, 1> Fnew_d;
  source_apply_no_source_fallback(2.0, Fstar_d, &chi, &Enew, &Fnew_d);
  require_close(chi, 1.0 / 3.0, 0.0,
                "no-source fallback uses Eddington closure");
  require_close(Enew, 2.0, 0.0, "no-source fallback preserves E");
  require_close(Fnew_d(1), Fstar_d(1), 0.0,
                "no-source fallback preserves Fx");
  require_close(Fnew_d(2), Fstar_d(2), 0.0,
                "no-source fallback preserves Fy");
  require_close(Fnew_d(3), Fstar_d(3), 0.0,
                "no-source fallback preserves Fz");

  FlatSourceState flat;
  tensor::generic<CCTK_REAL, 4, 1> Fold_d = Fstar_d;
  CCTK_REAL closure_status = 0.0;
  SourceUpdateContext ctx(
      nullptr, 0, 0, 0, 0, 1.0e-12, 128, true, 0.01, 1.0, flat.g_dd,
      flat.g_uu, flat.n_d, flat.n_u, flat.gamma_ud, flat.u_d, flat.u_u,
      flat.v_d, flat.v_u, flat.proj_ud, 1.0, 2.0, Fold_d, 2.0, Fstar_d, 0.0,
      0.0, 0.0, &closure_status);

  chi = 1.0 / 3.0;
  Enew = 2.0;
  Fnew_d = Fstar_d;
  const int status =
      source_update(ctx, CLOSURE_MINERBO, &chi, &Enew, &Fnew_d, 20.0, -1.0,
                    64, 1.0e-15, 1.0e-12);
  require(status != NUX_M1_SOURCE_FAIL,
          "source_update succeeds for zero-rate cells");
  require_close(Enew, 2.0, 1.0e-14,
                "source_update preserves E for zero-rate cells");
  require_close(Fnew_d(1), Fstar_d(1), 1.0e-14,
                "source_update preserves Fx for zero-rate cells");
  require_close(Fnew_d(2), Fstar_d(2), 1.0e-14,
                "source_update preserves Fy for zero-rate cells");
  require_close(Fnew_d(3), Fstar_d(3), 1.0e-14,
                "source_update preserves Fz for zero-rate cells");
}
