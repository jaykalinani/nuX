#include <cctk.h>

#include "../test.hxx"
#include <loop_device.hxx>
#include "../../nuX_Utils/src/nuX_roots.hxx"

#include <cmath>
#include <limits>
#include <random>

void nuX_M1_Tests::test_roots(std::mt19937_64 &engine, int repetitions) {
  using namespace nuX_Utils::roots;

  const auto parabola = [](CCTK_REAL x) { return x * x - 2.0; };
  brent_solver<CCTK_REAL> solver;
  int ierr = solver_set(&solver, parabola, CCTK_REAL(0.0), CCTK_REAL(2.0));
  require(ierr == code(status::success), "brent solver accepts a bracket");

  int interval_status = code(status::continue_iter);
  for (int n = 0; n < 128 && interval_status == code(status::continue_iter);
       ++n) {
    ierr = solver_iterate(&solver, parabola);
    require(ierr == code(status::success), "brent solver iteration succeeds");
    interval_status =
        test_interval(solver_x_lower(&solver), solver_x_upper(&solver),
                      CCTK_REAL(1.0e-13), CCTK_REAL(0.0));
  }
  require(interval_status == code(status::success),
          "brent solver converges on x^2-2");
  require_close(solver_root(&solver), std::sqrt(CCTK_REAL(2.0)), 1.0e-12,
                "brent solver root value");

  const auto endpoint = [](CCTK_REAL x) { return x; };
  ierr = solver_set(&solver, endpoint, CCTK_REAL(0.0), CCTK_REAL(2.0));
  require(ierr == code(status::success),
          "brent solver accepts an endpoint root");
  ierr = solver_iterate(&solver, endpoint);
  require(ierr == code(status::success),
          "brent solver endpoint iteration succeeds");
  require_close(solver_root(&solver), 0.0, 0.0,
                "brent solver returns endpoint root");

  const auto no_bracket = [](CCTK_REAL x) { return x * x + 1.0; };
  ierr = solver_set(&solver, no_bracket, CCTK_REAL(-1.0), CCTK_REAL(1.0));
  require(ierr == code(status::einval), "brent solver rejects no bracket");

  const auto bad_func = [](CCTK_REAL) {
    return std::numeric_limits<CCTK_REAL>::quiet_NaN();
  };
  CCTK_REAL y = 0.0;
  ierr = safe_func_call(bad_func, CCTK_REAL(0.0), &y);
  require(ierr == code(status::ebadfunc),
          "safe_func_call catches nonfinite values");

  using Vec2 = Arith::vec<CCTK_REAL, 2>;
  using Mat2 = Arith::mat<CCTK_REAL, 2>;
  auto residual = [](Vec2 x, Vec2 &f) {
    f(0) = x(0) - 3.0;
    f(1) = 2.0 * x(1) + 4.0;
    return code(status::success);
  };
  auto jacobian = [](Vec2, Mat2 &J) {
    J(0, 0) = 1.0;
    J(0, 1) = 0.0;
    J(1, 0) = 0.0;
    J(1, 1) = 2.0;
    return code(status::success);
  };
  auto residual_jacobian = [&](Vec2 x, Vec2 &f, Mat2 &J) {
    residual(x, f);
    return jacobian(x, J);
  };

  (void)engine;
  (void)repetitions;

  hybridsj_solver<CCTK_REAL, 2> hybrid;
  Vec2 x0{0.0, 0.0};
  ierr = hybridsj_set(&hybrid, residual_jacobian, x0);
  require(ierr == code(status::success), "hybridsj setup succeeds");

  int residual_status = test_residual(hybrid.f, CCTK_REAL(1.0e-13));
  for (int n = 0; n < 64 && residual_status == code(status::continue_iter);
       ++n) {
    ierr = hybridsj_iterate(&hybrid, residual, jacobian);
    require(ierr == code(status::success),
            "hybridsj iteration succeeds for a linear system");
    residual_status = test_residual(hybrid.f, CCTK_REAL(1.0e-13));
  }
  require(residual_status == code(status::success),
          "hybridsj converges on a linear system");
  require_close(hybrid.x(0), 3.0, 1.0e-12, "hybridsj x0 solution");
  require_close(hybrid.x(1), -2.0, 1.0e-12, "hybridsj x1 solution");
}
