#include <cctk.h>
#include <cctk_Arguments.h>
#include <loop_device.hxx>

#include <algorithm>
#include <array>
#include <cmath>

#include "m1_unit_tests.hxx"
#include "nuX_M1_weak_equil.hxx"

namespace nuX_M1 {
namespace {

struct mock_range_t {
  CCTK_REAL min;
  CCTK_REAL max;
};

struct mock_equilibrium_eos_t {
  mock_range_t rgtemp{0.1, 100.0};
  mock_range_t rgye{0.01, 0.6};

  CCTK_HOST CCTK_DEVICE CCTK_REAL mu_lepton_from_rho_temp_ye(
      CCTK_REAL, const CCTK_REAL temp, const CCTK_REAL ye) const {
    return temp * (0.4 + 2.0 * ye);
  }

  CCTK_HOST CCTK_DEVICE CCTK_REAL eps_from_rho_temp_ye(
      CCTK_REAL, const CCTK_REAL temp, const CCTK_REAL ye) const {
    return 0.02 * temp + 0.3 * ye;
  }
};

int weak_equilibrium_unit_tests() {
  constexpr CCTK_REAL rho = 1.0;
  constexpr CCTK_REAL nb = 0.05;
  constexpr CCTK_REAL particle_mass = 930.18478708;
  constexpr CCTK_REAL target[2] = {5.0, 0.2};
  const mock_equilibrium_eos_t eos;

  CCTK_REAL target_residual[2];
  if (func_eq_weak(rho, nb, particle_mass, 0.0, 0.0, target, target_residual,
                   &eos) != 0)
    return 1;
  const CCTK_REAL Yl = target_residual[0];
  const CCTK_REAL eps_total = target_residual[1];

  CCTK_REAL T_eq = 8.0;
  CCTK_REAL Y_eq = 0.3;
  int failures = BetaEquilibriumTrapped(rho, nb, particle_mass, eps_total, Yl,
                                        T_eq, Y_eq, T_eq, Y_eq, &eos);
  failures += std::abs(T_eq - target[0]) > 2.0e-7;
  failures += std::abs(Y_eq - target[1]) > 2.0e-9;

  const CCTK_REAL equilibrium_state[2] = {T_eq, Y_eq};
  CCTK_REAL residual[2];
  failures += func_eq_weak(rho, nb, particle_mass, eps_total, Yl,
                           equilibrium_state, residual, &eos) != 0;
  failures += weak_equilibrium_residual_norm(Yl, eps_total, residual) >
              nu_2DNR_residual_tol;

  CCTK_REAL analytic[2][2];
  failures += jacobi_eq_weak(rho, nb, particle_mass, target, analytic, &eos);
  for (int column = 0; column < 2; ++column) {
    const CCTK_REAL step = 1.0e-6 * std::max(CCTK_REAL(1), target[column]);
    CCTK_REAL plus[2] = {target[0], target[1]};
    CCTK_REAL minus[2] = {target[0], target[1]};
    plus[column] += step;
    minus[column] -= step;
    CCTK_REAL fplus[2], fminus[2];
    failures += func_eq_weak(rho, nb, particle_mass, eps_total, Yl, plus, fplus,
                             &eos) != 0;
    failures += func_eq_weak(rho, nb, particle_mass, eps_total, Yl, minus,
                             fminus, &eos) != 0;
    for (int row = 0; row < 2; ++row) {
      const CCTK_REAL numerical = (fplus[row] - fminus[row]) / (2.0 * step);
      const CCTK_REAL scale = std::max(
          {CCTK_REAL(1), std::abs(numerical), std::abs(analytic[row][column])});
      failures += std::abs(numerical - analytic[row][column]) > 2.0e-5 * scale;
    }
  }

  // Exercise inward Newton steps from both EOS boundaries.  The historical
  // projection discarded all steps at a bound, including valid inward ones.
  for (const auto &boundary_target : {std::array<CCTK_REAL, 2>{0.1005, 0.0105},
                                      std::array<CCTK_REAL, 2>{99.5, 0.599}}) {
    CCTK_REAL boundary_values[2];
    failures += func_eq_weak(rho, nb, particle_mass, 0.0, 0.0,
                             boundary_target.data(), boundary_values, &eos) !=
                0;
    CCTK_REAL boundary_T = boundary_target[0] < 1.0 ? -1.0 : 200.0;
    CCTK_REAL boundary_Y = boundary_target[1] < 0.1 ? -1.0 : 2.0;
    failures += BetaEquilibriumTrapped(
        rho, nb, particle_mass, boundary_values[1], boundary_values[0],
        boundary_T, boundary_Y, boundary_T, boundary_Y, &eos);
    failures += std::abs(boundary_T - boundary_target[0]) > 2.0e-5;
    failures += std::abs(boundary_Y - boundary_target[1]) > 2.0e-7;
  }

  CCTK_REAL invalid_T = target[0];
  CCTK_REAL invalid_Y = target[1];
  failures +=
      BetaEquilibriumTrapped(rho, 0.0, particle_mass, eps_total, Yl, invalid_T,
                             invalid_Y, invalid_T, invalid_Y, &eos) == 0;
  return failures;
}

} // namespace

extern "C" void nuX_M1_UnitTests(CCTK_ARGUMENTS) {
  const int source_failures = source_math_unit_tests(cctkGH);
  const int equilibrium_failures = weak_equilibrium_unit_tests();
  if (source_failures != 0)
    CCTK_VINFO("nuX M1 source-math self-test failed %d checks",
               source_failures);
  if (equilibrium_failures != 0)
    CCTK_VINFO("nuX M1 weak-equilibrium self-test failed %d checks",
               equilibrium_failures);
  const int failures = source_failures + equilibrium_failures;
  if (failures != 0)
    CCTK_VERROR("nuX M1 self-test failed %d checks", failures);
  CCTK_INFO("nuX M1 source/equilibrium self-test passed");
}

} // namespace nuX_M1
