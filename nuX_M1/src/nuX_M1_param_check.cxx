#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

namespace nuX_M1 {

extern "C" void nuX_M1_ParamCheck(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_M1_ParamCheck;
  DECLARE_CCTK_PARAMETERS

  if (optimize_prolongation) {
    if (cctk_nghostzones[0] < 4 || cctk_nghostzones[1] < 4 ||
        cctk_nghostzones[2] < 4) {
      CCTK_PARAMWARN("nuX_M1::optimize_prolongation requires at least "
                     "four ghost points");
    }
  }

  const bool is_imex =
      CCTK_Equals(method, "IMEX42L") || CCTK_Equals(method, "IMEX32L") ||
      CCTK_Equals(method, "IMEX122") || CCTK_Equals(method, "Implicit Euler");
  const bool is_semi_implicit = CCTK_Equals(method, "semi-implicit");
  const bool has_implicit_source = is_imex || is_semi_implicit;
  const bool supports_backreaction =
      CCTK_Equals(method, "IMEX42L") || CCTK_Equals(method, "IMEX32L") ||
      is_semi_implicit;
  if (!has_implicit_source) {
    CCTK_PARAMWARN(
        "nuX_M1 collision terms require an implicit-capable ODESolvers "
        "method (semi-implicit, IMEX42L, IMEX32L, IMEX122, or Implicit "
        "Euler); "
        "explicit-only methods omit the collision evolution");
  }
  if (backreact && !supports_backreaction) {
    CCTK_PARAMWARN(
        "nuX_M1::backreact requires ODESolvers::method=semi-implicit, "
        "IMEX42L, or IMEX32L");
  }

  if (is_imex && source_limiter >= 0) {
    CCTK_PARAMWARN("nuX_M1::source_limiter must be -1 with IMEX methods; "
                   "limiting individual diagonal source updates changes "
                   "the additive Runge-Kutta method");
  }

  if (!(rad_eps >= 0.0 && rad_eps < 1.0)) {
    CCTK_PARAMWARN("nuX_M1::rad_eps must satisfy 0 <= rad_eps < 1");
  }

  if (CCTK_Equals(rates_lib, "FakeRates")) {
    if (!CCTK_IsThornActive("nuX_FakeRates")) {
      CCTK_PARAMWARN("nuX_M1 requires the nuX_FakeRates thorn when "
                     "nuX_M1::rates_lib is FakeRates");
    }
    if (set_to_equilibrium || reset_to_equilibrium) {
      CCTK_PARAMWARN("nuX_M1::set_to_equilibrium and "
                     "nuX_M1::reset_to_equilibrium are not supported with "
                     "nuX_M1::rates_lib=FakeRates");
    }
  } else if (CCTK_Equals(rates_lib, "WeakRates")) {
    if (!CCTK_IsThornActive("nuX_WeakRates")) {
      CCTK_PARAMWARN("nuX_M1 requires the nuX_WeakRates thorn when "
                     "nuX_M1::rates_lib is WeakRates");
    }
  } else if (!CCTK_IsThornActive("nuX_NuRates")) {
    CCTK_PARAMWARN("nuX_M1 requires the nuX_NuRates thorn when "
                   "nuX_M1::rates_lib is NuRates");
  }
}

} // namespace nuX_M1
