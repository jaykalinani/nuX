#include <cassert>

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

  int method_type;
  const void *const method_p =
      CCTK_ParameterGet("method", "ODESolvers", &method_type);
  assert(method_p);
  assert(method_type == PARAMETER_KEYWORD);
  const char *const method =
      *static_cast<const char *const *>(method_p);
  const bool has_implicit_source =
      CCTK_Equals(method, "IMEX42L") || CCTK_Equals(method, "IMEX32L") ||
      CCTK_Equals(method, "IMEX122") ||
      CCTK_Equals(method, "Implicit Euler");
  if (!has_implicit_source) {
    CCTK_PARAMWARN(
        "nuX_M1 collision terms require an implicit-capable ODESolvers "
        "method (IMEX42L, IMEX32L, IMEX122, or Implicit Euler); "
        "explicit-only methods omit the collision evolution");
  }
  if ((CCTK_Equals(method, "IMEX42L") || CCTK_Equals(method, "IMEX32L")) &&
      source_limiter >= 0) {
    CCTK_PARAMWARN("nuX_M1::source_limiter must be -1 with IMEX42L and "
                   "IMEX32L; limiting individual diagonal source updates "
                   "destroys the L-stable stage cancellation");
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
