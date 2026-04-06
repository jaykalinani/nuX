#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "Symmetry.h"

namespace nuX_LeakageBase {

extern "C" void nuX_LeakageBase_SetSym(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_LeakageBase_SetSym;
  DECLARE_CCTK_PARAMETERS;

  int const sym_scal[3] = {1, 1, 1};
  SetCartSymGN(cctkGH, sym_scal, "nuX_LeakageBase::nuX_leakage_optd");
  SetCartSymGN(cctkGH, sym_scal, "nuX_LeakageBase::nuX_leakage_eff_rates");
}

} // namespace nuX_LeakageBase
