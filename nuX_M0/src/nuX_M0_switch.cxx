#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

namespace nuX_M0 {

extern "C" void nuX_M0_InitData(CCTK_ARGUMENTS);

extern "C" void nuX_M0_Switch(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_M0_Switch;
  DECLARE_CCTK_PARAMETERS;

  if ((cctk_iteration - 1) % compute_every != 0) {
    return;
  }

  if (verbose && CCTK_MyProc(cctkGH) == 0) {
    CCTK_INFO("nuX_M0_Switch");
  }

  if (wait_until_time > 0.0 && cctk_time < wait_until_time) {
    return;
  }

  bool const was_on = *nuX_M0_is_on;

  if (bns_sep_threshold > 0.0) {
    CCTK_REAL const *bns_sep_tot = static_cast<CCTK_REAL const *>(
        CCTK_VarDataPtr(cctkGH, 0, "BNSTrackerGen::bns_sep_tot"));
    if (bns_sep_tot && *bns_sep_tot < bns_sep_threshold) {
      *nuX_M0_is_on = 1;
    } else {
      *nuX_M0_is_on = 0;
    }
  } else {
    *nuX_M0_is_on = 1;
  }

  if (!was_on && *nuX_M0_is_on) {
    *nuX_M0_time = cctk_time;
    if (verbose && CCTK_MyProc(cctkGH) == 0) {
      CCTK_INFO("nuX_M0_Switch: M0 is on");
    }
  } else if (was_on && !*nuX_M0_is_on) {
    nuX_M0_InitData(CCTK_PASS_CTOC);
    if (verbose && CCTK_MyProc(cctkGH) == 0) {
      CCTK_INFO("nuX_M0_Switch: M0 is off");
    }
  }
}

} // namespace nuX_M0
