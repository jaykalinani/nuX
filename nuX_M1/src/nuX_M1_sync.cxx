#include <cstdio>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

namespace nuX_M1 {

// Sync conserved radiation fields (rN, rE, rF) before computing
// closures/fluxes.
extern "C" void nuX_M1_SyncConserved(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_M1_SyncConserved;
  DECLARE_CCTK_PARAMETERS;

  if (verbose) {
    CCTK_INFO("nuX_M1_SyncConserved");
  }

  // Do nothing, let SYNC handle BCs.
}

// Sync derived radiation fields (closure outputs/opacities) before fluxes.
extern "C" void nuX_M1_SyncDerived(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_M1_SyncDerived;
  DECLARE_CCTK_PARAMETERS;

  if (verbose) {
    CCTK_INFO("nuX_M1_SyncDerived");
  }

  // Do nothing, let SYNC handle BCs.
}

extern "C" void nuX_M1_Sync(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_M1_Sync;
  DECLARE_CCTK_PARAMETERS;

  if (verbose) {
    CCTK_INFO("nuX_M1_Sync");
  }

  // Do nothing, let SYNC handle BCs.
}

} // namespace nuX_M1
