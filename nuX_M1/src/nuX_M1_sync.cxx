#include <cstdio>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

namespace nuX_M1 {

extern "C" void nuX_M1_Sync(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS_nuX_M1_Sync;
    DECLARE_CCTK_PARAMETERS;

    if(verbose) {
        CCTK_INFO("nuX_M1_Sync");
    }

    // Do nothing, let SYNC handle BCs.
}

}
