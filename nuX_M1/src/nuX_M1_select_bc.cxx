#include <cstdio>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

namespace nuX_M1 {

extern "C" void nuX_M1_SelectBC(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  char vname[BUFSIZ];

  if (verbose) {
    CCTK_INFO("nuX_M1_SelectBC");
  }

  char const *mybc;
  if (CCTK_Equals(bc_type, "flat")) {
    mybc = "flat";
  } else {
    mybc = "none";
  }

  assert(cctk_nghostzones[0] == cctk_nghostzones[1]);
  assert(cctk_nghostzones[1] == cctk_nghostzones[2]);

  int ierr = 0;

#define nuX_SELECT_VAR_FOR_BC(VARIABLE)                                        \
  ierr += Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, cctk_nghostzones[0], \
                                  -1, VARIABLE, mybc)

  for (int ig = 0; ig < nspecies; ++ig) {
    std::snprintf(vname, BUFSIZ, "nuX_M1::rN[%d]", ig);
    nuX_SELECT_VAR_FOR_BC(vname);

    std::snprintf(vname, BUFSIZ, "nuX_M1::rE[%d]", ig);
    nuX_SELECT_VAR_FOR_BC(vname);

    std::snprintf(vname, BUFSIZ, "nuX_M1::rFx[%d]", ig);
    nuX_SELECT_VAR_FOR_BC(vname);

    std::snprintf(vname, BUFSIZ, "nuX_M1::rFy[%d]", ig);
    nuX_SELECT_VAR_FOR_BC(vname);

    std::snprintf(vname, BUFSIZ, "nuX_M1::rFz[%d]", ig);
    nuX_SELECT_VAR_FOR_BC(vname);
  }

#undef nuX_SELECT_VAR_FOR_BC

  if (ierr) {
    CCTK_ERROR("Failed to select the BCs");
  }
}

extern "C" void nuX_M1_ControlProlongation(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_M1_ControlProlongation;
  DECLARE_CCTK_PARAMETERS
  int ierr = 0;

  if (verbose) {
    CCTK_INFO("nuX_M1_ControlProlongation");
  }

  // Only prolongate at the end of the full timestep
  if (*TimeIntegratorStage == 1) {
    ierr += EnableProlongating(0);
  } else if (*TimeIntegratorStage == 0) {
    ierr += EnableProlongating(1);
  } else {
    CCTK_VERROR("Unepxected TimeIntegratorStage: %d expected 1 or 0",
                (int)*TimeIntegratorStage);
  }

  if (ierr) {
    CCTK_ERROR("Failed to set prolongation");
  }
}

} // namespace nuX_M1
