#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

namespace nuX_M1 {

extern "C" void nuX_M1_ParamCheck(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_M1_ParamCheck;
  DECLARE_CCTK_PARAMETERS;

  constexpr int max_groupspec_for_sources = 3;
  const int groupspec = nspecies * ngroups;

  if (groupspec > max_groupspec_for_sources) {
    CCTK_PARAMWARN("nuX_M1::nspecies * nuX_M1::ngroups exceeds the local "
                   "source-update accumulator size");
  }

  if (include_fluxes &&
      (cctk_nghostzones[0] < 2 || cctk_nghostzones[1] < 2 ||
       cctk_nghostzones[2] < 2)) {
    CCTK_PARAMWARN("nuX_M1 transport requires at least two ghost zones in "
                   "every direction");
  }

  if (include_GR_sources && grsource_spatial_order == 4 &&
      (cctk_nghostzones[0] < 3 || cctk_nghostzones[1] < 3 ||
       cctk_nghostzones[2] < 3)) {
    CCTK_PARAMWARN("nuX_M1::grsource_spatial_order=4 requires at least "
                   "three ghost zones in every direction");
  }

  if (optimize_prolongation) {
    if (cctk_nghostzones[0] < 4 || cctk_nghostzones[1] < 4 ||
        cctk_nghostzones[2] < 4) {
      CCTK_PARAMWARN("nuX_M1::optimize_prolongation requires at least "
                     "four ghost points");
    }
  }
}

} // namespace nuX_M1
