#include <loop_device.hxx>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

namespace nuX_M1 {

using namespace Loop;

extern "C" void nuX_M1_InitSemiImplicit(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_M1_InitSemiImplicit;

  *semi_implicit_stage = 2;
}

extern "C" void nuX_M1_SaveSemiImplicitBase(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_M1_SaveSemiImplicitBase;
  DECLARE_CCTK_PARAMETERS;

  const GridDescBaseDevice grid(cctkGH);
  const GF3D2layout layout_cc(cctkGH, {1, 1, 1});
  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        for (int ig = 0; ig < ngroups * nspecies; ++ig) {
          const int i4D = layout_cc.linear(p.i, p.j, p.k, ig);
          rN_base[i4D] = rN[i4D];
          rE_base[i4D] = rE[i4D];
          rFx_base[i4D] = rFx[i4D];
          rFy_base[i4D] = rFy[i4D];
          rFz_base[i4D] = rFz[i4D];
        }
      });
}

extern "C" void nuX_M1_PrepareSemiImplicitStage(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_M1_PrepareSemiImplicitStage;
  DECLARE_CCTK_PARAMETERS;

  const int stage = *semi_implicit_stage;
  if (stage != 1 && stage != 2)
    CCTK_VERROR("Invalid nuX semi-implicit stage %d", stage);
  const CCTK_REAL stage_dt = CCTK_DELTA_TIME / CCTK_REAL(stage);

  const GridDescBaseDevice grid(cctkGH);
  const GF3D2layout layout_cc(cctkGH, {1, 1, 1});
  grid.loop_int_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        for (int ig = 0; ig < ngroups * nspecies; ++ig) {
          const int i4D = layout_cc.linear(p.i, p.j, p.k, ig);
          rN[i4D] = rN_base[i4D] + stage_dt * rN_rhs[i4D];
          rE[i4D] = rE_base[i4D] + stage_dt * rE_rhs[i4D];
          rFx[i4D] = rFx_base[i4D] + stage_dt * rFx_rhs[i4D];
          rFy[i4D] = rFy_base[i4D] + stage_dt * rFy_rhs[i4D];
          rFz[i4D] = rFz_base[i4D] + stage_dt * rFz_rhs[i4D];
        }
      });
}

extern "C" void nuX_M1_AdvanceSemiImplicitStage(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_M1_AdvanceSemiImplicitStage;
  DECLARE_CCTK_PARAMETERS;

  if (*semi_implicit_stage <= 0)
    CCTK_VERROR("Invalid nuX semi-implicit stage %d",
                int(*semi_implicit_stage));
  --*semi_implicit_stage;
}

extern "C" void nuX_M1_FinalizeSemiImplicit(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_M1_FinalizeSemiImplicit;
  DECLARE_CCTK_PARAMETERS;

  if (*semi_implicit_stage != 0)
    CCTK_VERROR("Semi-implicit update ended at stage %d instead of zero",
                int(*semi_implicit_stage));
}

} // namespace nuX_M1
