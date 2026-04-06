#include <loop_device.hxx>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

namespace nuX_LeakageBase {

using namespace Loop;

extern "C" void nuX_LeakageBase_InitData(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_nuX_LeakageBase_InitData;
  DECLARE_CCTK_PARAMETERS;

  if (verbose && CCTK_MyProc(cctkGH) == 0) {
    CCTK_INFO("nuX_LeakageBase_InitData");
  }

  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        kappa_0_nue(p.I) = 0.0;
        kappa_0_nua(p.I) = 0.0;
        kappa_0_nux(p.I) = 0.0;
        kappa_1_nue(p.I) = 0.0;
        kappa_1_nua(p.I) = 0.0;
        kappa_1_nux(p.I) = 0.0;

        R_free_nue(p.I) = 0.0;
        R_free_nua(p.I) = 0.0;
        R_free_nux(p.I) = 0.0;
        Q_free_nue(p.I) = 0.0;
        Q_free_nua(p.I) = 0.0;
        Q_free_nux(p.I) = 0.0;

        R_eff_nue(p.I) = 0.0;
        R_eff_nua(p.I) = 0.0;
        R_eff_nux(p.I) = 0.0;
        Q_eff_nue(p.I) = 0.0;
        Q_eff_nua(p.I) = 0.0;
        Q_eff_nux(p.I) = 0.0;

        abs_number(p.I) = 0.0;
        abs_energy(p.I) = 0.0;
        luminosity_nue(p.I) = 0.0;
        luminosity_nua(p.I) = 0.0;
        luminosity_nux(p.I) = 0.0;

        optd_0_nue(p.I) = 0.0;
        optd_0_nua(p.I) = 0.0;
        optd_0_nux(p.I) = 0.0;
        optd_1_nue(p.I) = 0.0;
        optd_1_nua(p.I) = 0.0;
        optd_1_nux(p.I) = 0.0;

        optd_0_nue_p(p.I) = 0.0;
        optd_0_nua_p(p.I) = 0.0;
        optd_0_nux_p(p.I) = 0.0;
        optd_1_nue_p(p.I) = 0.0;
        optd_1_nua_p(p.I) = 0.0;
        optd_1_nux_p(p.I) = 0.0;

        optd_0_nue_p_p(p.I) = 0.0;
        optd_0_nua_p_p(p.I) = 0.0;
        optd_0_nux_p_p(p.I) = 0.0;
        optd_1_nue_p_p(p.I) = 0.0;
        optd_1_nua_p_p(p.I) = 0.0;
        optd_1_nux_p_p(p.I) = 0.0;
      });
}

} // namespace nuX_LeakageBase
