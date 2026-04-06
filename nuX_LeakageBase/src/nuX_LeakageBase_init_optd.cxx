#include <loop_device.hxx>

#include <cmath>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

namespace nuX_LeakageBase {

using namespace Loop;

extern "C" void nuX_LeakageBase_InitOpticalDepthSimple(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_nuX_LeakageBase_InitOpticalDepthSimple;
  DECLARE_CCTK_PARAMETERS;
  using std::isfinite;

  if (verbose && CCTK_MyProc(cctkGH) == 0) {
    CCTK_INFO("nuX_LeakageBase_InitOpticalDepthSimple");
  }

  grid.loop_int_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        CCTK_REAL rho_local = rho(p.I);
        if (!isfinite(rho_local) || rho_local < atmo_rho) {
          rho_local = atmo_rho;
        }
        CCTK_REAL const tau = rho_local / atmo_rho - 1.0;
        optd_0_nue(p.I) = tau;
        optd_0_nua(p.I) = tau;
        optd_0_nux(p.I) = tau;
        optd_1_nue(p.I) = tau;
        optd_1_nua(p.I) = tau;
        optd_1_nux(p.I) = tau;
        optd_0_nue_p(p.I) = tau;
        optd_0_nua_p(p.I) = tau;
        optd_0_nux_p(p.I) = tau;
        optd_1_nue_p(p.I) = tau;
        optd_1_nua_p(p.I) = tau;
        optd_1_nux_p(p.I) = tau;
        optd_0_nue_p_p(p.I) = tau;
        optd_0_nua_p_p(p.I) = tau;
        optd_0_nux_p_p(p.I) = tau;
        optd_1_nue_p_p(p.I) = tau;
        optd_1_nua_p_p(p.I) = tau;
        optd_1_nux_p_p(p.I) = tau;
      });
}

} // namespace nuX_LeakageBase
