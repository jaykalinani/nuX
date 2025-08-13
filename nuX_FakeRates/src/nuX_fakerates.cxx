#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include <loop_device.hxx>
#include <float.h>

#include "nuX_fakerates.hxx"

namespace nuX_FakeRates {
using namespace nuX_Rates;

CCTK_HOST void FakeRatesDef::init() {
  DECLARE_CCTK_PARAMETERS;

  kabs_nue = kappa_abs_nue;
  kabs_nua = kappa_abs_nua;
  kabs_nux = kappa_abs_nux;
  kabs_anux = kappa_abs_anux;

  kscat_nue = kappa_scat_nue;
  kscat_nua = kappa_scat_nua;
  kscat_nux = kappa_scat_nux;
  kscat_anux = kappa_scat_anux;

  et_nue = eta_nue;
  et_nua = eta_nua;
  et_nux = eta_nux;
  et_anux = eta_anux;
  return;
}
} // namespace

