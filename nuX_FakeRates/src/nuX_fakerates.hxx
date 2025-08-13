#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include <loop_device.hxx>
#include <float.h>

#include "m1_opacities.hpp"

namespace nuX_FakeRates {

class FakeRatesDef {
public:
  CCTK_REAL kabs_nue;
  CCTK_REAL kabs_nua;
  CCTK_REAL kabs_nux;
  CCTK_REAL kabs_anux;

  CCTK_REAL kscat_nue;
  CCTK_REAL kscat_nua;
  CCTK_REAL kscat_nux;
  CCTK_REAL kscat_anux;

  CCTK_REAL et_nue;
  CCTK_REAL et_nua;
  CCTK_REAL et_nux;
  CCTK_REAL et_anux;

  CCTK_HOST void init();

  CCTK_DEVICE inline M1Opacities ComputeFakeOpacities(const CCTK_REAL rho) {

    M1Opacities m1_opacities = {0};
    m1_opacities.eta_0[0] = rho * et_nue;
    m1_opacities.eta_0[1] = rho * et_nua;
    m1_opacities.eta_0[2] = rho * et_nux;
    m1_opacities.eta_0[3] = rho * et_anux;

    m1_opacities.kappa_0_a[0] = rho * kabs_nue;
    m1_opacities.kappa_0_a[1] = rho * kabs_nua;
    m1_opacities.kappa_0_a[2] = rho * kabs_nux;
    m1_opacities.kappa_0_a[3] = rho * kabs_anux;

    m1_opacities.eta[0] = rho * et_nue;
    m1_opacities.eta[1] = rho * et_nua;
    m1_opacities.eta[2] = rho * et_nux;
    m1_opacities.eta[3] = rho * et_anux;

    m1_opacities.kappa_a[0] = rho * kabs_nue;
    m1_opacities.kappa_a[1] = rho * kabs_nua;
    m1_opacities.kappa_a[2] = rho * kabs_nux;
    m1_opacities.kappa_a[3] = rho * kabs_anux;

    m1_opacities.kappa_s[0] = rho * kscat_nue;
    m1_opacities.kappa_s[1] = rho * kscat_nua;
    m1_opacities.kappa_s[2] = rho * kscat_nux;
    m1_opacities.kappa_s[3] = rho * kscat_anux;
    return m1_opacities;
  }

  CCTK_DEVICE CCTK_HOST inline 
  void FakeNeutrinoDens(CCTK_REAL rho, CCTK_REAL &num_nue, CCTK_REAL &num_nua,
                  CCTK_REAL &num_nux, CCTK_REAL &ene_nue, CCTK_REAL &ene_nua, CCTK_REAL &ene_nux) {

    if(rho*kabs_nue > FLT_EPSILON*et_nue) {
        num_nue = et_nue/(rho*kabs_nue);
        ene_nue = et_nue/(rho*kabs_nue);
    }
    else {
        num_nue = 1.0;
        ene_nue = 1.0;
    }

    if(rho*kabs_nua > FLT_EPSILON*et_nua) {
        num_nua = et_nua/(rho*kabs_nua);
        ene_nua = et_nua/(rho*kabs_nua);
    }
    else {
        num_nua = 1.0;
        ene_nua = 1.0;
    }

    num_nux = 1.0;
    ene_nux = 1.0;
  }
};

extern FakeRatesDef* global_fakerates;

} // namespace
