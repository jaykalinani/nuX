#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

namespace nuX_FakeRates {

CCTK_DEVICE CCTK_HOST NUX_ATTRIBUTE_NOINLINE
M1Opacities ComputeFakeOpacities(const MyQuadrature* quad_1d,
                               const MyQuadrature* quad_2d,
                               GreyOpacityParams* my_grey_opacity_params)
{
    DECLARE_CCTK_PARAMETERS;
    M1Opacities m1_opacities = {0};

    const BS_REAL rhoL = my_grey_opacity_params->eos_pars.nb; // this should just be plain fluid mass density in CU
    // Electron neutrinos
    m1_opacities.kappa_0_a[0] = rhoL * kappa_abs_nue;
    m1_opacities.kappa_a[0]   = rhoL * kappa_abs_nue;
    m1_opacities.kappa_s[0]   = rhoL * kappa_scat_nue;
    m1_opacities.eta_0[0]     = rhoL * eta_nue;
    m1_opacities.eta[0]       = rhoL * eta_nue;

    // Anti-electron neutrinos
    m1_opacities.kappa_0_a[1] = rhoL * kappa_abs_nua;
    m1_opacities.kappa_a[1]   = rhoL * kappa_abs_nua;
    m1_opacities.kappa_s[1]   = rhoL * kappa_scat_nua;
    m1_opacities.eta_0[1]     = rhoL * eta_nua;
    m1_opacities.eta[1]       = rhoL * eta_nua;

    // Heavy neutrinos
    m1_opacities.kappa_0_a[2] = rhoL * kappa_abs_nux;
    m1_opacities.kappa_a[2]   = rhoL * kappa_abs_nux;
    m1_opacities.kappa_s[2]   = rhoL * kappa_scat_nux;
    m1_opacities.eta_0[2]     = rhoL * eta_nux;
    m1_opacities.eta[2]       = rhoL * eta_nux;

    // Heavy neutrinos
    m1_opacities.kappa_0_a[3] = rhoL * kappa_abs_anux;
    m1_opacities.kappa_a[3]   = rhoL * kappa_abs_anux;
    m1_opacities.kappa_s[3]   = rhoL * kappa_scat_anux;
    m1_opacities.eta_0[3]     = rhoL * eta_anux;
    m1_opacities.eta[3]       = rhoL * eta_anux;

    return m1_opacities;
}

} // namespace
