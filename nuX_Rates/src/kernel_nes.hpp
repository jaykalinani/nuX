//=================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
//! \file  kernel_nes.hpp
//  \brief contains kernels for inelastic neutrino scattering
//         on electrons and positrons
//
// Computation of inelastic neutrino-electron and neutrino-positron
// scattering using Eq. (43) from Mezzacappa & Bruenn, ApJ v.410, p.740 (1993)
// https://ui.adsabs.harvard.edu/abs/1993ApJ...410..740M/abstract

#ifndef BNS_NURATES_INCLUDE_KERNEL_NEPS_HPP_
#define BNS_NURATES_INCLUDE_KERNEL_NEPS_HPP_

#include "bns_nurates.hpp"
#include "functions.hpp"
#include "constants.hpp"

// Numerical constants
//---------------------------------------------------------------------------------------------------------------------
constexpr BS_REAL kTaylorSeriesEpsilon = 1e-3;


// Physical constats
//---------------------------------------------------------------------------------------------------------------------

// Saves the expression that are needed for the unapproximated integral
CCTK_DEVICE CCTK_HOST NUX_ATTRIBUTE_NOINLINE
void ComputeFDIForInelastic(BS_REAL w, BS_REAL wp, BS_REAL eta,
                            BS_REAL* fdi_diff_w, BS_REAL* fdi_diff_abs)
{
    BS_REAL abs_val = fabs(w - wp);

    fdi_diff_w[0] = FDI_p1(eta - wp) - FDI_p1(eta - w);
    fdi_diff_w[1] = FDI_p2(eta - wp) - FDI_p2(eta - w);
    fdi_diff_w[2] = FDI_p3(eta - wp) - FDI_p3(eta - w);
    fdi_diff_w[3] = FDI_p4(eta - wp) - FDI_p4(eta - w);
    fdi_diff_w[4] = FDI_p5(eta - wp) - FDI_p5(eta - w);

    fdi_diff_abs[0] = FDI_p3(eta) - FDI_p3(eta - abs_val);
    fdi_diff_abs[1] = FDI_p4(eta) - FDI_p4(eta - abs_val);
    fdi_diff_abs[2] = FDI_p5(eta) - FDI_p5(eta - abs_val);
}

// Functions that calculates the kernel integral
//=========================================================================================================================================

// Not approximated out kernel integral
CCTK_DEVICE CCTK_HOST NUX_ATTRIBUTE_NOINLINE
BS_REAL MezzacappaIntOut(BS_REAL w, BS_REAL wp, BS_REAL x, BS_REAL y, int sign,
                         BS_REAL b1, BS_REAL b2, BS_REAL* fdi_diff_w,
                         BS_REAL* fdi_diff_abs)
{
    constexpr BS_REAL one_fifth = 0.2;
    constexpr BS_REAL two       = 2;
    constexpr BS_REAL three     = 3;
    constexpr BS_REAL six       = 6;

    return (((b1 + b2) * (sign * fdi_diff_abs[2] - fdi_diff_w[4]) *
                 one_fifth // All G5 terms

             - b1 * (w + wp) *
                   (

                       fdi_diff_w[3]

                       + two * ((w + wp) * fdi_diff_w[2] +
                                three * w * wp * fdi_diff_w[1])

                           ) // G4(eta - wp) + G3(eta - wp) term

             + sign * ((b1 * x - b2 * y) * fdi_diff_abs[1] +
                       two * (b1 * x * x + b2 * y * y) * fdi_diff_abs[0])) /
                (w * w * wp * wp)

            - six * b1 * fdi_diff_w[0]);
}

// Taylor expansion in the lowest energy of the function MezzacappaIntOut and
// MezzacappaIntIn
CCTK_DEVICE CCTK_HOST NUX_ATTRIBUTE_NOINLINE
BS_REAL MezzacappaIntOneEnergy(BS_REAL x, BS_REAL y, int sign, BS_REAL b1,
                               BS_REAL b2, const BS_REAL* fdis)
{
    constexpr BS_REAL one_fifth = 0.2;
    constexpr BS_REAL two       = 2;
    constexpr BS_REAL three     = 3;
    constexpr BS_REAL four      = 4;
    constexpr BS_REAL six       = 6;

    return -sign * y *
           (two * (b1 + b2) * fdis[0]

            + (b1 * (y + four * x) + three * b2 * y) * fdis[1]

            + (b1 * (three * y - four * x) + b2 * y) * fdis[2]

            + ((b1 + six * b2) * y * y * one_fifth + b1 * x * y +
               two * b1 * x * x) *
                  fdis[3]

            - ((six * b1 + b2) * y * y * one_fifth - three * b1 * x * y +
               two * b1 * x * x) *
                  fdis[4]) /
           (x * x);
}

// Taylor expansion in both energies of the function MezzacappaIntOut and
// MezzacappaIntIn
CCTK_DEVICE CCTK_HOST NUX_ATTRIBUTE_NOINLINE
BS_REAL MezzacappaIntTwoEnergies(BS_REAL w, BS_REAL wp, BS_REAL x, BS_REAL y,
                                 BS_REAL b1, BS_REAL b2, const BS_REAL* fdis)
{
    constexpr BS_REAL two        = 2;
    constexpr BS_REAL four       = 4;
    constexpr BS_REAL eleven     = 11;
    constexpr BS_REAL twenty     = 20;
    constexpr BS_REAL twentyfive = 25;
    constexpr BS_REAL thirthy    = 30;

    return y * (w - wp) *
           ((b1 + b2) *
                (four * fdis[2]

                 + fdis[0] *
                       (eleven * y * y - twentyfive * x * y + twenty * x * x) /
                       thirthy)

            + (b1 - b2) * (two * x - y) * fdis[1]) /
           (x * x);
}

// Functions for kernel calculation
//===================================================================================================

// Calculates and saves the neutrino electron scattering in and out kernel for
// every neutrino species
CCTK_DEVICE CCTK_HOST NUX_ATTRIBUTE_NOINLINE
MyKernelOutput NESKernels(InelasticScattKernelParams* kernel_params,
                          MyEOSParams* eos_params)
{
    const BS_REAL T                    = eos_params->temp;
    const BS_REAL w                    = kernel_params->omega / T;
    const BS_REAL wp                   = kernel_params->omega_prime / T;
    const BS_REAL x                    = fmax(w, wp);
    const BS_REAL y                    = fmin(w, wp);
    const BS_REAL eta_e                = eos_params->mu_e / T;
    const BS_REAL exp_factor           = NEPSExpFunc(wp - w);
    const BS_REAL exp_factor_exchanged = NEPSExpFunc(w - wp);

    MyKernelOutput output;

    constexpr BS_REAL zero = 0;
    constexpr BS_REAL one  = 1;
    constexpr BS_REAL six  = 6;

    if (y > eta_e * kTaylorSeriesEpsilon)
    {
        BS_REAL fdi_diff_abs[3], fdi_diff_w[5];

        const int sign = 2 * signbit(wp - w) - 1;

        ComputeFDIForInelastic(w, wp, eta_e, fdi_diff_w, fdi_diff_abs);

        output.abs[id_nue] =
            kBS_NEPS_Const * POW2(T) *
            MezzacappaIntOut(w, wp, x, y, sign, kBS_NEPS_BPlus, kBS_NEPS_BZero,
                             fdi_diff_w, fdi_diff_abs);
        output.abs[id_anue] =
            kBS_NEPS_Const * POW2(T) *
            MezzacappaIntOut(w, wp, x, y, sign, kBS_NEPS_BZero, kBS_NEPS_BPlus,
                             fdi_diff_w, fdi_diff_abs);
        output.abs[id_nux] =
            kBS_NEPS_Const * POW2(T) *
            MezzacappaIntOut(w, wp, x, y, sign, kBS_NEPS_BMinus, kBS_NEPS_BZero,
                             fdi_diff_w, fdi_diff_abs);
        output.abs[id_anux] =
            kBS_NEPS_Const * POW2(T) *
            MezzacappaIntOut(w, wp, x, y, sign, kBS_NEPS_BZero, kBS_NEPS_BMinus,
                             fdi_diff_w, fdi_diff_abs);
    }
    else if (x > eta_e * kTaylorSeriesEpsilon)
    {
        const int sign = 2 * signbit(wp - w) - 1;

        const BS_REAL fdis[5] = {
            FDI_p2(eta_e - x) - FDI_p2(eta_e), FDI_p1(eta_e - x), FDI_p1(eta_e),
            FDI_0(eta_e - x) + y * FermiDistr(zero, one, eta_e - x) / six,
            FDI_0(eta_e) - y * FermiDistr(zero, one, eta_e) / six};

        output.abs[id_nue] = kBS_NEPS_Const * POW2(T) *
                             MezzacappaIntOneEnergy(x, y, sign, kBS_NEPS_BPlus,
                                                    kBS_NEPS_BZero, fdis);
        output.abs[id_anue] = kBS_NEPS_Const * POW2(T) *
                              MezzacappaIntOneEnergy(x, y, sign, kBS_NEPS_BZero,
                                                     kBS_NEPS_BPlus, fdis);
        output.abs[id_nux] = kBS_NEPS_Const * POW2(T) *
                             MezzacappaIntOneEnergy(x, y, sign, kBS_NEPS_BMinus,
                                                    kBS_NEPS_BZero, fdis);
        output.abs[id_anux] = kBS_NEPS_Const * POW2(T) *
                              MezzacappaIntOneEnergy(x, y, sign, kBS_NEPS_BZero,
                                                     kBS_NEPS_BMinus, fdis);
    }
    else
    {
        const BS_REAL fdis[3] = {FermiDistr(0., 1., eta_e), FDI_0(eta_e),
                                 FDI_p1(eta_e)};

        output.abs[id_nue] =
            kBS_NEPS_Const * POW2(T) *
            MezzacappaIntTwoEnergies(w, wp, x, y, kBS_NEPS_BPlus,
                                     kBS_NEPS_BZero, fdis);
        output.abs[id_anue] =
            kBS_NEPS_Const * POW2(T) *
            MezzacappaIntTwoEnergies(w, wp, x, y, kBS_NEPS_BZero,
                                     kBS_NEPS_BPlus, fdis);
        output.abs[id_nux] =
            kBS_NEPS_Const * POW2(T) *
            MezzacappaIntTwoEnergies(w, wp, x, y, kBS_NEPS_BMinus,
                                     kBS_NEPS_BZero, fdis);
        output.abs[id_anux] =
            kBS_NEPS_Const * POW2(T) *
            MezzacappaIntTwoEnergies(w, wp, x, y, kBS_NEPS_BZero,
                                     kBS_NEPS_BMinus, fdis);
    }

    for (int idx = 0; idx < total_num_species; ++idx)
    {
        output.em[idx]  = output.abs[idx];
        output.abs[idx] = output.abs[idx] * exp_factor;
        output.em[idx]  = -output.em[idx] * exp_factor_exchanged;
    }

    return output;
}

// Calculates and saves the neutrino positron scattering in and out kernel for
// every neutrino species
CCTK_DEVICE CCTK_HOST NUX_ATTRIBUTE_NOINLINE
MyKernelOutput NPSKernels(InelasticScattKernelParams* kernel_params,
                          MyEOSParams* eos_params)
{
    const BS_REAL T                    = eos_params->temp;
    const BS_REAL w                    = kernel_params->omega / T;
    const BS_REAL wp                   = kernel_params->omega_prime / T;
    const BS_REAL x                    = fmax(w, wp);
    const BS_REAL y                    = fmin(w, wp);
    const BS_REAL eta_p                = -eos_params->mu_e / T;
    const BS_REAL exp_factor           = NEPSExpFunc(wp - w);
    const BS_REAL exp_factor_exchanged = NEPSExpFunc(w - wp);

    MyKernelOutput output;

    constexpr BS_REAL zero = 0;
    constexpr BS_REAL one  = 1;
    constexpr BS_REAL six  = 6;

    if (y > eta_p * kTaylorSeriesEpsilon)
    {
        BS_REAL fdi_diff_abs[3], fdi_diff_w[5];

        const int sign = 2 * signbit(wp - w) - 1;

        ComputeFDIForInelastic(w, wp, eta_p, fdi_diff_w, fdi_diff_abs);

        output.abs[id_nue] =
            kBS_NEPS_Const * POW2(T) *
            MezzacappaIntOut(w, wp, x, y, sign, kBS_NEPS_BZero, kBS_NEPS_BPlus,
                             fdi_diff_w, fdi_diff_abs);
        output.abs[id_anue] =
            kBS_NEPS_Const * POW2(T) *
            MezzacappaIntOut(w, wp, x, y, sign, kBS_NEPS_BPlus, kBS_NEPS_BZero,
                             fdi_diff_w, fdi_diff_abs);
        output.abs[id_nux] =
            kBS_NEPS_Const * POW2(T) *
            MezzacappaIntOut(w, wp, x, y, sign, kBS_NEPS_BZero, kBS_NEPS_BMinus,
                             fdi_diff_w, fdi_diff_abs);
        output.abs[id_anux] =
            kBS_NEPS_Const * POW2(T) *
            MezzacappaIntOut(w, wp, x, y, sign, kBS_NEPS_BMinus, kBS_NEPS_BZero,
                             fdi_diff_w, fdi_diff_abs);
    }
    else if (x > eta_p * kTaylorSeriesEpsilon)
    {
        const int sign = 2 * signbit(wp - w) - 1;

        const BS_REAL fdis[5] = {
            FDI_p2(eta_p - x) - FDI_p2(eta_p), FDI_p1(eta_p - x), FDI_p1(eta_p),
            FDI_0(eta_p - x) + y * FermiDistr(zero, one, eta_p - x) / six,
            FDI_0(eta_p) - y * FermiDistr(zero, one, eta_p) / six};

        output.abs[id_nue] = kBS_NEPS_Const * POW2(T) *
                             MezzacappaIntOneEnergy(x, y, sign, kBS_NEPS_BZero,
                                                    kBS_NEPS_BPlus, fdis);
        output.abs[id_anue] = kBS_NEPS_Const * POW2(T) *
                              MezzacappaIntOneEnergy(x, y, sign, kBS_NEPS_BPlus,
                                                     kBS_NEPS_BZero, fdis);
        output.abs[id_nux] = kBS_NEPS_Const * POW2(T) *
                             MezzacappaIntOneEnergy(x, y, sign, kBS_NEPS_BZero,
                                                    kBS_NEPS_BMinus, fdis);
        output.abs[id_anux] =
            kBS_NEPS_Const * POW2(T) *
            MezzacappaIntOneEnergy(x, y, sign, kBS_NEPS_BMinus, kBS_NEPS_BZero,
                                   fdis);
    }
    else
    {
        const BS_REAL fdis[3] = {FermiDistr(0., 1., eta_p), FDI_0(eta_p),
                                 FDI_p1(eta_p)};

        output.abs[id_nue] =
            kBS_NEPS_Const * POW2(T) *
            MezzacappaIntTwoEnergies(w, wp, x, y, kBS_NEPS_BZero,
                                     kBS_NEPS_BPlus, fdis);
        output.abs[id_anue] =
            kBS_NEPS_Const * POW2(T) *
            MezzacappaIntTwoEnergies(w, wp, x, y, kBS_NEPS_BPlus,
                                     kBS_NEPS_BZero, fdis);
        output.abs[id_nux] =
            kBS_NEPS_Const * POW2(T) *
            MezzacappaIntTwoEnergies(w, wp, x, y, kBS_NEPS_BZero,
                                     kBS_NEPS_BMinus, fdis);
        output.abs[id_anux] =
            kBS_NEPS_Const * POW2(T) *
            MezzacappaIntTwoEnergies(w, wp, x, y, kBS_NEPS_BMinus,
                                     kBS_NEPS_BZero, fdis);
    }

    for (int idx = 0; idx < total_num_species; ++idx)
    {
        output.em[idx]  = output.abs[idx];
        output.abs[idx] = output.abs[idx] * exp_factor;
        output.em[idx]  = -output.em[idx] * exp_factor_exchanged;
    }

    return output;
}

// Calculates the full in and out kernels
CCTK_DEVICE CCTK_HOST NUX_ATTRIBUTE_NOINLINE
MyKernelOutput InelasticScattKernels(InelasticScattKernelParams* kernel_params,
                                     MyEOSParams* eos_params)
{
    MyKernelOutput nes_kernel = NESKernels(kernel_params, eos_params);
    MyKernelOutput nps_kernel = NPSKernels(kernel_params, eos_params);

    MyKernelOutput tot_kernel = {0};

    for (int idx = 0; idx < total_num_species; ++idx)
    {
        tot_kernel.em[idx]  = nes_kernel.em[idx] + nps_kernel.em[idx];
        tot_kernel.abs[idx] = nes_kernel.abs[idx] + nps_kernel.abs[idx];
    }

    return tot_kernel;
}

CCTK_DEVICE CCTK_HOST NUX_ATTRIBUTE_NOINLINE
void InelasticKernelsTable(const int n, BS_REAL* nu_array,
                           GreyOpacityParams* grey_pars, M1MatrixKokkos2D* out)
{
    MyKernelOutput inel_1, inel_2;

    InelasticScattKernelParams inelastic_pars =
        grey_pars->kernel_pars.inelastic_kernel_params;
    for (int i = 0; i < n; ++i)
    {

        for (int j = i; j < n; ++j)
        {

            // compute the pair kernels
            inelastic_pars.omega       = nu_array[i];
            inelastic_pars.omega_prime = nu_array[j];
            inel_1 =
                InelasticScattKernels(&inelastic_pars, &grey_pars->eos_pars);

            inelastic_pars.omega       = nu_array[j];
            inelastic_pars.omega_prime = nu_array[i];
            inel_2 =
                InelasticScattKernels(&inelastic_pars, &grey_pars->eos_pars);


            for (int idx = 0; idx < total_num_species; ++idx)
            {
                out->m1_mat_em[idx][i][j] = inel_1.em[idx];
                out->m1_mat_em[idx][j][i] = inel_2.em[idx];

                out->m1_mat_ab[idx][i][j] = inel_1.abs[idx];
                out->m1_mat_ab[idx][j][i] = inel_2.abs[idx];
            }
        }
    }

    return;
}

#endif // BNS_NURATES_INCLUDE_KERNEL_NEPS_HPP_
