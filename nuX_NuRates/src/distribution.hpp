// ================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
//! \file distribution.hpp
//  \brief header file for distribution function reconstruction from M1
//  parameters
//         supports three different species: electron neutrino, electron
//         anti-neutrino and mu/tau neutrino/antineutrino

#ifndef BNS_NURATES_SRC_DISTRIBUTION_DISTRIBUTION_HPP_
#define BNS_NURATES_SRC_DISTRIBUTION_DISTRIBUTION_HPP_

#include "bns_nurates.hpp"
#include "constants.hpp"
#include "functions.hpp"

#define CONST_C_F 0.6

/* ===========================================================================
 * Functions for the optically thick distribution function
 * ===========================================================================
 */

/* Neutrino distribution function in optically thick regime: Fermi-Dirac
distribution
 *
 * omega:       neutrino energy
 * distr_pars:  uses temp_t (fluid temperature) and eta_t (degeneracy parameter)
 * species:     species of neutrino
 */
CCTK_HOST CCTK_DEVICE inline
BS_REAL NuFThick(const BS_REAL omega, const NuDistributionParams* distr_pars,
                 const int nuid)
{
    const BS_REAL arg = distr_pars->eta_t[nuid] -
                        omega / distr_pars->temp_t[nuid];
    // Both exponential arguments are nonpositive: no overflow guard is
    // needed. Allow true underflow instead of imposing SafeExp's occupation
    // floor, which would renormalize an extremely dilute distribution.
    const BS_REAL e = exp(-fabs(arg));
    return arg <= 0 ? e / (1 + e) : 1 / (1 + e);
}

// Check normalization, not just finiteness, against both original moments.
CCTK_HOST CCTK_DEVICE inline
bool ThickDistributionMatchesM1(const BS_REAL n, const BS_REAL J,
                                const BS_REAL T, const BS_REAL eta)
{
    if (!(n > 0 && J > 0 && T > 0) || !isfinite(n) || !isfinite(J) ||
        !isfinite(T) || !isfinite(eta))
        return false;

    // Logarithmic moments avoid underflow for very dilute populations.
    // Below eta = -20 the leading dilute terms are accurate to < 3e-10.
    const BS_REAL log_F2 = eta < -20 ? log(BS_REAL(2)) + eta
                                    : log(FDI_p2(eta));
    const BS_REAL log_F3 = eta < -20 ? log(BS_REAL(6)) + eta
                                    : log(FDI_p3(eta));
    const BS_REAL log_T = log(T);
    const BS_REAL log_n = log(kBS_FourPi_hc3) + 3 * log_T + log_F2;
    const BS_REAL log_J = log(kBS_FourPi_hc3) + 4 * log_T + log_F3;

    // Allow the existing inverse fit's measured 4.25e-4 relative error.
    // This tolerance does not clip or rescale either radiation moment.
    constexpr BS_REAL tolerance = 5e-4;
    return isfinite(log_n) && isfinite(log_J) &&
           fabs(log_n - log(n)) <= tolerance &&
           fabs(log_J - log(J)) <= tolerance;
}

// Positive-eta reflection identities, scaled by eta^3 and eta^4. In the
// high-degeneracy inverse use s=eta^-2: its Pauli deficit has a finite,
// nonzero derivative as s -> 0, unlike the ratio as a function of eta.
CCTK_HOST CCTK_DEVICE inline
BS_REAL DegenerateFDLogGap(const BS_REAL s, BS_REAL* derivative,
                           BS_REAL* log_ratio)
{
    constexpr BS_REAL pi2 = kBS_Pi * kBS_Pi;
    const BS_REAL root_s = sqrt(s);
    const BS_REAL eta = 1 / root_s;
    const BS_REAL fm1 = FDI_p1(-eta);
    const BS_REAL fm2 = FDI_p2(-eta);
    const BS_REAL fm3 = FDI_p3(-eta);
    const BS_REAL a = pi2 * s + 3 * fm2 * s * root_s;
    const BS_REAL b = 2 * pi2 * s + 7 * pi2 * pi2 * s * s / 15 -
                      4 * fm3 * s * s;
    const BS_REAL da = pi2 + 3 * fm1 + BS_REAL(4.5) * fm2 * root_s;
    const BS_REAL db = 2 * pi2 + 14 * pi2 * pi2 * s / 15 -
                       6 * fm2 * root_s - 8 * fm3 * s;
    const BS_REAL la = log1p(a), lb = log1p(b);
    *derivative = 3 * db / (1 + b) - 4 * da / (1 + a);
    *log_ratio = log(BS_REAL(4) / 3) + BS_REAL(0.5) * log(s) + la - lb;
    return 3 * lb - 4 * la;
}

CCTK_HOST CCTK_DEVICE inline
bool CalculateDegenerateThickParams(const BS_REAL log_n, const BS_REAL log_J,
                                    BS_REAL* temperature, BS_REAL* eta)
{
    const BS_REAL log_c = log(kBS_FourPi_hc3);
    const BS_REAL log_y_max = log(BS_REAL(64) / 81);
    const BS_REAL target = log_y_max - (4 * log_n - 3 * log_J - log_c);
    constexpr BS_REAL eps = std::numeric_limits<BS_REAL>::epsilon();
    // At the zero-T boundary these rounded moments cannot identify a finite
    // spectral width. Do not clip them to an invented maximum degeneracy.
    const BS_REAL uncertainty = 8 * eps *
        (4 * fabs(log_n) + 3 * fabs(log_J) + fabs(log_c) + 1);
    if (!isfinite(target) || target <= uncertainty)
        return false;

    // Y(eta=10) < 0.7, the existing high-fit matching point. Thus this
    // interval brackets every supported high-degeneracy root, without an
    // upper bound on eta. The lower endpoint's analytic gap is zero.
    BS_REAL lower = 0, upper = BS_REAL(0.01);
    BS_REAL derivative, log_ratio;
    const BS_REAL upper_gap = DegenerateFDLogGap(upper, &derivative, &log_ratio);
    if (!(target < upper_gap))
        return false;
    BS_REAL s = target / (2 * kBS_Pi * kBS_Pi);
    for (int iteration = 0; iteration < 64; ++iteration)
    {
        const BS_REAL gap = DegenerateFDLogGap(s, &derivative, &log_ratio);
        const BS_REAL residual = gap - target;
        if (!isfinite(gap) || !isfinite(derivative) || derivative <= 0)
            return false;
        const BS_REAL correction = residual / derivative;
        if (fabs(residual) <= 64 * eps * target &&
            fabs(correction) <= 64 * eps * s)
        {
            *eta = 1 / sqrt(s);
            *temperature = exp(log_J - log_n + log_ratio);
            return isfinite(*eta) && isfinite(*temperature) && *temperature > 0;
        }
        if (residual > 0)
            upper = s;
        else
            lower = s;
        const BS_REAL candidate = s - correction;
        const BS_REAL next = candidate > lower && candidate < upper
                                 ? candidate : BS_REAL(0.5) * (lower + upper);
        if (next == s)
            return false;
        s = next;
    }
    return false;
}

/* Recover distribution function parameters for optically thick regime from M1
quantities
 *
 * M1_params:       uses n (neutrino number density) and J (neutrino energy
density)
 * eos_pars:        uses fluid temperature
 * out_distr_pars:  computes trapped neutrino temperature and degeneracy
parameter
 */
CCTK_HOST CCTK_DEVICE inline
void CalculateThickParamsFromM1(const M1Quantities* M1_pars,
                                NuDistributionParams* out_distribution_pars)
{
    // constexpr BS_REAL zero         = 0;
    constexpr BS_REAL one          = 1;
    constexpr BS_REAL three_halves = 1.5;
    [[maybe_unused]] constexpr BS_REAL one_third    = 1. / 3.;

    constexpr BS_REAL thirty = 30;

    constexpr BS_REAL y1 = 0.005;
    constexpr BS_REAL y2 = 0.7;

    constexpr BS_REAL kBS_ThickDistr_num_1[6] = {
        0.19926987701997, 38865.0149220478, 14364.6331099737,
        5750.1878570758,  1120.71610972194, 1.60356108438235e-8};
    constexpr BS_REAL kBS_ThickDistr_num_2[7] = {
        41.3836568203438,    32.5515666786612, 157.774993512235,
        66.5726772426253,    14.4883415579211, 0.315360380575709,
        0.000660414331285249};

    constexpr BS_REAL kBS_ThickDistr_den_1[6] = {1.0,
                                                 38840.0743174942,
                                                 99.9009680656931,
                                                 171.874844843596,
                                                 75.7101579899442,
                                                 83.0160130941424};
    constexpr BS_REAL kBS_ThickDistr_den_2[6] = {
        1.8888797407042,   5.35488690539183,   1.94673781342617,
        0.483128792035557, 0.0113386564109086, 2.64160073447322e-5};

    for (int nuid = 0; nuid < total_num_species; ++nuid)
    {
        // Preserve every positive input, including very dilute populations.
        // Exact vacuum is assigned zero weights by CalculateDistrParamsFromM1.
        const BS_REAL n   = (M1_pars->n[nuid] > 0) ? M1_pars->n[nuid] : 1e-100; // [nm^-3]
        const BS_REAL J   = (M1_pars->J[nuid] > 0) ? M1_pars->J[nuid] : 1e-100; // [MeV nm^-3]
        const BS_REAL chi = M1_pars->chi[nuid];

        // BS_ASSERT(
        // n > zero,
        // "Neutrino (species=%d) number density is non-positive (n[%d]=%e).",
        // nuid, nuid, n);
        BS_ASSERT(chi >= one_third && chi <= one,
                  "Invalid Eddington factor (chi[%d]=%e).", nuid, chi);

        out_distribution_pars->w_t[nuid] = three_halves * (one - chi);
        if (out_distribution_pars->w_t[nuid] == 0)
        {
            // UNUSED component: a valid thin spectrum need not admit any FD
            // inverse. Keep finite placeholders because TotalNuF evaluates
            // both components before applying their weights.
            out_distribution_pars->eta_t[nuid] = 0;
            out_distribution_pars->temp_t[nuid] = 1;
            continue;
        }

        // y = F_2(eta)^4 / F_3(eta)^3. The rational fit below has a
        // nonzero constant term and must not be extrapolated to y -> 0.
        // Work in logs so dilute moments do not underflow before inversion.
        const BS_REAL log_n = log(n), log_J = log(J);
        const BS_REAL log_y = 4 * log_n - 3 * log_J - log(kBS_FourPi_hc3);
        if (log_y > log(y2))
        {
            if (!CalculateDegenerateThickParams(
                    log_n, log_J, &out_distribution_pars->temp_t[nuid],
                    &out_distribution_pars->eta_t[nuid]))
            {
                // Existing callers reject invalid trials through their finite
                // checks; never replace an unresolved/forbidden FD state.
                out_distribution_pars->temp_t[nuid] =
                    std::numeric_limits<BS_REAL>::quiet_NaN();
                out_distribution_pars->eta_t[nuid] =
                    std::numeric_limits<BS_REAL>::quiet_NaN();
            }
            continue;
        }
        if (log_y < log(BS_REAL(1e-4)))
        {
            // For fugacity z << 1, y = (2/27) z (1 - 5 z/16 + ...).
            // The first correction keeps the moment error below 4e-8 at
            // the matching point. Avoid ratios of underflowed Fermi integrals.
            const BS_REAL log_z = log(BS_REAL(13.5)) + log_y;
            const BS_REAL z = exp(log_z);
            out_distribution_pars->eta_t[nuid] = log_z + 5 * z / 16;
            out_distribution_pars->temp_t[nuid] =
                (J / n) * (1 - z / 16) / 3;
            continue;
        }
        const BS_REAL y = exp(log_y);

        if (y < y1)
        {
            out_distribution_pars->eta_t[nuid] =
                log((y * (y * (y * (y * (y * (y + kBS_ThickDistr_num_1[0]) +
                                         kBS_ThickDistr_num_1[1]) +
                                    kBS_ThickDistr_num_1[2]) +
                               kBS_ThickDistr_num_1[3]) +
                          kBS_ThickDistr_num_1[4]) +
                     kBS_ThickDistr_num_1[5]) /
                    (y * (y * (y * (y * (y * (y + kBS_ThickDistr_den_1[0]) +
                                         kBS_ThickDistr_den_1[1]) -
                                    kBS_ThickDistr_den_1[2]) -
                               kBS_ThickDistr_den_1[3]) +
                          kBS_ThickDistr_den_1[4]) +
                     kBS_ThickDistr_den_1[5]));
        }
        else
        {
            out_distribution_pars->eta_t[nuid] =
                (y * (y * (y * (y * (y * (kBS_ThickDistr_num_2[0] * y +
                                          kBS_ThickDistr_num_2[1]) -
                                     kBS_ThickDistr_num_2[2]) +
                                kBS_ThickDistr_num_2[3]) +
                           kBS_ThickDistr_num_2[4]) +
                      kBS_ThickDistr_num_2[5]) +
                 kBS_ThickDistr_num_2[6]) /
                    (y * (y * (y * (y * (y * (y + kBS_ThickDistr_den_2[0]) -
                                         kBS_ThickDistr_den_2[1]) +
                                    kBS_ThickDistr_den_2[2]) +
                               kBS_ThickDistr_den_2[3]) +
                          kBS_ThickDistr_den_2[4]) +
                     kBS_ThickDistr_den_2[5]) -
                thirty;
        }
        out_distribution_pars->temp_t[nuid] =
            FDI_p2(out_distribution_pars->eta_t[nuid]) * J /
            (FDI_p3(out_distribution_pars->eta_t[nuid]) * n);
    }
}

/* ===========================================================================
 * Functions for the optically thin distribution function
 * ===========================================================================
 */

/* Neutrino distribution function for optically thin regime: Maxwell-Boltzmann
 * distribution
 *
 * omega:       neutrino energy
 * distr_pars:  optically thin parameters
 * species:     species of neutrino
 */
CCTK_HOST CCTK_DEVICE inline
BS_REAL NuFThin(const BS_REAL omega, const NuDistributionParams* distr_pars,
                const int nuid)
{
    const BS_REAL T_f    = distr_pars->temp_f[nuid];
    const BS_REAL c_f    = distr_pars->c_f[nuid];
    const BS_REAL beta_f = distr_pars->beta_f[nuid];

    return beta_f * pow(omega, c_f) * exp(-omega / T_f);
}

/* Recover distribution function parameters for optically thin regime from M1
 * quantities
 *
 * M1_params:       uses n (neutrino number density) and J (neutrino energy
 * density) out_distr_pars:  computes free neutrino temperature and c_f
 */
CCTK_HOST CCTK_DEVICE inline
void CalculateThinParamsFromM1(const M1Quantities* M1_pars,
                               NuDistributionParams* out_distribution_pars)
{
    [[maybe_unused]] constexpr BS_REAL zero      = 0;
    constexpr BS_REAL one       = 1;
    constexpr BS_REAL three     = 3;
    constexpr BS_REAL half      = 0.5;
    [[maybe_unused]] constexpr BS_REAL one_third = 1. / 3.;

    constexpr BS_REAL c_f = CONST_C_F;
    // Integral of x^(c_f+2) exp(-x). Leading-order Stirling at 3.6
    // misnormalizes both reconstructed moments by about 2.34 percent.
    // c_f is fixed here; avoid a runtime special function on the GPU.
    static_assert(c_f == BS_REAL(0.6));
    constexpr BS_REAL gamma_c_f_plus_three =
        3.7170238530367914875762561372500703;

    for (int nuid = 0; nuid < total_num_species; ++nuid)
    {
        const BS_REAL n   = (M1_pars->n[nuid] > 0) ? M1_pars->n[nuid] : 1e-100; // [nm^-3]
        const BS_REAL J   = (M1_pars->J[nuid] > 0) ? M1_pars->J[nuid] : 1e-100; // [MeV nm^-3]
        const BS_REAL chi = M1_pars->chi[nuid];

        // BS_ASSERT(
        // n > zero,
        // "Neutrino (species=%d) number density is non-positive (n[%d]=%e).",
        // nuid, nuid, n);
        BS_ASSERT(chi >= one_third && chi <= one,
                  "Invalid Eddington factor (chi[%d]=%e).", nuid, chi);

        out_distribution_pars->w_f[nuid] = half * (three * chi - one);

        out_distribution_pars->c_f[nuid] = c_f;

        const BS_REAL Tnu                   = J / (n * (c_f + three));
        out_distribution_pars->temp_f[nuid] = Tnu;

        out_distribution_pars->beta_f[nuid] =
            n / (kBS_FourPi_hc3 * gamma_c_f_plus_three *
                 pow(Tnu, c_f + three));
    }
}

/* ===========================================================================
 * Functions for constructing the total distribution function
 * ===========================================================================
 */

/* Function for evaluating parameters of neutrino distribution function at
 * equilibrium
 */
CCTK_HOST CCTK_DEVICE inline
NuDistributionParams NuEquilibriumParams(const MyEOSParams* eos_pars)
{
    constexpr BS_REAL zero     = 0;
    constexpr BS_REAL one      = 1;
    [[maybe_unused]] constexpr BS_REAL thousand = 1e3;

    constexpr BS_REAL c_f = CONST_C_F;

    NuDistributionParams out;

    const BS_REAL T    = eos_pars->temp; // [MeV]
    const BS_REAL mu_e = eos_pars->mu_e; // [MeV]
    const BS_REAL mu_p = eos_pars->mu_p; // [MeV]
    const BS_REAL mu_n = eos_pars->mu_n; // [MeV]

    BS_ASSERT(T > zero && T < thousand,
              "Given temperature is either negative or beyond 1 GeV (T=%e).",
              T);

    for (int idx = 0; idx < total_num_species; ++idx)
    {
        out.w_t[idx]    = one;
        out.temp_t[idx] = T;

        out.w_f[idx]    = zero;
        out.c_f[idx]    = c_f;
        out.temp_f[idx] = T;
        out.beta_f[idx] = one;
    }

    out.eta_t[id_nue]  = (mu_e - mu_n + mu_p) / T;
    out.eta_t[id_anue] = -out.eta_t[id_nue];
    out.eta_t[id_nux]  = zero;
    out.eta_t[id_anux] = zero;

    return out;
}

/* Total neutrino distribution combining optically thick and thin regimes
 *
 * omega:       neutrino energy
 * distr_pars:  neutrino distribution parameters for thick and thin regimes
 * species:     species of neutrino
 */
CCTK_HOST CCTK_DEVICE inline
BS_REAL TotalNuF(const BS_REAL omega, const NuDistributionParams* distr_pars,
                 const int nuid)
{
    [[maybe_unused]] constexpr BS_REAL zero = 0;

    BS_ASSERT(omega >= zero, "Neutrino energy is negative.");
    BS_ASSERT(nuid >= 0 && nuid < total_num_species,
              "Invalid neutrino species ID.");

    const BS_REAL w_t = distr_pars->w_t[nuid];
    const BS_REAL w_f = distr_pars->w_f[nuid];

    if (w_t == 0 && w_f == 0)
        return 0;
    if (w_f == 0)
        return w_t * NuFThick(omega, distr_pars, nuid);
    if (w_t == 0)
        return w_f * NuFThin(omega, distr_pars, nuid);

    const BS_REAL f_thick = NuFThick(omega, distr_pars, nuid);
    const BS_REAL f_thin  = NuFThin(omega, distr_pars, nuid);

    return w_t * f_thick + w_f * f_thin;
}

// Preserve a representable hole probability even when the FD occupation
// rounds to one. The legacy mixture retains its original complement;
// this helper does not clip an inadmissible mixed spectrum.
CCTK_HOST CCTK_DEVICE inline
BS_REAL OneMinusTotalNuF(const BS_REAL omega,
                         const NuDistributionParams* parameters, const int sp)
{
    if (parameters->w_t[sp] == 1 && parameters->w_f[sp] == 0)
    {
        const BS_REAL arg = parameters->eta_t[sp] -
                            omega / parameters->temp_t[sp];
        const BS_REAL e = exp(-fabs(arg));
        return arg >= 0 ? e / (1 + e) : 1 / (1 + e);
    }
    return 1 - TotalNuF(omega, parameters, sp);
}

/* Calculate distribution function parameters in the thick and thin regime from
 * M1 quantities
 *
 * M1_params:   M1 quantities
 * eos_params:  parameters from EOS
 */
CCTK_HOST CCTK_DEVICE inline
NuDistributionParams CalculateDistrParamsFromM1(const M1Quantities* M1_pars,
                                                const MyEOSParams* eos_pars)
{
    NuDistributionParams out;

    CalculateThickParamsFromM1(M1_pars, &out);
    CalculateThinParamsFromM1(M1_pars, &out);

    // Exact vacuum has zero occupation, not the fictitious population implied
    // by the positive denominator guards in the reconstruction routines.
    for (int nuid = 0; nuid < total_num_species; ++nuid)
    {
        if (M1_pars->n[nuid] == 0 && M1_pars->J[nuid] == 0)
        {
            out.w_t[nuid] = out.w_f[nuid] = 0;
            out.temp_t[nuid] = out.temp_f[nuid] = eos_pars->temp;
            out.eta_t[nuid] = out.beta_f[nuid] = 0;
        }
    }

    return out;
}

// Scalar collision spectrum, independent of the TRANSPORT pressure closure.
// A positive FD/modified-MB mixture need not satisfy f<=1 in degenerate
// matter. Match the actual n,J with one FD spectrum instead; do not change
// either moment or replace a floor-sized state by a matter-temperature bath.
CCTK_HOST CCTK_DEVICE inline
bool CalculateFermiDistrParamsFromM1(const M1Quantities* moments,
                                    const MyEOSParams* matter,
                                    NuDistributionParams* out)
{
    if (!isfinite(matter->temp) || matter->temp <= 0)
        return false;
    M1Quantities isotropic = *moments;
    for (int sp = 0; sp < total_num_species; ++sp)
    {
        const BS_REAL n = moments->n[sp], J = moments->J[sp];
        if (!isfinite(n) || !isfinite(J) || n < 0 || J < 0 ||
            ((n == 0) != (J == 0)))
            return false;
        isotropic.chi[sp] = BS_REAL(1) / 3;
    }
    CalculateThickParamsFromM1(&isotropic, out);
    for (int sp = 0; sp < total_num_species; ++sp)
    {
        const BS_REAL n = moments->n[sp], J = moments->J[sp];
        out->w_f[sp] = 0;
        out->beta_f[sp] = 0;
        out->c_f[sp] = CONST_C_F;
        out->temp_f[sp] = 1; // finite, unused component
        if (n == 0)
        {
            out->w_t[sp] = 0;
            out->eta_t[sp] = 0;
            out->temp_t[sp] = matter->temp;
            continue;
        }
        BS_REAL eta = out->eta_t[sp];
        BS_REAL T = out->temp_t[sp];
        if (!isfinite(eta) || !isfinite(T) || T <= 0)
            return false;

        const BS_REAL log_n = log(n), log_J = log(J);
        const BS_REAL target = 4 * log_n - 3 * log_J - log(kBS_FourPi_hc3);
        // The ordinary rational inverse is an initial guess, not a new
        // normalization tolerance inside a tightly converged source solve.
        // eta=20 brackets target<=log(.7); the degenerate inverse above
        // handles larger ratios without an upper degeneracy cap.
        if (eta >= -20 && target <= log(BS_REAL(0.7)))
        {
            BS_REAL lower = -20, upper = 20;
            bool converged = false;
            for (int iteration = 0; iteration < 16; ++iteration)
            {
                const BS_REAL F1 = FDI_p1(eta);
                const BS_REAL F2 = FDI_p2(eta);
                const BS_REAL F3 = FDI_p3(eta);
                const BS_REAL residual = 4 * log(F2) - 3 * log(F3) - target;
                const BS_REAL tolerance = 32 * std::numeric_limits<BS_REAL>::epsilon() *
                                           (1 + fabs(target));
                if (fabs(residual) <= tolerance)
                {
                    converged = true;
                    break;
                }
                if (residual > 0)
                    upper = eta;
                else
                    lower = eta;
                const BS_REAL derivative = 8 * F1 / F2 - 9 * F2 / F3;
                const BS_REAL next = eta - residual / derivative;
                eta = isfinite(next) && next > lower && next < upper
                          ? next : (lower + upper) / 2;
            }
            if (!converged)
                return false;
            T = exp(log_J - log_n + log(FDI_p2(eta)) - log(FDI_p3(eta)));
        }
        if (!ThickDistributionMatchesM1(n, J, T, eta))
            return false;
        out->w_t[sp] = 1;
        out->temp_t[sp] = T;
        out->eta_t[sp] = eta;
    }
    return true;
}

/* Integrand for computing neutrino number density
 *
 * Computes this for three neutrino species
 */
CCTK_HOST CCTK_DEVICE inline
MyQuadratureIntegrand NuNumberIntegrand(BS_REAL* x, void* p)
{
    NuDistributionParams* distr_pars = (NuDistributionParams*)p;

    MyQuadratureIntegrand result;

    result.n = total_num_species;
    for (int idx = 0; idx < total_num_species; idx++)
    {
        result.integrand[idx] = x[0] * x[0] * TotalNuF(x[0], distr_pars, idx);
    }

    return result;
}

/* Compute neutrino number density
 *
 * Computes this for three neutrino species
 */
/*
inline MyQuadratureIntegrand NuNumber(NuDistributionParams* distr_pars)
{
    MyFunctionMultiD integrand;

    integrand.dim                       = 1;
    integrand.params                    = distr_pars;
    integrand.my_quadrature_integrand.n = total_num_species;

    MyQuadrature quad = quadrature_default;

    GaussLegendreMultiD(&quad);

    BS_REAL s[total_num_species];

    s[id_nue]  = fabs(distr_pars->temp_t[id_nue] * distr_pars->eta_t[id_nue]);
    s[id_anue] = fabs(distr_pars->temp_t[id_anue] * distr_pars->eta_t[id_anue]);
    s[id_nux]  = s[id_nue]; // @TODO: cannot be equal to zero
    s[id_anux] = s[id_nue]; // @TODO: cannot be equal to zero

    integrand.function = &NuNumberIntegrand;
    MyQuadratureIntegrand result =
        GaussLegendreIntegrate1D(&quad, &integrand, s);
    const BS_REAL result_factor = 4. * kBS_Pi / POW3(kBS_H * kBS_Clight);

    for (int species = 0; species < total_num_species; ++species)
    {
        result.integrand[species] = result.integrand[species] * result_factor;
    }

    return result;
}
*/

/* Integrand for computing the neutrino energy density
 *
 * Computes this for three neutrino species
 */
inline MyQuadratureIntegrand NuEnergyIntegrand(BS_REAL* x, void* p)
{
    MyQuadratureIntegrand result = NuNumberIntegrand(x, p);

    for (int nuid = 0; nuid < total_num_species; ++nuid)
    {
        result.integrand[nuid] = x[0] * result.integrand[nuid];
    }

    return result;
}

/* Compute neutrino energy density
 *
 * Computes this for three neutrino species
 */
/*
inline MyQuadratureIntegrand NuEnergy(NuDistributionParams* distr_pars)
{
    MyFunctionMultiD integrand;

    integrand.dim                       = 1;
    integrand.params                    = distr_pars;
    integrand.my_quadrature_integrand.n = total_num_species;

    MyQuadrature quad = quadrature_default;

    GaussLegendreMultiD(&quad);

    BS_REAL s[total_num_species];

    s[id_nue]  = fabs(distr_pars->temp_t[id_nue] * distr_pars->eta_t[id_nue]);
    s[id_anue] = fabs(distr_pars->temp_t[id_anue] * distr_pars->eta_t[id_anue]);
    s[id_nux]  = s[id_nue]; // @TODO: cannot be equal to zero
    s[id_anux] = s[id_nue]; // @TODO: cannot be equal to zero

    integrand.function = &NuEnergyIntegrand;
    MyQuadratureIntegrand result =
        GaussLegendreIntegrate1D(&quad, &integrand, s);
    const BS_REAL result_factor = 4. * kBS_Pi / POW3(kBS_H * kBS_Clight);

    for (int species = 0; species < total_num_species; ++species)
    {
        result.integrand[species] = result.integrand[species] * result_factor;
    }

    return result;
}
*/

CCTK_HOST CCTK_DEVICE inline
void ComputeM1DensitiesEq(const MyEOSParams* eos_pars,
                          const NuDistributionParams* nu_distribution_params,
                          M1Quantities* m1_pars)
{
    const BS_REAL n_prefactor = kBS_FourPi_hc3 * POW3(eos_pars->temp);
    const BS_REAL j_prefactor = n_prefactor * eos_pars->temp;

    for (int nuid = 0; nuid < total_num_species; ++nuid)
    {
        m1_pars->n[nuid] =
            n_prefactor * FDI_p2(nu_distribution_params->eta_t[nuid]);
        m1_pars->J[nuid] =
            j_prefactor * FDI_p3(nu_distribution_params->eta_t[nuid]);

        BS_ASSERT(isfinite(m1_pars->n[nuid]),
                  "Neutrino (species=%d) number density computed assuming "
                  "equilibrium is not finite (n=%e nm^-3).",
                  nuid, m1_pars->n[nuid]);
        BS_ASSERT(isfinite(m1_pars->J[nuid]),
                  "Neutrino (species=%d) energy density computed assuming "
                  "equilibrium is not finite (J=%e MeV nm^-3).",
                  nuid, m1_pars->J[nuid]);
    }

    return;
}

#endif // BNS_NURATES_SRC_DISTRIBUTION_DISTRIBUTION_HPP_
