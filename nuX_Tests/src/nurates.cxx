#include <cctk.h>
#include <cctk_Arguments.h>
#include <loop_device.hxx>

#include <algorithm>
#include <array>
#include <cmath>

#include "kernels.hpp"
#include "m1_opacities.hpp"

namespace nuX_NuRates {
namespace {

bool close_rel(const BS_REAL x, const BS_REAL y,
               const BS_REAL tolerance = 2.0e-12) {
  const BS_REAL scale = std::max(std::abs(x), std::abs(y));
  return scale == 0.0 || std::abs(x - y) <= tolerance * scale;
}

int check_kernel(const MyKernelOutput &kernel, const BS_REAL balance,
                 const bool zeroth_legendre) {
  int failures = 0;
  for (int species = 0; species < total_num_species; ++species) {
    const BS_REAL emission = kernel.em[species];
    const BS_REAL absorption = kernel.abs[species];
    failures += !std::isfinite(emission) || !std::isfinite(absorption);
    if (zeroth_legendre) {
      failures += emission < 0.0 || absorption < 0.0;
    }
    if (std::isfinite(emission) && std::isfinite(absorption)) {
      failures += !close_rel(absorption, balance * emission);
    }
  }
  return failures;
}

int check_zero_kernel(const MyKernelOutput &kernel) {
  int failures = 0;
  for (int species = 0; species < total_num_species; ++species) {
    failures += kernel.em[species] != 0.0;
    failures += kernel.abs[species] != 0.0;
  }
  return failures;
}

int check_brt06_units() {
  constexpr BS_REAL nb = 1.0e17;
  constexpr BS_REAL temperature = 10.0;
  constexpr BS_REAL xn = 0.9;
  constexpr BS_REAL xp = 0.1;
  constexpr BS_REAL coefficient = 2.0778e23;
  constexpr BS_REAL half = 0.5;
  constexpr BS_REAL mixed_factor = 28.0 / 3.0;

  const BS_REAL rho = nb * kBS_Mb;
  const BS_REAL expected =
      coefficient * half * erg_to_MeV *
      (xn * xn + xp * xp + mixed_factor * xn * xp) * rho * rho *
      std::pow(temperature, 5.5);
  return !close_rel(QBrem_BRT06(nb, temperature, xn, xp), expected);
}

int check_direct_beta_absorptivity() {
  constexpr BS_REAL omega = 3000.0;
  OpacityParams opacity = opacity_params_default_none;
  opacity.use_beta_low_density_fallback = false;

  MyEOSParams eos{};
  eos.nb = 1.0e17;
  eos.temp = 1.0;
  eos.yp = 0.1;
  eos.yn = 0.9;
  eos.mu_e = 180.0;
  eos.mu_p = 1000.0;
  eos.mu_n = 1200.0;

  BS_REAL rates[4]{};
  AbsOpacitySingleLep(omega, &opacity, &eos, kBS_Me, eos.mu_e, rates);

  const BS_REAL mu_hat = eos.mu_n - eos.mu_p - kBS_Q;
  const BS_REAL etanp =
      EtaNP(eos.nb * eos.yn, eos.nb * eos.yp, mu_hat, eos.temp);
  const BS_REAL electron_energy = omega + kBS_Q;
  const BS_REAL phase_space =
      electron_energy * electron_energy *
      std::sqrt(1.0 - kBS_Me * kBS_Me / (electron_energy * electron_energy));
  const BS_REAL expected =
      kBS_Beta_Const * etanp * phase_space *
      OneMinusFermiDistr(electron_energy, eos.temp, eos.mu_e);

  int failures = 0;
  failures += !std::isfinite(rates[0]) || rates[0] <= 0.0;
  failures += !close_rel(rates[0], expected);

  // This state exercises both numerical pathologies fixed by PR 20: the
  // literal complement rounds to zero, while the old detailed-balance
  // construction contains separately clamped tiny and huge exponentials.
  failures += 1.0 - FermiDistr(0.0, 1.0, 1000.0) != 0.0;
  failures += !(OneMinusFermiDistr(0.0, 1.0, 1000.0) > 0.0);
  const BS_REAL old_detailed_balance =
      rates[1] *
      SafeExp((omega - (eos.mu_p + eos.mu_e - eos.mu_n)) / eos.temp);
  failures += close_rel(rates[0], old_detailed_balance, 1.0e-6);
  return failures;
}

int check_nonthermal_separation() {
  MyQuadrature quadrature{};
  quadrature.type = kGauleg;
  quadrature.dim = 1;
  quadrature.nx = 6;
  quadrature.ny = 1;
  quadrature.nz = 1;
  quadrature.x1 = 0.0;
  quadrature.x2 = 1.0;
  GaussLegendre(&quadrature);

  GreyOpacityParams params{};
  params.opacity_flags = opacity_flags_default_none;
  params.opacity_flags.use_inelastic_scatt = 1;
  params.opacity_pars = opacity_params_default_none;
  params.eos_pars.nb = 4.208366627847035e17;
  params.eos_pars.temp = 12.406403541564941;
  params.eos_pars.ye = 0.07158458232879639;
  params.eos_pars.yp = params.eos_pars.ye;
  params.eos_pars.yn = 1.0 - params.eos_pars.yp;
  params.eos_pars.mu_e = 187.1814489040245;
  params.eos_pars.mu_p = 1011.017977368873;
  params.eos_pars.mu_n = 1221.590136808168;
  params.distr_pars = NuEquilibriumParams(&params.eos_pars);
  ComputeM1DensitiesEq(&params.eos_pars, &params.distr_pars,
                       &params.m1_pars);
  for (int species = 0; species < total_num_species; ++species)
    params.m1_pars.chi[species] = 1.0 / 3.0;

  const auto rates = ComputeM1OpacitiesNonThermalSeparated(
      &quadrature, &quadrature, &params);
  int failures = 0;
  for (int species = 0; species < total_num_species; ++species) {
    // With NEPS as the only enabled interaction, number and thermal-energy
    // coefficients must be structurally zero, while energy redistribution is
    // retained in the separate nonthermal pair.
    failures += rates.eta_0[species] != 0.0;
    failures += rates.kappa_0_a[species] != 0.0;
    failures += rates.eta_th[species] != 0.0;
    failures += rates.kappa_a_th[species] != 0.0;
    failures += rates.kappa_s[species] != 0.0;
    failures += !std::isfinite(rates.eta_non_th[species]) ||
                rates.eta_non_th[species] <= 0.0;
    failures += !std::isfinite(rates.kappa_a_non_th[species]) ||
                rates.kappa_a_non_th[species] <= 0.0;
  }
  return failures;
}

} // namespace

extern "C" void nuX_NuRates_TestKernels(CCTK_ARGUMENTS) {
  constexpr BS_REAL omega = 5.0;
  constexpr BS_REAL omega_prime = 7.0;

  MyEOSParams eos{};
  eos.nb = 1.0e17; // 0.1 fm^-3
  eos.temp = 10.0;
  eos.ye = 0.1;
  eos.yp = eos.ye;
  eos.yn = 1.0 - eos.yp;
  eos.mu_e = 20.0;

  PairKernelParams pair{};
  pair.omega = omega;
  pair.omega_prime = omega_prime;
  pair.cos_theta = 1.0;
  pair.mu = 1.0;
  pair.mu_prime = 1.0;

  BremKernelParams brem{};
  brem.omega = omega;
  brem.omega_prime = omega_prime;
  brem.l = 0;
  brem.use_NN_medium_corr = false;

  const BS_REAL balance = SafeExp((omega + omega_prime) / eos.temp);
  int failures = 0;
  failures += check_brt06_units();
  failures += check_direct_beta_absorptivity();
  failures += check_nonthermal_separation();
  failures += check_kernel(PairKernels(&eos, &pair), balance, true);
  failures += check_kernel(BremKernelsLegCoeff(&brem, &eos), balance, true);
  failures += check_kernel(BremKernelsBRT06(&brem, &eos), balance, true);
  failures += check_kernel(BremKernelAbsGP19(&brem, &eos), balance, true);

  // The first Legendre coefficient is allowed to be negative, but it must
  // remain finite and obey the same detailed-balance relation.
  brem.l = 1;
  failures += check_kernel(BremKernelsLegCoeff(&brem, &eos), balance, false);
  brem.l = 0;

  // Exercise the extrapolation paths which previously admitted negative or
  // nonfinite GP19 and HR98 rates.
  for (const auto state : std::array<std::array<BS_REAL, 3>, 3>{{
           {1.0e12, 1.0, 0.1},
           {1.0e17, 100.0, 0.1},
           {1.0e12, 100.0, 0.1},
       }}) {
    eos.nb = state[0];
    eos.temp = state[1];
    eos.ye = state[2];
    eos.yp = eos.ye;
    eos.yn = 1.0 - eos.yp;
    const BS_REAL local_balance = SafeExp((omega + omega_prime) / eos.temp);
    failures +=
        check_kernel(BremKernelAbsGP19(&brem, &eos), local_balance, true);
    failures +=
        check_kernel(BremKernelsLegCoeff(&brem, &eos), local_balance, true);
  }

  // Pure neutron and pure proton matter contain absent NN channels.  Those
  // channels must contribute zero without contaminating the valid channel.
  eos.nb = 1.0e17;
  eos.temp = 10.0;
  for (const BS_REAL ye : std::array<BS_REAL, 2>{0.0, 1.0}) {
    eos.ye = ye;
    eos.yp = ye;
    eos.yn = 1.0 - ye;
    const BS_REAL local_balance = SafeExp((omega + omega_prime) / eos.temp);
    const MyKernelOutput hr98 = BremKernelsLegCoeff(&brem, &eos);
    const MyKernelOutput brt06 = BremKernelsBRT06(&brem, &eos);
    failures += check_kernel(hr98, local_balance, true);
    failures += check_kernel(brt06, local_balance, true);
    failures += !(hr98.abs[id_nue] > 0.0);
    failures += !(brt06.abs[id_nue] > 0.0);
  }

  // Invalid thermodynamic states must be rejected before any divisions,
  // exponentials, or fitted fractional powers are evaluated.
  eos.temp = 0.0;
  failures += check_zero_kernel(PairKernels(&eos, &pair));
  failures += check_zero_kernel(BremKernelsLegCoeff(&brem, &eos));
  failures += check_zero_kernel(BremKernelsBRT06(&brem, &eos));
  failures += check_zero_kernel(BremKernelAbsGP19(&brem, &eos));

  eos.temp = 10.0;
  eos.nb = 0.0;
  failures += check_zero_kernel(BremKernelsLegCoeff(&brem, &eos));
  failures += check_zero_kernel(BremKernelsBRT06(&brem, &eos));
  failures += check_zero_kernel(BremKernelAbsGP19(&brem, &eos));

  eos.nb = 1.0e17;
  brem.omega = 0.0;
  brem.omega_prime = 0.0;
  failures += check_zero_kernel(BremKernelsLegCoeff(&brem, &eos));
  failures += check_zero_kernel(BremKernelsBRT06(&brem, &eos));
  failures += check_zero_kernel(BremKernelAbsGP19(&brem, &eos));
  brem.omega = omega;
  brem.omega_prime = omega_prime;

  pair.omega = 0.0;
  pair.omega_prime = 0.0;
  failures += check_zero_kernel(PairKernels(&eos, &pair));
  pair.omega = omega;
  failures += check_zero_kernel(PairKernels(&eos, &pair));
  pair.omega = 0.0;
  pair.omega_prime = omega_prime;
  failures += check_zero_kernel(PairKernels(&eos, &pair));

  if (failures != 0) {
    CCTK_VERROR("NuRates microphysics self-test failed %d checks",
               failures);
  }
  CCTK_INFO("NuRates microphysics self-test passed");
}

} // namespace nuX_NuRates
