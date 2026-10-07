#include <cctk.h>
#include <cctk_Arguments.h>
#include <loop_device.hxx>

#include <algorithm>
#include <array>
#include <cmath>

#include "kernels.hpp"

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
    CCTK_VERROR("NuRates pair/bremsstrahlung self-test failed %d checks",
                failures);
  }
  CCTK_INFO("NuRates pair/bremsstrahlung self-test passed");
}

} // namespace nuX_NuRates
