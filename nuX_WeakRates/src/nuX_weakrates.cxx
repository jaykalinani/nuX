#include <cctk.h>
#include <cctk_Parameters.h>

#include "nuX_weakrates.hxx"

namespace nuX_WeakRates {

CCTK_HOST void WeakRates::init() {
  DECLARE_CCTK_PARAMETERS;

  include_beta = use_beta;
  include_pair = use_pair_procs;
  include_plasmon = use_plasmon;
  include_nonbeta_absorption = use_nonbeta_detailed_balance;
  include_bremsstrahlung = use_bremsstrahlung;
  include_elastic = use_elastic;
  beta_low_density_threshold = beta_low_density_rho_threshold;
}

} // namespace nuX_WeakRates
