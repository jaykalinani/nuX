#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include "bns_nurates.hpp"

namespace nuX_Rates {

OpacityFlags global_opac_flags;
OpacityParams global_opac_params;

extern "C" void nuX_Rates_Setup(CCTK_ARGUMENTS) {
  DECLARE_CCTK_PARAMETERS;

  // reaction flags
  global_opac_flags.use_abs_em = use_beta;
  global_opac_flags.use_brem = use_bremsstrahlung;
  global_opac_flags.use_pair = use_pair_procs;
  global_opac_flags.use_iso = use_elastic;
  global_opac_flags.use_inelastic_scatt = use_inelastic;
  
  // other flags
  global_opac_params.use_WM_ab = corr_WM_abs;
  global_opac_params.use_WM_sc = corr_WM_scat;
  global_opac_params.use_dU = corr_dU;
  global_opac_params.use_dm_eff = corr_dmeff;
  global_opac_params.use_NN_medium_corr = corr_NN_medium;
  global_opac_params.neglect_blocking = neglect_blocking;
  global_opac_params.use_decay = corr_decay;
  global_opac_params.use_BRT_brem = corr_BRT_brem;

}

} // namespace
