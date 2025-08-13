#ifndef BNS_NURATES_SRC_AVG_BARYON_MASS_H_
#define BNS_NURATES_SRC_AVG_BARYON_MASS_H_

#include "constants.hpp"

CCTK_DEVICE CCTK_HOST NUX_ATTRIBUTE_NOINLINE
CCTK_REAL AverageBaryonMass(CCTK_REAL mev_mass) {
  // This factor sets rho / m_b to units of fm^-3, which is used for number density throughout the code
  const CCTK_REAL inv_ndens_CU_to_fm_m3 = nuX_ndens_conv / nuX_CUndens_conv;
  return mev_mass * kBS_MeV / (kBS_Clight * kBS_Clight) / nuX_mass_conv * inv_ndens_CU_to_fm_m3;
}


#endif
