#ifndef NUX_BARYON_MASS_HXX
#define NUX_BARYON_MASS_HXX

#include <cctk.h>

#include "nuX_rate_units.hxx"

namespace nuX_Utils {

CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
AverageBaryonMass(CCTK_REAL const mev_mass) {
  const CCTK_REAL inv_number_density_code_to_transport =
      rate_units::transport_number_density_to_nm3 /
      rate_units::code_number_density_to_nm3;
  return mev_mass * rate_units::mev_to_nm2_g_s2 /
         (rate_units::speed_of_light_nm_s * rate_units::speed_of_light_nm_s) /
         rate_units::code_mass_to_g * inv_number_density_code_to_transport;
}

} // namespace nuX_Utils

#endif // NUX_BARYON_MASS_HXX
