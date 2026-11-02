#ifndef NUX_M1_SPECTRAL_TYPES_HXX
#define NUX_M1_SPECTRAL_TYPES_HXX

// Reuse the established FD integrals/inverse and per-member moment units.
// This header does not initialize or call the NuRates reaction provider.
#include "m1_opacities.hpp"

namespace nuX_M1 {
struct spectral_opacities_t {
  CCTK_REAL abs_0, abs_1, scat_1, eta_0, eta_1;
};
} // namespace nuX_M1

#endif
