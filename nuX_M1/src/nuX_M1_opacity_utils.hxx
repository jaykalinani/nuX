#ifndef NUX_M1_OPACITY_UTILS_HXX
#define NUX_M1_OPACITY_UTILS_HXX

#include <cctk.h>

#include <cmath>

namespace nuX_M1 {

using std::fmax;
using std::fmin;
using std::isfinite;

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline bool
equilibrium_moments_are_valid(const CCTK_REAL number,
                              const CCTK_REAL energy) {
  if (!isfinite(number) || !isfinite(energy) || number < CCTK_REAL(0) ||
      energy < CCTK_REAL(0)) {
    return false;
  }

  // A vanishing equilibrium population has no energy. Conversely, a finite
  // temperature population must have both a number and an energy moment.
  return (number == CCTK_REAL(0)) == (energy == CCTK_REAL(0));
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
fallback_equilibrium_moments(CCTK_REAL &number, CCTK_REAL &energy,
                             const CCTK_REAL fallback_number,
                             const CCTK_REAL fallback_energy) {
  if (equilibrium_moments_are_valid(number, energy)) {
    return;
  }

  if (equilibrium_moments_are_valid(fallback_number, fallback_energy)) {
    number = fallback_number;
    energy = fallback_energy;
  } else {
    number = CCTK_REAL(0);
    energy = CCTK_REAL(0);
  }
}

// Opacities scale approximately with the square of the neutrino mean energy.
// Return unity whenever either mean energy is undefined; applying an arbitrary
// extreme correction to a zero or invalid population is less physical than
// retaining the uncorrected coefficient.
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
opacity_mean_energy_correction(const CCTK_REAL radiation_number,
                               const CCTK_REAL radiation_energy,
                               const CCTK_REAL equilibrium_number,
                               const CCTK_REAL equilibrium_energy,
                               const CCTK_REAL correction_max) {
  if (!isfinite(radiation_number) || !isfinite(radiation_energy) ||
      !isfinite(equilibrium_number) || !isfinite(equilibrium_energy) ||
      !isfinite(correction_max) || radiation_number <= CCTK_REAL(0) ||
      radiation_energy <= CCTK_REAL(0) ||
      equilibrium_number <= CCTK_REAL(0) ||
      equilibrium_energy <= CCTK_REAL(0) ||
      correction_max < CCTK_REAL(1)) {
    return CCTK_REAL(1);
  }

  const CCTK_REAL ratio =
      (radiation_energy * equilibrium_number) /
      (radiation_number * equilibrium_energy);
  if (!isfinite(ratio) || ratio <= CCTK_REAL(0)) {
    return CCTK_REAL(1);
  }

  const CCTK_REAL correction = ratio * ratio;
  if (!isfinite(correction)) {
    return correction_max;
  }
  return fmax(CCTK_REAL(1) / correction_max,
              fmin(correction, correction_max));
}

} // namespace nuX_M1

#endif // NUX_M1_OPACITY_UTILS_HXX
