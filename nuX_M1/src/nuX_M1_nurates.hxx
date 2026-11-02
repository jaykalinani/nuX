#ifndef NUX_M1_NURATES_HXX
#define NUX_M1_NURATES_HXX

#include "nuX_M1_opacity_utils.hxx"
#include "m1_opacities.hpp"
#include <limits>

namespace nuX_M1 {

// Keep thermal and number-conserving energy exchange separate at the API
// boundary. These are local values, not GFs.
struct nurates_raw_opacities_t {
  CCTK_REAL abs_0;
  CCTK_REAL abs_1_thermal;
  CCTK_REAL abs_1_nonthermal;
  CCTK_REAL scat_1;
  CCTK_REAL eta_0;
  CCTK_REAL eta_1_thermal;
  CCTK_REAL eta_1_nonthermal;
};

struct nurates_opacities_t {
  CCTK_REAL abs_0;
  CCTK_REAL abs_1;
  CCTK_REAL scat_1;
  CCTK_REAL eta_0;
  CCTK_REAL eta_1;
};

CCTK_HOST CCTK_DEVICE inline nurates_raw_opacities_t
convert_nurates_opacities(const M1OpacitiesNonThermalSeparated &rates,
                          const int species, const CCTK_REAL species_factor) {
  if (species == id_nux && species_factor == CCTK_REAL(4)) {
    // The evolved heavy group contains two neutrino and two antineutrino
    // flavors with equal moment weights. Their NEPS kernels need not agree
    // in electron-asymmetric matter, even for identical input spectra.
    return {
        (rates.kappa_0_a[id_nux] + rates.kappa_0_a[id_anux]) *
            (CCTK_REAL(0.5) * nuX_length_conv),
        (rates.kappa_a_th[id_nux] + rates.kappa_a_th[id_anux]) *
            (CCTK_REAL(0.5) * nuX_length_conv),
        (rates.kappa_a_non_th[id_nux] + rates.kappa_a_non_th[id_anux]) *
            (CCTK_REAL(0.5) * nuX_length_conv),
        (rates.kappa_s[id_nux] + rates.kappa_s[id_anux]) *
            (CCTK_REAL(0.5) * nuX_length_conv),
        (rates.eta_0[id_nux] + rates.eta_0[id_anux]) /
            nuX_ndens_conv * nuX_time_conv * 2,
        (rates.eta_th[id_nux] + rates.eta_th[id_anux]) /
            nuX_edens_conv * nuX_time_conv * 2,
        (rates.eta_non_th[id_nux] + rates.eta_non_th[id_anux]) /
            nuX_edens_conv * nuX_time_conv * 2};
  }
  return {rates.kappa_0_a[species] * nuX_length_conv,
          rates.kappa_a_th[species] * nuX_length_conv,
          rates.kappa_a_non_th[species] * nuX_length_conv,
          rates.kappa_s[species] * nuX_length_conv,
          rates.eta_0[species] / nuX_ndens_conv * nuX_time_conv * species_factor,
          rates.eta_th[species] / nuX_edens_conv * nuX_time_conv * species_factor,
          rates.eta_non_th[species] / nuX_edens_conv * nuX_time_conv * species_factor};
}

// The native API already performs the spectral averaging. Preserve its
// independent emissivity and stimulated absorption; do not impose a second
// grey Kirchhoff prescription or multiply by a mean-energy correction.
CCTK_HOST CCTK_DEVICE inline nurates_opacities_t
native_nurates_opacities(const nurates_raw_opacities_t &raw) {
  return {raw.abs_0, raw.abs_1_thermal + raw.abs_1_nonthermal, raw.scat_1,
          raw.eta_0, raw.eta_1_thermal + raw.eta_1_nonthermal};
}

} // namespace nuX_M1

#endif // NUX_M1_NURATES_HXX
