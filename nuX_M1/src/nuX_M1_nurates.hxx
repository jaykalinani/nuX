#ifndef NUX_M1_NURATES_HXX
#define NUX_M1_NURATES_HXX

#include "nuX_M1_opacity_utils.hxx"
#include "m1_opacities.hpp"
#include <limits>

namespace nuX_M1 {

// Keep thermal and number-conserving energy exchange separate until after
// the legacy grey equilibrium correction. These are local values, not GFs.
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

// Independent numerical N/E floors do not specify a physical spectrum.
// Define the temporary reference in the SAME comoving units on both the
// scheduled-opacity and source-residual paths. Never modify evolved moments
// or the unmodified moments used for the equilibrium conservation budget.
CCTK_HOST CCTK_DEVICE inline void nurates_rate_reference(
    M1Quantities &radiation, const CCTK_REAL volume,
    const CCTK_REAL matter_temperature, const CCTK_REAL number_floor,
    const CCTK_REAL energy_floor, bool floor_only[3]) {
  constexpr CCTK_REAL slack =
      1 + 64 * std::numeric_limits<CCTK_REAL>::epsilon();
  for (int sp = 0; sp < 3; ++sp) {
    const CCTK_REAL flavors = sp == 2 ? 4 : 1;
    const CCTK_REAL number_limit =
        number_floor / (volume * flavors) * nuX_ndens_conv;
    const CCTK_REAL energy_limit =
        energy_floor / (volume * flavors) * nuX_edens_conv;
    floor_only[sp] = radiation.n[sp] <= slack * number_limit &&
                     radiation.J[sp] <= slack * energy_limit;
    if (floor_only[sp]) {
      const CCTK_REAL reference_number =
          fmin(radiation.n[sp], radiation.J[sp] / (3 * matter_temperature));
      radiation.n[sp] = reference_number;
      radiation.J[sp] = 3 * matter_temperature * reference_number;
    }
  }
  radiation.n[id_anux] = radiation.n[id_nux];
  radiation.J[id_anux] = radiation.J[id_nux];
}

CCTK_HOST CCTK_DEVICE inline nurates_opacities_t
correct_nurates_opacities(const nurates_raw_opacities_t &raw,
                          const bool heavy_species, const bool floor_only,
                          const CCTK_REAL number, const CCTK_REAL energy,
                          const bool number_normalization_defined,
                          const bool energy_normalization_defined,
                          const CCTK_REAL equilibrium_number,
                          const CCTK_REAL equilibrium_energy,
                          const CCTK_REAL correction_max,
                          const CCTK_REAL number_floor,
                          const CCTK_REAL energy_floor) {
  const CCTK_REAL correction =
      floor_only ? CCTK_REAL(1)
                 : opacity_mean_energy_correction(
                       number, energy, equilibrium_number, equilibrium_energy,
                       correction_max);
  nurates_opacities_t result;
  result.scat_1 = correction * raw.scat_1;

  CCTK_REAL absorption_thermal, emission_thermal;
  if (heavy_species) {
    result.eta_0 = correction * raw.eta_0;
    emission_thermal = correction * raw.eta_1_thermal;
    result.abs_0 = equilibrium_number > number_floor
                       ? result.eta_0 / equilibrium_number
                       : CCTK_REAL(0);
    absorption_thermal = equilibrium_energy > energy_floor
                             ? emission_thermal / equilibrium_energy
                             : CCTK_REAL(0);
  } else {
    result.abs_0 = correction * raw.abs_0;
    absorption_thermal = correction * raw.abs_1_thermal;
    // A zero grey denominator below the library threshold does not mean
    // zero spontaneous emission. Preserve the integrated rate there.
    result.eta_0 = number_normalization_defined
                       ? result.abs_0 * equilibrium_number
                       : correction * raw.eta_0;
    emission_thermal = energy_normalization_defined
                           ? absorption_thermal * equilibrium_energy
                           : correction * raw.eta_1_thermal;
  }
  result.abs_1 = absorption_thermal + raw.abs_1_nonthermal;
  result.eta_1 = emission_thermal + raw.eta_1_nonthermal;
  return result;
}

} // namespace nuX_M1

#endif // NUX_M1_NURATES_HXX
