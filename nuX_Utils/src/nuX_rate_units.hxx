#ifndef NUX_RATE_UNITS_HXX
#define NUX_RATE_UNITS_HXX

#include <cctk.h>

namespace nuX_Utils::rate_units {

inline constexpr CCTK_REAL mev_to_nm2_g_s2 = 1.6021766341182763e+8;
inline constexpr CCTK_REAL speed_of_light_nm_s = 2.99792458e+17;
inline constexpr CCTK_REAL mev_mass_to_g =
    mev_to_nm2_g_s2 / (speed_of_light_nm_s * speed_of_light_nm_s);

// nuX evolves number density in a legacy transport unit. One transport unit
// corresponds to 1e9 nm^-3. This is not fm^-3: 1 fm^-3 is 1e18 nm^-3.
inline constexpr CCTK_REAL transport_number_density_to_nm3 = 1.0e9;
inline constexpr CCTK_REAL physical_number_density_fm3_to_nm3 = 1.0e18;
inline constexpr CCTK_REAL code_number_density_to_nm3 = 3.1059132074685616e-37;
inline constexpr CCTK_REAL code_mass_to_g = 1.9889199999999999e+33;
inline constexpr CCTK_REAL code_energy_to_mev = 1.1154161350360074e+60;
inline constexpr CCTK_REAL code_density_to_g_nm3 =
    code_number_density_to_nm3 * code_mass_to_g;
inline constexpr CCTK_REAL code_energy_density_to_mev_nm3 =
    code_number_density_to_nm3 * code_energy_to_mev;
inline constexpr CCTK_REAL code_time_to_s = 4.9257949707731345e-06;
inline constexpr CCTK_REAL code_length_to_nm = 1476625038050.1248;

inline constexpr CCTK_REAL per_cm3_to_per_nm3 = 1.0e-21;
inline constexpr CCTK_REAL per_cm_to_per_nm = 1.0e-7;

} // namespace nuX_Utils::rate_units

#endif // NUX_RATE_UNITS_HXX
