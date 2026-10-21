# nuX_NuRates

`nuX_NuRates` vendors the microphysics from the standalone `bns_nurates`
library and adapts it to the Cactus/CarpetX thorn build.

## bns_nurates Snapshot

This thorn was updated from:

- Repository: `git@github.com:RelNucAs/bns_nurates.git`
- Branch: `main`
- Commit: `c31f79cc0ae61c65c2d1cf012c344a9d5c64c8f7`
- Commit date: `2026-06-03`
- Commit subject: `Separate non-thermal processes (#18)`

The thorn also includes the following reviewed upstream fixes on top of that
snapshot:

- PR 19 (`61dcaedf`): restore the BRT06 bremsstrahlung unit conversion.
- PR 20 (`cb8f210`): evaluate beta-process inverse absorptivities directly.

The imported code keeps the `nuX_NuRates` Cactus style:

- `CCTK_HOST CCTK_DEVICE inline` annotations are used instead of Kokkos
  annotations.
- Plain math calls are used instead of `Kokkos::` math wrappers.
- Thorn-specific Cactus unit conversions are retained in `src/constants.hpp`.
- Thorn setup still initializes global opacity flags in
  `src/init_nuX_nurates.cxx`.

## Local Thorn Option

`nuX_NuRates` keeps one local option that is not part of the standalone
`bns_nurates` snapshot:

- `nuX_NuRates::beta_low_density_fallback`

When enabled, beta-process absorption uses the THC/WeakRates non-degenerate
low-density fallback below
`nuX_NuRates::beta_low_density_rho_threshold`. The default threshold is
`2.0e11 g/cm^3`.

## Bremsstrahlung Selection

`nuX_NuRates::brem_implementation` selects the standalone-library bremsstrahlung
implementation:

- `HR98`: Hannestad-Raffelt 1998
- `BRT06`: Burrows et al. 2006
- `GP19`: Guo-Martinez-Pinedo 2019

The GP19 lookup table is imported from `bns_nurates`, but the upstream GP19
path keeps the table host-only and provides a zero device stub. Use `HR98` or
`BRT06` for GPU production runs until the GP19 table is moved to device-safe
storage.

## Nonthermal Processes

M1 uses the separated opacity API: neutrino-electron inelastic scattering
contributes to energy exchange but not to neutrino-number emission or
absorption.  Thermal equilibrium corrections are applied only to beta, pair,
and bremsstrahlung coefficients before the thermal and nonthermal energy
sources are combined.

M0 and LeakageBase currently omit inelastic neutrino-electron scattering.  A
leakage free-emission rate cannot represent redistribution of energy between
matter and an existing neutrino distribution; that process requires a
dedicated energy-exchange closure.  Elastic scattering remains part of their
transport opacity.

## Dilute Radiation and Numerical Floors

The dilute Fermi-Dirac reconstruction uses its small-fugacity expansion rather
than extrapolating the rational fit to zero density. Numerical floors for
number and energy are in different units and must not define a physical
mean energy. In M1, only when both evolved moments are at their configured
floors, rate evaluation uses temporary reference moments
`n_ref = min(n, J / (3 * T_matter))` and `J_ref = 3 * T_matter * n_ref`.
Neither reference moment exceeds its input. The reconstructed FD component
must reproduce both reference moments within the inverse fit's accuracy;
nonrepresentable moments are reported, not silently accepted at the eta cap.
Evolved moments are unchanged, and populated spectra retain their measured
mean energy. Exact vacuum has zero occupation and can still receive
spontaneous emission.

Invalid inputs or coefficients terminate the NuRates evaluation with an error,
including in release GPU builds. They no longer silently zero all interactions
for the affected species. The existing `opacity_rate_max` limit is unchanged.

## Future Updates

When `bns_nurates` is updated, refresh `nuX_NuRates/src` from the new standalone
snapshot, preserve the Cactus-specific unit conversions and setup glue, and
update the snapshot information in this README.
