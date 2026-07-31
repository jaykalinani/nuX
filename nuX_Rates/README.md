# nuX_Rates

`nuX_Rates` vendors the microphysics from the standalone `bns_nurates`
library and adapts it to the Cactus/CarpetX thorn build.

## bns_nurates Snapshot

This thorn was updated from:

- Repository: `git@github.com:RelNucAs/bns_nurates.git`
- Branch: `main`
- Commit: `c31f79cc0ae61c65c2d1cf012c344a9d5c64c8f7`
- Commit date: `2026-06-03`
- Commit subject: `Separate non-thermal processes (#18)`

The imported code keeps the `nuX_Rates` Cactus style:

- `CCTK_HOST CCTK_DEVICE inline` annotations are used instead of Kokkos
  annotations.
- Plain math calls are used instead of `Kokkos::` math wrappers.
- Thorn-specific Cactus unit conversions are retained in `src/constants.hpp`.
- Thorn setup still initializes global opacity flags in
  `src/init_nuX_rates.cxx`.

## Local Thorn Option

`nuX_Rates` keeps one local option that is not part of the standalone
`bns_nurates` snapshot:

- `nuX_Rates::beta_low_density_fallback`

When enabled, beta-process absorption uses the THC/WeakRates non-degenerate
low-density fallback below
`nuX_Rates::beta_low_density_rho_threshold`. The default threshold is
`2.0e11 g/cm^3`.

## Bremsstrahlung Selection

`nuX_Rates::brem_implementation` selects the standalone-library bremsstrahlung
implementation:

- `HR98`: Hannestad-Raffelt 1998
- `BRT06`: Burrows et al. 2006
- `GP19`: Guo-Martinez-Pinedo 2019

The GP19 lookup table is imported from `bns_nurates`, but the upstream GP19
path keeps the table host-only and provides a zero device stub. Use `HR98` or
`BRT06` for GPU production runs until the GP19 table is moved to device-safe
storage.

## Future Updates

When `bns_nurates` is updated, refresh `nuX_Rates/src` from the new standalone
snapshot, preserve the Cactus-specific unit conversions and setup glue, and
update the snapshot information in this README.
