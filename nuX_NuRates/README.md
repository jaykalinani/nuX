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
mean energy. In this `nurates` feature branch, M1 classifies a state as
floor-only when both comoving moments are below the configured number and
energy floors converted with the cell volume and species multiplicity.
This differs from the original evolved-moment criterion for moving matter;
the two criteria are not algebraically equivalent. Scheduled opacities and
nonlinear source trials currently share the comoving criterion. Its behavior
in the full-star atmosphere is not yet validated.

For a classified floor-only state, rate evaluation uses temporary moments
`n_ref = min(n, J / (3 * T_matter))` and `J_ref = 3 * T_matter * n_ref`.
Neither reference moment exceeds its input. The reconstructed FD component
must reproduce both reference moments within the inverse's accuracy;
nonrepresentable moments are reported, not silently replaced by another spectrum.
Evolved moments are unchanged, and spectra outside the floor-only reference
retain their measured mean energy. Exact vacuum has zero occupation and can
still receive spontaneous emission.

Invalid inputs or coefficients terminate the NuRates evaluation with an error,
including in release GPU builds. They no longer silently zero all interactions
for the affected species. The existing `opacity_rate_max` limit is unchanged.

## Highly Degenerate FD Reconstruction

The high-degeneracy inverse no longer caps eta at 20 or extrapolates the
last rational fit. It solves the forward FD moment ratio in `s=eta^-2`, using
scaled reflection identities and a safeguarded scalar iteration. The physical
upper ratio is `n^4/(C*J^3)=64/81`, where `C=4*pi/(hc)^3`; it is not the old
fit endpoint. Both the gap residual and relative scalar correction must
converge. No moments, rates, source tolerances or radiation floors are rescaled.

The ordinary and dilute fits are retained. At the existing ratio matching
point 0.7 the ordinary fit still has a small mismatch (about 2.3e-6 in moment
normalization); the inverse is therefore not claimed to be globally smooth.
Near the filled-sphere limit, input roundoff can make finite temperature
unidentifiable. The routine reports failure there or outside the FD domain,
rather than choosing a maximum eta. Exactly unused thick components retain
finite placeholders, and exact vacuum retains zero occupation.
The coupled solver requests actual FD coordinates for its baseline guess.
If the provisional moments themselves have no FD inverse, a finite
mean-energy-matched trial supplies an initial guess only; it never replaces
the physical moments, rates, or final residual verification. A valid relaxed
root is not rejected solely because its provisional state had no inverse.

The optional test thorn checks forward/inverse round trips, eta=33/860
production-callback reconstruction, the matching point, dilute states,
vacuum, unused components and Pauli-boundary rejection on host and the active
backend. These checks and the thin/provisional source-guess tests passed on
the host and CUDA backend in validation job 3471324 on 2026-10-27. The
configured EOS/matter-source test and historical cold-cell replay also passed;
cold accepted host/device moments agreed to about 7.1e-15 on the declared
radiation scales. These tests use no diagnostic grid functions.

This inverse fix does not resolve the full-star atmosphere failure. The
separately captured warm cell from job 3470487 still rejects its LTE trial
and fails the full coupled solve on both host and CUDA in job 3471324, despite
valid baseline rates. The later rejection requires separate localization;
neither the unit checks nor cold replay establish full TOV stability.

## Future Updates

When `bns_nurates` is updated, refresh `nuX_NuRates/src` from the new standalone
snapshot, preserve the Cactus-specific unit conversions and setup glue, and
update the snapshot information in this README.
