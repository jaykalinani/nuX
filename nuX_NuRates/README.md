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
absorption. The integrated emissivities and stimulated absorptivities are
used directly. No additional LTE Kirchhoff reconstruction or mean-energy
correction is applied to this actual-spectrum calculation. Those are
alternative grey prescriptions, not corrections to apply a second time.

M0 and LeakageBase currently omit inelastic neutrino-electron scattering.  A
leakage free-emission rate cannot represent redistribution of energy between
matter and an existing neutrino distribution; that process requires a
dedicated energy-exchange closure.  Elastic scattering remains part of their
transport opacity.

## M1 Collision Spectrum and Numerical Floors

M1 uses one moment-matched Fermi-Dirac spectrum for its scalar collision
integrals at every transport flux factor. This is a declared spectral-model
choice: it replaces the FD/modified-Maxwell-Boltzmann mixture in the M1
adapter, whose occupation can exceed one in degenerate matter. The original
vendor mixture API is retained separately. The Minerbo transport closure is
unchanged; a bounded scalar occupation is not a proof of a full fermionic
angular closure.

Scheduled opacities, source trials and final source verification call the
same evaluator on the actual comoving number and energy moments. There is
no floor-triggered substitution of a matter-temperature spectrum. Such a
substitution previously caused discontinuous rate jumps near the numerical
floors. Exact vacuum has zero occupation but can receive spontaneous
emission. Every positive population is normalized without an absolute
population cutoff. Transport floors and realizability repair remain
separate from the physical collision operator.

Independent number/energy floors can still imply an unphysical mean energy;
removing the rate discontinuity does not make those initial conditions, the
kernel energy range, or finite-quadrature errors physically validated.

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

The vendor mixture retains the ordinary and dilute fits. In the M1 collision
constructor the ordinary fit supplies an initial guess for a safeguarded
FD moment-ratio refinement, so its normalization error is not imported into
the nonlinear source equation. The dilute branch uses its small-fugacity
expansion rather than extrapolating the ordinary fit to zero population.
Near the filled-sphere limit, input roundoff can make finite temperature
unidentifiable. The routine reports failure there or outside the FD domain,
rather than choosing a maximum eta. Exactly unused thick components retain
finite placeholders, and exact vacuum retains zero occupation.
The coupled solver requests actual FD coordinates for its baseline guess.
If the provisional moments themselves have no FD inverse, a finite
mean-energy-matched trial supplies an initial guess only; it never replaces
the physical moments, rates, or final residual verification. A valid relaxed
root is not rejected solely because its provisional state had no inverse.

The optional test thorn covers FD moment recovery across transport closure
factors, occupation/blocking bounds, eta=33/860, dilute/vacuum states,
heavy-flavor conversion and native coefficient preservation. It uses no
diagnostic grid functions. CPU/CUDA cold/warm source replays have passed.
The current number-current transport and test changes have compiled with CUDA
and completed a four-rank/four-GPU two-step poisoned-grid smoke with finite
norms. Full-resolution profiling and longer TOV validation are separate runs.
Earlier results from different collision-spectrum models are not validation
of this version. None of the local or small-grid checks establishes full TOV
stability or finite-quadrature convergence. A subsequent full-resolution TOV
profiling attempt failed its coupled source solve from floor-sized radiation
moments, followed by the device-abort/CUDA error. The run did not complete its
iteration window and must not be used as a throughput measurement. Its root
cause and the physical consistency of the initialized number/energy floors
remain under investigation; no parameter or rate correction is implied by
the earlier successful small-grid checks.

## Future Updates

When `bns_nurates` is updated, refresh `nuX_NuRates/src` from the new standalone
snapshot, preserve the Cactus-specific unit conversions and setup glue, and
update the snapshot information in this README.
