# Optional WeakLibInterp provider

`nuX_WeakLib` is a real, opt-in grey M1 provider for WeakLib's main EmAb
corrected-absorption arrays and both elastic-scattering Legendre moments.
EOSX remains the hydrodynamic EOS and supplies chemical potentials. No WLI EOS
table is loaded; no EOS energy or rest-energy shift is changed. Existing
NuRates/WeakRates/FakeRates configurations do not acquire a WLI dependency.

This is a **partial WeakLib coupling**, not an implementation of all WeakLib
interactions. Pair, NES, Brem, plasmon and the separate nuclear `EC_table` are
not yet coupled. The Pair/NES/Brem/EC reserved switches fail parameter checking
AND table setup when enabled. A file containing `EC_table` can still supply its
main EmAb arrays; startup reports that the separate nuclear subtable is not used.

## Opt-in integration

Add the existing upstream wrapper thorn
`WeakLibInterp/cactus/thorns/WeakLibInterp` to your build's arrangements and add
that thorn plus `nuX/nuX_WeakLib` to its thornlist. Do not copy or rebuild a
second AMReX: the wrapper already requires the ET AMReX, HDF5 and MPI providers
and builds/locates `wli_lib` against the AMReX CarpetX uses. To consume an
existing compatible installation, set `WEAKLIBINTERP_DIR` in the option list.
No repository thornlist or local machine configuration is changed here.

Activate `nuX_WeakLib`, select `nuX_M1::rates_lib="WeakLib"` before setup, and
adapt [par/weaklib-provider.par](par/weaklib-provider.par), a configuration
fragment rather than a standalone simulation. The other rates thorns need not
be active. The historically steerable `rates_lib` parameter does not support
switching *into* WeakLib after table setup; an uninitialized provider aborts.

Use opacity tables generated with matter physics consistent with the EOSX EOS.
This integration cannot infer that physics-model compatibility from filenames.
`opacity_tau_trap=-1` and `opacity_corr_fac_max=1` are required: the other
providers' trapped-equilibrium blend and mean-energy-squared correction are
not part of this spectral-average closure. Standard radiation floors and M1
source-solver settings still apply. `opacity_rate_max` is not used by this
provider: invalid/non-finite coefficients abort instead of being zeroed.

## Conventions and grey closure

The actual public API and production consumption were checked against
[WeakLibInterp](https://github.com/jaykalinani/WeakLibInterp), its local
`specs/{cactus-integration,table-format-and-io,opacity-emab-iso,amrex-device-interface}.md`,
and the production [thornado opacity assembly](https://github.com/endeve/thornado/blob/master/Modules/Opacities/NeutrinoOpacitiesComputationModule.F90)
and [matter solver](https://github.com/endeve/thornado/blob/master/Modules/TwoMoment/TwoMoment_NeutrinoMatterSolverModule.F90).
WeakLib's [Bruenn85 table creator](https://github.com/starkiller-astro/weaklib/blob/5836de983d2ea1305af31bf971099702759b1889/TableCreator/OpacitySource/Bruenn85/wlCreateOpacityTable.f90)
records explicit Iso kernel units. The native FullWeakPhysics creator uses a
different label for the same raw-kernel normalization; see the explicit host-only
format declaration below.

- WLI value interpolation receives the exact INTEGER energy-node index,
  `log10(rho/[g cm^-3])`, `log10(T/K)`, and raw Ye. EOSX temperature is MeV
  and is converted with `k_B=8.617333262145e-11 MeV/K`.
  The HDF5 energy/thermodynamic/kernel units,
  uniform matter-coordinate axes, shape, species count, declared log/linear
  conventions and offsets are validated before use. Unsupported units/layouts
  fail loudly.
  WLI's index formula is not a general nonuniform-grid search: log rho/T and
  linear Ye must remain uniformly spaced. Energy is queried only at stored
  nodes: a resident `{0,...,ne-1}` coordinate array REPLACES the old log-energy
  buffer, and value-only wrappers select the exact integer node. This supports
  finite positive increasing nonuniform physical energy grids without physical
  energy interpolation; enabled channels still need identical physical grids.
  Actual energies still determine FD spectra, physical-E trapezoidal weights,
  E^2/E^3 factors and tail endpoints. The declared `LogInterp=1` file schema is
  retained, but physical log(E) need not be uniformly spaced. The node-index
  adapter has NO off-node physical-energy or d/dE contract; do not use WLI
  energy-derivative routines with the artificial index coordinates.
- With `chemical_potentials="total"` (default), use
  `mu_nue=mu_p+mu_e-mu_n` and `mu_anue=-mu_nue`; heavy-lepton chemical
  potentials vanish. This also applies when the SAME reference has been
  subtracted from both nucleon potentials: that common shift cancels.
  Negative individual potentials alone do not justify another correction.
  Only if each nucleon's OWN rest mass has been separately subtracted, the
  `chemical_potentials="rest-mass-subtracted"` option restores
  `m_n-m_p=1.29333236 MeV` in `mu_n-mu_p`. Electron chemical potential must
  include electron rest energy. EOSX's internal-energy shift is unrelated and
  is never used to guess this convention. For tables containing `munu`,
  compare it with the direct combination before choosing the option; that
  check does not establish overall EOS/opacity-model compatibility.
- Define `I_q=integral E^q f_FD(E;T,mu) dE`. Number and energy equilibrium
  moments are `4*pi/(hc)^3 I_2` and `4*pi/(hc)^3 I_3`, with one helicity per
  species, `hc=1.2398419843320026e-10 MeV cm`. Heavy-lepton moments include
  all four species (`nu_mu`, `anti-nu_mu`, `nu_tau`, `anti-nu_tau`). Their
  neutral-current scattering coefficient is the mean of the neutrino and
  antineutrino tables, NOT four times a single-species coefficient.
- Main EmAb entries are corrected absorption `chi(E)` in cm^-1. The two grey
  coefficients are `integral chi E^2 f / I_2` and
  `integral chi E^3 f / I_3`; emissivities are each coefficient times its
  equilibrium moment. No evolved number/energy density is divided into an
  emissivity, so emission is finite in an initially empty radiation field.
  Heavy-lepton EmAb is zero, not borrowed from electron-flavor data.
  Negative-chemical-potential spectra share a factored fugacity when forming
  ratios; densities and emissivities are restored in logarithms. A physically
  negligible absolute population may round to zero without destroying its
  finite mean energy or attenuation. No occupation floor or eta cap is added.
- Accepted Iso entries are interpreted as raw `cm^-1 MeV^-2` kernels. They are Legendre expansion coefficients,
  not already-normalized angular moments: WeakLib's table creator multiplies
  the l-th moment by `(2*l+1)/2`. The physical transport opacity is therefore
  `4*pi*E^2*(phi_0-phi_1/3)`, as in thornado's `Sigma_Iso_P` assembly, averaged
  with `E^3 f`. Signed `phi_1` is valid; neither it nor the physical coefficient
  is artificially clipped.
- Composite trapezoidal weights integrate in **physical E**, not log E,
  between the first and last tabulated energy nodes. All enabled channels
  must have identical energy nodes. Equilibrium moments, grey averages,
  emission and equilibrium initialization use the same quadrature, so the
  discrete grey collision source vanishes at equilibrium. These are
  equilibrium-spectrum (Planck) grey averages, not a reconstruction of a
  nonequilibrium radiation spectrum.
- Low/high omitted FD number and energy tails have conservative bounds.
  Both must be below `spectral_tail_tolerance` times the in-band integral at
  every unmasked point. This checks spectral *coverage*, not quadrature error
  or opacity-weighted tail error; production use needs independent spectral
  quadrature/table-extent convergence tests. No opacity is extrapolated into
  tails. The finite-band approximation is explicit, not a silent full-band
  accuracy claim.
- Conversions reuse `nuX_rate_units.hxx`, including nuX's legacy number-density
  unit (`1 transport unit = 1e9 nm^-3`, not fm^-3). `nueave` is the ratio of
  the converted energy and number moments, in the unit expected by M1.

Only masked cells have zero rates. An unmasked state outside either the EOSX
or any enabled opacity table, a failed coverage test, a negative physical
rate, or non-finite arithmetic aborts. WLI deliberately extrapolates by
default; this consumer never calls its interpolator outside validated ranges.
There is no atmosphere-density cutoff, EOS boundary clamp, arbitrary opacity
cap, or silent invalid-rate fallback.

## Iso file convention (host-only, non-steerable)

`iso_table_convention="explicit-kernel-units"` is the strict default: both
species must carry one of the supported explicit `cm^-1 MeV^-2` labels. It does
not accept `Per Centimeter` or a literal `1/cm` Iso opacity. For a table whose
provenance you have checked against WeakLib's native FullWeakPhysics raw-kernel
convention, explicitly select `"weaklib-fullweakphysics-5836de98"`. This names
[generator revision 5836de983d2ea1305af31bf971099702759b1889](https://github.com/starkiller-astro/weaklib/blob/5836de983d2ea1305af31bf971099702759b1889/TableCreator/OpacitySource/FullWeakPhysics/wlCreateOpacityTable.F90),
which writes `Per Centimeter` but stores the raw Legendre expansion coefficients
with the `(2*l+1)/2` factor. Its [kernel source](https://github.com/starkiller-astro/weaklib/blob/5836de983d2ea1305af31bf971099702759b1889/TableCreator/OpacitySource/FullWeakPhysics/Chimera/scatical_weaklib.f90)
and thornado's consumer establish the numerical interpretation; the legacy
label alone does not. The native option accepts only that spelling (ignoring
case/space padding), not arbitrary inverse-length aliases. Cross-format or
mixed-species labels fail; no format is inferred from a filename or magnitudes.

Both selections pass identical unscaled arrays to the existing
`4*pi*E^2*(phi_0-phi_1/3)` contraction. Neither changes angular physics, an EOS,
or any reaction activation. The selector is not copied into device views.
The host validator checks actual HDF types, ranks and extents, including
two-dimensional `{2,2}` offsets, two matching species arrays, and integer
species/moment counts and coordinate flags. The native option additionally
requires the generator's `{1}` integer `weak_magnetism_corr`, `ion_ion_corr`,
`many_body_corr` (each 0 or 1), and finite double `ga_strange`. Their values
are logged, not changed. Existing numeric-energy, uniform-matter-coordinate,
table-range and physical-kernel checks remain in force. The convention selector
does not relax them; exact-node access separately permits nonuniform physical
energy grids, not nonuniform matter axes.

Only the AMReX I/O root reads HDF metadata. Failure status and bounded diagnostic
text reach every rank before a WLI table-load broadcast. Startup records the
declaration, original unit labels and native corrections. This is a user-declared
convention plus schema validation, **not authenticated generator provenance or
verified opacity/EOSX-model compatibility**. Check those separately; no WLI
matter EOS is loaded. A successful native Iso load is not full WeakLib coupling.

## Lifecycle, tests and remaining work

Tables load collectively at `CCTK_WRAGH` after `EOSX_Setup_EOS`. CarpetX has
already initialized AMReX in `SetupGH` at that point; `Driver_Startup` is too
early. WLI performs root-read/broadcast. Unit/schema checks use the same I/O
root and broadcast failure status before entering the loader. Arrays upload
once to resident WLI tables/device axis vectors; temporary host arrays are
released. Cleanup synchronizes and destroys all owners at shutdown before
`Driver_Shutdown` finalizes AMReX. No new grid functions are allocated. One
transient integer error flag per grid call provides strict host-side errors
even when release GPU builds compile out device assertions.

`run_self_tests=yes` runs [src/test_weaklib.cxx](src/test_weaklib.cxx) on the
active AMReX backend during table setup, plus opt-in host-only
[actual-HDF Iso fixtures](src/test_weaklib_io.cxx). The latter create a private
`/tmp/nux-weaklib-iso-XXXXXX` directory only on the I/O root, exclusively create
small files (768 bytes of raw species data per fixture), and remove only those
owned paths after readers finish. They exercise both declarations through the
real metadata validator and public WLI reader, node/interior recovery, unchanged
angular normalization, malformed units/counts/types/shapes/native flags, and
collective missing-file/root-creation failures. No production table is rewritten.
The native-shaped fixture set also includes a nonuniform physical energy grid.
Fixture data need no EOS, but their opt-in startup entry still follows the
provider's normal EOSX/production-table setup. Synthetic device tests cover
distinct E^2 and E^3 opacity averages, signed first scattering moment, angular/E^2
normalization, four-heavy-species moments, detailed balance, positive vacuum
emission, elastic-only behavior, EOS-independent table bounds, invalid inputs
and inadequate energy coverage. The new host interoperability code and fixture
sources have only been statically checked: **no compilation or fixture execution
was performed for this change**; their presence is not a passing test result.
The opt-in [energy-node tests](src/test_weaklib_nodes.cxx) additionally cover
every node (including both endpoints and ne=2), nonlinear energy-plane values,
genuine varying rho/T/Ye interpolation, an old-locator failure sensitivity
control, and full `evaluate()` on a broad 512-node nonuniform physical grid
against independent physical-E sums at two matter states. Uniform-grid and
previous degenerate/underflow regressions remain. Bounded temporary test-table
allocations are opt-in only; the production resident buffer count is unchanged.
Validate CPU, GPU, multi-rank missing/malformed files, restart cleanup and
source/equilibrium evolution before scientific use.

WLI also exposes two-energy Pair/NES/Brem interpolation. Those APIs alone are
not grey collision models. Pair requires lower-triangle JI/JII interpolation,
upper-triangle energy AND component swapping, flavor C1/C2 assembly and a
coupled evolved-antineutrino spectral closure that gives identical net nue
and anti-nue pair-number sources. An equilibrium-partner bath applied
independently to each grey species violates electron-lepton conservation away
from equilibrium, so it is deliberately not substituted here. NES must
conserve number while redistributing energy; Brem additionally needs the
`rho*xp`, `rho*xn`, `rho*sqrt(xp*xn)` density decomposition with weights
`[1,1,28/3]`, its own normalization and coupled pair closure. Nuclear capture
needs a loader for `EC_table`, which the current public WLI reader does not
expose. These processes remain unfinished work, explicitly rejected rather
than silently approximated.
Plasmon also remains unimplemented and has no native channel in WLI's five-table
API; the format adapter does not resolve that required full-reaction gap.
