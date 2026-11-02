# Optional WeakLibInterp provider

`nuX_WeakLib` is an opt-in grey M1 provider for WeakLib interaction tables.
It includes the main EmAb corrected-absorption arrays, both elastic-scattering
Legendre coefficients, and an optional coupled collision model.
EOSX remains the hydrodynamic EOS and supplies chemical potentials. No WLI EOS
table is loaded; no EOS energy or rest-energy shift is changed. Existing
NuRates/WeakRates/FakeRates configurations do not acquire a WLI dependency.

The original equilibrium-spectrum EmAb/Iso model remains the default. The
opt-in `collision_model="monopole"` adds current-spectrum electron/positron
Pair, NES and native HR98 bremsstrahlung, with an optional analytic plasmon
supplement. This is an angular-monopole approximation, not full first-order
angular WeakLib collision physics. The separate nuclear `EC_table` remains
unsupported and fails setup if selected. A file containing it may still supply
its main EmAb arrays; those arrays can already contain nuclear capture, so
the separate table must never be added blindly.

| Interaction | Coupling | Enabled by default |
|---|---|---|
| Charged-current emission/absorption | Native EmAb arrays | Yes |
| Elastic scattering | Native zeroth/first Legendre coefficients | Yes |
| Electron-positron pair annihilation | Native Pair kernels, angular monopole | No |
| Inelastic electron/positron scattering | Native NES kernels; zero number source | No |
| Nucleon bremsstrahlung | Native HR98 response tables | No |
| Plasmon decay | Analytic supplement with partner blocking and inverse balance | No |
| Separate nuclear electron capture | Not coupled | No |

The coupled source, table loading and optional tests compile with CUDA.
CPU checks and earlier GPU table/source/two-step smokes have passed. The
latest clean executable's GPU startup gate and full TOV validation remain
pending. Compilation and local checks do not establish stellar accuracy.
The physics, comparison with WeakRates/NuRates, and table preparation are
described in [the LaTeX note](../Docs/nuX_WeakLib.tex).

## Opt-in integration

Add the existing upstream wrapper thorn
`WeakLibInterp/cactus/thorns/WeakLibInterp` to your build's arrangements and add
that thorn plus `nuX/nuX_WeakLib` to its thornlist. Do not copy or rebuild a
second AMReX: the wrapper already requires the ET AMReX, HDF5 and MPI providers
and builds/locates `wli_lib` against the AMReX CarpetX uses. To consume an
existing compatible installation, set `WEAKLIBINTERP_DIR` in the option list.
No repository thornlist or local machine configuration is changed here.
The reusable FD math and channel-only plasmon helper require the NuRates and
WeakRates thorns in the build, but neither provider needs to be active at
runtime. Their EOS adapters and complete reaction providers are not called.

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
  E^2/E^3 factors and tail endpoints. Declared `LogInterp=0` and `1` energy
  grids are accepted, including native geometric grids; physical log(E) need
  not be uniformly spaced. The node-index
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

## Coupled monopole reactions

Set `collision_model="monopole"`, enable the desired `include_pair`,
`include_nes`, `include_brem` and `include_plasmon` switches, and supply
`pair_table`, `nes_table` and `brem_table` paths for the native channels.
All remain disabled by default. See the commented additions in
[par/weaklib-provider.par](par/weaklib-provider.par). Keep Minerbo (or the
already selected transport pressure closure): choosing an FD scalar collision
spectrum does not change angular transport or force transport chi to 1/3.

- EOSX supplies fixed matter T, rho, Ye, chemical potentials and free-nucleon
  mass fractions. No WeakLib EOS is read or called. The supplied opacity family
  must still match the EOSX model; an SFHo file is not an APRLDP substitute.
- A bounded FD scalar spectrum is reconstructed from the **current trial**
  comoving number and energy of each species. Heavy occupations are per member
  of the evolved four-species sum. The existing FD integrals/inverse are reused,
  with an explicit conversion between their h*c convention and this provider's
  phase-space measure. No MB mixing, occupation clipping or floor-temperature
  replacement is used. Positive moments outside the FD domain reject a trial.
- Pair and HR98 Brem integrate one physical neutrino/antineutrino event once.
  Both receive the same net particle-number exchange. Heavy sources sum two
  mu/tau pairs, not four copies of a partner population. Pair upper triangles
  transpose energies **and** I/II components/offsets; NES uses its separate
  detailed-balance transpose, without a Pair-style coupling swap.
  Exact Pair crossing also requires matching physical I0/II0 diagonals;
  supplied production tables still need that generation/reciprocity check.
  The operator never averages or clamps an inconsistent table diagonal.
- NES leaves its integrated number coefficients identically zero while
  exchanging energy. Its stimulated A+B attenuation includes diagonal angular
  damping. The optional grey angular ansatz is h(E)=g(E)H/J, not a full l=1
  partner-flux/pressure contraction or a fermionic-angular realizability proof.
- Pair/NES require cm^-1 MeV^-3 kernels, Kelvin matter grids, positive tabulated
  electron eta=mu_e/T and identical physical energy nodes. NES requires NPS=1
  (electron **and** positron scattering). No absolute-value eta substitution,
  extrapolation, silent channel switch or new thermodynamic cutoff is added.
- Brem consumes S_sigma in MeV^-1, not an opacity. It sums the three effective
  densities rho*Xp, rho*Xn, rho*sqrt(Xp*Xn) with weights 1,1,28/3. Its native
  prefactor follows the pinned [thornado consumer](https://github.com/endeve/thornado/blob/1c0ed66b702656d9f7d300317c5debc1bd3ec55a/Modules/Opacities/NeutrinoOpacitiesComputationModule.F90)
  and [physical constants](https://github.com/endeve/thornado/blob/1c0ed66b702656d9f7d300317c5debc1bd3ec55a/Modules/Units/PhysicalConstantsModule.f90),
  including the reference-density normalization already implicit in S_sigma.
  It does not multiply by the local baryon density again. Exactly absent
  constituents contribute zero; positive effective densities outside the table
  fail explicitly. Native HR98 is not the NuRates BRT06 model.
- Plasmon has no native table family. Its separately selectable supplement
  reuses **only** the unblocked WeakRates plasmon kinematics/prefactor helper,
  with a named monoenergetic approximation and a conservative inverse reaction.
  Actual partner blocking and exp[-(E+Eprime)/T] balance replace independent
  per-species Kirchhoff estimates. No beta, pair or elastic WeakRates provider
  is overlaid. Stable FD hole probabilities avoid loss of dilute partner
  emissivity when the occupied spectrum rounds to one.
- The same callback feeds scheduled rates and every all-species nonlinear
  source trial. A Cactus function alias passes existing M1 stage workspace to
  the optional provider, so M1 builds without WLI remain independent. The shared
  six-coordinate source/inner solver is reused; no partner coefficients are
  frozen. The accepted root is stored directly when source_limiter<0 to avoid
  base+(root-base) cancellation. The normal M1 application/backreaction logic
  stays shared; IMEX does not receive a final unweighted source projection.
  One optional host-side alias query selects physical number-current Gamma
  consistently in closure, transport, source RHS and matter backreaction.
  These collision/temporary-current calculations do not replace sub-floor
  moments by configured positive floors. The limiter/application baseline
  uses the same physical moments; any configured final numerical repair is
  applied only after the collision increment. The separately scheduled numerical
  `RepairState` policy is unchanged and remains distinct from collisions.
- Equilibrium initialization uses continuous FD moments consistent with this
  inverse, not coarse finite-band moment sums followed by a different inverse.
  Emissivity and attenuation share the same actual energy quadrature, preserving
  discrete LTE balance. Tail checks bound omitted number/energy populations,
  not opacity-weighted tails or coarse-grid quadrature error: both still require
  physical resolution/energy-band convergence studies.

All two-energy contractions stream scalar accumulators through the resident
tables. No private energy-pair matrices or new diagnostic GFs are allocated.
The large evaluator/source entries are not forcibly inlined. The CUDA build
has completed, but WeakLib throughput and register use have not yet been
measured in the profiling campaign.

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
and inadequate energy coverage. The host interoperability and fixture sources
have compiled; earlier GPU startup checks exercised the actual-HDF fixtures
and provider self-tests. This does not establish exhaustive multi-rank failure
coverage, restart correctness or full-resolution evolution stability.
The opt-in [energy-node tests](src/test_weaklib_nodes.cxx) additionally cover
every node (including both endpoints and ne=2), nonlinear energy-plane values,
genuine varying rho/T/Ye interpolation, an old-locator failure sensitivity
control, and full `evaluate()` on a broad 512-node nonuniform physical grid
against independent physical-E sums at two matter states. Uniform-grid and
previous degenerate/underflow regressions remain. Bounded temporary test-table
allocations are opt-in only; the production resident buffer count is unchanged.
Validate CPU, GPU, multi-rank missing/malformed files, restart cleanup and
source/equilibrium evolution before scientific use.

The new opt-in [collision tests](src/test_weaklib_collision.cxx) add native
Pair/Brem and analytic-plasmon LTE balance, strongly degenerate opposite
chemical potentials (+/-40 and +/-80), finite vacuum emission, off-equilibrium
pair lepton conservation, heavy multiplicity, NES zero-number redistribution,
heating/cooling signs, continuous-FD initialization/inverse consistency and
sub-floor physical number-current consistency, plus explicit domain rejection.
CPU collision/source checks and earlier GPU self-tests and short physical
evolutions have passed. APRLDP candidate tables have been generated and
sample-checked; spectral resolution, table extent and microscopic-model
compatibility still need convergence studies. The latest clean build removes
temporary failure capture without changing the collision equations; its GPU
startup gate precedes the queued full-resolution TOV run. A previous combined
gate ended nonzero at its final provenance-path check after the physical
checks had completed; that job is not reported as an entirely passing gate.
Full TOV stability, broad CPU/GPU comparisons, multi-rank failures and restart
coverage remain to be demonstrated. Exercise the real
scheduled limiter/application path with a sub-floor state,
`provider_update=1` and a nontrivial limiter fraction; the direct stored-root
path alone does not test that baseline contract. Run CST/alias binding and
compare native HDF-family results against an independent absolute-unit oracle,
not only the synthetic-kernel conservation tests. Unit-level native monopole
invariants do not establish angular completeness,
stellar accuracy, general IMEX stability or a stable TOV run.
