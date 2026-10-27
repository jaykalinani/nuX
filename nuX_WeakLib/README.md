# Optional WeakLibInterp provider

`nuX_WeakLib` is a real, opt-in grey M1 provider for WeakLib's main EmAb
corrected-absorption arrays and both elastic-scattering Legendre moments.
EOSX remains the hydrodynamic EOS and supplies chemical potentials. No WLI EOS
table is loaded; no EOS energy or rest-energy shift is changed. Existing
NuRates/WeakRates/FakeRates configurations do not acquire a WLI dependency.

This is a **partial WeakLib coupling**, not an implementation of all WeakLib
interactions. Pair, NES, Brem and the separate nuclear `EC_table` are not yet
coupled. Their reserved switches fail parameter checking AND table setup when
enabled. A file containing `EC_table` can still supply its main EmAb arrays;
startup explicitly reports that the separate nuclear subtable is not used.

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
WeakLib's [table creator](https://github.com/starkiller-astro/weaklib/blob/5836de983d2ea1305af31bf971099702759b1889/TableCreator/OpacitySource/Bruenn85/wlCreateOpacityTable.f90)
records the kernel unit strings.

- WLI interpolation receives `log10(E/MeV)`, `log10(rho/[g cm^-3])`,
  `log10(T/K)`, and raw Ye. EOSX temperature is MeV and is converted with
  `k_B=8.617333262145e-11 MeV/K`. The HDF5 energy/thermodynamic/kernel units,
  uniform interpolation-coordinate axes, shape, species count, log/linear conventions and offsets are
  validated before use. Unsupported units/layouts fail loudly.
  WLI's index formula is not a general nonuniform-grid search: log E/rho/T
  and linear Ye must be uniformly spaced. Zoomed/nonuniform tables are not
  silently interpolated using the wrong cell.
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
- Iso entries are in `cm^-1 MeV^-2`. They are Legendre expansion coefficients,
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
active AMReX backend during table setup. Synthetic tests cover distinct E^2
and E^3 opacity averages, signed first scattering moment, angular/E^2
normalization, four-heavy-species moments, detailed balance, positive vacuum
emission, elastic-only behavior, EOS-independent table bounds, invalid inputs
and inadequate energy coverage. This implementation was only statically
validated: **no C++ compilation, library rebuild, or simulation was run**.
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
