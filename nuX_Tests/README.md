# Numerical checks

`test/unit_tests.par` runs the M1, NuRates and WeakRates checks through the
existing Cactus test workflow. It activates CarpetX and the required thorns;
activating only `nuX_Tests` is not a supported standalone parameter file.
The checks use local values, not diagnostic grid functions.

## Replay an exact source input

In a working tabulated-EOSX/NuRates parameter file, activate `nuX_Tests` and set:

```text
nuX_Tests::run_nurates_source_replay = yes
nuX_Tests::nurates_source_replay_file = "/absolute/path/to/nuX/nuX_Tests/data/warm_source_3470487.txt"
```

Use the normal MPI launcher with `-P` to execute these startup checks without
evolution. An empty fixture path retains the original built-in cold-cell case.
The optional replay is off by default and is not an additional automatic CI
simulation. It calls the production source on host and the active backend;
`host-parallel` output from a CPU build is not GPU evidence.

The text format has ordered, labeled rows: version, cell, coordinates,
scalars, matter, solver, six four-vectors, four four-by-four tensors, three
species, and rate options. Values are finite round-trip doubles; integer
fields must be integers. No geometry or moments are reconstructed or repaired
while loading. Reactions, corrections and EOS come from the active parameter
file, not from the fixture. Both source coordinate time and macro proper time
are stored separately and must not be interchanged.
`N`, `E` and `F_i` are densitized evolved-group code-unit moments; the third
group sums four heavy flavors. Matter `n_b` is in nm^-3 and temperature and
chemical potentials are in MeV. Both stored time intervals are in code units.

`data/warm_source_3470487.txt` preserves the captured cell `(66,65,60)` at
`(-1,-3,-13)` from job 3470487, byte-for-byte (SHA256
`1a28d518275979507bcc0850683e3b85a603c6f81112aa3603a845d8f215a96a`).
It was captured with all five NuRates channels, BRT06, blocking retained,
the low-density beta fallback, and the APRLDP tabulated EOS. The file alone
does not pin that EOS or the runtime reaction settings: retain/check those
separately when reproducing it.

In validation job 3471324 the historical cold replay passed on host and CUDA;
this warm case still failed on both. It is retained as a known failing input,
not a claimed TOV fix. A successful source replay alone does not validate
stellar evolution, and a failed LTE probe does not prove a failed Jacobian:
that derivative is skipped when the probe cannot produce a residual.
