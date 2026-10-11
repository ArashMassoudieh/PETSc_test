# BTC investigation — 2026-10-01

## Findings

The archive is `/mnt/2nd/Projects_H/PETSc_UpscalingResults`. These are 2026
runs, produced by the legacy `sim_runner.cpp` / `SolveTransport` pipeline.
The supplied patches affect `script2_executor.cpp` / `SolveCopulaTransport`,
a different pipeline. They cannot retroactively explain or correct the archive.

In the two selected archived runs, `BTC_mean.csv` contains derivatives (PDFs),
not cumulative concentration. Compare it with
`upscaled_mean/upscaled_BTC_Upscaled_derivative.csv`. Generic BTC filenames are
not reliable indicators of representation. Current script2 `save_btc` also
writes derivatives. Particle-arrival PDFs and concentration derivatives should
not be treated as interchangeable without checking injection/observation weights.

## Measured archive mismatch

Relative L2 norms are integrated over the common time interval, using 20,001
uniform comparison points, without extrapolation or normalization. The archived
mean files were used as supplied, not rebuilt from every realization.

| Archive run | x=0.5 | x=1.5 | x=2.5 |
|---|---:|---:|---:|
| 100Realizations_20260309_175930_std2_D0_aniso | 17.86% | 21.25% | 29.98% |
| 100Realizations_20260202_003241_std2_D0.1_aniso1&0.1 | 17.41% | 15.00% | 16.01% |

Both `mean_params.txt` files report nReal=100. Other directories bearing
“100Realizations” report nReal=20. The March D=0.01 neighbor has an empty
`BTC_mean.csv`; it was excluded. The April run reports only two realizations.

Reproduce from the repository root (NumPy required):

```sh
python3 analysis/compare_archive_btc.py \
  '/mnt/2nd/Projects_H/PETSc_UpscalingResults/Finished Runs/100Realizations_20260309_175930_std2_D0_aniso' \
  '/mnt/2nd/Projects_H/PETSc_UpscalingResults/Finished Runs/100Realizations_20260202_003241_std2_D0.1_aniso1&0.1' \
  --fine-pdf BTC_mean.csv --output analysis/archive_btc_metrics.csv
```

## Confirmed new-solver defect and fix

The copula solver added an inlet diffusion source to the RHS but omitted the
corresponding inlet loss from the matrix diagonal for NX>1. Uniform inlet
concentration therefore was not a stationary solution.

Use the half-cell distance from inlet face to first cell centre: inlet diffusion
source and loss are both 2D/dx². Add D/dx² for each interior face and zero for
the zero-gradient outlet. This also handles NX=1. Advective coefficients are
unchanged. This defect does not account for the archived D=0 discrepancy.

A real PETSc uniform-velocity run with nx=8, Lx=1, v=1, D=0.1, dt=0.01,
c_left=1, t_end=5 gave reconstructed final concentrations:

| Solver | x=0.1 | x=0.5 | x=0.9 |
|---|---:|---:|---:|
| Original grid.cpp | 1.799992 | 1.799936 | 1.799860 |
| Corrected grid.cpp | 1.000000 | 0.999993 | 0.999983 |

The executable stores interval derivatives, so concentration was reconstructed
by cumulative sum times dt, not trapezoidal integration of midpoint samples.
The production-assembly regression failed four positive-D cases before the
fix; all nine NX/D combinations pass afterward.

## Correlation output and supplied patches

`save_corr` was parsed but never written. It now writes the existing sampled-field
Pearson diagnostic, with checked output errors, bare/nested path support, and
collective status propagation. Rank-zero directory creation avoids calling the
barrier-containing `ensureDir` inside a rank-zero-only block. Two-process tests
cover valid paths and a blocked output path.

The normal-score transformation and u_cap=0.99 patches are deliberately not
adopted as established physical fixes. Gaussian-copula latent correlation,
spatial correlation length, and copula sampling distance are distinct concepts.
The local Massoudieh & Dentz (2020) paper, Eq. (6), explicitly distinguishes
sampling distance epsilon from longitudinal correlation length lc. Its general
copula equation is not restricted to a Gaussian copula. The script currently
pairs lambda with extraction distance; changing lambda alone changes the rate
without re-extracting a compatible transition kernel.

The sampler also uses Gaussian-distributed spatial offsets and clamps at domain
boundaries. It is not a fixed-distance streamline transition sampler. Fitting
its diagnostic as a fixed-lag exponential needs a separate justification.

## Validation environment

The saved executable requires newer glibc and unavailable libraries. A temporary
build used Ubuntu 22.04 PETSc 3.15.5/VTK 9.1 dependencies extracted under /tmp
(no system installation), a PetscCallAbort compatibility macro, and the clean
sibling `/mnt/3rd900/Projects/Utilities` checkout at c70c6ab. The embedded Utilities
checkout is older and has pre-existing uncommitted changes; it was not altered
or included in this commit. The temporary build is not proof that the default
qmake project or pinned submodule builds unchanged on this machine.

Run assembly tests with:

```sh
python3 -m unittest discover -s tests -p 'test_copula_boundary.py'
```

Run the executable checks with a configured runtime library path and MPI access:

```sh
PETSC_TEST_EXECUTABLE=/path/to/PETSc_test python3 -m unittest discover -s tests
```

## Next physical-model experiment

First reproduce one archived D=0 case with its original solver/configuration.
Then compare the mean of per-realization upscaled solutions with the solution
using averaged parameters and averaged inverse CDFs. Those operations need not
commute. This separates ensemble averaging error from the closure approximation.
Only afterward vary the velocity cap, sampling method, and diffusion closure one
at a time with matching fine/upscaled settings and held-out BTC locations.

The archived residual mismatch is still unresolved. Neither the boundary fix
nor these diagnostic repairs are claimed to close it.
