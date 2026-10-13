# BTC averaging-order comparison

Compare the same 100 realizations in three ways: the fine-scale mean, the mean of separately upscaled realizations, and one upscaled solution using arithmetic mean inputs. These scripts reproduce the D=0 benchmark; they are not a general launcher for arbitrary diffusion models or geometries.

## Build

Requirements: C++17 compiler, CMake 3.16+, MPI, PETSc, GSL, Armadillo, FFTW3, VTK (CommonCore, CommonDataModel, CommonExecutionModel, IOXML), and a compatible Utilities checkout. The target writes CSVs, but the existing production solver requires VTK to compile. It does not change the application's build.

```sh
cmake -S analysis/averaging -B analysis/averaging/build -DCMAKE_BUILD_TYPE=Release
cmake --build analysis/averaging/build -j 2
python3 -m pip install -r analysis/averaging/requirements.txt
```

PETSc is discovered through pkg-config (`PETSc` or `petsc`); set `PKG_CONFIG_PATH` for a custom installation. Alternatively supply `-DPETSC_INCLUDE_DIR=... -DPETSC_LIBRARY=...`. For a nonstandard VTK installation use `VTK_DIR`, or explicit `VTK_INCLUDE_DIR` and `VTK_LIBRARIES`. Supply `-DUTILITIES_DIR=/path/to/Utilities` to use a different Utilities checkout; the default is the repository submodule. This code was verified with Utilities revision `c70c6ab`, PETSc 3.15.5, and the conditional `PetscCallAbort` compatibility definition in `compat.h`. Shared libraries must be discoverable through the normal runtime loader configuration. No temporary build directory or personal filesystem path is embedded in these scripts.

## Supported archive

Required files:

- `fine_params_all.csv`: 100 sequential realization rows with lc, lambda_x, lambda_y.
- `qx_inverse_cdfs.txt`: 100 rows, 100 paired rank/velocity columns; ranks 0–0.99.
- `mean_params.txt`, `mean_qx_inverse_cdf.txt`, and `BTC_mean.csv`.
- `upscaled_mean/upscaled_BTC_Upscaled_derivative.csv`.
- For each `fine_rNNNN/`: `rNNNN_meta.txt`, `rNNNN_qx_inverse_cdf.txt`, `rNNNN_BTC_FineScaled.csv`, and `rNNNN_BTC_FineScaled_derivative.csv`.

Every individual metadata file must confirm D=0, nx=300, ny=100, nu=100, Lx=3, Ly=1, and matching fitted lengths. Each stack column pair is checked against its named individual distribution; header numbers alone are not trusted. Source hashes and the executable hash are recorded. Reusing an output directory with changed inputs or executable is rejected; completed outputs are checked before reuse. Keep different experiments in different output directories, outside the input archive. Run only one launcher against an output directory at a time.

## Run

Replace the three paths below. Use a small worker count initially; each worker runs a separate single-rank solver. Do not wrap the driver or launcher in `mpirun -n >1`.

```sh
archive=/path/to/verified/archive
output=/path/to/new/experiment
executable="$PWD/analysis/averaging/build/ensemble_probe"

python3 analysis/averaging/run_ensemble.py --archive "$archive" --output "$output" --executable "$executable" --mode prepare
python3 analysis/averaging/run_ensemble.py --archive "$archive" --output "$output" --executable "$executable" --mode pilot --workers 2
python3 analysis/averaging/run_ensemble.py --archive "$archive" --output "$output" --executable "$executable" --mode baseline --workers 2
python3 analysis/averaging/run_ensemble.py --archive "$archive" --output "$output" --executable "$executable" --mode time-step --workers 2
python3 analysis/averaging/run_ensemble.py --archive "$archive" --output "$output" --executable "$executable" --mode grid --workers 2
python3 analysis/averaging/run_ensemble.py --archive "$archive" --output "$output" --executable "$executable" --mode solver

python3 analysis/averaging/rebuild_fine_mean.py --output "$output"
python3 analysis/averaging/validate_results.py --output "$output"
python3 analysis/averaging/analyze_ensemble.py --output "$output"
```

The pilot is reused by the baseline batch. Baseline runs comprise 100 individual realizations plus arithmetic-mean and saved-mean controls. Time-step checks cover the two mean controls and realizations 1, 50, 100. Grid checks independently double the spatial and rank resolution for arithmetic-mean inputs. The solver check reruns realization 3 with tighter tolerances. Full analysis expects baseline, time-step, and grid runs to be complete; the solver check is optional for the validator. Computation can take tens of minutes depending on hardware and worker count.

## Method and limits

The driver calls production `Grid2D::SolveTransport` with a 300×100 grid, domain 3×1, unit inlet, zero initial concentration, porosity one, D=0, dt=0.005, and end time 20. It observes x=0.5, 1.5, 2.5. At D=0 the mixing rate is v/lc; transverse fitted lengths do not affect this rate. Velocity ranks use j/(nu−1) with constant extension beyond the supplied rank 0.99, matching the legacy setup. The arithmetic-mean control averages the same 100 parameter rows and quantile columns.

Fine-scale means are rebuilt from archived concentration curves and their derivatives with equal weights. The comparison uses 20,001 points on 0.005–19.995 with linear interpolation and no extrapolation, normalization, or time shifting. PDF means concentration derivative, not a particle-arrival histogram. Reported moments are conditional on that finite interval. Reusing fine-scale output does not independently verify the historical fine solver or its full boundary configuration.

Validation records small unit-inlet upper-bound violations separately instead of clipping them or claiming the bound passed. Inspect `validation.json`, logs, and numerical-sensitivity tables before interpreting any error reduction. Time/grid checks are not a proof of full convergence for all realizations. A finding for this D=0 benchmark does not establish the same effect for nonzero diffusion or anisotropic cases.

Generated curves, logs, manifests, figures, and `results.md` belong in the chosen output folder, not in the source commit. The local `build/` and `output/` directories are ignored.
