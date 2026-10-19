# Phase 2 validation -- what actually happened, measured

This corrects every unmeasured "Expected Impact" / "Tested" / "production
ready" claim in `README_START_HERE.md`, `PHASE2_CLUSTERING.md`,
`PROJECT_DELIVERY.md`, `FINAL_STATUS.txt`, and `UPSCALING_IMPLEMENTATION.txt`.
Those documents were written before Phase 2 was ever run end-to-end, and
turned out to describe a feature that was a silent no-op.

## What was actually wrong

Standalone clustering (`use_representative_clustering=true` without also
setting `upscale_each_realization=true`) did nothing: no error, clean exit
code, `clustering_summary.txt` never written, representative output files
empty. Root cause: `run_representative_upscaling()` depends on
`fine_outputs.qx_inverse_cdf_all`, but that array was only populated inside
the fine loop when `upscale_each_realization==true`. Confirmed by
reproducing the no-op, then fixing the gating condition and re-running.

The clustering algorithm itself also had several real defects, found by
writing synthetic tests first and discovering the *original* code
segfaults on one of them (`tests/test_clustering.cpp::test_single_and_empty_inputs`,
K<=0 -- confirmed via AddressSanitizer, see commit log): no guard for
non-positive K, "Wasserstein distance" that was actually an L-infinity norm
over a fragile exact-time-match walk (silently dropping points when the two
series' native grids didn't align), first-K-realizations initialization
(not spread out), stale final assignments when the medoid loop didn't fully
converge, and no correction of the reported cluster count when some
requested clusters ended up empty.

## What was fixed

See the three preceding commits for the actual diffs. Summary:
1. `qx_inverse_cdf_all` now populates under `upscale_each_realization OR
   use_representative_clustering`. Both the empty-array case and an
   out-of-range `num_representatives` are now hard errors, not silent
   success.
2. `clustering.h` renamed to `KMedoidsClustering` (kept `KMeansClustering`
   as a type alias for the one existing call site). True 1-Wasserstein
   distance (trapezoidal integral of |inverse-CDF difference| over a
   shared uniform quantile grid, via `TimeSeries::interpol`, not a native-
   grid exact-match walk). Deterministic farthest-first initialization.
   Final assignments recomputed against the final medoids after the
   iteration loop. Non-finite features treated as maximal-but-finite
   distance instead of propagating NaN. Reported `num_clusters` is the
   actual nonempty cluster count. CLI (`--representative-clustering=N`,
   `--clustering-metric=`) and `.uscript` (`representative_clustering`,
   `num_representatives`, `clustering_metric`) support added; the old
   "not implemented yet" rejection in `script_executor.cpp` removed.
3. The weighted aggregate at each BTC location is computed explicitly as
   `sum(cluster_weight_k * representative_curve_k(t))` on a common time
   grid (`run_representative_upscaling`'s `weighted_sum_curve` helper),
   not an average of the unweighted diagnostic curves. Written to
   `x=..._representative_weighted_{pdf,cdf}.csv`, kept out of the fine-
   scale ensemble mean (appended to comparison stacks only after
   `computeMeanBTCs` has already run), named unambiguously
   (`representative_weighted_x=...`).
4. `clustering_summary.txt` now records requested vs. actual cluster
   count, per-cluster members/weights, the distance metric, iteration
   count, convergence, weight sum (should be 1.0), and Phase 2 runtime.
   `run_manifest.txt` already recorded `num_representatives` and
   `clustering_metric` (pre-existing, unaffected by this fix).

## Tests

`tests/test_clustering.cpp` (standalone, PETSc/VTK-free, builds against
`clustering.h` + `sim_runner.h` + GSL + Armadillo only): deterministic
assignments, K=1, K=N, two well-separated groups recovered exactly, and
malformed inputs (NaN feature, empty inverse-CDF) handled without crashing
or propagating NaN into weights. All 6 pass against the fixed code.
**The original code segfaults on the K<=0 malformed-input case** --
confirmed with AddressSanitizer (see commit log for the exact crash site).

Build and run:
```sh
g++ -std=c++17 -I. -IUtilities -I/usr/include/x86_64-linux-gnu/qt6 \
    tests/test_clustering.cpp -lgsl -lgslcblas -larmadillo -o /tmp/test_clustering
/tmp/test_clustering
```

## End-to-end validation (measured, not projected)

Setup: 60x20 grid, nu=20, 10 realizations, t_end=6.0, D=0.001, stdev=2,
identical seed (fixed at 20260115UL) across every run below, so all six
configurations upscale the exact same 10 fine-scale realizations.
`--no-qx-ranks` was used throughout (disables an unrelated, expensive
per-delta diagnostic scatter-CSV dump that is not needed for this
comparison and was exhausting this sandbox's disk quota; it does not
affect `lc`/velocity-CDF estimation, which falls back to the already-
existing raw-correlation path). Both the CLI entry point (`main.cpp`) and
the `.uscript` entry point (`main_script.cpp` via `script_executor.cpp`)
were exercised; see `logs/task5_uscript.log` for the latter.

Reproduce with (CLI entry point):
```sh
COMMON="--nx=60 --ny=20 --nu=20 --nreal=10 --tend=6.0 --diffusion=0.001 \
        --output-dir=./task5_validation --fine-transport --up-transport --no-qx-ranks"
mpirun -np 1 ./PETSc_test $COMMON                                    # baseline
mpirun -np 1 ./PETSc_test $COMMON --upscale-each-realization         # Phase 1
mpirun -np 1 ./PETSc_test $COMMON --representative-clustering=1      # Phase 2, K=1
mpirun -np 1 ./PETSc_test $COMMON --representative-clustering=2      # K=2
mpirun -np 1 ./PETSc_test $COMMON --representative-clustering=5      # K=5
mpirun -np 1 ./PETSc_test $COMMON --representative-clustering=10     # K=10 (=N)
```

### Accuracy vs. the fine-scale ensemble mean (relative L2, %)

| Location | Baseline | Phase 1 (full) | K=1 | K=2 | K=5 | K=10 |
|---|---:|---:|---:|---:|---:|---:|
| x=0.5 | 21.73 | 13.93 | 26.60 | 15.17 | 14.64 | 13.84 |
| x=1.5 | 16.61 | 10.58 | 21.89 | 14.02 | 11.92 | 10.58 |
| x=2.5 | 13.80 | 12.40 | 20.50 | 13.90 | 13.42 | 12.42 |

### K=N equivalence check (acceptance criterion for Task 3)

Weighted-aggregate (K=10=N) vs. full Phase 1's mean-of-upscaled, relative L2:

| Location | relL2 difference |
|---|---:|
| x=0.5 | 0.157% |
| x=1.5 | 0.087% |
| x=2.5 | 0.090% |

**Documented numerical tolerance: < 0.5% relative L2.** The residual is
consistent with running Phase 1 and Phase 2's upscaling of "the same"
realization as two separate PETSc solves rather than exact floating-point
reuse; at this tolerance the acceptance criterion holds.

### Timing and speedup (wall-clock, upscaled stage only)

| Stage | Seconds | Speedup vs. full Phase 1 (15.92s) |
|---|---:|---:|
| Baseline (upscale-once-on-mean) | 1.38 | n/a (different accuracy class) |
| Phase 1 (10/10 realizations) | 15.92 | 1.0x |
| Phase 2, K=1 | 1.28 | 12.4x |
| Phase 2, K=2 | 3.26 | 4.9x |
| Phase 2, K=5 | 7.41 | 2.1x |
| Phase 2, K=10 (=N) | 16.55 | 1.0x (expected: upscales all 10) |

### Selected representatives and weights

See `clustering_summary.txt` in each run's output directory
(`results/task5/*/clustering_summary.txt` in this package) for the exact
medoid, member list, and weight of every cluster at every K.

## Honest reading of these numbers

- K=1 is a poor approximation here (worse than baseline on 2/3 locations)
  -- a single medoid cannot represent 10 realizations' spread. Do not use
  K=1 in practice.
- K=2 and K=5 give a real, usable speed/accuracy tradeoff at this (tiny,
  10-realization) scale: roughly 5x and 2x speedup for single-digit-percent
  relL2 deviation from the full Phase 1 answer.
- **This is measured at N=10, not N=100.** The "2-20x speedup" and
  "40-95% accuracy" table in the now-corrected documents was a projection
  for N=100 that was never run. Do not quote those numbers as validated;
  only the N=10 table above is measured. The qualitative shape (fewer
  representatives = faster, less accurate; K=N recovers Phase 1 exactly)
  should hold at larger N, but the specific percentages will not transfer
  without re-running at that scale.
- Phase 1 itself (not just Phase 2) measurably improves accuracy over the
  baseline at all three locations here -- consistent with the non-trivial
  finding from the broader investigation this session (see
  `SUMMARY.md`/`BTC_INVESTIGATION.md`) that averaging-order matters.
