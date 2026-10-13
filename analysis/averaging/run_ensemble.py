"""Reproduce the 100-realization, D=0 averaging-order experiment.

The supported archive layout and fixed physical settings are documented in README.md.
Inputs and output directory are always explicit; no files are written into the archive.
"""
import argparse
import concurrent.futures
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time

import numpy as np


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def parameters(path):
    return {k: float(v) for k, v in (
        line.split('=', 1) for line in path.read_text().splitlines()
        if '=' in line and not line.lstrip().startswith('#'))}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def prepare(archive, output, executable):
    require(output != archive and archive not in output.parents,
            'Output must be outside the input archive.')
    require(output != Path(__file__).parent.resolve(),
            'Use a separate output directory, not the scripts directory.')
    require(executable.is_file(), 'Executable does not exist.')
    stack = np.loadtxt(archive / 'qx_inverse_cdfs.txt', delimiter=',', skiprows=1)
    rows = np.loadtxt(archive / 'fine_params_all.csv', delimiter=',', skiprows=1)
    require(stack.shape == (100, 200) and rows.shape == (100, 7),
            'This benchmark requires 100 realizations with 100 quantile rows each.')
    require(np.isfinite(stack).all() and np.isfinite(rows[:, :4]).all(), 'Nonfinite input.')
    mean_q = np.loadtxt(archive / 'mean_qx_inverse_cdf.txt', delimiter=',')
    means = parameters(archive / 'mean_params.txt')
    sources = {}
    cases = []
    inputs = {}

    def record(path):
        sources[str(path.relative_to(archive))] = sha(path)

    for i, row in enumerate(rows, 1):
        require(row[0] == i and np.all(row[1:4] > 0), 'Invalid realization IDs or lengths.')
        name = f'r{i:04d}'
        individual = archive / f'fine_{name}' / f'{name}_qx_inverse_cdf.txt'
        meta_path = archive / f'fine_{name}' / f'{name}_meta.txt'
        meta = parameters(meta_path)
        expected = dict(realization=i, D=0, nx=300, ny=100, nu=100, Lx=3, Ly=1,
                        lc=row[1], fitted_velocity_lambda_x=row[2], fitted_velocity_lambda_y=row[3])
        require(all(meta.get(k) == v for k, v in expected.items()),
                f'{name}: metadata differs from the supported D=0 configuration or parameter row.')
        pair = stack[:, 2*(i-1):2*i]
        require(np.allclose(np.loadtxt(individual, delimiter=','), pair, rtol=0, atol=1e-5),
                f'{name}: inverse-CDF column pairing failed.')
        require(np.all(np.diff(pair[:, 0]) > 0) and np.all(np.diff(pair[:, 1]) >= 0)
                and np.all(pair[:, 1] >= 0) and np.allclose(pair[:, 0], np.arange(100)/100),
                f'{name}: invalid quantile grid.')
        inputs[name] = pair
        cases.append(dict(name=name, lc=row[1], lx=row[2], ly=row[3]))
        for path in [individual, meta_path,
                     archive / f'fine_{name}' / f'{name}_BTC_FineScaled.csv',
                     archive / f'fine_{name}' / f'{name}_BTC_FineScaled_derivative.csv']:
            record(path)
    inputs['mean'] = np.column_stack([stack[:, 0], stack[:, 1::2].mean(axis=1)])
    inputs['archived_mean'] = mean_q
    cases.append(dict(zip(['name', 'lc', 'lx', 'ly'], ['mean', *rows[:, 1:4].mean(axis=0)])))
    cases.append(dict(name='archived_mean', lc=means['lc_mean'],
                      lx=means['lambda_x_mean'], ly=means['lambda_y_mean']))
    for name in ['fine_params_all.csv', 'qx_inverse_cdfs.txt', 'mean_params.txt',
                 'mean_qx_inverse_cdf.txt', 'BTC_mean.csv',
                 'upscaled_mean/upscaled_BTC_Upscaled_derivative.csv']:
        record(archive / name)
    manifest = dict(archive=str(archive), binary=str(executable), binary_sha256=sha(executable),
                    source_files=sources, cases=cases,
                    settings=dict(D=0, nx=300, nu=100, Lx=3, Lu=1, dt=.005, end=20,
                                  inlet=1, initial=0, porosity=1, locations=[.5, 1.5, 2.5]))
    existing = output / 'manifest.json'
    if existing.exists():
        require(json.loads(existing.read_text()) == manifest,
                'Inputs, executable, or settings changed. Use a new output directory.')
    else:
        require(not output.exists() or not any(output.iterdir()),
                'Output must be empty or contain a matching experiment manifest.')
    (output / 'inputs').mkdir(parents=True, exist_ok=True)
    (output / 'runs').mkdir(exist_ok=True)
    for name, values in inputs.items():
        np.savetxt(output / 'inputs' / f'{name}.csv', values, delimiter=',', fmt='%.17g')
    existing.write_text(json.dumps(manifest, indent=2) + '\n')
    return cases


def run(case, output, executable, dt=.005, suffix='', nx=300, nu=100, options=()):
    name = case['name'] + suffix
    prefix = output / 'runs' / name
    input_path = output / 'inputs' / (case['name'] + '.csv')
    cmd = [str(executable), str(input_path), '0', str(dt), str(nx), str(nu), '20',
           str(prefix), *(str(case[k]) for k in ['lc', 'lx', 'ly']), *options]
    done = Path(str(prefix) + '_done.json')
    if done.exists():
        previous = json.loads(done.read_text())
        require(previous['command'] == cmd and previous['input_sha256'] == sha(input_path),
                f'{name}: cached settings differ; use a new output directory.')
        require(all(sha(output / file) == digest for file, digest in previous['outputs'].items()),
                f'{name}: cached output was modified.')
        return name, 'cached'
    env = os.environ.copy()
    env['OPENBLAS_NUM_THREADS'] = '1'
    env['OMP_NUM_THREADS'] = '1'
    start = time.monotonic()
    with Path(str(prefix) + '.log').open('w') as log:
        result = subprocess.run(cmd, cwd=output, env=env, stdout=log, stderr=subprocess.STDOUT)
    require(result.returncode == 0, f'{name}: solver failed; inspect its log.')
    hashes = {}
    for kind in ['cdf', 'pdf']:
        file = Path(str(prefix) + f'_{kind}.csv')
        a = np.loadtxt(file, delimiter=',', skiprows=1)
        require(a.ndim == 2 and a.shape[1] == 6 and np.isfinite(a).all()
                and a[-1, ::2].min() >= 19.995, f'{name}: incomplete or nonfinite output.')
        hashes[str(file.relative_to(output))] = sha(file)
    seconds = time.monotonic() - start
    done.write_text(json.dumps(dict(command=cmd, seconds=seconds, returncode=0,
                                   input_sha256=sha(input_path), outputs=hashes), indent=2) + '\n')
    return name, round(seconds, 1)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--archive', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--executable', type=Path, required=True)
    parser.add_argument('--workers', type=int, default=2)
    parser.add_argument('--mode', choices=['prepare', 'pilot', 'baseline', 'time-step', 'grid', 'solver'],
                        default='baseline')
    args = parser.parse_args()
    require(args.workers > 0, '--workers must be positive.')
    archive, output, executable = [p.resolve() for p in (args.archive, args.output, args.executable)]
    cases = prepare(archive, output, executable)
    if args.mode == 'prepare':
        print('Verified 100 paired realizations and wrote the input manifest.')
        return
    jobs = [(c, {}) for c in cases]
    if args.mode == 'pilot':
        jobs = [(c, {}) for c in [cases[0], cases[-2]]]
    elif args.mode == 'time-step':
        jobs = [(c, dict(dt=.0025, suffix='_refined'))
                for c in [cases[-2], cases[-1], cases[0], cases[49], cases[99]]]
    elif args.mode == 'grid':
        jobs = [(cases[-2], dict(suffix='_space', nx=600)),
                (cases[-2], dict(suffix='_rank', nu=200))]
    elif args.mode == 'solver':
        jobs = [(cases[2], dict(suffix='_tight', options=(
            '-probe_ksp_rtol', '1e-10', '-probe_ksp_atol', '1e-12', '-probe_ksp_error_if_not_converged')))]
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = [pool.submit(run, c, output, executable, **settings) for c, settings in jobs]
        for i, future in enumerate(concurrent.futures.as_completed(futures), 1):
            print(i, len(jobs), future.result(), flush=True)


if __name__ == '__main__':
    main()
