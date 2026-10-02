"""Compare explicitly selected transport PDFs on their shared time interval.

Legacy BTC_mean.csv may contain a PDF or a CDF depending on the producing
version. The caller must identify its representation; this tool never infers it
from the filename. No archived data is modified.
"""
import argparse
import csv
from pathlib import Path
import numpy as np


def read_pairs(path):
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError(f'Missing or empty input: {path}')
    data = np.genfromtxt(path, delimiter=',', skip_header=1, ndmin=2)
    if not data.size or data.shape[1] % 2:
        raise ValueError(f'Expected paired time/value columns: {path}')
    pairs = []
    for j in range(0, data.shape[1], 2):
        p = data[:, j:j+2]
        p = p[np.isfinite(p).all(axis=1)]
        if len(p) < 2 or np.any(np.diff(p[:, 0]) <= 0):
            raise ValueError(f'Insufficient or non-increasing times: {path}, pair {j//2}')
        pairs.append(p)
    return pairs


def compare(fine, upscaled, samples):
    left, right = max(fine[0, 0], upscaled[0, 0]), min(fine[-1, 0], upscaled[-1, 0])
    if right <= left:
        raise ValueError('Inputs have no shared time interval')
    t = np.linspace(left, right, samples)
    f = np.interp(t, fine[:, 0], fine[:, 1])
    u = np.interp(t, upscaled[:, 0], upscaled[:, 1])
    denom = np.trapz(f*f, t)
    if denom <= 0:
        raise ValueError('Fine-scale reference has zero norm')
    return dict(t_start=left, t_end=right,
                relative_l2_percent=100*np.sqrt(np.trapz((u-f)**2,t)/denom),
                fine_integral=np.trapz(f,t), upscaled_integral=np.trapz(u,t),
                fine_peak_time=t[f.argmax()], upscaled_peak_time=t[u.argmax()])


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('runs', nargs='+', type=Path)
    p.add_argument('--fine-pdf', required=True, help='Explicitly verified PDF filename relative to each run')
    p.add_argument('--upscaled-pdf', default='upscaled_mean/upscaled_BTC_Upscaled_derivative.csv')
    p.add_argument('--samples', type=int, default=20001)
    p.add_argument('--output', required=True, type=Path)
    args = p.parse_args()
    if args.samples < 2: p.error('--samples must be at least 2')
    rows=[]
    for run in args.runs:
        f,u=read_pairs(run/args.fine_pdf),read_pairs(run/args.upscaled_pdf)
        if len(f)!=len(u): raise ValueError(f'Different column counts in {run}')
        for i,(a,b) in enumerate(zip(f,u)):
            rows.append(dict(run=run.name, pair_index=i, **compare(a,b,args.samples)))
    with args.output.open('w',newline='') as out:
        writer=csv.DictWriter(out,fieldnames=list(rows[0]),lineterminator="\n")
        writer.writeheader();writer.writerows(rows)
    print(f'Wrote {len(rows)} comparisons to {args.output}')

if __name__=='__main__': main()
