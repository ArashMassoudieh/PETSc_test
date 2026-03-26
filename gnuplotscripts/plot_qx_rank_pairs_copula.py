#!/usr/bin/env python3
"""Plot Fig.3-style qx rank-pair panels with Gaussian copula contours.

Usage:
  python gnuplotscripts/plot_qx_rank_pairs_copula.py \
      --input-dir /path/to/fine_scale_realization \
      --output qx_rank_pairs_fig3.png
"""

from __future__ import annotations

import argparse
import csv
import math
import re
from pathlib import Path
from typing import List, Tuple

import matplotlib.pyplot as plt
import numpy as np


def norminv_approx(p: np.ndarray) -> np.ndarray:
    # Acklam approximation (vectorized)
    a = np.array([
        -3.969683028665376e01,
        2.209460984245205e02,
        -2.759285104469687e02,
        1.383577518672690e02,
        -3.066479806614716e01,
        2.506628277459239e00,
    ])
    b = np.array([
        -5.447609879822406e01,
        1.615858368580409e02,
        -1.556989798598866e02,
        6.680131188771972e01,
        -1.328068155288572e01,
    ])
    c = np.array([
        -7.784894002430293e-03,
        -3.223964580411365e-01,
        -2.400758277161838e00,
        -2.549732539343734e00,
        4.374664141464968e00,
        2.938163982698783e00,
    ])
    d = np.array([
        7.784695709041462e-03,
        3.224671290700398e-01,
        2.445134137142996e00,
        3.754408661907416e00,
    ])

    plow = 0.02425
    phigh = 1.0 - plow
    p = np.clip(p, 1e-12, 1 - 1e-12)
    x = np.zeros_like(p)

    low = p < plow
    high = p > phigh
    mid = ~(low | high)

    if np.any(low):
        q = np.sqrt(-2.0 * np.log(p[low]))
        x[low] = (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / \
                 ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0)

    if np.any(high):
        q = np.sqrt(-2.0 * np.log(1.0 - p[high]))
        x[high] = -(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / \
                  ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0)

    if np.any(mid):
        q = p[mid] - 0.5
        r = q * q
        x[mid] = (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q / \
                 (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1.0)

    return x


def gaussian_copula_density(u1: np.ndarray, u2: np.ndarray, rho: float) -> np.ndarray:
    z1 = norminv_approx(u1)
    z2 = norminv_approx(u2)
    den = max(1e-12, 1.0 - rho * rho)
    expo = (2.0 * rho * z1 * z2 - rho * rho * (z1 * z1 + z2 * z2)) / (2.0 * den)
    return np.exp(expo) / math.sqrt(den)


def read_rank_pair_csv(path: Path, max_points: int) -> Tuple[np.ndarray, np.ndarray]:
    u1: List[float] = []
    u2: List[float] = []
    with path.open("r", newline="") as f:
        r = csv.DictReader(f)
        for i, row in enumerate(r):
            if i >= max_points:
                break
            u1.append(float(row["u1"]))
            u2.append(float(row["u2"]))
    return np.array(u1), np.array(u2)


def parse_delta_from_name(path: Path) -> float:
    m = re.search(r"qx_rank_pairs_dx_([0-9]+\.[0-9]+)\.csv$", path.name)
    if not m:
        raise ValueError(f"Cannot parse delta_x from: {path.name}")
    return float(m.group(1))


def choose_four_indices(n: int) -> List[int]:
    if n <= 4:
        return list(range(n))
    return [0, n // 3, (2 * n) // 3, n - 1]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input-dir", required=True, type=Path)
    ap.add_argument("--output", required=True, type=Path)
    ap.add_argument("--max-points", type=int, default=15000)
    ap.add_argument("--point-alpha", type=float, default=0.15)
    args = ap.parse_args()

    files = sorted(args.input_dir.glob("*qx_rank_pairs_dx_*.csv"))
    if not files:
        raise SystemExit(f"No rank-pair files found in {args.input_dir}")

    files_with_dx = sorted([(parse_delta_from_name(p), p) for p in files], key=lambda x: x[0])
    chosen = choose_four_indices(len(files_with_dx))

    fig, axs = plt.subplots(2, 2, figsize=(11, 9), constrained_layout=True)
    axs = axs.ravel()

    for ax_i, idx in enumerate(chosen):
        dx, path = files_with_dx[idx]
        u1, u2 = read_rank_pair_csv(path, args.max_points)
        if u1.size < 5:
            continue

        z1 = norminv_approx(np.clip(u1, 1e-8, 1 - 1e-8))
        z2 = norminv_approx(np.clip(u2, 1e-8, 1 - 1e-8))
        rho = float(np.corrcoef(z1, z2)[0, 1])

        ax = axs[ax_i]
        ax.scatter(u1, u2, s=3, c="red", alpha=args.point_alpha, edgecolors="none")

        grid = np.linspace(0.01, 0.99, 140)
        U1, U2 = np.meshgrid(grid, grid)
        C = gaussian_copula_density(U1, U2, rho)
        ax.contour(U1, U2, C, levels=10, colors="black", linewidths=0.8, alpha=0.55)

        ax.set_xlim(0.0, 1.0)
        ax.set_ylim(0.0, 1.0)
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlabel(r"$u_1$")
        ax.set_ylabel(r"$u_2$")
        ax.set_title(f"Δx={dx:.4g}, ρ={rho:.3f}, n={u1.size}")

    fig.suptitle("Sampled velocity-rank pairs with fitted Gaussian-copula contours", fontsize=13)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=180)
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
