#!/usr/bin/env python3
"""Plot the two refined ground-branch crossings."""

import argparse
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/xdiag_tri_hubbard_matplotlib")

import matplotlib.pyplot as plt
import numpy as np


def lowest_branch(data, nx, ny, parity):
    mask = (
        (data["nx"] == nx)
        & (data["ny"] == ny)
        & (data["C2_parity"] == parity)
        & (data["level_index"] == 1)
    )
    rows = data[mask]
    return rows[np.argsort(rows["U_over_t"])]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("table", type=Path)
    parser.add_argument("output_base", type=Path)
    args = parser.parse_args()

    data = np.genfromtxt(args.table, delimiter=",", names=True)
    gamma = lowest_branch(data, 0, 0, 1)
    corner = lowest_branch(data, 1, 0, 0)
    if not np.array_equal(gamma["U_over_t"], corner["U_over_t"]):
        raise ValueError("refined branches do not share the same U/t grid")

    u_values = gamma["U_over_t"]
    delta = corner["energy"] - gamma["energy"]
    gaps = np.diff(u_values)
    split = int(np.argmax(gaps)) + 1 if len(gaps) and np.max(gaps) > 1.5 * np.min(gaps) else len(u_values)
    windows = [slice(0, split), slice(split, len(u_values))]

    figure, axes = plt.subplots(1, 2, figsize=(8.4, 3.7), constrained_layout=True)
    for axis, window in zip(axes, windows):
        axis.axhline(0.0, color="#555555", linewidth=0.9)
        axis.plot(u_values[window], delta[window], "o-", color="#7A3E9D", linewidth=1.5)
        axis.set_xlabel(r"$U/t$")
        axis.grid(alpha=0.22, linewidth=0.6)
    axes[0].set_ylabel(r"$[E_{(1,0)}-E_{(0,0),C_2=+}]/t$")
    figure.suptitle("3x4 triangular Hubbard: refined ground-branch crossings")
    args.output_base.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(args.output_base.with_suffix(".pdf"))
    figure.savefig(args.output_base.with_suffix(".png"), dpi=220)
    plt.close(figure)


if __name__ == "__main__":
    main()
