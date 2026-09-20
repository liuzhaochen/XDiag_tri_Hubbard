#!/usr/bin/env python3
"""Plot the validated long-form spectrum table produced by summarize_scan.jl."""

import argparse
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/xdiag_tri_hubbard_matplotlib")

import matplotlib.pyplot as plt
import numpy as np


COLORS = {
    (0, 0, 1): "#0072B2",
    (0, 0, -1): "#56B4E9",
    (0, 1, 0): "#009E73",
    (0, 2, 1): "#D55E00",
    (0, 2, -1): "#E69F00",
    (1, 0, 0): "#CC79A7",
    (1, 1, 0): "#000000",
    (1, 2, 0): "#7A3E9D",
    (1, 3, 0): "#777777",
}


def sector_label(key):
    nx, ny, parity = key
    suffix = "" if parity == 0 else f", C2={'+' if parity > 0 else '-'}"
    return f"({nx},{ny}{suffix})"


def sector_keys(data):
    return sorted(
        {(int(row["nx"]), int(row["ny"]), int(row["C2_parity"])) for row in data},
        key=lambda key: (key[0], key[1], -key[2]),
    )


def select(data, key, level):
    nx, ny, parity = key
    mask = (
        (data["nx"] == nx)
        & (data["ny"] == ny)
        & (data["C2_parity"] == parity)
        & (data["level_index"] == level)
    )
    rows = data[mask]
    return rows[np.argsort(rows["U_over_t"])]


def plot_panel(data, keys, output_base, title, all_levels=True):
    figure, axis = plt.subplots(figsize=(7.2, 4.8), constrained_layout=True)
    levels = sorted({int(value) for value in data["level_index"]}) if all_levels else [1]
    for key in keys:
        for level in levels:
            rows = select(data, key, level)
            if len(rows) == 0:
                continue
            axis.plot(
                rows["U_over_t"],
                rows["energy_minus_ground"],
                color=COLORS[key],
                linewidth=1.55 if level == 1 else 0.9,
                alpha=1.0 if level == 1 else 0.58,
                marker="o",
                markersize=2.8 if level == 1 else 1.8,
                label=sector_label(key) if level == 1 else None,
            )
    axis.set_xlabel(r"$U/t$")
    axis.set_ylabel(r"$(E-E_0)/t$")
    axis.set_title(title)
    axis.grid(alpha=0.22, linewidth=0.6)
    axis.legend(fontsize=8, ncol=2, frameon=False)
    figure.savefig(output_base.with_suffix(".pdf"))
    figure.savefig(output_base.with_suffix(".png"), dpi=220)
    plt.close(figure)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("table", type=Path)
    parser.add_argument("output_directory", type=Path)
    args = parser.parse_args()
    args.output_directory.mkdir(parents=True, exist_ok=True)
    data = np.genfromtxt(args.table, delimiter=",", names=True)
    keys = sector_keys(data)

    plot_panel(
        data,
        keys,
        args.output_directory / "complete_inequivalent_spectrum",
        "3x4 triangular Hubbard: inequivalent spectrum",
    )
    plot_panel(
        data,
        [(0, 0, 1), (0, 0, -1), (1, 0, 0)],
        args.output_directory / "Ky_0_spectrum",
        r"$K_y=0$",
    )
    plot_panel(
        data,
        [(0, 2, 1), (0, 2, -1), (1, 2, 0)],
        args.output_directory / "Ky_pi_spectrum",
        r"$K_y=\pi$",
    )
    plot_panel(
        data,
        [
            (0, 0, 1),
            (0, 0, -1),
            (1, 0, 0),
            (0, 2, 1),
            (0, 2, -1),
            (1, 2, 0),
        ],
        args.output_directory / "Ky_0_pi_lowest_branches",
        r"Lowest branches on $K_y=0,\pi$",
        all_levels=False,
    )


if __name__ == "__main__":
    main()
