#!/usr/bin/env python3
"""Plot the validated long-form spectrum table produced by summarize_scan.jl."""

import argparse
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/xdiag_tri_hubbard_matplotlib")

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import AutoMinorLocator

from plot_style import apply_plot_style


apply_plot_style()


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

LINESTYLES = {
    (0, 0, 1): "-",
    (0, 0, -1): "--",
    (0, 1, 0): "-",
    (0, 2, 1): "-",
    (0, 2, -1): "--",
    (1, 0, 0): "-",
    (1, 1, 0): "-",
    (1, 2, 0): "-",
    (1, 3, 0): "-",
}


def sector_label(key):
    nx, ny, parity = key
    if parity != 0:
        irrep = "A" if parity > 0 else "B"
        momentum = r"\Gamma" if (nx, ny) == (0, 0) else r"M_2"
        return rf"${momentum}({irrep})$"
    return rf"$({nx},{ny})$"


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


def draw_panel(axis, data, keys, panel_label, all_levels=True):
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
                linestyle=LINESTYLES[key],
                linewidth=1.35 if level == 1 else 0.60,
                alpha=1.0 if level == 1 else 0.45,
                marker="o",
                markersize=2.25 if level == 1 else 1.10,
                markeredgewidth=0.0,
                label=sector_label(key) if level == 1 else None,
            )
    axis.set_xlim(float(np.min(data["U_over_t"])), float(np.max(data["U_over_t"])))
    axis.set_ylim(-0.015, 1.0)
    axis.set_xticks(np.arange(6, 15, 2))
    axis.set_yticks(np.arange(0.0, 1.01, 0.2))
    axis.xaxis.set_minor_locator(AutoMinorLocator(2))
    axis.yaxis.set_minor_locator(AutoMinorLocator(2))
    axis.tick_params(axis="both", which="major", direction="in",
                     top=True, right=True, width=0.8, length=4.0,
                     labelsize=8.3, pad=2.5)
    axis.tick_params(axis="both", which="minor", direction="in",
                     top=True, right=True, width=0.65, length=2.2)
    for spine in axis.spines.values():
        spine.set_linewidth(0.85)
    axis.text(0.0, 1.035, panel_label, transform=axis.transAxes,
              ha="left", va="bottom", fontsize=9.8,
              fontweight="bold", clip_on=False)


def plot_overview(data, keys, output_base):
    figure, axes = plt.subplots(
        1, 3, figsize=(7.1, 3.05), sharex=True, sharey=True,
    )
    draw_panel(axes[0], data, keys, "(a)", all_levels=False)
    draw_panel(
        axes[1], data,
        [(0, 0, 1), (0, 0, -1), (1, 0, 0)],
        r"(b)  $K_y=0$",
    )
    draw_panel(
        axes[2], data,
        [(0, 2, 1), (0, 2, -1), (1, 2, 0)],
        r"(c)  $K_y=\pi$",
    )

    axes[0].set_ylabel(r"$(E-E_0)/t$", fontsize=9.6)
    for axis in axes:
        axis.set_xlabel(r"$U/t$", fontsize=9.6)

    handles, labels = axes[0].get_legend_handles_labels()
    figure.legend(
        handles, labels, fontsize=7.8, ncol=5, frameon=False,
        loc="upper center", bbox_to_anchor=(0.54, 0.985),
        handlelength=1.45, handletextpad=0.35, columnspacing=0.85,
        labelspacing=0.35,
    )
    figure.subplots_adjust(
        left=0.085, right=0.992, bottom=0.17, top=0.73, wspace=0.10,
    )

    figure.savefig(output_base.with_suffix(".pdf"), bbox_inches="tight", pad_inches=0.03)
    figure.savefig(output_base.with_suffix(".png"), dpi=300,
                   bbox_inches="tight", pad_inches=0.03)
    plt.close(figure)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("table", type=Path)
    parser.add_argument("output_directory", type=Path)
    args = parser.parse_args()
    args.output_directory.mkdir(parents=True, exist_ok=True)
    data = np.genfromtxt(args.table, delimiter=",", names=True)
    keys = sector_keys(data)

    plot_overview(data, keys, args.output_directory / "spectrum_overview")


if __name__ == "__main__":
    main()
