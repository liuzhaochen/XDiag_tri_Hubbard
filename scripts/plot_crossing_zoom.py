#!/usr/bin/env python3
"""Plot the two refined ground-branch crossings."""

import argparse
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/xdiag_tri_hubbard_matplotlib")

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import AutoMinorLocator

from plot_style import apply_plot_style


apply_plot_style()


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

    figure, axes = plt.subplots(1, 2, figsize=(7.1, 2.75))
    for panel_index, (axis, window) in enumerate(zip(axes, windows)):
        local_u = u_values[window]
        local_delta = delta[window]
        axis.axhline(0.0, color="#555555", linewidth=0.8, zorder=0)
        sign_changes = np.flatnonzero(local_delta[:-1] * local_delta[1:] <= 0.0)
        for index in sign_changes:
            axis.axvspan(local_u[index], local_u[index + 1],
                         color="#BBBBBB", alpha=0.28, linewidth=0, zorder=0)
        axis.plot(
            local_u, local_delta, "o-", color="#7A3E9D",
            linewidth=1.45, markersize=4.1, markeredgewidth=0.0,
        )
        axis.set_xlabel(r"$U/t$", fontsize=10)
        axis.xaxis.set_minor_locator(AutoMinorLocator(2))
        axis.yaxis.set_minor_locator(AutoMinorLocator(2))
        axis.tick_params(axis="both", which="major", direction="in",
                         top=True, right=True, width=0.8, length=4.0,
                         labelsize=8.6, pad=2.5)
        axis.tick_params(axis="both", which="minor", direction="in",
                         top=True, right=True, width=0.65, length=2.2)
        for spine in axis.spines.values():
            spine.set_linewidth(0.85)
        axis.text(-0.16, 1.0, "(a)" if panel_index == 0 else "(b)",
                  transform=axis.transAxes, ha="right", va="top",
                  fontsize=10.5, fontweight="bold", clip_on=False)
    axes[0].set_ylabel(r"$[E_{(1,0)}-E_{\Gamma(A)}]/t$", fontsize=10)
    figure.subplots_adjust(left=0.12, right=0.985, bottom=0.19, top=0.96,
                           wspace=0.42)
    args.output_base.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(args.output_base.with_suffix(".pdf"),
                   bbox_inches="tight", pad_inches=0.03)
    figure.savefig(args.output_base.with_suffix(".png"), dpi=300,
                   bbox_inches="tight", pad_inches=0.03)
    plt.close(figure)


if __name__ == "__main__":
    main()
