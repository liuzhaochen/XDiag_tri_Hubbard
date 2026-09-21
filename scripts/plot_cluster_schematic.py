#!/usr/bin/env python3
"""Plot the real- and reciprocal-space geometry of the 3x4 torus."""

import argparse
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/xdiag_tri_hubbard_matplotlib")

import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects
import numpy as np

from plot_style import apply_plot_style


apply_plot_style()


A1 = np.array([1.0, 0.0])
A2 = np.array([0.5, np.sqrt(3) / 2])
B1 = np.array([2 * np.pi, -2 * np.pi / np.sqrt(3)])
B2 = np.array([0.0, 4 * np.pi / np.sqrt(3)])
BASIS = np.column_stack([B1, B2])


def site_position(x, y):
    return x * A1 + y * A2


def folded_momentum(nx, ny):
    candidates = []
    for m in range(-2, 3):
        for n in range(-2, 3):
            cartesian = BASIS @ np.array([nx / 3 + m, ny / 4 + n])
            candidates.append((np.dot(cartesian, cartesian), m, n, cartesian))
    return min(candidates, key=lambda item: (round(item[0], 12), item[1], item[2]))[3]


def displayed_momentum(nx, ny):
    """Choose display representatives that keep the lower BZ region clear."""
    momentum = folded_momentum(nx, ny)
    if (nx, ny) == (0, 2):
        momentum = momentum + B2
    elif (nx, ny) == (2, 2):
        # Show this reciprocal image to the right of the first BZ, away from
        # the compact legend above it.
        momentum = momentum + B1 + B2
    return momentum


def plot_real_space(axis):
    directions = ((1, 0), (0, 1), (-1, 1))

    # One shell of periodic images makes the boundary-crossing bonds explicit.
    extended_sites = {
        (x, y) for x in range(-1, 4) for y in range(-1, 5)
    }
    primary_sites = {
        (x, y) for x in range(3) for y in range(4)
    }
    for x, y in sorted(extended_sites):
        start = site_position(x, y)
        for dx, dy in directions:
            destination = (x + dx, y + dy)
            if destination in extended_sites:
                end = site_position(*destination)
                axis.plot([start[0], end[0]], [start[1], end[1]],
                          color="#D8D8D8", linewidth=0.75, zorder=0)

    image_positions = np.asarray([
        site_position(x, y)
        for x, y in sorted(extended_sites - primary_sites)
    ])
    axis.scatter(image_positions[:, 0], image_positions[:, 1], s=20,
                 facecolors="white", edgecolors="#B8B8B8", linewidths=0.7,
                 zorder=1)

    for x in range(3):
        for y in range(4):
            start = site_position(x, y)
            for dx, dy in directions:
                destination = (x + dx, y + dy)
                if 0 <= destination[0] < 3 and 0 <= destination[1] < 4:
                    end = site_position(*destination)
                    axis.plot([start[0], end[0]], [start[1], end[1]],
                              color="#8A8A8A", linewidth=1.25, zorder=1)

    positions = []
    for x in range(3):
        for y in range(4):
            position = site_position(x, y)
            positions.append(position)
            site = 1 + y + 4 * x
            if x == 0:
                label_offset = np.array([0.0, -0.16])
                horizontal_alignment = "center"
                vertical_alignment = "top"
            else:
                label_offset = np.array([0.055, 0.07])
                horizontal_alignment = "left"
                vertical_alignment = "baseline"
            axis.text(*(position + label_offset), str(site),
                      fontsize=10, color="#222222", zorder=3,
                      ha=horizontal_alignment, va=vertical_alignment)
    positions = np.asarray(positions)
    axis.scatter(positions[:, 0], positions[:, 1], s=42, color="#0072B2",
                 edgecolors="white", linewidths=0.7, zorder=2)

    corners = np.array([
        [0.0, 0.0], 3 * A1, 3 * A1 + 4 * A2, 4 * A2, [0.0, 0.0]
    ])
    axis.plot(corners[:, 0], corners[:, 1], color="#292929", linewidth=1.8,
              zorder=3)
    axis.annotate("", xy=3 * A1, xytext=np.zeros(2),
                  arrowprops=dict(arrowstyle="->", color="#D55E00", lw=1.5))
    axis.annotate("", xy=4 * A2, xytext=np.zeros(2),
                  arrowprops=dict(arrowstyle="->", color="#D55E00", lw=1.5))
    t1_label = axis.text(
        1.50, -0.42, r"$\mathbf{T}_1=3\mathbf{a}_1$",
        color="#D55E00", ha="center", va="top", fontsize=11.5, zorder=5,
    )
    t1_label.set_path_effects([
        path_effects.withStroke(linewidth=2.8, foreground="white")
    ])
    axis.text(0.48, 1.88, r"$\mathbf{T}_2=4\mathbf{a}_2$", color="#D55E00",
              rotation=60, ha="center", fontsize=11.5)
    axis.set_aspect("equal")
    axis.set_xlim(-1.8, 5.75)
    axis.set_ylim(-1.15, 4.35)
    axis.axis("off")
    axis.text(0.02, 0.98, "(a)", transform=axis.transAxes,
              ha="left", va="top", fontsize=14, fontweight="bold")


def plot_reciprocal_space(axis):
    vertices = np.array([
        [4 * np.pi / 3, 0.0],
        [2 * np.pi / 3, 2 * np.pi / np.sqrt(3)],
        [-2 * np.pi / 3, 2 * np.pi / np.sqrt(3)],
        [-4 * np.pi / 3, 0.0],
        [-2 * np.pi / 3, -2 * np.pi / np.sqrt(3)],
        [2 * np.pi / 3, -2 * np.pi / np.sqrt(3)],
        [4 * np.pi / 3, 0.0],
    ])
    axis.plot(vertices[:, 0], vertices[:, 1], color="#333333", linewidth=1.6)

    momenta = []
    labels = []
    for nx in range(3):
        for ny in range(4):
            momenta.append(displayed_momentum(nx, ny))
            labels.append(f"({nx},{ny})")
    momenta = np.asarray(momenta)
    axis.scatter(momenta[:, 0], momenta[:, 1], s=47, color="#0072B2",
                 edgecolors="white", linewidths=0.7, zorder=3)

    calculated_sectors = {
        (0, 0), (0, 1), (0, 2),
        (1, 0), (1, 1), (1, 2), (1, 3),
    }
    calculated_momenta = np.asarray([
        displayed_momentum(nx, ny) for nx, ny in sorted(calculated_sectors)
    ])
    axis.scatter(calculated_momenta[:, 0], calculated_momenta[:, 1],
                 marker="o", s=92, facecolors="none", edgecolors="#E69F00",
                 linewidths=1.45, zorder=3.5, label="ED representatives")

    for position, label in zip(momenta, labels):
        if label == "(2,2)":
            offset_x, offset_y = 0.10, 0.0
            alignment, vertical_alignment = "left", "center"
        else:
            offset_x = 0.08 if position[0] <= 0 else -0.08
            alignment = "left" if position[0] <= 0 else "right"
            offset_y = 0.09 if position[1] <= 0 else -0.09
            vertical_alignment = "bottom" if position[1] <= 0 else "top"
        axis.text(position[0] + offset_x, position[1] + offset_y, label,
                  fontsize=9.2, ha=alignment,
                  va=vertical_alignment)

    # Use the BZ images adjacent to the folded representatives highlighted below.
    exact_k = np.array([-2 * np.pi / 3, 2 * np.pi / np.sqrt(3)])
    exact_m = displayed_momentum(0, 2)
    k_cluster = folded_momentum(2, 1)
    axis.scatter(*exact_k, marker="*", s=135, color="#D55E00", zorder=4,
                 label=r"$K$")
    axis.scatter(*exact_m, marker="s", s=58, facecolors="none",
                 edgecolors="#009E73", linewidths=1.6, zorder=4,
                 label=r"$M_2$")
    axis.scatter(*k_cluster, marker="o", s=125, facecolors="none",
                 edgecolors="#CC79A7", linewidths=2.0, zorder=4,
                 label=r"$K_{\rm cl}$")

    axis.set_aspect("equal")
    axis.set_xlim(-4.65, 4.85)
    axis.set_ylim(-4.1, 5.15)
    axis.set_xticks([])
    axis.set_yticks([])
    for spine in axis.spines.values():
        spine.set_visible(False)
    axis.legend(loc="upper right", ncol=2, frameon=False, fontsize=9.4,
                handletextpad=0.4, labelspacing=0.25, columnspacing=0.8,
                borderaxespad=0.15)
    axis.text(0.045, 0.94, "(b)", transform=axis.transAxes,
              ha="left", va="top", fontsize=14, fontweight="bold")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("output_directory", type=Path, nargs="?", default=Path("plots/N12"))
    args = parser.parse_args()
    args.output_directory.mkdir(parents=True, exist_ok=True)

    # Render at the intended full-width Supplementary-page size so that text
    # remains at its designed point size when included with width=\textwidth.
    figure, axes = plt.subplots(1, 2, figsize=(7.1, 3.35), constrained_layout=True)
    plot_real_space(axes[0])
    plot_reciprocal_space(axes[1])
    for extension in ("pdf", "png"):
        kwargs = {"dpi": 300} if extension == "png" else {}
        kwargs["bbox_inches"] = "tight"
        kwargs["pad_inches"] = 0.04
        figure.savefig(args.output_directory / f"cluster_schematic.{extension}", **kwargs)
    plt.close(figure)


if __name__ == "__main__":
    main()
