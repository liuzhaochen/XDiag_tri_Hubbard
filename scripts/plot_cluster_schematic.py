#!/usr/bin/env python3
"""Plot the real- and reciprocal-space geometry of the 3x4 torus."""

import argparse
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/xdiag_tri_hubbard_matplotlib")

import matplotlib.pyplot as plt
import numpy as np


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


def plot_real_space(axis):
    directions = ((1, 0), (0, 1), (-1, 1))
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
            axis.text(position[0] + 0.055, position[1] + 0.07, str(site),
                      fontsize=8, color="#222222", zorder=3)
    positions = np.asarray(positions)
    axis.scatter(positions[:, 0], positions[:, 1], s=42, color="#0072B2",
                 edgecolors="white", linewidths=0.7, zorder=2)

    corners = np.array([
        [0.0, 0.0], 3 * A1, 3 * A1 + 4 * A2, 4 * A2, [0.0, 0.0]
    ])
    axis.plot(corners[:, 0], corners[:, 1], color="#D55E00", linewidth=2.0)
    axis.annotate("", xy=3 * A1, xytext=np.zeros(2),
                  arrowprops=dict(arrowstyle="->", color="#D55E00", lw=1.5))
    axis.annotate("", xy=4 * A2, xytext=np.zeros(2),
                  arrowprops=dict(arrowstyle="->", color="#D55E00", lw=1.5))
    axis.text(1.45, -0.25, r"$3\mathbf{a}_1$", color="#D55E00", ha="center")
    axis.text(0.62, 1.88, r"$4\mathbf{a}_2$", color="#D55E00",
              rotation=60, ha="center")
    axis.text(2.5, 3.6, "opposite edges identified", fontsize=8.5,
              color="#555555", ha="center")
    axis.set_aspect("equal")
    axis.set_xlim(-0.35, 5.35)
    axis.set_ylim(-0.42, 3.92)
    axis.set_title(r"Real space: $3\times4$ triangular torus")
    axis.set_xlabel(r"$x$")
    axis.set_ylabel(r"$y$")


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
            momenta.append(folded_momentum(nx, ny))
            labels.append(f"({nx},{ny})")
    momenta = np.asarray(momenta)
    axis.scatter(momenta[:, 0], momenta[:, 1], s=47, color="#0072B2",
                 edgecolors="white", linewidths=0.7, zorder=3,
                 label="allowed cluster momenta")
    for position, label in zip(momenta, labels):
        offset_x = 0.08 if position[0] <= 0 else -0.08
        alignment = "left" if position[0] <= 0 else "right"
        offset_y = 0.09 if position[1] <= 0 else -0.09
        axis.text(position[0] + offset_x, position[1] + offset_y, label,
                  fontsize=7.3, ha=alignment,
                  va="bottom" if position[1] <= 0 else "top")

    # Use the BZ images adjacent to the folded representatives highlighted below.
    exact_k = np.array([-2 * np.pi / 3, 2 * np.pi / np.sqrt(3)])
    exact_m = folded_momentum(0, 2)
    k_cluster = folded_momentum(2, 1)
    axis.scatter(*exact_k, marker="*", s=135, color="#D55E00", zorder=4,
                 label=r"exact $K$")
    axis.scatter(*exact_m, marker="s", s=58, facecolors="none",
                 edgecolors="#009E73", linewidths=1.6, zorder=4,
                 label=r"exact allowed $M_2$")
    axis.scatter(*k_cluster, marker="o", s=125, facecolors="none",
                 edgecolors="#CC79A7", linewidths=2.0, zorder=4,
                 label=r"selected $K_{\rm cl}=(2,1)$")

    axis.set_aspect("equal")
    axis.set_xlim(-4.65, 4.65)
    axis.set_ylim(-4.1, 4.1)
    axis.set_title("Reciprocal space: first Brillouin zone")
    axis.set_xlabel(r"$q_x$")
    axis.set_ylabel(r"$q_y$")
    axis.legend(loc="upper right", frameon=False, fontsize=7.6)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("output_directory", type=Path, nargs="?", default=Path("plots/N12"))
    args = parser.parse_args()
    args.output_directory.mkdir(parents=True, exist_ok=True)

    figure, axes = plt.subplots(1, 2, figsize=(11.0, 4.8), constrained_layout=True)
    plot_real_space(axes[0])
    plot_reciprocal_space(axes[1])
    figure.suptitle("3x4 triangular-lattice cluster and its 12 allowed momenta")
    for extension in ("pdf", "png"):
        kwargs = {"dpi": 220} if extension == "png" else {}
        figure.savefig(args.output_directory / f"cluster_schematic.{extension}", **kwargs)
    plt.close(figure)


if __name__ == "__main__":
    main()
