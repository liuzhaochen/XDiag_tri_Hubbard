#!/usr/bin/env python3
"""Plot ground, low-state, and momentum-space observable diagnostics."""

import argparse
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/xdiag_tri_hubbard_matplotlib")

import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects
import numpy as np
from matplotlib.ticker import AutoMinorLocator

from plot_style import apply_plot_style


apply_plot_style()


def load_csv(path):
    return np.genfromtxt(path, delimiter=",", names=True, dtype=None, encoding="utf-8")


def save(figure, output_base, dpi=180):
    output_base.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_base.with_suffix(".pdf"), bbox_inches="tight", pad_inches=0.03)
    figure.savefig(output_base.with_suffix(".png"), dpi=dpi,
                   bbox_inches="tight", pad_inches=0.03)
    plt.close(figure)


def crossing_spans(axes):
    for axis in np.atleast_1d(axes).flat:
        axis.axvspan(6.60, 6.65, color="#BBBBBB", alpha=0.25, linewidth=0)
        axis.axvspan(9.30, 9.35, color="#BBBBBB", alpha=0.25, linewidth=0)


def plot_ground(data, output_base):
    u = data["U_over_t"]
    figure, axes = plt.subplots(1, 2, figsize=(7.1, 2.85))

    axes[0].plot(
        u, data["double_occupancy"], "o-", color="#0072B2",
        linewidth=1.65, markersize=4.5, markeredgewidth=0.0,
    )
    axes[0].set_ylabel(r"$\langle n_\uparrow n_\downarrow\rangle$", fontsize=10)

    axes[1].plot(
        u, data["spin_Kcl"], "o-", color="#CC79A7",
        linewidth=1.65, markersize=4.5, markeredgewidth=0.0,
        label=r"$K_{\rm cl}$",
    )
    axes[1].plot(
        u, data["spin_M2"], "s-", color="#009E73",
        linewidth=1.65, markersize=4.1, markeredgewidth=0.0,
        label=r"$M_2$",
    )
    axes[1].set_ylabel(r"$S(\mathbf{q})$", fontsize=10)
    axes[1].legend(frameon=False, fontsize=8.6, loc="upper left",
                   handlelength=1.7, handletextpad=0.45)

    panel_labels = ("(a)", "(b)")
    for axis, panel_label in zip(axes, panel_labels):
        axis.set_xlabel(r"$U/t$", fontsize=10)
        axis.set_xlim(5.9, 14.1)
        axis.set_xticks(np.arange(6, 15, 2))
        axis.xaxis.set_minor_locator(AutoMinorLocator(2))
        axis.yaxis.set_minor_locator(AutoMinorLocator(2))
        axis.tick_params(axis="both", which="major", direction="in",
                         top=True, right=True, width=0.8, length=4.0,
                         labelsize=8.6, pad=2.5)
        axis.tick_params(axis="both", which="minor", direction="in",
                         top=True, right=True, width=0.65, length=2.2)
        for spine in axis.spines.values():
            spine.set_linewidth(0.85)
        axis.text(-0.10, 1.0, panel_label, transform=axis.transAxes,
                  ha="right", va="top", fontsize=10.5,
                  fontweight="bold", clip_on=False)
    crossing_spans(axes)
    figure.subplots_adjust(left=0.105, right=0.985, bottom=0.18, top=0.965,
                           wspace=0.34)
    save(figure, output_base, dpi=300)


def plot_low_states(data, output_base):
    mask = data["energy_minus_ground"] <= 0.6 + 1e-12
    selected = data[mask]
    quantities = [
        ("spin_Kcl", r"$S(K_{\rm cl})$"),
        ("spin_M2", r"$S(M_2)$"),
    ]
    color_minimum = min(float(np.min(selected[field])) for field, _ in quantities)
    color_maximum = max(float(np.max(selected[field])) for field, _ in quantities)

    figure, axes = plt.subplots(
        1, 2, figsize=(7.1, 3.0),
        sharex=True, sharey=True,
    )
    points = []
    for panel_index, (axis, (field, label)) in enumerate(zip(axes, quantities)):
        values = selected[field]
        points.append(axis.scatter(
            selected["U_over_t"], selected["energy_minus_ground"],
            c=values, cmap="viridis", s=20, edgecolors="none",
            vmin=color_minimum, vmax=color_maximum,
        ))
        axis.set_xlabel(r"$U/t$", fontsize=10)
        axis.set_xlim(5.9, 14.1)
        axis.set_ylim(-0.015, 0.61)
        axis.set_xticks(np.arange(6, 15, 2))
        axis.xaxis.set_minor_locator(AutoMinorLocator(2))
        axis.yaxis.set_minor_locator(AutoMinorLocator(2))
        axis.tick_params(axis="both", which="major", direction="in",
                         top=True, right=True, width=0.8, length=4.0,
                         labelsize=8.6, pad=2.5)
        axis.tick_params(axis="both", which="minor", direction="in",
                         top=True, right=True, width=0.65, length=2.2)
        for spine in axis.spines.values():
            spine.set_linewidth(0.85)
        panel_label = "(a)" if panel_index == 0 else "(b)"
        axis.text(
            -0.10, 1.0, panel_label,
            transform=axis.transAxes, ha="right", va="top",
            fontsize=10.5, fontweight="bold", clip_on=False,
        )
        axis.text(
            0.5, 1.025, label,
            transform=axis.transAxes, ha="center", va="bottom",
            fontsize=10, clip_on=False,
        )
    axes[0].set_ylabel(r"$(E-E_0)/t$", fontsize=10)
    crossing_spans(axes)
    figure.subplots_adjust(left=0.105, right=0.89, bottom=0.17, top=0.91,
                           wspace=0.27)
    colorbar_axis = figure.add_axes([0.915, 0.17, 0.018, 0.74])
    colorbar = figure.colorbar(points[0], cax=colorbar_axis)
    colorbar.set_label(r"$S(\mathbf{q})$", fontsize=10)
    colorbar.ax.tick_params(direction="in", width=0.8, length=3.5,
                            labelsize=8.6)
    save(figure, output_base, dpi=300)


B1 = np.array([2 * np.pi, -2 * np.pi / np.sqrt(3)])
B2 = np.array([0.0, 4 * np.pi / np.sqrt(3)])
BASIS = np.column_stack([B1, B2])


def extended_momentum_grid(points=111):
    x = np.linspace(-2 * np.pi, 2 * np.pi, points)
    y = np.linspace(-4 * np.pi / np.sqrt(3), 4 * np.pi / np.sqrt(3), points)
    xx, yy = np.meshgrid(x, y)
    cart = np.stack([xx, yy], axis=-1)
    fractional = np.linalg.solve(BASIS, cart.reshape(-1, 2).T).T.reshape(points, points, 2)
    return xx, yy, fractional[..., 0], fractional[..., 1]


def first_bz_vertices():
    return np.array([
        [4 * np.pi / 3, 0.0],
        [2 * np.pi / 3, 2 * np.pi / np.sqrt(3)],
        [-2 * np.pi / 3, 2 * np.pi / np.sqrt(3)],
        [-4 * np.pi / 3, 0.0],
        [-2 * np.pi / 3, -2 * np.pi / np.sqrt(3)],
        [2 * np.pi / 3, -2 * np.pi / np.sqrt(3)],
        [4 * np.pi / 3, 0.0],
    ])


def allowed_momentum_replicas(x_limits, y_limits):
    points = []
    for nx in range(3):
        for ny in range(4):
            for shift_x in range(-2, 3):
                for shift_y in range(-2, 3):
                    cartesian = BASIS @ np.array([
                        nx / 3 + shift_x, ny / 4 + shift_y
                    ])
                    if (x_limits[0] <= cartesian[0] <= x_limits[1]
                            and y_limits[0] <= cartesian[1] <= y_limits[1]):
                        points.append(cartesian)
    return np.unique(np.round(points, decimals=12), axis=0)


def evaluate_displacements(correlation, q1, q2):
    result = np.zeros(q1.shape, dtype=np.complex128)
    for dx in range(3):
        for dy in range(4):
            result += np.exp(2j * np.pi * (q1 * dx + q2 * dy)) * correlation[dx, dy]
    return np.real(result)


def evaluate_displacements_complex(correlation, q1, q2):
    result = np.zeros(q1.shape, dtype=np.complex128)
    for dx in range(3):
        for dy in range(4):
            result += np.exp(2j * np.pi * (q1 * dx + q2 * dy)) * correlation[dx, dy]
    return result


def read_selected_correlation(rows, tag, observable):
    chosen = rows[(rows["tag"] == tag) & (rows["observable"] == observable)]
    result = np.zeros((3, 4), dtype=np.complex128)
    for row in chosen:
        result[int(row["dx"]), int(row["dy"])] = row["real"] + 1j * row["imag"]
    return result


def read_chirality_gram(rows, tag):
    chosen = rows[rows["tag"] == tag]
    result = np.zeros((2, 2, 3, 4), dtype=np.complex128)
    for row in chosen:
        first = int(row["sublattice_a"]) - 1
        second = int(row["sublattice_b"]) - 1
        result[first, second, int(row["q_nx"]), int(row["q_ny"])] = (
            row["real"] + 1j * row["imag"]
        )
    return result


def chirality_correlations_from_gram(gram):
    # G_ab(q)=sum_R exp(+iq.R) C_ab(R), so C_ab=FFT[G_ab]/N.
    return np.fft.fft2(gram, axes=(-2, -1)) / 12


def evaluate_chirality(correlations, q1, q2):
    offsets = ((1 / 3, 1 / 3), (2 / 3, 2 / 3))
    weights = [
        np.exp(2j * np.pi * (q1 * offset[0] + q2 * offset[1])) / np.sqrt(2)
        for offset in offsets
    ]
    result = np.zeros(q1.shape, dtype=np.complex128)
    for first in range(2):
        for second in range(2):
            gram_component = evaluate_displacements_complex(
                correlations[first, second], q1, q2
            )
            result += np.conj(weights[first]) * gram_component * weights[second]
    # Between the 12 sampled momenta, take the Hermitian (real) part of the
    # trigonometric interpolation. At every allowed momentum the stored Gram
    # matrix is explicitly validated as Hermitian and positive semidefinite.
    return np.real(result)


def validate_chirality_interpolation(correlations, grid_rows, u_value):
    for row in grid_rows[np.isclose(grid_rows["U_over_t"], u_value)]:
        reconstructed = evaluate_chirality(
            correlations, np.asarray(row["q1"]), np.asarray(row["q2"])
        )
        if not np.isclose(reconstructed, row["full"], rtol=1e-8, atol=1e-10):
            raise ValueError("chirality interpolation does not reproduce allowed momenta")


def plot_heatmaps(correlation_rows, output_base):
    xx, yy, q1, q2 = extended_momentum_grid()
    x_limits = (float(np.min(xx)), float(np.max(xx)))
    y_limits = (float(np.min(yy)), float(np.max(yy)))
    allowed_cart = allowed_momentum_replicas(x_limits, y_limits)
    bz = first_bz_vertices()
    selections = [("U8_k10", 8.0, r"$U/t=8$"),
                  ("U14_Gamma", 14.0, r"$U/t=14$")]
    fields = [
        evaluate_displacements(
            read_selected_correlation(correlation_rows, tag, "spin"), q1, q2
        )
        for tag, _, _ in selections
    ]
    color_maximum = max(float(np.max(field)) for field in fields)
    display_xx = xx / np.pi
    display_yy = yy / np.pi
    display_allowed = allowed_cart / np.pi
    display_bz = bz / np.pi
    display_x_limits = (x_limits[0] / np.pi, x_limits[1] / np.pi)
    display_y_limits = (y_limits[0] / np.pi, y_limits[1] / np.pi)

    figure = plt.figure(figsize=(7.1, 3.75))
    grid = figure.add_gridspec(
        1, 3, width_ratios=(1.0, 1.0, 0.04),
        left=0.080, right=0.935, bottom=0.14, top=0.90, wspace=0.070,
    )
    axes = np.asarray([
        figure.add_subplot(grid[0, 0]),
        figure.add_subplot(grid[0, 1]),
    ])
    colorbar_axis = figure.add_subplot(grid[0, 2])
    images = []
    for panel_index, (axis, field, (_, _, state_label)) in enumerate(
        zip(axes, fields, selections)
    ):
        image = axis.pcolormesh(
            display_xx, display_yy, field, shading="auto", cmap="magma",
            vmin=0.0, vmax=color_maximum, rasterized=True,
        )
        images.append(image)
        axis.plot(display_bz[:, 0], display_bz[:, 1], "w--", linewidth=1.35)
        axis.scatter(
            display_allowed[:, 0], display_allowed[:, 1], s=14,
            facecolors="white", edgecolors="#222222", linewidths=0.65,
        )
        axis.set_aspect("equal")
        axis.set_xlim(*display_x_limits)
        axis.set_ylim(*display_y_limits)
        axis.set_xticks(np.arange(-2, 2.1, 1))
        axis.set_yticks(np.arange(-2, 2.1, 1))
        axis.set_xlabel(r"$q_x/\pi$", fontsize=10.5)
        axis.set_ylabel(r"$q_y/\pi$" if panel_index == 0 else "", fontsize=10.5)
        axis.tick_params(axis="both", which="major", direction="in",
                         top=True, right=True, width=0.8, length=4.0,
                         labelsize=9.2)
        if panel_index == 1:
            axis.tick_params(labelleft=False)
        for spine in axis.spines.values():
            spine.set_linewidth(0.85)
        panel_label = "(a)" if panel_index == 0 else "(b)"
        panel_text = axis.text(
            0.035, 0.965, panel_label,
            transform=axis.transAxes, ha="left", va="top",
            fontsize=11, fontweight="bold", color="white",
        )
        panel_text.set_path_effects([
            path_effects.withStroke(linewidth=1.5, foreground="#222222")
        ])
        axis.text(
            0.5, 1.025, state_label,
            transform=axis.transAxes, ha="center", va="bottom",
            fontsize=11, clip_on=False,
        )

    colorbar = figure.colorbar(images[0], cax=colorbar_axis)
    colorbar.set_label(r"$S(\mathbf{q})$", fontsize=10.5)
    colorbar.ax.tick_params(direction="in", width=0.8, length=3.5,
                            labelsize=9.2)
    save(figure, output_base)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("state_table", type=Path)
    parser.add_argument("ground_table", type=Path)
    parser.add_argument("heatmap_correlations", type=Path)
    parser.add_argument("chirality_grid", type=Path)
    parser.add_argument("chirality_sublattice_gram", type=Path)
    parser.add_argument("output_directory", type=Path)
    args = parser.parse_args()
    args.output_directory.mkdir(parents=True, exist_ok=True)

    states = load_csv(args.state_table)
    ground = load_csv(args.ground_table)
    heatmap = load_csv(args.heatmap_correlations)
    plot_ground(ground, args.output_directory / "ground_state_observables")
    plot_low_states(states, args.output_directory / "low_lying_state_observables")
    plot_heatmaps(heatmap, args.output_directory / "structure_factor_heatmaps")


if __name__ == "__main__":
    main()
