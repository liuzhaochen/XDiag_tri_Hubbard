#!/usr/bin/env python3
"""Plot ground, low-state, and momentum-space observable diagnostics."""

import argparse
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/xdiag_tri_hubbard_matplotlib")

import matplotlib.pyplot as plt
import numpy as np


def load_csv(path):
    return np.genfromtxt(path, delimiter=",", names=True, dtype=None, encoding="utf-8")


def save(figure, output_base):
    output_base.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_base.with_suffix(".pdf"))
    figure.savefig(output_base.with_suffix(".png"), dpi=220)
    plt.close(figure)


def crossing_spans(axes):
    for axis in np.atleast_1d(axes).flat:
        axis.axvspan(6.60, 6.65, color="#BBBBBB", alpha=0.25, linewidth=0)
        axis.axvspan(9.30, 9.35, color="#BBBBBB", alpha=0.25, linewidth=0)


def plot_ground(data, output_base):
    u = data["U_over_t"]
    figure, axes = plt.subplots(2, 2, figsize=(10.2, 7.2), constrained_layout=True)

    axes[0, 0].plot(u, data["double_occupancy"], "o-", color="#0072B2")
    axes[0, 0].set_ylabel(r"$\langle n_\uparrow n_\downarrow\rangle$")

    spin_m = np.vstack([data["spin_M1cl"], data["spin_M2"], data["spin_M3cl"]])
    axes[0, 1].plot(u, data["spin_Gamma"], "o-", label=r"$\Gamma$")
    axes[0, 1].plot(u, data["spin_Kcl"], "o-", label=r"$K_{\rm cl}$")
    axes[0, 1].plot(u, np.mean(spin_m, axis=0), "o-", label=r"mean $M_{\rm cl}$")
    axes[0, 1].fill_between(u, np.min(spin_m, axis=0), np.max(spin_m, axis=0),
                            color="#999999", alpha=0.18, label=r"$M_{\rm cl}$ range")
    axes[0, 1].set_ylabel(r"$S(\mathbf{q})$")
    axes[0, 1].legend(frameon=False, fontsize=8)

    charge_m = np.vstack([data["charge_M1cl"], data["charge_M2"], data["charge_M3cl"]])
    axes[1, 0].plot(u, data["charge_Gamma"], "o-", label=r"$\Gamma$")
    axes[1, 0].plot(u, data["charge_Kcl"], "o-", label=r"$K_{\rm cl}$")
    axes[1, 0].plot(u, np.mean(charge_m, axis=0), "o-", label=r"mean $M_{\rm cl}$")
    axes[1, 0].fill_between(u, np.min(charge_m, axis=0), np.max(charge_m, axis=0),
                            color="#999999", alpha=0.18)
    axes[1, 0].set_ylabel(r"$N(\mathbf{q})$")
    axes[1, 0].legend(frameon=False, fontsize=8)

    axes[1, 1].plot(u, data["chirality_projected_abs"], "o-",
                    label=r"projected $|\bar\chi|$")
    axes[1, 1].plot(u, data["chirality_sf_Gamma"], "o-",
                    label=r"$C_\chi(\Gamma)$")
    axes[1, 1].plot(u, data["chirality_sf_Kcl"], "o-",
                    label=r"$C_\chi(K_{\rm cl})$")
    axes[1, 1].plot(u, data["chirality_sf_max"], "o-",
                    label=r"max allowed $C_\chi$")
    axes[1, 1].set_ylabel("chirality diagnostics")
    axes[1, 1].legend(frameon=False, fontsize=8)

    for axis in axes.flat:
        axis.set_xlabel(r"$U/t$")
        axis.grid(alpha=0.22, linewidth=0.6)
    crossing_spans(axes)
    figure.suptitle("3x4 triangular Hubbard: ground-state observables")
    save(figure, output_base)


def plot_low_states(data, output_base):
    mask = data["energy_minus_ground"] <= 0.6 + 1e-12
    selected = data[mask]
    quantities = [
        ("double_occupancy", r"$\langle n_\uparrow n_\downarrow\rangle$"),
        ("spin_Kcl", r"$S(K_{\rm cl})$"),
        ("spin_M2", r"$S(M_2)$"),
        ("chirality_mean", r"$|\bar\chi|$"),
    ]
    figure, axes = plt.subplots(2, 2, figsize=(10.4, 7.4), constrained_layout=True,
                                sharex=True, sharey=True)
    for axis, (field, label) in zip(axes.flat, quantities):
        values = np.abs(selected[field]) if field == "chirality_mean" else selected[field]
        points = axis.scatter(selected["U_over_t"], selected["energy_minus_ground"],
                              c=values, cmap="viridis", s=24, edgecolors="none")
        figure.colorbar(points, ax=axis, label=label, shrink=0.88)
        axis.set_xlabel(r"$U/t$")
        axis.set_ylabel(r"$(E-E_0)/t$")
        axis.grid(alpha=0.18, linewidth=0.5)
    crossing_spans(axes)
    figure.suptitle(r"Low-lying states colored by observables ($E-E_0\leq0.6t$)")
    save(figure, output_base)


B1 = np.array([2 * np.pi, -2 * np.pi / np.sqrt(3)])
B2 = np.array([0.0, 4 * np.pi / np.sqrt(3)])
BASIS = np.column_stack([B1, B2])


def first_bz_grid(points=241):
    x = np.linspace(-4.35, 4.35, points)
    y = np.linspace(-3.75, 3.75, points)
    xx, yy = np.meshgrid(x, y)
    cart = np.stack([xx, yy], axis=-1)
    norm0 = np.sum(cart**2, axis=-1)
    mask = np.ones(xx.shape, dtype=bool)
    for m in range(-1, 2):
        for n in range(-1, 2):
            if m == 0 and n == 0:
                continue
            reciprocal = m * B1 + n * B2
            mask &= norm0 <= np.sum((cart - reciprocal) ** 2, axis=-1) + 1e-10
    fractional = np.linalg.solve(BASIS, cart.reshape(-1, 2).T).T.reshape(points, points, 2)
    return xx, yy, fractional[..., 0], fractional[..., 1], mask


def fold_to_first_bz(q1, q2):
    best = None
    for m in range(-1, 2):
        for n in range(-1, 2):
            cart = BASIS @ np.array([q1 + m, q2 + n])
            if best is None or np.dot(cart, cart) < np.dot(best, best):
                best = cart
    return best


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


def plot_heatmaps(correlation_rows, chirality_rows, chirality_gram_rows, output_base):
    xx, yy, q1, q2, mask = first_bz_grid()
    allowed_cart = np.array([
        fold_to_first_bz(nx / 3, ny / 4) for nx in range(3) for ny in range(4)
    ])
    selections = [("U8_k10", 8.0, r"$U/t=8$, $(1,0)$"),
                  ("U14_Gamma", 14.0, r"$U/t=14$, $\Gamma,C_2=+$")]
    figure, axes = plt.subplots(2, 3, figsize=(15.2, 8.2), constrained_layout=True)
    figure.suptitle(
        "Continuous Fourier visualization; dots are the 12 allowed momenta",
        fontsize=14,
    )
    for row_index, (tag, u_value, state_label) in enumerate(selections):
        spin = read_selected_correlation(correlation_rows, tag, "spin")
        charge = read_selected_correlation(correlation_rows, tag, "charge")
        chi_gram = read_chirality_gram(chirality_gram_rows, tag)
        chi_correlations = chirality_correlations_from_gram(chi_gram)
        validate_chirality_interpolation(chi_correlations, chirality_rows, u_value)
        fields = [
            (evaluate_displacements(spin, q1, q2), r"$S(\mathbf{q})$"),
            (evaluate_displacements(charge, q1, q2), r"$N(\mathbf{q})$"),
            (evaluate_chirality(chi_correlations, q1, q2), r"$C_\chi(\mathbf{q})$"),
        ]
        for column, (field, title) in enumerate(fields):
            axis = axes[row_index, column]
            shown = np.where(mask, field, np.nan)
            image = axis.pcolormesh(xx, yy, shown, shading="auto", cmap="magma")
            axis.scatter(allowed_cart[:, 0], allowed_cart[:, 1], s=17,
                         facecolors="white", edgecolors="black", linewidths=0.6)
            axis.set_aspect("equal")
            axis.set_xlabel(r"$q_x$")
            axis.set_ylabel(r"$q_y$")
            axis.set_title(f"{state_label}: {title}", fontsize=11)
            figure.colorbar(image, ax=axis, shrink=0.82)
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
    chirality = load_csv(args.chirality_grid)
    chirality_gram = load_csv(args.chirality_sublattice_gram)
    plot_ground(ground, args.output_directory / "ground_state_observables")
    plot_low_states(states, args.output_directory / "low_lying_state_observables")
    plot_heatmaps(
        heatmap, chirality, chirality_gram,
        args.output_directory / "structure_factor_heatmaps",
    )


if __name__ == "__main__":
    main()
