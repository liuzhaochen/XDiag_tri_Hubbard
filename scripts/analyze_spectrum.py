#!/usr/bin/env python3
"""Text analysis of a validated long-form spectrum table."""

import argparse
from pathlib import Path

import numpy as np


def key_from_row(row):
    return (int(row["nx"]), int(row["ny"]), int(row["C2_parity"]))


def label(key):
    nx, ny, parity = key
    if parity != 0:
        irrep = "A" if parity > 0 else "B"
        momentum = "Gamma" if (nx, ny) == (0, 0) else "M2"
        return f"{momentum}({irrep})"
    return f"({nx},{ny})"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("table", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    data = np.genfromtxt(args.table, delimiter=",", names=True)
    u_values = np.unique(data["U_over_t"])
    spacings = np.diff(u_values)
    nominal_spacing = float(np.min(spacings)) if len(spacings) else np.inf
    contiguous = spacings <= 1.5 * nominal_spacing
    keys = sorted({key_from_row(row) for row in data}, key=lambda x: (x[0], x[1], -x[2]))

    level_one = {}
    for key in keys:
        mask = (
            (data["nx"] == key[0])
            & (data["ny"] == key[1])
            & (data["C2_parity"] == key[2])
            & (data["level_index"] == 1)
        )
        rows = data[mask]
        rows = rows[np.argsort(rows["U_over_t"])]
        level_one[key] = rows["energy"]

    ground_keys = []
    ground_energies = []
    first_gaps = []
    for u_value in u_values:
        rows = data[data["U_over_t"] == u_value]
        order = np.argsort(rows["energy"])
        ground = rows[order[0]]
        ground_keys.append(key_from_row(ground))
        ground_energies.append(float(ground["energy"]))
        first_gaps.append(float(rows[order[1]]["energy"] - ground["energy"]))

    intervals = []
    start = 0
    for index in range(1, len(u_values) + 1):
        separated = index < len(u_values) and not contiguous[index - 1]
        if index == len(u_values) or separated or ground_keys[index] != ground_keys[start]:
            intervals.append((u_values[start], u_values[index - 1], ground_keys[start]))
            start = index

    crossing_brackets = []
    for first_index, first_key in enumerate(keys):
        for second_key in keys[first_index + 1 :]:
            delta = level_one[first_key] - level_one[second_key]
            for index in range(len(u_values) - 1):
                if not contiguous[index]:
                    continue
                if delta[index] == 0 or delta[index] * delta[index + 1] < 0:
                    crossing_brackets.append(
                        (u_values[index], u_values[index + 1], first_key, second_key)
                    )

    lines = [
        "3x4 triangular Hubbard spectrum analysis",
        f"U/t grid span: {u_values[0]:g} ... {u_values[-1]:g} ({len(u_values)} points)",
        f"symmetry blocks: {len(keys)}; levels per block: {int(np.max(data['level_index']))}",
        f"maximum reported residual: {np.max(data['residual_norm']):.3e}",
        f"maximum explicit residual: {np.max(data['explicit_residual_norm']):.3e}",
        f"maximum orthogonality error: {np.max(data['orthogonality_error']):.3e}",
        "",
        "Lowest sectors among the included blocks:",
    ]
    for first_u, last_u, key in intervals:
        lines.append(f"  U/t={first_u:g} ... {last_u:g}: {label(key)}")
    gap_index = int(np.argmin(first_gaps))
    lines.extend(
        [
            "",
            f"Smallest sampled many-body gap: {first_gaps[gap_index]:.12g} t "
            f"at U/t={u_values[gap_index]:g}",
            "",
            "Crossings of block-lowest branches (sampled brackets):",
        ]
    )
    if crossing_brackets:
        for left, right, first_key, second_key in crossing_brackets:
            lines.append(
                f"  U/t={left:g} ... {right:g}: {label(first_key)} <-> {label(second_key)}"
            )
    else:
        lines.append("  none")

    text = "\n".join(lines) + "\n"
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(text)
    print(text, end="")


if __name__ == "__main__":
    main()
