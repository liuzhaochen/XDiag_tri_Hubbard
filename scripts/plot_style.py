"""Shared publication plotting style for the project."""

from pathlib import Path

import matplotlib as mpl
from matplotlib import font_manager


FONT_PATH = (
    Path(__file__).resolve().parents[1]
    / "assets"
    / "fonts"
    / "firamath"
    / "FiraMath-Regular.otf"
)


def apply_plot_style():
    """Use the vendored Fira Math font for text and math in every figure."""
    if not FONT_PATH.is_file():
        raise FileNotFoundError(f"required plot font is missing: {FONT_PATH}")

    font_manager.fontManager.addfont(FONT_PATH)
    family = font_manager.FontProperties(fname=FONT_PATH).get_name()
    mpl.rcParams.update(
        {
            "font.family": family,
            "mathtext.fontset": "custom",
            "mathtext.rm": family,
            "mathtext.it": family,
            "mathtext.bf": family,
            "mathtext.sf": family,
            "mathtext.tt": family,
            "mathtext.cal": family,
            "mathtext.fallback": "stixsans",
            "axes.unicode_minus": False,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )
