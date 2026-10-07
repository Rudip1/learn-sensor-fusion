"""Shared matplotlib style for the figures and notebooks (presentation only, no algorithms)."""

from __future__ import annotations

import matplotlib as mpl
import matplotlib.pyplot as plt

#: Categorical colours in fixed order: blue, orange, aqua, yellow, magenta, green, violet, red.
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948"]
INK = "#0b0b0b"
INK_SECONDARY = "#52514e"
GRID = "#d9d8d4"
SURFACE = "#fcfcfb"


def use_style() -> None:
    """Thin marks, recessive axes and grid, fixed colour order."""
    mpl.rcParams.update(
        {
            "figure.facecolor": SURFACE,
            "axes.facecolor": SURFACE,
            "savefig.facecolor": SURFACE,
            "figure.dpi": 110,
            "savefig.dpi": 130,
            "axes.prop_cycle": mpl.cycler(color=SERIES),
            "axes.edgecolor": GRID,
            "axes.labelcolor": INK_SECONDARY,
            "axes.titlesize": 11,
            "axes.titlecolor": INK,
            "axes.grid": True,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "grid.color": GRID,
            "grid.linewidth": 0.6,
            "xtick.color": INK_SECONDARY,
            "ytick.color": INK_SECONDARY,
            "lines.linewidth": 1.6,
            "lines.markersize": 4,
            "legend.frameon": False,
            "font.size": 9.5,
        }
    )


def figure(nrows: int = 1, ncols: int = 1, width: float = 7.0, height: float = 3.2, **kw):
    use_style()
    return plt.subplots(nrows, ncols, figsize=(width, height), constrained_layout=True, **kw)
