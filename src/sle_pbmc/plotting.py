"""Figure style shared by all phase-1 plots.

Colours come from a validated categorical palette (checked for colour-blind
separation). Two-group comparisons use slots 1-2. Anything with more than two
groups on a scatter (batches, cell types) is drawn as small multiples: one
panel per group, highlighted against the rest in grey, which stays readable
for every reader and never needs more than one hue.
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"
BACKGROUND_CELLS = "#dcdbd5"
SLOT = ["#2a78d6", "#eb6834", "#1baf7a"]  # blue, orange, aqua

# fixed colour per disease label, whatever order the data arrive in
DISEASE_COLORS = {"normal": SLOT[0], "systemic lupus erythematosus": SLOT[1]}
DISEASE_SHORT = {"normal": "Healthy", "systemic lupus erythematosus": "SLE"}


def apply_style() -> None:
    plt.rcParams.update(
        {
            "figure.facecolor": SURFACE,
            "axes.facecolor": SURFACE,
            "savefig.facecolor": SURFACE,
            "font.family": "sans-serif",
            "font.size": 9,
            "axes.edgecolor": AXIS,
            "axes.labelcolor": INK_2,
            "axes.titlesize": 10,
            "axes.titlecolor": INK,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "xtick.color": MUTED,
            "ytick.color": MUTED,
            "axes.grid": False,
            "grid.color": GRID,
            "grid.linewidth": 0.6,
            "legend.frameon": False,
            "savefig.dpi": 150,
            "savefig.bbox": "tight",
        }
    )


def _clean_axis(ax) -> None:
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)


def umap_highlight(ax, xy: np.ndarray, mask: np.ndarray, title: str, color: str = SLOT[0], s: float = 0.6) -> None:
    """All cells in grey, the highlighted group on top in one colour."""
    ax.scatter(xy[~mask, 0], xy[~mask, 1], s=s, c=BACKGROUND_CELLS, linewidths=0, rasterized=True)
    ax.scatter(xy[mask, 0], xy[mask, 1], s=s, c=color, linewidths=0, rasterized=True)
    ax.set_title(f"{title}\n{int(mask.sum()):,} cells", fontsize=8.5, color=INK)
    _clean_axis(ax)


def umap_two_groups(ax, xy: np.ndarray, labels: np.ndarray, title: str, seed: int = 0, s: float = 0.6,
                    legend: bool = False) -> None:
    """Disease on one panel: points drawn in random order so neither group hides the other."""
    order = np.random.default_rng(seed).permutation(len(labels))
    colors = np.array([DISEASE_COLORS.get(str(v), MUTED) for v in labels])
    ax.scatter(xy[order, 0], xy[order, 1], s=s, c=colors[order], linewidths=0, rasterized=True)
    ax.set_title(title, fontsize=8.5, color=INK)
    _clean_axis(ax)
    if legend:
        add_disease_legend(ax.figure)


def add_disease_legend(fig) -> None:
    """One legend for the whole figure, top right, clear of the data."""
    from matplotlib.lines import Line2D

    handles = [Line2D([], [], marker="o", ls="", ms=6, color=c, label=DISEASE_SHORT[n])
               for n, c in DISEASE_COLORS.items()]
    fig.legend(handles=handles, loc="upper right", ncol=2, fontsize=8.5, frameon=False,
               bbox_to_anchor=(0.99, 1.0))
