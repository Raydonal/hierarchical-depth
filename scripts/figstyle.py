"""One figure style for the whole manuscript (Physica A, elsarticle
``preprint,12pt``: text width 390 pt = 5.40 in).

Figures are produced at their final physical size, so LaTeX includes them at
``width=\\textwidth`` (full) or ``0.485\\textwidth`` (half) without rescaling
the fonts: 8 pt tick labels, 9 pt axis labels, legends and bold panel labels.

Colour: the Okabe-Ito palette for categories and a truncated viridis map for
the ordered depth ``N`` -- both safe under colour-vision deficiency. No
distinction is carried by colour alone: every series also gets its own line
style and/or marker (``depth_style``, ``LINESTYLES``, ``MARKERS``).
"""

from __future__ import annotations

import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib import font_manager  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parents[1]
FIGURES = REPO_ROOT / "figures"

FULL_WIDTH = 5.40     # in  (= 390 pt text width)
HALF_WIDTH = 2.62     # in  (= 0.485 text width)

OKABE_ITO = {
    "black": "#000000",
    "orange": "#E69F00",
    "sky": "#56B4E9",
    "green": "#009E73",
    "yellow": "#F0E442",     # avoid on white
    "blue": "#0072B2",
    "vermillion": "#D55E00",
    "purple": "#CC79A7",
}
# Category order (yellow left out: unreadable on white).
CATEGORY = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#000000"]
MARKERS = ["o", "s", "^", "D", "v", "P", "X"]
LINESTYLES = ["-", "--", "-.", ":", (0, (5, 1.5, 1, 1.5, 1, 1.5)), (0, (8, 2))]
GREY = "0.45"


def _prefer_text_size_latin_modern():
    """Latin Modern ships one file per optical size under a single family
    name; keep only the 10 pt design so matplotlib's choice is deterministic."""
    fm = font_manager.fontManager
    fm.ttflist = [
        f for f in fm.ttflist
        if not (f.name == "Latin Modern Roman" and "lmroman10-" not in Path(f.fname).name)
    ]


def apply():
    _prefer_text_size_latin_modern()
    # Missing candidate families (CMU Serif is rarely installed) are not an error.
    import logging
    logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)
    warnings.filterwarnings("ignore", message="findfont")
    plt.rcParams.update({
        "font.family": "serif",
        # "Latin Modern Roman" is the Unicode build of Computer Modern shipped
        # with TeX Live; it stands in when CMU Serif is not installed.
        "font.serif": ["CMU Serif", "Computer Modern Roman", "Latin Modern Roman",
                       "DejaVu Serif"],
        "mathtext.fontset": "cm",
        "axes.unicode_minus": True,
        "font.size": 9,
        "axes.labelsize": 9,
        "axes.titlesize": 9,
        "legend.fontsize": 9,
        "legend.title_fontsize": 9,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "lines.linewidth": 1.2,
        "lines.markersize": 4,
        "axes.linewidth": 0.6,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "xtick.minor.width": 0.4,
        "ytick.minor.width": 0.4,
        "xtick.major.size": 3,
        "ytick.major.size": 3,
        "xtick.minor.size": 1.8,
        "ytick.minor.size": 1.8,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "legend.frameon": False,
        "legend.handlelength": 2.4,
        "legend.borderaxespad": 0.3,
        "legend.labelspacing": 0.3,
        "figure.constrained_layout.use": True,
        "figure.constrained_layout.w_pad": 0.04,
        "figure.constrained_layout.h_pad": 0.04,
        "figure.constrained_layout.wspace": 0.08,
        "figure.constrained_layout.hspace": 0.08,
        "figure.dpi": 150,
        "savefig.dpi": 300,
        "savefig.bbox": None,
        "pdf.fonttype": 3,
    })


def depth_colors(n: int):
    """``n`` colours for the ordered depths ``N = 0..n-1``: viridis truncated to
    0.0-0.85 (no near-white yellow)."""
    return [plt.cm.viridis(x) for x in np.linspace(0.0, 0.85, n)]


def depth_style(N: int, n: int):
    """Colour + line style + marker for depth ``N`` out of ``n`` depths."""
    return {"color": depth_colors(n)[N], "linestyle": LINESTYLES[N % len(LINESTYLES)],
            "marker": MARKERS[N % len(MARKERS)]}


def panel_label(ax, label: str, x: float = 0.0, y: float = 1.0, dx: float = 0.0,
                dy: float = 3.0):
    """Bold 9 pt panel label, "(a)", anchored at the top-left corner of the
    axes box and placed just above it."""
    ax.annotate(label, xy=(x, y), xycoords="axes fraction", xytext=(dx, dy),
                textcoords="offset points", ha="left", va="bottom",
                fontsize=9, fontweight="bold")


def save(fig, stem: str):
    """PDF (vector) + PNG (300 dpi preview) at the exact physical size."""
    FIGURES.mkdir(exist_ok=True)
    fig.savefig(FIGURES / f"{stem}.pdf", bbox_inches=None)
    fig.savefig(FIGURES / f"{stem}.png", bbox_inches=None, dpi=300)
    plt.close(fig)
    w, h = fig.get_size_inches()
    print(f"wrote figures/{stem}.pdf (+ .png)  [{w:.2f} x {h:.2f} in]")
