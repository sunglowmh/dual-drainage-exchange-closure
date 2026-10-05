# -*- coding: utf-8 -*-
"""Shared Nature/JOH figure house style for the paper-1 figure set.

Everything the figures need that is *style*, not data:
  * rcParams (Arial, 7 pt ticks, no top/right spines, editable vector text)
  * an Okabe-Ito (colour-blind safe) palette with a fixed event mapping
  * helpers: panel labels, scale bar, north arrow, column widths, save()

Run any fig*.py with the project interpreter:
  C:/Users/DELL/.workbuddy/binaries/python/envs/swmm/Scripts/python.exe
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrow, Polygon as MPoly
import numpy as np
ROOT = os.environ.get("EMS_SWMM_BASE",
                   os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# --------------------------------------------------------------------- paths
ROOT = ROOT
WORK = ROOT + r"/work"
OUT = ROOT + r"/output"
FIGDIR = OUT + r"/figs"
os.makedirs(FIGDIR, exist_ok=True)

# ------------------------------------------------------------------ palette
# Okabe & Ito (2008) -- safe for all common colour-vision deficiencies.
OI = {
    "blue": "#0072B2",
    "orange": "#E69F00",
    "green": "#009E73",
    "vermillion": "#D55E00",
    "purple": "#CC79A7",
    "sky": "#56B4E9",
    "yellow": "#F0E442",
    "black": "#000000",
    "grey": "#7F7F7F",
    "lgrey": "#D9D9D9",
}

# Fixed event identity across ALL panels of ALL figures.
EVENTS = ["20250627", "20250730", "20260719", "20260809"]
EVCOL = {"20250627": OI["blue"], "20250730": OI["green"],
         "20260719": OI["vermillion"], "20260809": OI["purple"]}
EVMK = {"20250627": "o", "20250730": "s", "20260719": "^", "20260809": "D"}
EVLBL = {
    "20250627": "2025-06-27  63.6 mm / 0.8 h",
    # depths are trapezoidal integrals over the irregular timestamps in the
    # INP (scripts/_M6_rain.py), NOT "assumed step x sum"
    "20250730": "2025-07-30  125.3 mm / 13.2 h",
    "20260719": "2026-07-19  178.6 mm / 5.2 h",
    "20260809": "2026-08-09  242.5 mm / 36 h",
}
EVSHORT = {"20250627": "2025-06-27", "20250730": "2025-07-30",
           "20260719": "2026-07-19", "20260809": "2026-08-09"}

# --------------------------------------------------------------- rcParams
MM = 1.0 / 25.4
W1, W15, W2 = 89 * MM, 120 * MM, 183 * MM       # single / 1.5 / double column

# Global font scale.  Every explicit `fontsize=` in the fig*.py scripts is
# wrapped in FS(), so a figure's type size can be changed without touching its
# layout code.  The scale a figure is PROOFED at is declared in that figure's
# own `use_style(scale=...)` call rather than left to whoever set the
# environment at render time -- otherwise re-rendering silently changes the art:
#     Fig. 1, Fig. 2   scale 1.0   (design sizes 7.0-7.5 pt; already >= floor)
#     Fig. 3 - Fig. 7  scale 1.15  (design sizes reach down to 5.3 pt)
FONT_SCALE = float(os.environ.get("FIGFONT", "1.0"))
#
# FS() also holds every size at or above FLOOR: Elsevier's artwork floor is 7 pt
# (line art / axis text), and a few annotations were designed at 5.3-6.0 pt to
# fit.  Enforcing the floor here rather than at 30-odd call sites means a size
# can never silently slip back below it.
FLOOR = float(os.environ.get("FIGFLOOR", "7.0"))


def FS(pt):
    """Scale one font size (pt) by FIGFONT, never returning less than FLOOR."""
    return max(FLOOR, pt * FONT_SCALE)


def _flr(v):
    """Floor a size that is NOT routed through FIGFONT (the rcParams defaults)."""
    return max(FLOOR, v)


def use_style(scale=None):
    """Install the house rcParams.  `scale` pins FONT_SCALE for this figure."""
    global FONT_SCALE
    if scale is not None:
        FONT_SCALE = float(scale)
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": _flr(7),
        "axes.labelsize": _flr(7.5),
        "axes.titlesize": _flr(7.5),
        "xtick.labelsize": _flr(7),
        "ytick.labelsize": _flr(7),
        "legend.fontsize": _flr(7),
        "axes.linewidth": 0.6,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "xtick.minor.width": 0.4,
        "ytick.minor.width": 0.4,
        "xtick.major.size": 2.6,
        "ytick.major.size": 2.6,
        "xtick.minor.size": 1.5,
        "ytick.minor.size": 1.5,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "legend.frameon": False,
        "legend.handlelength": 1.6,
        "legend.labelspacing": 0.3,
        "legend.borderpad": 0.2,
        "figure.dpi": 200,
        "savefig.dpi": 400,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.02,
        "pdf.fonttype": 42,       # TrueType -> editable text in the PDF
        "ps.fonttype": 42,
        "svg.fonttype": "none",   # keep text as text in the SVG
        # mathtext must use the SAME face as the rest of the figure.  The
        # default fontset is DejaVu, so every $...$ (Greek letters, symbols)
        # silently switched typeface and the artwork ended up mixing Arial
        # with DejaVu -- DejaVu is not an Elsevier-recommended font.
        "mathtext.fontset": "custom",
        "mathtext.rm": "Arial",
        "mathtext.it": "Arial:italic",
        "mathtext.bf": "Arial:bold",
        "mathtext.default": "it",
        "axes.unicode_minus": False,
        "grid.linewidth": 0.4,
        "grid.color": OI["lgrey"],
        "grid.alpha": 0.8,
    })


# ------------------------------------------------------------------ helpers
def panel(ax, label, dx=-0.13, dy=1.06, fs=None):
    """Lowercase bold panel label in the top-left corner.

    `fs` is a FINAL point size; leave it None to take the house size (>= 7 pt).
    """
    ax.text(dx, dy, label, transform=ax.transAxes,
            fontsize=FS(7) if fs is None else _flr(fs),
            fontweight="bold", va="bottom", ha="left")


def despine(ax, left=False, bottom=False):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    if left:
        ax.spines["left"].set_visible(False)
    if bottom:
        ax.spines["bottom"].set_visible(False)


def scale_bar(ax, length_m, label=None, y=0.045, x=0.045, color="k", lw=1.6,
              fs=None, unit="m"):
    """Horizontal scale bar drawn in axes fraction, length in map units.

    `fs` is a FINAL point size (pass FS(7) for the house size); None -> FS(7).
    """
    x0, x1 = ax.get_xlim()
    y0, y1 = ax.get_ylim()
    frac = length_m / (x1 - x0)
    X0, Y0 = x0 + x * (x1 - x0), y0 + y * (y1 - y0)
    ax.plot([X0, X0 + frac * (x1 - x0)], [Y0, Y0], color=color, lw=lw,
            solid_capstyle="butt", zorder=30, clip_on=False)
    ax.text(X0 + 0.5 * frac * (x1 - x0), Y0 + 0.012 * (y1 - y0),
            label or ("%g %s" % (length_m, unit)), ha="center", va="bottom",
            fontsize=FS(7) if fs is None else _flr(fs), color=color, zorder=30)


def north_arrow(ax, x=0.955, y=0.915, size=0.055, color="k", fs=None):
    """Small north arrow in axes fraction (`fs` = final pt; None -> FS(7))."""
    ax.annotate("", xy=(x, y + size), xytext=(x, y),
                xycoords="axes fraction",
                arrowprops=dict(arrowstyle="-|>", color=color, lw=1.0,
                                mutation_scale=8), zorder=30)
    ax.text(x, y + size + 0.012, "N", ha="center", va="bottom",
            fontsize=FS(7) if fs is None else _flr(fs),
            color=color, zorder=30, transform=ax.transAxes)


def nice_step(x):
    """A 1/2/5 x 10^k step near x -- for ticks that must not print on top of
    each other.  A short axis (e.g. a 14 mm rainfall strip) cannot carry the
    hand-picked 0/25/50 labels of a 200 mm/h range: pick the step from the range.
    """
    if x <= 0:
        return 1.0
    e = float(np.floor(np.log10(x)))
    f = x / 10.0 ** e
    m = 1.0 if f <= 1.5 else (2.0 if f <= 3.5 else (5.0 if f <= 7.5 else 10.0))
    return m * 10.0 ** e


def gaussian_arrow(ax, xy, xytext, color="k", lw=0.7, ls="-", rad=0.0,
                   head=6, shrinkA=1.5, shrinkB=1.5, zorder=10, alpha=1.0):
    ax.annotate("", xy=xy, xytext=xytext,
                arrowprops=dict(arrowstyle="-|>", color=color, lw=lw,
                                linestyle=ls, mutation_scale=head,
                                shrinkA=shrinkA, shrinkB=shrinkB,
                                connectionstyle="arc3,rad=%.3f" % rad,
                                joinstyle="miter"), zorder=zorder, alpha=alpha)


# ------------------------------------------------------------- log axis labels
# matplotlib's own log formatter emits mathtext ($10^{4}$), and mathtext sets the
# exponent at 0.7x the base size -- 4.9 pt on a 7 pt base, under the 7 pt floor.
# Arial has every superscript DIGIT (U+2070, U+00B9, U+00B2, U+00B3, U+2074-2079)
# but NOT the superscript minus (U+207B), so negative powers use decimals.
_SUP = {"0": "\u2070", "1": "\u00b9", "2": "\u00b2", "3": "\u00b3",
        "4": "\u2074", "5": "\u2075", "6": "\u2076", "7": "\u2077",
        "8": "\u2078", "9": "\u2079"}


def log_label(v, pos=None):
    """Plain-text label for one log-scale tick: 10**4 -> '10' + SUPERSCRIPT 4."""
    if v <= 0:
        return ""
    e = int(round(float(np.log10(v))))
    if abs(10.0 ** e - v) > 1e-6 * v:
        return "%g" % v          # not a power of ten
    if e >= 0:
        return "10" + "".join(_SUP[c] for c in str(e))
    return "%g" % v              # 0.1, 0.01 ... (no superscript minus in Arial)


def log_ticks(ax, which="both"):
    """Force plain-text labels on a log-scaled axis (see log_label)."""
    from matplotlib.ticker import FuncFormatter, NullFormatter
    f = FuncFormatter(log_label)
    if which in ("x", "both"):
        ax.xaxis.set_major_formatter(f)
        ax.xaxis.set_minor_formatter(NullFormatter())
    if which in ("y", "both"):
        ax.yaxis.set_major_formatter(f)
        ax.yaxis.set_minor_formatter(NullFormatter())


def light_rain_axis(ax, t, rain_mm_h, color="#9EC9E2", width_frac=1.0,
                    ylim=None, fs=None):
    """Inverted rainfall bars on a twin axis (standard hydrological layout).

    Units are written `mm/h`, not `mm h^-1`: Arial has no SUPERSCRIPT MINUS
    (U+207B), so that form would silently fall back to another face.
    """
    fs = FS(7) if fs is None else _flr(fs)
    ax2 = ax.twinx()
    ax2.spines["top"].set_visible(False)
    ax2.spines["left"].set_visible(False)
    ax2.spines["right"].set_color(color)
    ax2.tick_params(axis="y", colors=color, labelsize=fs)
    ax2.bar(t, rain_mm_h, width=width_frac, color=color, lw=0, zorder=0)
    ax2.set_ylim(ylim if ylim else (max(rain_mm_h) * 3.2, 0))
    ax2.set_ylabel("Rainfall (mm/h)", color=color, fontsize=fs)
    return ax2


def save(fig, name, formats=("pdf", "png"), tight=True):
    out = []
    for f in formats:
        p = "%s/%s.%s" % (FIGDIR, name, f)
        if tight:
            fig.savefig(p, format=f, facecolor="white",
                        bbox_inches="tight", pad_inches=0.02)
        else:
            fig.savefig(p, format=f, facecolor="white", bbox_inches=None)
        out.append(p)
    plt.close(fig)
    print("saved:", ", ".join(out))
    return out
