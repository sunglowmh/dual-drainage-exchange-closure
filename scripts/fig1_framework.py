# -*- coding: utf-8 -*-
"""Figure 1 -- the four-layer dual-drainage modelling chain and the three
knowledge gaps this paper closes.

Pure vector schematic (no axes).  Width = double column (183 mm).

Content is tied to the manuscript by a fixed mapping (Table 2, Fig. 1 caption):
    G1 transferability        <- E1        -> Fig. 4
    G2 identifiability AND
       admissibility          <- E2 + E3   -> Figs. 5, 6, 7
    G3 decision consequence   <- E4        -> Table 6

Every type size is >= 7 pt (Elsevier: 7 pt for normal text, >= 6 pt for
sub/superscripts, measured at the FINAL PRINTED size -- these figures are
drawn at their printed width, so design size == printed size).
"""
import sys
import os

ROOT = os.environ.get("EMS_SWMM_BASE",
                   os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

import figstyle as F

F.use_style()

W, H = 183.0, 137.0
fig = plt.figure(figsize=(W * F.MM, H * F.MM))
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, W)
ax.set_ylim(0, H)
ax.axis("off")

TINT = {"L1": "#E8F1F8", "L2": "#E7F4EF", "L3": "#F3EFF7", "L4": "#FBEADF"}
EDGE = {"L1": "#2E6E9E", "L2": "#1F7A5C", "L3": "#7A5C9E", "L4": F.OI["vermillion"]}
ACCENT = {"gap": "#F5F2E9", "gapedge": "#8A7B4F"}

LX0, LX1 = 3.0, 104.0
RX0, RX1 = 110.0, 180.0
TOP = 128.0
# The four layer boxes need a gap wide enough for an ARROW, not just for a
# head: the box outline carries `pad=0.6`, so the clear space is GAPY - 1.2.
GAPY = 4.0
LS = 1.25                      # line spacing, in units of the font size
PT = 7.0                       # the body size; everything else is >= this

# =============================================================== LEFT: layers
layers = [
    ("L1", "OBSERVATION LAYER", "known \u2192 usable as a constraint", 28.0,
     "1 tipping-bucket rain gauge (1 min)  |  12 ultrasonic manhole sensors\n"
     "5 buried road-depth sensors  |  2 m bare-earth DEM, sewers, land use",
     "error model: telemetry freeze, full-shaft truncation, sensor zero offset;\n"
     "threshold alarms \u2192 missingness correlates with the signal (informative censoring)"),
    ("L2", "SEWER LAYER (1D)", "SWMM 5.2  DYNWAVE", 28.0,
     "3 277 nodes / 3 271 conduits \u2192 SWMM 3 283 nodes / 3 278 links;\n"
     "10 outlet structures (DWG); 2 pumping stations (wet well + PUMP3)",
     "outputs: node flooding rate (from .out), node head, outlet export"),
    ("L3", "SURFACE LAYER (2D)", "in-house solver", 24.0,
     "2 m diffusive-wave grid; buildings, greens and river as no-flow cells;\n"
     "open campus gates; lake stage\u2013area\u2013volume curve with gate and pump",
     "outputs: depth field, wetted area, lake stage"),
    ("L4", "EXCHANGE LAYER", "the only prescribed layer \u2014 this paper", 32.0,
     "(a) no exchange (one-way)\n"
     "(b) uniform ponding area A_ponded   (current practice; E1)\n"
     "(c) two-way return-flux closures   (three candidates; E3)",
     "L1\u2013L3 are each testable against their own observation class,\n"
     "so the effect of L4 can be ISOLATED"),
]

yy = TOP
gap_centres = []
for tag, name, sub, BH, body, foot in layers:
    y1 = yy - BH
    core = tag == "L4"
    ax.add_patch(FancyBboxPatch((LX0, y1), LX1 - LX0, BH,
                                boxstyle="round,pad=0.6,rounding_size=1.6",
                                linewidth=1.6 if core else 0.8,
                                edgecolor=EDGE[tag], facecolor=TINT[tag], zorder=3))
    ax.add_patch(FancyBboxPatch((LX0 + 1.4, yy - 7.4), 12.6, 5.6,
                                boxstyle="round,pad=0.2,rounding_size=0.8",
                                linewidth=0, facecolor=EDGE[tag], zorder=5))
    ax.text(LX0 + 7.7, yy - 4.6, tag, color="white", fontsize=F.FS(7.5),
            fontweight="bold", ha="center", va="center", zorder=6)
    ax.text(LX0 + 15.8, yy - 3.4, name, fontsize=F.FS(7.5), fontweight="bold",
            color=EDGE[tag], ha="left", va="center", zorder=6)
    ax.text(LX1 - 2.2, yy - 3.4, sub, fontsize=F.FS(PT), color="#555555",
            ha="right", va="center", style="italic", zorder=6)
    ax.plot([LX0 + 1.4, LX1 - 1.4], [yy - 8.6, yy - 8.6], color=EDGE[tag],
            lw=0.4, alpha=0.45, zorder=4)
    ax.text(LX0 + 2.2, yy - 10.2, body, fontsize=F.FS(PT), color="#1A1A1A",
            ha="left", va="top", linespacing=LS, zorder=6)
    ax.text(LX0 + 2.2, y1 + 2.6, foot, fontsize=F.FS(PT), color="#5A5A5A",
            ha="left", va="bottom", linespacing=LS, zorder=6)
    yy = y1 - GAPY
    gap_centres.append(y1 - GAPY / 2.0)      # centre of the gap BELOW this box
    if core:
        l4_centre = (y1 + BH) / 2.0

# flow arrows between the layers (top -> bottom).  Draw across the WHOLE clear
# gap with no shrink: the previous version spanned a fixed +/-0.9 units, i.e.
# 1.8 mm, of which 2 x 1.5 pt was shrink -- leaving less shaft than the head is
# long, so the arrow rendered as a bare triangle.
for yc in gap_centres[:-1]:
    F.gaussian_arrow(ax, (LX0 + 50.5, yc - GAPY / 2.0 + 0.75),
                     (LX0 + 50.5, yc + GAPY / 2.0 - 0.75),
                     color="#9A9A9A", lw=0.8, head=4.0, shrinkA=0.0,
                     shrinkB=0.0, zorder=2)

# ============================================================== RIGHT: gaps
ax.text(RX0, 134.5, "KNOWLEDGE GAPS  \u2192  EXPERIMENTS", fontsize=F.FS(7.5),
        fontweight="bold", color="#4D4D4D", ha="left", va="center")
ax.text(RX0, 131.4, "all three act on L4, the only prescribed layer",
        fontsize=F.FS(PT), color="#888888", ha="left", va="center", style="italic")

PAD_L = 2.4                      # card left padding
ROWT = 3.6                       # first row centre, below the panel top
ROWP = 3.08                      # line pitch at 7 pt with linespacing 1.25

gaps = [
    ("G1", "Transferability",
     "Can one calibrated ponding area be\ncarried from event to event?",
     [("E1", "Fig. 4",
       "0 / 100 / \u2265400 / \u2265400 m\u00b2\n\u03C1 = +1.00 with event size")],
     35.4, 15.8),
    ("G2", "Identifiability & admissibility",
     "Can the closure be inverted from a\nhandful of sensors, and is it admissible?",
     [("E2", "Fig. 5", "2\u20135 stations suffice; sampling binds (flips at 20\u201360 min)"),
      ("E3", "Figs. 6, 7", "\u03BA = 0.089, \u03C6k = 0.83: fixed point unreachable")],
     41.6, 22.0),
    ("G3", "Decision consequence",
     "Does the closure change which flood\nmitigation works, not just the depths?",
     [("E4", "Table 6",
       "gate closure \u2192 export 0; lake max 1.24 \u2192 2.50 m,\ncampus surface depth unchanged")],
     35.4, 15.8),
]

gy = TOP
for tag, name, q, rows, GH, PH in gaps:
    y1 = gy - GH
    ax.add_patch(FancyBboxPatch((RX0, y1), RX1 - RX0, GH,
                                boxstyle="round,pad=0.6,rounding_size=1.6",
                                linewidth=0.8, edgecolor=ACCENT["gapedge"],
                                facecolor=ACCENT["gap"], zorder=3))
    ax.text(RX0 + PAD_L, gy - 3.8, tag, fontsize=F.FS(8.0), fontweight="bold",
            color=ACCENT["gapedge"], ha="left", va="center")
    ax.text(RX0 + 9.0, gy - 3.8, name, fontsize=F.FS(7.5), fontweight="bold",
            color="#3D3D3D", ha="left", va="center")
    ax.text(RX0 + PAD_L, gy - 8.6, q, fontsize=F.FS(PT), color="#1A1A1A",
            ha="left", va="top", linespacing=LS)

    ptop = gy - GH + 3.0 + PH          # panel is anchored to the card bottom
    ax.add_patch(FancyBboxPatch((RX0 + PAD_L, ptop - PH), RX1 - RX0 - 2 * PAD_L, PH,
                                boxstyle="round,pad=0.3,rounding_size=0.8",
                                linewidth=0, facecolor="white", zorder=4))
    ycur = ptop - ROWT
    for k, (exp, figno, res) in enumerate(rows):
        ax.text(RX0 + PAD_L + 1.4, ycur, exp, fontsize=F.FS(PT),
                fontweight="bold", color="#1A1A1A", ha="left", va="center",
                zorder=5)
        ax.text(RX1 - PAD_L - 1.4, ycur, figno, fontsize=F.FS(PT),
                color="#888888", ha="right", va="center", zorder=5,
                style="italic")
        ax.text(RX0 + PAD_L + 1.4, ycur - 3.4, res, fontsize=F.FS(PT),
                color="#1A1A1A", ha="left", va="top", zorder=5,
                linespacing=LS)
        nlines = res.count("\n") + 1
        ycur -= 3.4 + nlines * ROWP + 2.0
    gy = y1 - 4.6

# arrow from the exchange layer into the gap column
F.gaussian_arrow(ax, (RX0 - 1.2, l4_centre), (LX1 + 1.2, l4_centre),
                 color=F.OI["vermillion"], lw=1.0, head=7, zorder=8)

F.save(fig, "Fig1_framework")
