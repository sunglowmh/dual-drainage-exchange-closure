# -*- coding: utf-8 -*-
"""Figure 5 -- C3: is the head-driven two-way closure admissible here?

(a) the three candidate closures return very different amounts of water in the
    1.5 h peak window, and the two extremes are both artefacts
(b) overflow budget over the 13.7 h window: every cubic metre returned adds
    0.91 m3 of overflow -- the return flow damps rather than amplifies
(c) export budget: only 9 % of the returned water actually leaves the system
(d) the node-level capacity constraint is a false premise -- summing node-level
    spare capacity overestimates the network's acceptance by ~770x
"""
import json
import sys
import os

ROOT = os.environ.get("EMS_SWMM_BASE",
                   os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
import matplotlib.pyplot as plt
import numpy as np

import figstyle as F

F.use_style(scale=1.15)   # design scale for this figure (see figstyle)
C = json.load(open(F.WORK + "/fig_c3.json", encoding="utf-8"))
EC, CP = C["export_check"], C["capacity"]

fig = plt.figure(figsize=(F.W2, 132 * F.MM))
gs = fig.add_gridspec(2, 2, left=0.075, right=0.985, bottom=0.075, top=0.945,
                      wspace=0.30, hspace=0.52)
axa, axb = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])
axc, axd = fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1])

# ------------------------------------------------- (a) return volume by closure
cl = C["closures"]
xs = np.arange(3)
vals = [c["return_peak_m3"] for c in cl]
cols = ["#C8CDD2", "#E3B8AE", F.OI["blue"]]
axa.bar(xs, vals, 0.58, color=cols, edgecolor="#5A5A5A", linewidth=0.5, zorder=3)
for x_, c in zip(xs, cl):
    axa.text(x_, c["return_peak_m3"] + 250,
             "%.0f m³" % c["return_peak_m3"], ha="center", va="bottom",
             fontsize=F.FS(6.4), fontweight="bold")
    axa.text(x_, c["return_peak_m3"] * 0.52, "%.0f %%\nin peak"
             % (100 * c["peak_share"]), ha="center", va="center", fontsize=F.FS(6.0),
             color="white" if c["key"] == "c" else "#404040",
             fontweight="bold", linespacing=1.3)
    axa.text(x_, 15.6e3, c["verdict"], ha="center", va="center", fontsize=F.FS(6.1),
             color=F.OI["vermillion"] if c["key"] != "c" else "#333333",
             style="italic")
axa.set_xticks(xs)
axa.set_xticklabels(["(" + c["key"] + ") " + c["name"] for c in cl], fontsize=F.FS(6.3))
axa.set_xlim(-0.62, 2.62)
axa.set_ylim(0, 17.2e3)
axa.set_yticks([0, 4000, 8000, 12000, 16000])
axa.set_ylabel("Returned into the sewer within\nthe 1.5 h peak window (m³)")
axa.grid(axis="y")
F.panel(axa, "a", dx=-0.22, dy=1.02)

# ----------------------------------------------------------------- budgets
def budget(ax, tag, base, delta, after, ylab, note, color_delta):
    z = ["baseline", "injected return", "after injection"]
    cols = [F.OI["grey"], "#3D3D3D"]
    ax.bar([0, 2], [base, after], 0.58, color=cols,
           edgecolor="white", linewidth=0.5, zorder=3)
    ax.bar([1], [base], 0.58, color=F.OI["grey"], edgecolor="white",
           linewidth=0.5, zorder=3)
    ax.bar([1], [delta], 0.58, bottom=[base], color=color_delta,
           edgecolor="white", linewidth=0.5, zorder=3)
    ax.plot([0.29, 0.71], [base, base], color="#333333", lw=0.5, ls=":", zorder=4)
    ax.plot([0.71, 1.29], [base + delta, base + delta], color="#333333", lw=0.5,
            ls=":", zorder=4)
    for xi, v in zip([0, 1, 2], [base, delta, after]):
        ax.text(xi, max(base, base + delta, after) * 1.045 if xi == 1 else v * 1.045,
                ("+%.0f" % v) if xi == 1 else "%.0f" % v, ha="center",
                va="bottom", fontsize=F.FS(6.3), fontweight="bold",
                color=color_delta if xi == 1 else "#333333")
    ax.set_xticks([0, 1, 2])
    ax.set_xticklabels(z, fontsize=F.FS(6.3))
    ax.set_xlim(-0.62, 2.62)
    ax.set_ylim(0, max(base + delta, after) * 1.30)
    ax.set_ylabel(ylab)
    ax.grid(axis="y")
    ax.text(0.02, 0.975, note, transform=ax.transAxes, fontsize=F.FS(6.0),
            va="top", ha="left", color="#333333", linespacing=1.5)
    F.panel(ax, tag, dx=-0.22, dy=1.02)


budget(axb, "b", EC["flood_base_m3"], EC["d_flood_m3"], EC["flood_after_m3"],
       "Node overflow over the\n13.7 h window (m³)",
       "response k=Δoverflow/inflow= %.2f\n(damps, does not amplify)"
       % EC["amplification"],
       "#D9836E")
budget(axc, "c", EC["export_base_m3"], EC["d_export_m3"], EC["export_after_m3"],
       "Export at the 7 outlets over\nthe 13.7 h window (m³)",
       # Kept to two SHORT lines on purpose: a long second line runs under the
       # +delta label above the middle bar.  The rest of the argument ("fixed
       # point exists but needs 3.8x the event water") is in the caption.
       "κ = Δexport/inflow = %.2f\n"
       "φk = %.2f < 1"
       % (EC["kappa"], EC["phi_k"]),
       F.OI["vermillion"])

# --------------------------------------------- (d) local vs system capacity
names = ["Node-level\nsum (p50)", "Node-level\nsum (max)", "Observed\nexport"]
vals = [CP["local_sum_p50_cms"], CP["local_sum_max_cms"], CP["system_export_cms"]]
axd.bar([0, 1, 2], vals, 0.56, color=["#C8CDD2", "#C8CDD2", F.OI["blue"]],
        edgecolor="#5A5A5A", linewidth=0.5, zorder=3)
for xi, v in zip([0, 1, 2], vals):
    axd.text(xi, v * 1.30, "%.2f" % v, ha="center", va="bottom", fontsize=F.FS(6.3),
             fontweight="bold")
axd.set_yscale("log")
F.log_ticks(axd, "y")          # plain-text 10^n labels, never mathtext
axd.set_ylim(0.05, 900)
axd.set_xticks([0, 1, 2])
axd.set_xticklabels(names, fontsize=F.FS(6.3))
axd.set_xlim(-0.62, 2.62)
axd.set_ylabel("Acceptance / export\ncapacity (m³/s)")
axd.grid(axis="y")
axd.axhline(CP["flood_peak_cms"], color=F.OI["vermillion"], lw=0.9, ls="--",
            zorder=4)
axd.text(2.45, CP["flood_peak_cms"] * 1.18,
         "peak overflow rate\n%.1f m³/s" % CP["flood_peak_cms"],
         fontsize=F.FS(5.9), color=F.OI["vermillion"], va="bottom", ha="right",
         linespacing=1.4)
axd.annotate("", xy=(2.0, 2.0), xytext=(0.0, 2.0),
             arrowprops=dict(arrowstyle="<->", lw=0.7, color="#333333"))
axd.text(1.0, 2.6, "× %d overestimate" % CP["overestimate"],
         fontsize=F.FS(6.3), va="bottom", ha="center", color="#333333")
F.panel(axd, "d", dx=-0.22, dy=1.02)

F.save(fig, "Fig6_E3_closures")
