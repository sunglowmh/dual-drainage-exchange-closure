# -*- coding: utf-8 -*-
"""Figure 3 -- E1: is one uniform ponding area transferable across events?

(a) bias-corrected RMSE vs A_ponded, four events
(b) the equal-weight composite score actually used to select the value
(c) selected A_ponded against the event size (m3 of overflow per flooding node)
(d) the node-level equivalent area V_node / max(hbar_wet, 0.10 m) inverted from
    the corrected one-way two-dimensional solution of 2026-07-19 (log abscissa)
"""
import json
import sys
import os

ROOT = os.environ.get("EMS_SWMM_BASE",
                   os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

import figstyle as F

F.use_style(scale=1.15)   # design scale for this figure (see figstyle)
D = json.load(open(F.WORK + "/fig_e1.json", encoding="utf-8"))
EFF = np.load(F.WORK + "/fig_e1_eff.npy")

fig = plt.figure(figsize=(F.W2, 160 * F.MM))
gs = fig.add_gridspec(2, 2, left=0.065, right=0.985, bottom=0.145, top=0.955,
                      wspace=0.30, hspace=0.34)
axa, axc = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])
axb, axd = fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1])

AREAS = [0, 50, 100, 200, 300, 400]        # labelled major ticks
AREAS_MI = [25, 75, 150, 250]              # unlabelled (the sweep is 9 levels)
LABA = ["0", "50", "100", "200", "300", "400"]

# ------------------------------------------------------------ (a) RMSE curves
for ev in F.EVENTS:
    d = D[ev]["cfg"]
    xs = [d[c]["area"] for c in d]
    ys = [d[c]["rmse"] for c in d]
    o = np.argsort(xs)
    axa.plot(np.array(xs)[o], np.array(ys)[o], "-", lw=0.9,
             color=F.EVCOL[ev], marker=F.EVMK[ev], ms=3.4,
             markeredgecolor="white", markeredgewidth=0.35,
             label=F.EVLBL[ev], zorder=3)
    b = D[ev]["best"]
    axa.plot([d[b]["area"]], [d[b]["rmse"]], marker="*", ms=10,
             markerfacecolor="none", markeredgecolor=F.EVCOL[ev],
             markeredgewidth=0.8, zorder=4)
axa.set_xticks(AREAS)
axa.set_xticklabels(LABA)
axa.set_xticks(AREAS_MI, minor=True)
axa.tick_params(which="minor", length=2.0)
axa.set_xlim(-30, 440)
# 0.80, not 0.72: the 2026-07-19 run at A_ponded = 25 has RMSE 0.761, which the
# old ceiling cut off -- the curve then entered and left through the top, and
# the reader saw a break instead of a peak.
axa.set_ylim(0.15, 0.80)
axa.set_xlabel("Prescribed ponding area  A_ponded  (m²)")
axa.set_ylabel("Bias-corrected RMSE (m)")
axa.grid(axis="y")
axa.annotate("RMSE alone cannot choose:\nat 2026-07-19 the 0 and 400 m²\n"
             "runs are 0.06 m apart",
             # 0.32, not 0.20: the 2026-08-09 curve runs at 0.28-0.30 through
             # this x-range and was drawn straight across the first line.
             xy=(400, 0.51), xytext=(150, 0.32), fontsize=F.FS(6.0),
             color="#444444", linespacing=1.45,
             arrowprops=dict(arrowstyle="-", lw=0.5, color="#888888"))
F.panel(axa, "a")
axa.text(0.02, 0.985, "optimal value marked with an open star",
         transform=axa.transAxes, fontsize=F.FS(6.0), color="#777777", va="top")

# --------------------------------------------------------- (b) score curves
for ev in F.EVENTS:
    d = D[ev]["cfg"]
    xs = [d[c]["area"] for c in d]
    ys = [d[c]["score"] for c in d]
    o = np.argsort(xs)
    axb.plot(np.array(xs)[o], np.array(ys)[o], "-", lw=0.9,
             color=F.EVCOL[ev], marker=F.EVMK[ev], ms=3.4,
             markeredgecolor="white", markeredgewidth=0.35, zorder=3)
    b = D[ev]["best"]
    axb.plot([d[b]["area"]], [d[b]["score"]], marker="*", ms=10,
             color=F.EVCOL[ev], zorder=4)
axb.axhline(1.0, color="#BBBBBB", lw=0.5, ls="--", zorder=1)
axb.set_xticks(AREAS)
axb.set_xticklabels(LABA)
axb.set_xticks(AREAS_MI, minor=True)
axb.tick_params(which="minor", length=2.0)
axb.set_xlim(-30, 440)
axb.set_ylim(0.36, 1.10)
axb.set_xlabel("Prescribed ponding area  A_ponded  (m²)")
axb.set_ylabel("Composite score  S  (\u2013)")
axb.grid(axis="y")
axb.text(0.98, 0.96, "S = 0.5 RMSE/RMSE_max"
         " + 0.5 |ΔT|/|ΔT|_max",
         transform=axb.transAxes, fontsize=F.FS(6.0), color="#555555",
         ha="right", va="top")
F.panel(axb, "b")

# --------------------------------------------------- (c) size vs optimum area
vpn = [D[ev]["v_per_node"] for ev in F.EVENTS]
opt = [D[ev]["cfg"][D[ev]["best"]]["area"] for ev in F.EVENTS]
for ev, x, y in zip(F.EVENTS, vpn, opt):
    axc.scatter([x], [y], s=42, c=F.EVCOL[ev], marker=F.EVMK[ev],
                edgecolors="white", linewidths=0.5, zorder=4)
    dx, dy = {"20250627": (6, -13), "20250730": (-8, 9),
              "20260719": (6, 7), "20260809": (2, -17)}[ev]
    axc.annotate(F.EVSHORT[ev], (x, y), textcoords="offset points",
                 xytext=(dx, dy), fontsize=F.FS(6.2), color=F.EVCOL[ev],
                 ha="right" if dx < 0 else "left")
xs = np.array([vpn[i] for i in range(4) if opt[i] > 0])
ys = np.array([opt[i] for i in range(4) if opt[i] > 0])
k, c = np.polyfit(xs, np.log(ys), 1)
xx = np.linspace(xs.min() - 2, xs.max() + 3, 40)
axc.plot(xx, np.exp(c) * np.exp(k * xx), "--", lw=0.7, color="#888888",
         zorder=2)
# R^2 of the log-linear fit, computed rather than hard-coded (only n = 3
# points, so it is reported with its sample size, as a description only)
_r2 = 1.0 - ((np.log(ys) - (c + k * xs)) ** 2).sum() / \
      ((np.log(ys) - np.log(ys).mean()) ** 2).sum()
axc.text(6, 285, "descriptive only (n = %d):\nA ≈ %.1f exp(%.3f V_n)"
         "\nR²_log = %.2f" % (len(xs), np.exp(c), k, _r2),
         fontsize=F.FS(6.0), color="#666666", linespacing=1.45)
axc.annotate("Spearman ρ = +1.00", xy=(0.97, 0.06),
             xycoords="axes fraction", ha="right", va="bottom", fontsize=F.FS(6.6),
             color="#333333")
axc.axhline(0, color="#CCCCCC", lw=0.5, ls=":")
axc.set_xlim(5, 58)
axc.set_ylim(-60, 470)
axc.set_yticks([0, 100, 200, 300, 400])
axc.set_xlabel("Event size  V_n  (m³ of overflow per flooding node)")
axc.set_ylabel("Selected A_ponded (m²)")
axc.grid(axis="y")
F.panel(axc, "c")

# -------------------------------------------- (d) node-level equivalent area
# Estimator: A_eff = V_node / max(hbar_wet, 0.10 m), hbar_wet = mean depth over
# wetted cells in a 5x5 (10 m) window (Section 5.4, item 7); the two-dimensional field is the
# corrected one-way 2026-07-19 solution, not the retired compressed-release run.
# The distribution spans four decades, so the abscissa is logarithmic.
XLO, XHI = 3.0, 3.0e4
bins = np.logspace(np.log10(XLO), np.log10(XHI), 29)
axd.hist(np.clip(EFF, XLO, XHI), bins=bins, color="#B9CDDD",
         edgecolor="#7FA0BC", linewidth=0.35, zorder=3)
axd.set_xscale("log")
_cnt, _ = np.histogram(np.clip(EFF, XLO, XHI), bins=bins)
YTOP = float(_cnt.max()) * 1.34
axd.set_xlim(XLO, XHI)
axd.set_ylim(0, YTOP)
axd.set_xticks([10, 100, 1000, 10000])
axd.set_xticklabels(["10", "100", "1000", "10⁴"])
med = np.median(EFF)
p25, p75 = np.percentile(EFF, [25, 75])
axd.axvline(med, color=F.OI["orange"], lw=1.3, zorder=5)
axd.axvline(100, color=F.OI["vermillion"], lw=1.1, ls="--", zorder=5)
axd.text(med * 1.28, YTOP * 0.78, "median\n%.0f m²" % med,
         fontsize=F.FS(6.2), color=F.OI["orange"], va="top", ha="left",
         linespacing=1.4)
axd.text(100 * 0.78, YTOP * 0.94, "prescribed\n100 m²", fontsize=F.FS(6.2),
         color=F.OI["vermillion"], va="top", ha="right", linespacing=1.4)
axd.annotate("", xy=(p25, YTOP * 0.20), xytext=(p75, YTOP * 0.20),
             arrowprops=dict(arrowstyle="<->", lw=0.6, color="#333333"))
axd.text(np.sqrt(p25 * p75), YTOP * 0.235,
         "IQR %.0f\u2013%.0f m²  \u00B7  p75/p25 %.0f\u00D7"
         % (p25, p75, p75 / p25), fontsize=F.FS(6.0), ha="center", va="bottom",
         color="#333333")
axd.set_xlabel("Node-level equivalent area  V_node/h_wet  (m²)")
axd.set_ylabel("Number of nodes")
axd.grid(axis="y")
axd.text(0.985, 0.985, "2026-07-19  \u00B7  corrected two-dimensional field\nn = %d flooded nodes"
         % EFF.size, transform=axd.transAxes, ha="right", va="top",
         fontsize=F.FS(6.0), color="#555555", linespacing=1.4)
F.panel(axd, "d")

hs = [Line2D([], [], color=F.EVCOL[e], marker=F.EVMK[e], ms=3.6, lw=0.9,
             markeredgecolor="white", markeredgewidth=0.35, label=F.EVLBL[e])
      for e in F.EVENTS]
hs.append(Line2D([], [], color="#333333", marker="*", ms=9, lw=0,
                 markerfacecolor="none", label="selected A_ponded (a) / minimum (b)"))
fig.legend(handles=hs, loc="lower center", bbox_to_anchor=(0.5, 0.002),
           ncol=3, frameon=False, fontsize=F.FS(6.3), handletextpad=0.5,
           columnspacing=2.0, labelspacing=0.4)

F.save(fig, "Fig4_E1_calibration")
