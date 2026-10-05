# -*- coding: utf-8 -*-
"""Figure 4 -- E2: how much observation does the ponding area actually need?

(a) deployed sensor records vs usable observations inside the evaluation window
(b) agreement with the full-set selection as a function of the station count K
(c) degradation under observation thinning (selection flip vs sampling interval)
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
from matplotlib.patches import Patch

import figstyle as F

F.use_style(scale=1.15)   # design scale for this figure (see figstyle)
E = json.load(open(F.WORK + "/e2_obs_sufficiency.json", encoding="utf-8"))

fig = plt.figure(figsize=(F.W2, 78 * F.MM))
gs = fig.add_gridspec(1, 3, left=0.052, right=0.988, bottom=0.30, top=0.90,
                      wspace=0.46)
axa, axb, axc = (fig.add_subplot(gs[0, i]) for i in range(3))

# ------------------------------------------- (a) deployed vs usable stations
DEPLOY = {"20250627": 11, "20250730": 8, "20260719": 12, "20260809": 12}
USE = {k: len(v) for k, v in E["usable"].items()}
x = np.arange(4)
w = 0.34
dep = [DEPLOY[e] for e in F.EVENTS]
use = [USE[e] for e in F.EVENTS]
axa.bar(x - w / 2, dep, w, facecolor="white", edgecolor="#8C8C8C",
        linewidth=0.6, hatch="///", label="records in the dataset", zorder=3)
axa.bar(x + w / 2, use, w, facecolor=F.OI["blue"], edgecolor="white",
        linewidth=0.4, label="usable in the evaluation window", zorder=3)
for xi, (a_, b_) in enumerate(zip(dep, use)):
    axa.text(xi - w / 2, a_ + 0.25, "%d" % a_, ha="center", va="bottom",
             fontsize=F.FS(6.0), color="#666666")
    axa.text(xi + w / 2, b_ + 0.25, "%d" % b_, ha="center", va="bottom",
             fontsize=F.FS(6.0), fontweight="bold", color=F.OI["blue"])
axa.set_xticks(x)
axa.set_xticklabels([F.EVSHORT[e][5:] for e in F.EVENTS], fontsize=F.FS(6.2))
axa.set_ylim(0, 15.2)
axa.set_yticks([0, 4, 8, 12])
axa.set_ylabel("Ultrasonic sensors")
axa.grid(axis="y")
F.panel(axa, "a", dx=-0.22)
axa.annotate("", xy=(2.17, 3.2), xytext=(2.17, 11.4),
             arrowprops=dict(arrowstyle="<->", lw=0.6, color=F.OI["vermillion"]))
axa.text(2.28, 7.3, "\u221275%", fontsize=F.FS(6.0), color=F.OI["vermillion"],
         ha="left", va="center")

# --------------------------------------------- (b) station-count sufficiency
for ev in F.EVENTS:
    d = E["thin_k"][ev]
    ks = [r[0] for r in d]
    fr = [r[1] for r in d]
    axb.plot(ks, fr, "-", lw=0.9, color=F.EVCOL[ev], marker=F.EVMK[ev],
             ms=3.2, markeredgecolor="white", markeredgewidth=0.3,
             label=F.EVSHORT[ev], zorder=3)
axb.axhline(80, color="#8C8C8C", lw=0.6, ls="--", zorder=1)
axb.text(10.6, 82, "80 % agreement", fontsize=F.FS(5.9), color="#666666",
         ha="right", va="bottom")
axb.set_xlim(0.4, 11.6)
axb.set_ylim(25, 108)
axb.set_xticks([1, 2, 3, 4, 5, 6, 8, 10, 11])
axb.set_yticks([25, 50, 75, 100])
axb.set_xlabel("Stations used, K (\u2013)")
axb.set_ylabel("Agreement with the full-set\nselection (%)")
axb.grid(axis="y")
axb.legend(loc="lower right", ncol=1, fontsize=F.FS(5.9), labelspacing=0.3,
           columnspacing=1.0)
F.panel(axb, "b", dx=-0.20)

# ------------------------------------------------------ (c) observation thinning
for ev in F.EVENTS:
    d = E["thin_t"][ev]
    tt = [r[0] for r in d]
    rr = [r[1] for r in d]
    ok = [r[3] for r in d]
    axc.plot(tt, rr, "-", lw=0.9, color=F.EVCOL[ev], zorder=3)
    for t_, r_, o_ in zip(tt, rr, ok):
        axc.plot([t_], [r_], marker=F.EVMK[ev] if o_ else "X",
                 ms=4.4 if o_ else 5.2, color=F.EVCOL[ev],
                 markerfacecolor=F.EVCOL[ev] if o_ else "white",
                 markeredgecolor=F.EVCOL[ev], markeredgewidth=0.6,
                 zorder=4, linestyle="none")
axc.axhline(1.0, color="#CCCCCC", lw=0.5, ls=":", zorder=1)
axc.set_xscale("log")
axc.set_xticks([10, 20, 30, 60])
axc.set_xticklabels(["10", "20", "30", "60"])
axc.set_xlim(9, 68)
axc.xaxis.set_minor_locator(plt.NullLocator())
axc.xaxis.set_major_formatter(plt.FuncFormatter(lambda v, p: "%g" % v))
axc.set_ylim(0.72, 1.34)
axc.set_xlabel("Observation interval (min)")
axc.set_ylabel("RMSE relative to the 10 min series (\u2013)")
axc.grid(axis="y")
axc.text(0.03, 0.17, "filled: A_ponded unchanged\ncross: A_ponded flips",
         transform=axc.transAxes, fontsize=F.FS(5.9), color="#555555",
         ha="left", va="top", linespacing=1.45,
         # Bottom band, left: nothing is plotted below 0.84 (the lowest point on
         # any of the four curves), so the key sits in clear space instead of on
         # the 2025-06-27 curve, which leaves the panel through the top-right
         # corner.  Two SHORT lines so the block fits a 33 mm panel; the wording
         # matches the caption ("unchanged" / "a flip").
         bbox=dict(boxstyle="round,pad=0.18", facecolor="white",
                   edgecolor="none", alpha=0.86))
F.panel(axc, "c", dx=-0.20)

hs = [Line2D([], [], color=F.EVCOL[e], marker=F.EVMK[e], ms=3.4, lw=0.9,
             markeredgecolor="white", markeredgewidth=0.3,
             label=F.EVLBL[e]) for e in F.EVENTS]
hs.append(Patch(facecolor="white", edgecolor="#8C8C8C", hatch="///",
                label="sensor records in the dataset"))
hs.append(Patch(facecolor=F.OI["blue"], edgecolor="white",
                label="usable in the window"))
fig.legend(handles=hs, loc="lower center", bbox_to_anchor=(0.5, 0.005),
           ncol=3, frameon=False, fontsize=F.FS(6.3), handletextpad=0.6,
           columnspacing=2.0, labelspacing=0.45)

F.save(fig, "Fig5_E2_identifiability")
