# -*- coding: utf-8 -*-
"""Fig 7 -- the outlet-stress number and the boundary of admissibility.

(a) export efficiency eta and marginal recovery kappa against the dimensionless
    outlet-stress number L = V_in / V_onset (controlled design-storm family);
    the two-sided transition sits at the flooding onset, L = 1.
(b) peak concurrent export at the seven outfalls against the load: it grows by
    25x (the ratio is computed at draw time, so it cannot drift) and never
    plateaus -> the network is resistance-limited, not
    capacity-limited; the municipal nameplate ratings bracket, but do not
    control, the response.
(c) kappa for the two perturbation directions from the same state (5 a,
    120 min): intensity (more rain in the same 2 h) vs duration (same return
    period, longer storm).  kappa is a directional derivative, not a state
    function.
"""
import json
import os
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = os.environ.get("EMS_SWMM_BASE",
                   os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
import figstyle as F

F.use_style(scale=1.15)   # design scale for this figure (see figstyle)
OI = F.OI
WORK = os.path.join(ROOT, "work", "_M6")

D = json.load(open(WORK + "/final.json", encoding="utf-8"))
A = [r for r in D["A"] if r["tag"] != "M6_P2000a"]
onset = D["onset"]
V0, q0 = onset["V"], onset["q"]


def kappa_of(rows):
    out = {}
    for x, y in zip(rows[:-1], rows[1:]):
        dv = y["ww"] - x["ww"]
        out[x["lbl"]] = (1 - (y["flood"] - x["flood"]) / dv) if dv > 0 else None
    return out


kA = kappa_of(D["A"])

# peak export, parsed from peak.log
peaks = {}
for l in open(WORK + "/peak.log", errors="replace"):
    m = re.match(r"\s*(\S+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)", l)
    if m:
        peaks[m.group(1)] = float(m.group(2))

fig = plt.figure(figsize=(F.W2, F.W2 * 0.345))
gs = fig.add_gridspec(1, 3, wspace=0.42, left=0.055, right=0.985,
                      bottom=0.185, top=0.93)
axa, axb, axc = (fig.add_subplot(gs[0, i]) for i in range(3))

# ---------------------------------------------------------------- panel a
L = np.array([r["ww"] / V0 for r in A], float)
eta = np.array([r["eta"] for r in A], float)
kap = np.array([kA.get(r["lbl"], np.nan) for r in A], float)

axa.axvspan(0.4, 1.0, color=OI["green"], alpha=0.10, lw=0)
axa.axvline(1.0, color="k", lw=0.9, ls="--")
axa.text(1.06, 0.055, "flooding\nonset", fontsize=F.FS(6.0), ha="left", va="bottom",
         color="k", linespacing=1.05)
axa.text(0.47, 0.90, "admissible", fontsize=F.FS(6.2), ha="left", va="center",
         color=OI["green"])
axa.plot(L, eta, "-o", color=OI["blue"], ms=2.7, lw=1.0, mfc="white", mew=0.8,
         label="export efficiency η")
axa.plot(L, kap, "-s", color=OI["vermillion"], ms=2.5, lw=1.0, mfc="white", mew=0.8,
         label="marginal recovery κ")
# observed events
OBS = {"0711": (7763, 0.756), "0627": (34422, 0.258), "0730": (57110, 0.419),
       "0719": (95540, 0.240), "0809": (117982, 0.386)}
for k, (v, e) in OBS.items():
    axa.plot(v / V0, e, marker="D", ms=3.0, color="k", mfc="none", mew=0.7, zorder=5)
axa.plot([], [], "D", ms=3.0, color="k", mfc="none", mew=0.7, label="gauged events (\u03b7)")
axa.set_xscale("log")
F.log_ticks(axa, "x")          # plain-text 10^n labels, never mathtext
axa.set_xlim(0.42, 30)
axa.set_ylim(-0.03, 1.09)
axa.set_xlabel("outlet-stress number  Λ = V_in/V_onset")
axa.set_ylabel("dimensionless response  (-)")
axa.legend(loc="upper right", frameon=False, fontsize=F.FS(5.9), handlelength=1.5,
           borderaxespad=0.2, labelspacing=0.25)
F.panel(axa, "a", dx=-0.19)

# ---------------------------------------------------------------- panel b
Vi = np.array([r["ww"] for r in A], float)
pk = np.array([peaks.get(r["lbl"], np.nan) for r in A], float)
m = ~np.isnan(pk)
axb.plot(Vi[m], pk[m], "-o", color=OI["purple"], ms=2.7, lw=1.0, mfc="white", mew=0.8)
q1, q10 = 5 * 0.293 + 0.789, 5 * 0.926 + 0.789
axb.axhline(q1, color=OI["grey"], lw=0.9, ls=":")
axb.axhline(q10, color=OI["grey"], lw=0.9, ls=":")
axb.text(2450, q10 * 1.10, "municipal capacity\nof the 5 gravity outlets\n(dn700, 1-10 permille)",
         fontsize=F.FS(5.4), color=OI["grey"], va="bottom", ha="left", linespacing=1.15)
axb.set_xscale("log")
axb.set_yscale("log")
axb.set_xlabel("event load  V_in  (m³)")
axb.set_ylabel("peak campus export  (m³/s)")
axb.set_xlim(2200, 2.0e5)
axb.set_ylim(0.16, 30.0)
axb.set_yticks([0.2, 0.5, 1, 2, 5, 10])
axb.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
# matplotlib's own log formatter would set the exponent in mathtext at 0.7x the
# base size (4.9 pt at a 7 pt base) -- below the journal floor.
F.log_ticks(axb, "x")
# the ratio is computed, not hard-coded, so it cannot drift when the
# probe members change
_ratio = pk[m][-1] / pk[m][0]
axb.annotate("%.0f×, no plateau" % _ratio,
             xy=(Vi[m][-1] * 0.72, pk[m][-1] * 0.86), xytext=(3.2e3, 0.42),
             fontsize=F.FS(6.0), color=OI["purple"],
             arrowprops=dict(arrowstyle="-|>", lw=0.7, color=OI["purple"],
                             shrinkA=1, shrinkB=2,
                             connectionstyle="arc3,rad=0.18"))
F.panel(axb, "b", dx=-0.19)

# ---------------------------------------------------------------- panel c
seg_int = [("3→5", kA.get("3a")), ("5→10", kA.get("5a")),
           ("10→20", kA.get("10a"))]
kB = kappa_of(D["B"])
seg_dur = [("30→60", kB.get("5a/30min")), ("60→120", kB.get("5a/60min")),
           ("120→240", kB.get("5a/120min")), ("240→360", kB.get("5a/240min"))]
pos_d = [0, 1.0, 2.0, 3.0]
pos_i = [4.4, 5.4, 6.4]
axc.bar(pos_d, [v for _, v in seg_dur], width=0.68, color=OI["sky"],
        edgecolor="k", lw=0.4)
axc.bar(pos_i, [v for _, v in seg_int], width=0.68, color=OI["vermillion"],
        edgecolor="k", lw=0.4)
axc.set_xticks(pos_d + pos_i)
# 90 degrees, not 40: at 7 pt a label like "120→240" is ~8.6 mm long, while the
# categories are 6.6 mm apart, so a slanted label overlaps its neighbour.
axc.set_xticklabels([l for l, _ in seg_dur] + [l for l, _ in seg_int],
                    rotation=90, ha="center", va="top", fontsize=F.FS(6.1))
axc.set_ylabel("marginal recovery  κ  (-)")
axc.set_ylim(0, 0.80)
axc.set_xlim(-0.75, 7.15)
axc.axvline(3.7, color=F.OI["lgrey"], lw=0.8)
axc.text(1.5, 0.66, "longer storm", fontsize=F.FS(5.9), color=OI["sky"],
         ha="center", va="bottom")
axc.text(1.5, 0.615, "T↑, same 5 a", fontsize=F.FS(5.4), color=OI["sky"],
         ha="center", va="bottom")
axc.text(5.4, 0.66, "heavier storm", fontsize=F.FS(5.9), color=OI["vermillion"],
         ha="center", va="bottom")
axc.text(5.4, 0.615, "T = 2 h, a↑", fontsize=F.FS(5.4), color=OI["vermillion"],
         ha="center", va="bottom")
F.panel(axc, "c", dx=-0.24)

F.save(fig, "Fig7_outlet_authority")
