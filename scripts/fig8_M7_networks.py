# -*- coding: utf-8 -*-
"""Fig 8 -- the kappa(Lambda) structure reproduces on three conceptual networks.

(a) marginal export recovery kappa against the outlet-stress number
    Lambda = V_in/V_onset for the gauged campus and for three rule-generated
    conceptual networks (grid / random tree / sparse flat), each driven by the
    same 16+2-member Chicago load family and analysed with the same extraction;
    the transition sits at the flooding onset (Lambda ~ 1) in every system.
(b-d) the three conceptual layouts (junction dots, conduits, single outfall).
"""
import json
import sys
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.environ.get("EMS_SWMM_BASE",
                   os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
import figstyle as F
import _M7_nets

F.use_style(scale=1.15)
OI = F.OI

M6 = json.load(open(os.path.join(ROOT, "work", "_M6", "final.json"), encoding="utf-8"))
M7 = json.load(open(os.path.join(ROOT, "work", "_M7", "final.json"), encoding="utf-8"))

V0_CAMPUS = M6["onset"]["V"]


def seg_curve(rows, V0):
    """(Lambda_mid, kappa) segment list from a family's rows."""
    out = []
    for x, y in zip(rows[:-1], rows[1:]):
        dv = y["ww"] - x["ww"]
        if dv <= 0:
            continue
        kap = 1.0 - (y["flood"] - x["flood"]) / dv
        out.append((0.5 * (x["ww"] + y["ww"]) / V0, kap))
    return out


campus = seg_curve(M6["A"], V0_CAMPUS)
NETCOL = {"N1": OI["blue"], "N2": OI["vermillion"], "N3": OI["green"]}
NETLBL = {"N1": "N1 grid", "N2": "N2 random tree", "N3": "N3 sparse flat"}

fig = plt.figure(figsize=(F.W2, F.W2 * 0.34))
gs = fig.add_gridspec(1, 2, width_ratios=[1.45, 1.0], wspace=0.34,
                     left=0.065, right=0.985, bottom=0.20, top=0.92)
axa = fig.add_subplot(gs[0, 0])
gsr = gs[0, 1].subgridspec(3, 1, hspace=0.62)

# ---------------------------------------------------------------- (a) kappa-Lambda
axa.axvline(1.0, color=OI["black"], lw=0.7, ls="--", zorder=1)
axa.axhspan(0, 1, xmin=0, xmax=1, color=OI["lgrey"], alpha=0.0, zorder=0)
cx, cy = zip(*campus)
axa.plot(cx, cy, "-", color=OI["grey"], lw=1.0, marker="o", ms=2.4,
         mfc="white", mew=0.7, label="Campus (3 277 nodes)", zorder=3)
for nm in ("N1", "N2", "N3"):
    d = M7[nm]
    pts = seg_curve(d["rows"], d["onset"]["V"])
    x, y = zip(*pts)
    axa.plot(x, y, "-", color=NETCOL[nm], lw=1.0, marker="o", ms=2.4,
             mfc="white", mew=0.7, label="%s (%d nodes)"
             % (NETLBL[nm], d["spec"]["junctions"]), zorder=4)
axa.set_xscale("log")
F.log_ticks(axa, "x")
axa.set_xlim(0.08, 30)
axa.set_ylim(0, 1.02)
axa.set_xlabel("Outlet-stress number Λ = V_in/V_onset")
axa.set_ylabel("Marginal export recovery κ")
axa.annotate("flooding onset (Λ = 1)", xy=(0.995, 0.35), xytext=(0.90, 0.35),
             ha="right", va="center", fontsize=F.FS(7), color=OI["black"],
             arrowprops=dict(arrowstyle="-", lw=0.6, color=OI["black"]))
axa.legend(loc="upper right", frameon=False, handlelength=1.6, borderaxespad=0.2)
F.panel(axa, "a", dx=-0.155)

# ---------------------------------------------------------------- (b-d) layouts
for i, nm in enumerate(("N1", "N2", "N3")):
    ax = fig.add_subplot(gsr[i])
    net = {"N1": _M7_nets.build_N1, "N2": _M7_nets.build_N2,
           "N3": _M7_nets.build_N3}[nm]()
    for _cn, a, b, _l in net.edges:
        xa, ya = net.xy[a]
        xb, yb = net.xy[b] if b in net.xy else (0.0, 0.0)
        ax.plot([xa, xb], [ya, yb], "-", color=OI["grey"], lw=0.4, zorder=2)
    xs = [net.xy[j][0] for j in net.junc]
    ys = [net.xy[j][1] for j in net.junc]
    ax.plot(xs, ys, "o", color=NETCOL[nm], ms=1.7, mew=0, zorder=3)
    ax.plot([0.0], [0.0], marker="v", color=OI["black"], ms=3.2, zorder=4)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlim(min(xs) - 60, max(xs) + 60)
    ax.set_ylim(min(ys) - 60, max(ys) + 60)
    for s in ax.spines.values():
        s.set_linewidth(0.5)
    s = M7[nm]["spec"]
    ax.set_title("%s: %d nodes, %.1f ha\nV_onset = %s m³"
                 % (NETLBL[nm], s["junctions"], s["area_ha"],
                    "{:,}".format(round(M7[nm]["onset"]["V"])).replace(",", " ")),
                 fontsize=F.FS(6.6), loc="left", pad=2.5)
    F.panel(ax, "bcd"[i], dx=-0.14, dy=1.16)

F.save(fig, "Fig8_M7_networks")
