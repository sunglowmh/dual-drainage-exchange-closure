# -*- coding: utf-8 -*-
"""Figure 6 -- the 2026-07-19 storm: observation vs model at the three
stations that survived inside the evaluation window.

Draws, per station, the bias-corrected observation (valid samples solid,
frozen / full-shaft-truncated samples open and excluded from the metrics),
the model with and without the ponding area, the rain hyetograph and the
alarm-induced telemetry gaps (informative censoring).
"""
import datetime as dt
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
from metrics_lib import (ME, frozen_mask, load_sim, metrics, event_window,
                         obs_manhole_series)

F.use_style(scale=1.15)   # design scale for this figure (see figstyle)

EV = "20260719"
CODES = ["867948079409621", "867948079405843", "867948079403574"]
TITLE = {"867948079409621": "S11  Tongwen Building (gate)",
         "867948079405843": "S02  Structural Laboratory (east)",
         "867948079403574": "S01  Chifeng Rd 200 (west)"}
GAPS = {"867948079403574": ("13:56", "20:06"), "867948079405843": ("15:03", "21:11")}

obs = json.load(open("%s/obs_%s.json" % (ME, EV), encoding="utf-8"))
simA = load_sim(EV, "A100")
sim0 = load_sim(EV, "base")
t_lo, t_hi = event_window(obs, sim0)

fig = plt.figure(figsize=(F.W2, 138 * F.MM))
gs = fig.add_gridspec(4, 1, left=0.095, right=0.975, bottom=0.135, top=0.885,
                      hspace=0.18, height_ratios=[1.35, 2.5, 2.5, 2.5])
axr = fig.add_subplot(gs[0])
axes = [fig.add_subplot(gs[i], sharex=axr) for i in (1, 2, 3)]

X0 = dt.datetime(2026, 7, 19, 10, 0)
X1 = dt.datetime(2026, 7, 20, 6, 0)

# ------------------------------------------------------------------- rainfall
rt = [dt.datetime.fromisoformat(t) for t, _ in obs["rain"]]
rv = np.array([v * 60.0 / obs["dt_min"] for _, v in obs["rain"]])
axr.bar(rt, rv, width=1.0 / (24 * 60), color="#8FBEDA", linewidth=0, zorder=3)
axr.set_ylim(max(rv) * 1.15, 0)
# The rainfall strip is only ~14 mm tall and the 1-min peak is 192 mm/h, so a
# hand-picked 0/25/50 puts "25" and "50" 1.6 mm apart -- on top of each other.
# Take the step from the range instead.
_rtop = float(max(rv)) * 1.15
axr.set_yticks(np.arange(0.0, _rtop + 1e-9, F.nice_step(_rtop / 2.5)))
axr.set_ylabel("Rainfall\n(mm/h)", fontsize=F.FS(6.6))
axr.tick_params(labelbottom=False)
axr.grid(axis="y")
F.despine(axr)
axr.text(0.012, 0.86, "Campus rain gauge (1 min)", transform=axr.transAxes,
         fontsize=F.FS(6.2), color="#3E6E8C", va="top")
F.panel(axr, "a", dx=-0.095, dy=1.02)

# ---------------------------------------------------------------- stations
osr = obs_manhole_series(obs, t_lo, t_hi)
for k, (ax, code) in enumerate(zip(axes, CODES)):
    s = osr[code]
    ot = [t for t, _ in s]
    ov = np.array([v for _, v in s])
    keep = np.array(frozen_mask([v for _, v in s]))
    mt = [dt.datetime.fromisoformat(p[0]) for p in simA["nodes"][code]]
    mv = np.array([p[1] for p in simA["nodes"][code]])
    m0 = np.array([p[1] for p in sim0["nodes"][code]])
    mm = metrics(s, list(zip(mt, mv)), t_lo, t_hi)
    off = mm["offset"]      # metrics() scores (model - off) against the raw observation,
    oc = ov + off           # i.e. the observation is re-referenced to the model datum

    # alarm-driven telemetry gaps
    if code in GAPS:
        g0 = dt.datetime.fromisoformat("2026-07-19T" + GAPS[code][0])
        g1 = dt.datetime.fromisoformat("2026-07-19T" + GAPS[code][1])
        ax.axvspan(g0, g1, facecolor="#F2DEDA", edgecolor="none", zorder=0)
        ax.text((g0 + (g1 - g0) * 0.5), 0.93, "alarm gap",
                transform=ax.get_xaxis_transform(), ha="center", va="top",
                fontsize=F.FS(5.9), color="#B4573F")

    ax.plot(mt, m0, lw=0.6, color="#A9A9A9", ls="--", zorder=3)
    ax.plot(mt, mv, "-", lw=0.9, color=F.OI["vermillion"], zorder=4)
    ax.plot(ot, oc, "-", lw=0.7, color="#B9CFE0", zorder=2)
    ax.plot([t for t, k_ in zip(ot, keep) if not k_],
            oc[~keep], "o", ms=3.0, mfc="white", mec=F.OI["blue"],
            mew=0.5, zorder=6)
    ax.plot([t for t, k_ in zip(ot, keep) if k_], oc[keep], "o", ms=3.2,
            color=F.OI["blue"], zorder=6)

    ax.set_ylabel("Water depth above\nthe well invert (m)", fontsize=F.FS(6.6))
    ax.grid(axis="y")
    ax.set_xlim(X0, X1)
    ymax = max(np.nanmax(mv), np.nanmax(oc)) * 1.30
    ax.set_ylim(0, ymax)
    ax.text(0.012, 0.94, TITLE[code], transform=ax.transAxes, fontsize=F.FS(6.8),
            fontweight="bold", va="top")
    ax.text(0.985, 0.94,
            "r = %.2f   RMSE = %.2f m" % (mm["r"], mm["rmse"]),
            transform=ax.transAxes, fontsize=F.FS(6.3), va="top", ha="right",
            color="#555555")
    if code == "867948079403574":
        ax.annotate("post-gap plateau:\nthe observed peak is a\nlower bound",
                    xy=(dt.datetime(2026, 7, 19, 21, 0), 2.45),
                    xytext=(dt.datetime(2026, 7, 20, 0, 45), 1.45),
                    fontsize=F.FS(5.9), color="#7A4A2A", linespacing=1.45,
                    arrowprops=dict(arrowstyle="-|>", lw=0.6, color="#7A4A2A",
                                    shrinkA=2, shrinkB=2))
    F.panel(ax, "bcd"[k], dx=-0.095, dy=1.03)

axes[-1].set_xlabel("Time (local, 2026-07-19 to 07-20)")
for ax in axes[:-1]:
    ax.tick_params(labelbottom=False)

hs = [Line2D([], [], color=F.OI["blue"], marker="o", ms=3.2, lw=0.7,
             label="Observed (re-referenced)"),
      Line2D([], [], color=F.OI["blue"], marker="o", ms=3.2, lw=0,
             mfc="white", label="Frozen / truncated (\u2265 4 identical values)"),
      Line2D([], [], color=F.OI["vermillion"], lw=0.9,
             label="A_ponded = 100 m² (selected)"),
      Line2D([], [], color="#A9A9A9", lw=0.7, ls="--",
             label="no ponding area"),
      Patch(facecolor="#F2DEDA", label="Alarm-induced telemetry gap"),
      Patch(facecolor="#8FBEDA", label="Rainfall (right-hand bars)")]
fig.legend(handles=hs, loc="lower center", bbox_to_anchor=(0.5, 0.002),
           ncol=3, frameon=False, fontsize=F.FS(6.3), handletextpad=0.6,
           columnspacing=2.0, labelspacing=0.45)

F.save(fig, "Fig3_hydrographs_0719")
