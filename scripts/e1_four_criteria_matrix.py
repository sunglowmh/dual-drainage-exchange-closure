# -*- coding: utf-8 -*-
"""P2-4 -- do the M2 claims survive the densified grid?

Reads work/fig_e1_dense.json (9 Aponded levels x 4 events) and re-checks:

  Claim A  no single value is transferable -> the transfer-cost matrix
           (apply one event's optimum to another), for the composite score and
           for the three component metrics;
  Claim B  the mid-segment optimum is not identifiable -> the width of the
           near-optimal set at a declared tolerance, with the number of
           in-window stations beside it;
  Claim C  identifiability is high at the two ends and low in the middle.

Also reports, per event, the optimum under each of four equally defensible
criteria, because a "shift" that appears under only one criterion is a property
of the criterion, not of the system.
"""
import json

import numpy as np

WORK = r"E:/work/tongji/swmm_model/work"
D = json.load(open(WORK + "/fig_e1_dense.json", encoding="utf-8"))
EVL = {"20250627": "0627 63.6 mm/0.8 h", "20250730": "0730 125.3 mm/13.2 h",
       "20260719": "0719 178.6 mm/5.2 h", "20260809": "0809 242.5 mm/36 h"}
SHORT = {k: v.split()[0] for k, v in EVL.items()}
COMPLETE = [e for e in EVL if len(D.get(e, {}).get("xs", [])) == 9]
print("完整 9 档的事件：%s" % ", ".join(SHORT[e] for e in COMPLETE))
if len(COMPLETE) < 4:
    print("⚠ 只有 %d/4 场事件齐了 9 档，下面只对完整的事件作结论" % len(COMPLETE))

print()
print("=" * 104)
print("Per-event optimum under four equally defensible criteria")
print("=" * 104)
print("%-22s %6s %8s %10s %10s %8s %10s %10s" %
      ("event", "n_sta", "V/node", "composite", "min rmse", "max r", "min |dT|", "#distinct"))
opt = {}
for ev in COMPLETE:
    r = D[ev]
    xs, per = r["xs"], r["per"]
    keys = dict(composite=[per[str(x)]["score"] for x in xs],
                rmse=[per[str(x)]["rmse"] for x in xs],
                r=[-per[str(x)]["r"] for x in xs],
                dt=[per[str(x)]["dt"] for x in xs])
    o = {k: xs[int(np.argmin(v))] for k, v in keys.items()}
    opt[ev] = o
    print("%-22s %6d %8.1f %10d %10d %8d %10d %10d"
          % (SHORT[ev], r["n_station"], r["v_per_node"], o["composite"], o["rmse"],
             o["r"], o["dt"], len(set(o.values()))))

print()
print("=" * 104)
# The tolerance is an ABSOLUTE increment of the normalised composite score
# (which lies in [0,1] within each event), i.e. "+0.05" is 5 % of the SCORE
# SCALE -- the "+5 %" quoted in the manuscript -- and NOT "5 % of the best
# value".  The earlier column headers said "+2 %"/"+5 %" without saying which.
print("Claim B / C -- near-optimal SET width at a declared score tolerance")
print("=" * 104)
print("%-22s %6s | %-30s | %-30s"
      % ("event", "n_sta", "within +0.02 (2 % of scale)",
         "within +0.05 (5 % of scale, = paper)"))
for ev in COMPLETE:
    r = D[ev]
    xs, per = r["xs"], r["per"]
    sc = np.array([per[str(x)]["score"] for x in xs])
    best = sc.min()
    cells = []
    for tol in (0.02, 0.05):
        sel = [x for x in xs if per[str(x)]["score"] <= best + tol]
        cells.append("%s (=%d 档, 宽 %d m²)" % (str(sel), len(sel), max(sel) - min(sel)))
    print("%-22s %6d | %-30s | %-30s" % (SHORT[ev], r["n_station"], cells[0], cells[1]))

print()
print("=" * 104)
print("Claim A -- transfer cost (row: whose optimum is applied; col: event scored)")
print("=" * 104)
for key in ("score", "rmse", "r", "dt"):
    print()
    print("  %s：" % {"score": "composite score", "rmse": "bias-corrected RMSE",
                   "r": "correlation r", "dt": "peak-timing |dT|"}[key])
    print("    %-8s%s" % ("apply\\", "".join("%14s" % SHORT[e] for e in COMPLETE)))
    for a in COMPLETE:
        row = "    %-8s" % SHORT[a]
        for b in COMPLETE:
            xa = opt[a]["composite"] if key == "score" else opt[a][
                {"rmse": "rmse", "r": "r", "dt": "dt"}[key]]
            kk = {"score": "score", "rmse": "rmse", "r": "r", "dt": "dt"}[key]
            vb = D[b]["per"][str(xa)][kk]
            vbest = D[b]["per"][str(opt[b]["composite"] if key == "score" else opt[b][kk])][kk]
            if key == "r":
                row += "%14s" % ("%+.3f" % (vb - vbest))
            else:
                row += "%14s" % ("%+.0f %%" % (100 * (vb / vbest - 1)) if abs(vbest) > 1e-9 else "-")
        print(row)

print()
print("=" * 104)
print("Claim C -- the two ends against the middle")
print("=" * 104)
order = sorted(COMPLETE, key=lambda e: D[e]["v_per_node"])
for ev in order:
    r = D[ev]
    xs, per = r["xs"], r["per"]
    sc = np.array([per[str(x)]["score"] for x in xs])
    best = sc.min()
    sel = [x for x in xs if per[str(x)]["score"] <= best + 0.05]
    print("  %-6s V/node %6.1f | n_sta %2d | 最优 %3d (%s) | +5 %% 近优集 %s"
          % (SHORT[ev], r["v_per_node"], r["n_station"], r["grid_opt"], r["kind"],
             "%s (宽 %d)" % (sel, max(sel) - min(sel))))
