# -*- coding: utf-8 -*-
"""M6 -- final consolidation: the outlet-stress number against the response.

Three families, all on the SAME network (byte-identical JUNCTIONS / CONDUITS /
XSECTIONS / SUBCATCHMENTS / SUBAREAS / PUMPS / STORAGE / OUTFALLS sections):

  A  intensity family : Chicago 120 min, return period 0.1 - 2000 a
  B  duration family  : return period 5 a, duration 30 - 360 min
  C  observed events  : the five gauged storms used in the paper
  D  outlet variants  : the same event with the outlet conduits / pumps rescaled

Reported for each: V_in, V_flood, V_export, the export efficiency
eta = V_export/(V_export+V_flood), the marginal flooding share k = dV_flood/dV_in
and the marginal recovery kappa = 1 - k = dV_export/dV_in, plus two candidate
dimensionless numbers
    L_rate = (V_in / T) / q_onset          (rate ratio)
    L_vol  = V_in / V_onset                (volume ratio)
with the onset taken from family A.
"""
import json
import os
import re
import sys

import numpy as np

ROOT = os.environ.get("EMS_SWMM_BASE",
                   os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from _M6_rain import rain_spec

OUT = os.path.join(ROOT, "output")
WORK = os.path.join(ROOT, "work", "_M6")


def block(t, title):
    i = t.find(title)
    if i < 0:
        return ""
    j = t.find("\n  *", i + 200)
    return t[i:j if j > 0 else len(t)]


def gv(b, lbl):
    m = re.search(re.escape(lbl) + r"[. ]*(-?\d+\.\d+)\s+(-?\d+\.\d+)", b)
    return float(m.group(2)) * 1000.0 if m else None


def _pick(tag, ext, min_size=0):
    for p in ("%s/swmm_%s.%s" % (OUT, tag, ext), "%s/%s.%s" % (WORK, tag, ext),
              "%s/%s.%s" % (WORK, tag.replace("M6_", ""), ext)):
        if os.path.exists(p) and os.path.getsize(p) >= min_size:
            return p
    return None


def rpt_path(inp_tag):
    return _pick(inp_tag, "rpt", 300000)


def inp_path(inp_tag):
    return _pick(inp_tag, "inp")


A = [("0.1a", "M6_P0.1a"), ("0.15a", "M6_P0.15a"), ("0.2a", "M6_P0.2a"),
     ("0.25a", "M6_P0.25a"), ("0.5a", "M6_P0.5a"), ("1a", "M6_P1a"),
     ("2a", "v7_2a"), ("3a", "M6_P3a"), ("5a", "v7_5a"), ("10a", "M6_P10a"),
     ("20a", "M6_P20a"), ("50a", "v7_50a"), ("100a", "M6_P100a"),
     ("x1.2", "M6_L120"), ("x1.4", "M6_L140"), ("x1.6", "M6_L160")]
B = [("5a/30min", "M6_D30"), ("5a/60min", "M6_D60"), ("5a/120min", "v7_5a"),
     ("5a/240min", "M6_D240"), ("5a/360min", "M6_D360")]
C = [("0711", "20260711_base"), ("0627", "20250627_base"), ("0730", "20250730_base"),
     ("0719", "20260719_base"), ("0809", "20260809_base")]
D = [("0719 out x0.5", "20260719_cap050_base"), ("0719 base", "20260719_base"),
     ("0719 out x2", "20260719_cap200_base"), ("0719 pump x0.5", "20260719_ps050"),
     ("0719 pump x1.5", "20260719_ps150"),
     ("0719 tw0.30", "20260719_tw030_base"), ("0719 tw0.60", "20260719_tw060_base")]


def load(tag):
    r = rpt_path(tag)
    if not r:
        return None
    t = open(r, errors="replace").read()
    b = block(t, "Flow Routing Continuity")
    if not b:
        return None
    d = dict(tag=tag, ww=gv(b, "Wet Weather Inflow") or 0.0,
             flood=gv(b, "Flooding Loss") or 0.0,
             export=gv(b, "External Outflow") or 0.0)
    i = inp_path(tag)
    rs = rain_spec(i) if i else None
    d.update(depth=(rs or {}).get("depth_mm", float("nan")),
             dur_h=(rs or {}).get("dur_h", float("nan")),
             peak=(rs or {}).get("peak", float("nan")))
    d["eta"] = d["export"] / (d["export"] + d["flood"]) if (d["export"] + d["flood"]) > 0 else float("nan")
    d["qbar"] = d["ww"] / (d["dur_h"] * 3600) if d["dur_h"] else float("nan")
    return d


def seg(rows):
    """marginal shares between consecutive members"""
    out = []
    for x, y in zip(rows[:-1], rows[1:]):
        dv = y["ww"] - x["ww"]
        if dv <= 0:
            out.append((x["lbl"], None, None))
            continue
        k = (y["flood"] - x["flood"]) / dv
        out.append((x["lbl"] + "->" + y["lbl"], k, 1 - k))
    return out


def build(fam, title):
    rows = []
    for lbl, tag in fam:
        d = load(tag)
        if d:
            d["lbl"] = lbl
            rows.append(d)
    print("### %s" % title)
    print()
    print("| %-12s | depth mm | T h | V_in m3 | V_flood m3 | V_export m3 | eta | q_in m3/s |"
          % "member")
    print("|---|---|---|---|---|---|---|---|")
    for d in rows:
        print("| %-12s | %8.1f | %5.2f | %8.0f | %9.0f | %9.0f | %.3f | %8.3f |"
              % (d["lbl"], d["depth"], d["dur_h"], d["ww"], d["flood"], d["export"],
                 d["eta"], d["qbar"]))
    print()
    print("| segment | dV_in | dV_flood | dV_export | k = dF/dV | kappa = 1-k | sum |")
    print("|---|---|---|---|---|---|---|")
    for i, (x, y) in enumerate(zip(rows[:-1], rows[1:])):
        dv = y["ww"] - x["ww"]
        df = y["flood"] - x["flood"]
        de = y["export"] - x["export"]
        if dv <= 0:
            continue
        print("| %s->%s | %.0f | %.0f | %.0f | %.3f | %.3f | %.3f |"
              % (x["lbl"], y["lbl"], dv, df, de, df / dv, 1 - df / dv,
                 (df + de) / dv))
    print()
    return rows


print("# M6 -- response tables")
print()
rA = build(A, "A. Intensity family (Chicago 120 min, return period 0.1-2000 a)")
rB = build(B, "B. Duration family (5 a, duration 30-360 min)")
rC = build(C, "C. Observed events")
rD = build(D, "D. Outlet / pump / tailwater variants at 0719")

# ---- onset from family A --------------------------------------------------
Vi = np.array([d["ww"] for d in rA], float)
Ff = np.array([d["flood"] for d in rA], float)
fr = Ff / Vi
onset = None
for i in range(len(Vi) - 1):
    if fr[i] < 0.01 <= fr[i + 1]:
        w = (0.01 - fr[i]) / (fr[i + 1] - fr[i])
        onset = dict(V=float(Vi[i] + w * (Vi[i + 1] - Vi[i])),
                     q=float((Vi[i] + w * (Vi[i + 1] - Vi[i])) / (rA[i]["dur_h"] * 3600)))
        break
print("## Surcharge onset (flooding = 1 %% of V_in), interpolated in family A")
print()
print("  V_onset = %.0f m3 over a 120 min storm   ->   q_onset = %.3f m3/s"
      % (onset["V"], onset["q"]))
print()

q0, V0 = onset["q"], onset["V"]
print("## Dimensionless numbers")
print()
print("| family | member | L_rate = q_in/q_onset | L_vol = V_in/V_onset | eta | kappa |")
print("|---|---|---|---|---|---|")
for rows, nm in ((rA, "A"), (rB, "B"), (rC, "C"), (rD, "D")):
    ks = {}
    for i, (x, y) in enumerate(zip(rows[:-1], rows[1:])):
        dv = y["ww"] - x["ww"]
        ks[x["lbl"]] = (1 - (y["flood"] - x["flood"]) / dv) if dv > 0 else None
    for d in rows:
        print("| %s | %-12s | %10.3f | %10.3f | %.3f | %s |"
              % (nm, d["lbl"], d["qbar"] / q0 if q0 else float("nan"),
                 d["ww"] / V0 if V0 else float("nan"), d["eta"],
                 ("%.3f" % ks[d["lbl"]]) if ks.get(d["lbl"]) is not None else "-"))
print()

json.dump(dict(onset=onset, A=rA, B=rB, C=rC, D=rD),
          open(WORK + "/final.json", "w", encoding="utf-8"), indent=1)
print("-> work/_M6/final.json")
