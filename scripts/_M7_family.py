# -*- coding: utf-8 -*-
"""M7 -- run the 16-member Chicago load family on the three conceptual networks.

Same load family and same extraction conventions as the campus family
(_M6_build_family / _M6_final): Chicago 120 min, r = 0.4, 1-min steps, return
periods 0.1-100 a plus the three load probes (the 100-a hyetograph x1.2/1.4/1.6);
V_in / flooding / export read from the SWMM "Flow Routing Continuity" table;
kappa = 1 - dFlood/dV_in between consecutive members; V_onset = interpolated load
where flooding reaches 1 % of V_in; Lambda = V_in / V_onset.

Each network's diameter scale is calibrated first so that its onset falls inside
the family (between the 0.5-a and 5-a members).

    python scripts/_M7_family.py            -> work/_M7/final.json + summary
"""
import importlib.util
import json
import os
import re
import sys
import time

ROOT = os.environ.get("EMS_SWMM_BASE",
                   os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
import _M7_nets

WORK = os.path.join(ROOT, "work", "_M7")
os.makedirs(WORK, exist_ok=True)

spec = importlib.util.spec_from_file_location(
    "gen", os.path.join(ROOT, "scripts", "generate_swmm.py"))
gen = importlib.util.module_from_spec(spec)
spec.loader.exec_module(gen)

P_MEMBERS = [(0.1, "P0.1a"), (0.15, "P0.15a"), (0.2, "P0.2a"), (0.25, "P0.25a"),
             (0.5, "P0.5a"), (1.0, "P1a"), (2.0, "P2a"), (3.0, "P3a"), (5.0, "P5a"),
             (10.0, "P10a"), (20.0, "P20a"), (50.0, "P50a"), (100.0, "P100a")]
PROBES = [(1.2, "L120"), (1.4, "L140"), (1.6, "L160"), (2.0, "L200"), (3.0, "L300")]
BASE100 = [v for _, v in gen.chicago(100.0, 120, 1, 0.4)]


def rain_of(tag):
    if tag.startswith("L"):
        f = float(tag[1:]) / 100.0
        return [(t, f * v) for t, v in zip(range(len(BASE100)), BASE100)]
    P = float(tag[1:-1])
    return gen.chicago(P, 120, 1, 0.4)


# ------------------------------------------------------------------ run + parse
def block(t, title):
    i = t.find(title)
    if i < 0:
        return ""
    j = t.find("\n  *", i + 200)
    return t[i:j if j > 0 else len(t)]


def gv(b, lbl):
    m = re.search(re.escape(lbl) + r"[. ]*(-?\d+\.\d+)\s+(-?\d+\.\d+)", b)
    return float(m.group(2)) * 1000.0 if m else None


def run_one(net_builder, scale, tag, tmp=False):
    """Build one member INP, run it with pyswmm, parse the continuity table.

    tmp=True overwrites a fixed CAL_ name (calibration iterations); no file is
    ever deleted (the host sandbox blocks os.remove).
    """
    nm = net_builder.__name__[-2:]
    inp = os.path.join(WORK, ("CAL_%s.inp" % nm) if tmp else ("%s_%s.inp" % (nm, tag)))
    net = net_builder()
    txt = net.inp(rain_of(tag), scale, net.label)
    with open(inp, "w", encoding="utf-8") as f:
        f.write(txt)
    rpt, out = inp[:-4] + ".rpt", inp[:-4] + ".out"
    from pyswmm import Simulation
    t0 = time.time()
    with Simulation(inp, rpt, out) as sim:
        sim.step_advance(3600)
        for _ in sim:
            pass
    dt = time.time() - t0
    t = open(rpt, errors="replace").read()
    b = block(t, "Flow Routing Continuity")
    assert b, "no continuity block in " + rpt
    row = dict(tag=tag, ww=gv(b, "Wet Weather Inflow") or 0.0,
               flood=gv(b, "Flooding Loss") or 0.0,
               export=gv(b, "External Outflow") or 0.0)
    m = re.search(r"Continuity Error \(%\) \.+\s*(-?[\d.]+)", t)
    row["cont_err"] = float(m.group(1)) if m else None
    row["sec"] = round(dt, 1)
    row["flood_frac"] = row["flood"] / row["ww"] if row["ww"] else 0.0
    row["eta"] = (row["export"] / (row["export"] + row["flood"])
                  if (row["export"] + row["flood"]) > 0 else float("nan"))
    return row


# ------------------------------------------------------------------ calibration
def calibrate(builder):
    """Bisect the diameter scale until the flooding onset lies between 0.5 a and 5 a.

    Larger scale = larger pipes = onset at larger return period.  lo marks a
    scale whose onset is below 0.5 a (pipes too small), hi one whose onset is
    above 5 a (too large); the acceptable window lies between them.
    """
    scale = 1.0
    lo = hi = None
    hist = []
    for _ in range(14):
        r05 = run_one(builder, scale, "P0.5a", tmp=True)
        r5 = run_one(builder, scale, "P5a", tmp=True)
        f05, f5 = r05["flood_frac"], r5["flood_frac"]
        hist.append((round(scale, 4), f05, f5))
        print("    scale %.3f  floodfrac(0.5a)=%.4f  (5a)=%.4f" % (scale, f05, f5), flush=True)
        if f05 <= 0.01 <= f5:
            return scale, hist
        if f5 < 0.01:                       # onset above 5 a: pipes too large
            hi = scale
        else:                               # onset below 0.5 a: too small
            lo = scale
        if lo is not None and hi is not None:
            if hi - lo < 0.02:              # window may be empty; widen it
                break
            scale = 0.5 * (lo + hi)
        elif hi is not None:
            scale = max(0.15, scale * 0.6)
        else:
            scale = min(5.0, scale * 1.7)
    # fallback: accept any scale whose onset falls anywhere inside the family
    for s in sorted(set(h[0] for h in hist)):
        r02 = run_one(builder, s, "P0.2a", tmp=True)
        r20 = run_one(builder, s, "P20a", tmp=True)
        if r02["flood_frac"] <= 0.01 <= r20["flood_frac"]:
            print("    fallback window (0.2a, 20a) accepted at scale %.3f" % s, flush=True)
            return s, hist
    return None, hist


# ------------------------------------------------------------------ family
def family(builder, scale):
    rows = []
    for tag in [tag for _P, tag in P_MEMBERS] + [tag for _f, tag in PROBES]:
        r = run_one(builder, scale, tag)
        r["lbl"] = tag
        rows.append(r)
        print("    %-8s V_in %8.0f  flood %9.0f  export %9.0f  eta %.3f  cont %.2f%%  (%.0f s)"
              % (tag, r["ww"], r["flood"], r["export"], r["eta"], r["cont_err"], r["sec"]))
    return rows


def onset_and_kappa(rows):
    Vi = [r["ww"] for r in rows]
    onset = None
    for i in range(len(rows) - 1):
        fi, fj = rows[i]["flood_frac"], rows[i + 1]["flood_frac"]
        if fi < 0.01 <= fj:
            w = (0.01 - fi) / (fj - fi)
            onset = Vi[i] + w * (Vi[i + 1] - Vi[i])
            break
    segs = []
    for x, y in zip(rows[:-1], rows[1:]):
        dv = y["ww"] - x["ww"]
        if dv <= 0:
            segs.append(None)
            continue
        k = (y["flood"] - x["flood"]) / dv
        segs.append(dict(lbl=x["lbl"] + "->" + y["lbl"], k=k, kappa=1.0 - k,
                         Lam_mid=0.5 * (x["ww"] + y["ww"]) / onset if onset else None,
                         Lam_x=x["ww"] / onset if onset else None,
                         Lam_y=y["ww"] / onset if onset else None))
    return onset, segs


BUILDERS = [(_M7_nets.build_N1, "N1", "grid"),
            (_M7_nets.build_N2, "N2", "random tree"),
            (_M7_nets.build_N3, "N3", "sparse flat")]


def main():
    result = {}
    fp = os.path.join(WORK, "final.json")
    if os.path.exists(fp):                       # resume: keep completed networks
        result = json.load(open(fp, encoding="utf-8"))
        print("resume: %s already done" % ", ".join(result), flush=True)
    for builder, nm, label in BUILDERS:
        if nm in result:
            print("### %s (%s) -- already done, skip" % (nm, label), flush=True)
            continue
        print("### %s (%s)" % (nm, label), flush=True)
        print("  calibrating diameter scale ...", flush=True)
        scale, hist = calibrate(builder)
        assert scale, "calibration failed for %s" % nm
        print("  chosen scale %.3f" % scale, flush=True)
        rows = family(builder, scale)
        onset, segs = onset_and_kappa(rows)
        net = builder()
        result[nm] = dict(label=label, scale=scale,
                          spec=_M7_nets.summary(net),
                          onset=dict(V=onset, q=onset / (120 * 60)),
                          rows=rows, segments=segs, calib=hist)
        print("  V_onset = %.0f m3  (q_onset %.3f m3/s)" % (onset, onset / 7200.0), flush=True)
        print()
        json.dump(result, open(fp, "w", encoding="utf-8"),
                  indent=1)                       # incremental: keep what is done

    json.dump(result, open(os.path.join(WORK, "final.json"), "w", encoding="utf-8"), indent=1)

    print("# M7 summary: kappa-Lambda across three conceptual networks")
    print()
    for nm, d in result.items():
        s = d["spec"]
        print("## %s -- %s: %d junctions, %d conduits, %.1f ha, scale %.2f, V_onset %.0f m3"
              % (nm, d["label"], s["junctions"], s["conduits"], s["area_ha"],
                 d["scale"], d["onset"]["V"]))
        print("| segment | Lambda(mid) | kappa |")
        print("|---|---|---|")
        for g in d["segments"]:
            if g:
                print("| %s | %8.2f | %.3f |" % (g["lbl"], g["Lam_mid"], g["kappa"]))
        print()
    print("-> work/_M7/final.json")


if __name__ == "__main__":
    main()
