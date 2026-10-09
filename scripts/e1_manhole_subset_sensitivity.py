# -*- coding: utf-8 -*-
"""v3 review F4c -- manhole-only ponding-area sensitivity (M-family).  rev 3.

Question (v2 M4 / v3 F4): E1 imposes A_ponded on the mixed node set (3 277
nodes).  Does the event-dependent optimum survive if the area is applied ONLY
to the nodes that sit on an inspection-manhole symbol (circle_large, n = 1 414)?

rev 3 changes (after two lock-related batch deaths):
  * build() is NO LONGER CALLED when the A-family INP already exists (all 32
    do).  The A-INP is READ and only the ponded-area column is rewritten -- this
    removes every write to a pre-existing file, which is what the transient
    locks were denying.  It also makes the M-family bit-identical to the
    A-family INPs that produced the published E1 results (perfect comparability).
  * run window metadata comes from the A-family sim JSON (same windows as E1).
  * M-INP written to a temp file + os.replace (atomic), with retries.
  * optional --events filter so several workers can run in parallel on
    disjoint event sets (each worker writes only its own M-INP/RPT/OUT/JSON).
"""
import datetime as dt
import glob
import json
import os
import re
import sys
import time

from pyswmm import Simulation, Nodes, Links

sys.path.insert(0, r"E:/work/tongji/swmm_model/scripts")
from metrics_lib import (ME, WORK, coverage, event_window, load_sim, metrics,
                         obs_manhole_series)
from run_multi_event import NODE_NAME, NODES, PUMPS, build

OUT = r"E:/work/tongji/swmm_model/output"
EVENTS = ["20250627", "20250730", "20260719", "20260809"]
AREAS = [25, 50, 75, 100, 150, 200, 300, 400]

SYM = json.load(open(WORK + "/node_symbol_class.json", encoding="utf-8"))
MANHOLE = {int(k) for k, v in SYM.items() if v == "circle_large"}
print("manhole-symbol nodes: %d / %d" % (len(MANHOLE), len(SYM)), flush=True)

J_LINE = re.compile(r"^(N(\d+))\s+(-?[\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+(\d+)\s*$",
                    flags=re.M)


def robust_write(path, text, tries=6):
    """Write with lock-proof fallback: direct writes only, never os.replace.

    Locks in this environment hit existing files written shortly before
    (AV/indexer scanning).  New unique filenames have never failed, so: try
    the canonical path if it does not exist yet, then fall back to unique
    suffixes -- SWMM only needs a readable INP at some path.
    """
    stem, ext = os.path.splitext(path)
    candidates = [path] if not os.path.exists(path) else []
    candidates += ["%s_u%d_%d%s" % (stem, os.getpid(), k, ext) for k in range(tries)]
    for target in candidates:
        try:
            with open(target, "w", encoding="utf-8") as f:
                f.write(text)
            return target
        except PermissionError:
            print("  [write-blocked] %s -> next candidate" % target, flush=True)
            time.sleep(15.0)
    raise PermissionError("all write candidates failed for " + path)


def subset_from_text(txt, tag_stats):
    def repl(m):
        tag_stats["n_junc"] += 1
        nid = int(m.group(2))
        if nid in MANHOLE:
            tag_stats["n_man"] += 1
            return m.group(0)
        tag_stats["n_zeroed"] += 1
        return m.group(1) + "\t" + "\t".join(m.group(g) for g in (3, 4, 5, 6)) + "\t0"

    return J_LINE.sub(repl, txt)


def make_subset_inp(ev, area):
    """Read the existing A-family INP, zero non-manhole ponded areas.

    Falls back to build() only if the A-INP or its sim metadata are missing.
    Returns (dst_path, t0, t1, eval_end).
    """
    a_inp = "%s/swmm_%s_A%d.inp" % (OUT, ev, area)
    # A100 is the legacy 'B' alias in run_multi_event: its sim metadata may be
    # stored as sim_{ev}_B.json (e.g. sim_20260719_A100.json does not exist).
    sim_candidates = ["%s/sim_%s_A%d.json" % (ME, ev, area)]
    if area == 100:
        sim_candidates.append("%s/sim_%s_B.json" % (ME, ev))
    a_sim = next((c for c in sim_candidates if os.path.exists(c)), None)
    if os.path.exists(a_inp) and a_sim:
        txt = open(a_inp, encoding="utf-8").read()
        meta = json.load(open(a_sim, encoding="utf-8"))
        win = meta.get("window", [])
        t0 = dt.datetime.fromisoformat(win[0]) if win else None
        t1 = dt.datetime.fromisoformat(win[1]) if len(win) > 1 else t0
        ev_end = t1        # provenance only; collect() derives the metric
                           # window from obs + base sim, as E1 did
        src_label = os.path.basename(a_inp) + " (read-only reuse)"
    else:
        print("  [fallback] A-INP/sim missing for %s A%d -> build()" % (ev, area), flush=True)
        src, t0, t1, ev_end = build(ev, "A%d" % area)
        txt = open(src, encoding="utf-8").read()
        src_label = os.path.basename(src) + " (rebuilt)"

    stats = {"n_junc": 0, "n_man": 0, "n_zeroed": 0}
    txt = subset_from_text(txt, stats)
    assert stats["n_junc"] == stats["n_man"] + stats["n_zeroed"], stats
    assert stats["n_man"] == 1414, stats          # the manhole set is fixed

    dst = "%s/swmm_%s_M%d.inp" % (OUT, ev, area)
    dst = robust_write(dst, txt)
    print("  INP %s via %s: junctions %d, manholes kept %d, zeroed %d"
          % (os.path.basename(dst), src_label,
             stats["n_junc"], stats["n_man"], stats["n_zeroed"]), flush=True)
    return dst, t0, t1, ev_end


def sim_paths(ev, tag):
    """All result files for (ev, tag): canonical + unique-fallback names."""
    return (glob.glob("%s/sim_%s_%s.json" % (ME, ev, tag))
            + glob.glob("%s/sim_%s_%s_u*.json" % (ME, ev, tag)))


def load_m(ev, tag):
    paths = sim_paths(ev, tag)
    return json.load(open(paths[0], encoding="utf-8")) if paths else None


def run_one(ev, area):
    tag = "M%d" % area
    out_json = "%s/sim_%s_%s.json" % (ME, ev, tag)
    if sim_paths(ev, tag):
        print("[skip] %s/%s exists" % (ev, tag), flush=True)
        return
    inp, t0, t1, ev_end = make_subset_inp(ev, area)
    hist = {n: [] for n in NODES}
    pumps = {p: [] for p in PUMPS}
    t_start = time.time()
    with Simulation(inp, inp.replace(".inp", ".rpt"), inp.replace(".inp", ".out")) as sim:
        sim.step_advance(60)
        nodes, links = Nodes(sim), Links(sim)
        for _ in sim:
            ct = sim.current_time
            for n in NODES:
                hist[n].append([ct.isoformat(), round(nodes["N%d" % n].depth, 4)])
            for p in PUMPS:
                pumps[p].append([ct.isoformat(), round(links[p].flow, 4)])
    el = time.time() - t_start
    payload = {"event": ev, "cfg": tag,
               "window": [t0.isoformat(), t1.isoformat()],
               "eval_window": [t0.isoformat(), ev_end.isoformat()],
               "nodes": {NODE_NAME[n]: hist[n] for n in NODES if n in NODE_NAME},
               "pumps": pumps}
    robust_write(out_json, json.dumps(payload))
    print("[%s/%s] %.0f s done" % (ev, tag, el), flush=True)


def collect():
    """Metrics for the M-family (area 0 via the base run) -> fig_e1_M_subset.json."""
    res = {}
    for ev in EVENTS:
        obs = json.load(open("%s/obs_%s.json" % (ME, ev), encoding="utf-8"))
        sim0 = load_sim(ev, "base")
        t_lo, t_hi = event_window(obs, sim0)
        osr = obs_manhole_series(obs, t_lo, t_hi)
        per = {}
        for area, cfg in [(0, "base")] + [(a, "M%d" % a) for a in AREAS]:
            sim = load_sim(ev, cfg) if area == 0 else load_m(ev, cfg)
            if sim is None:
                print("  ! %s/%s missing" % (ev, cfg), flush=True)
                continue
            vals = []
            for code, s in sorted(osr.items()):
                if code not in sim["nodes"]:
                    continue
                mt = [dt.datetime.fromisoformat(p[0]) for p in sim["nodes"][code]]
                mv = [p[1] for p in sim["nodes"][code]]
                mm = metrics(s, list(zip(mt, mv)), t_lo, t_hi)
                if mm is not None:
                    vals.append(mm)
            if not vals:
                continue
            per[area] = dict(n=len(vals),
                             r=float(sum(v["r"] for v in vals if v["r"] is not None))
                             / max(1, sum(1 for v in vals if v["r"] is not None)),
                             rmse=float(sum(v["rmse"] for v in vals) / len(vals)),
                             dt=float(sum(abs(v["dt_h"]) for v in vals) / len(vals)))
        xs = sorted(per)
        rm_max = max(per[a]["rmse"] for a in xs) or 1.0
        dt_max = max(per[a]["dt"] for a in xs) or 1.0
        for a in xs:
            per[a]["score"] = 0.5 * per[a]["rmse"] / rm_max + 0.5 * per[a]["dt"] / dt_max
        opt = min(xs, key=lambda a: per[a]["score"])
        res[ev] = dict(xs=xs, per={str(a): per[a] for a in xs}, grid_opt=opt,
                       rmse_range=(max(per[a]["rmse"] for a in xs) - per[opt]["rmse"]))
        print("%s  M-optimum %d m2  (score curve: %s)"
              % (ev, opt, " ".join("%d:%.2f" % (a, per[a]["score"]) for a in xs)), flush=True)
    json.dump(res, open(WORK + "/fig_e1_M_subset.json", "w", encoding="utf-8"))
    print("saved work/fig_e1_M_subset.json", flush=True)


if __name__ == "__main__":
    if "--collect" in sys.argv:
        collect()
        sys.exit(0)
    events = list(EVENTS)
    if "--events" in sys.argv:
        events = sys.argv[sys.argv.index("--events") + 1].split(",")
    for rnd in (1, 2):          # second pass retries anything the first pass lost
        for ev in events:
            for a in AREAS:
                if sim_paths(ev, "M%d" % a):
                    continue
                try:
                    run_one(ev, a)
                except Exception as exc:   # one bad run must not kill the batch
                    print("[FAIL] %s/M%d: %r" % (ev, a, exc), flush=True)
        left = [(ev, a) for ev in events for a in AREAS
                if not sim_paths(ev, "M%d" % a)]
        if not left:
            break
        if rnd == 2:
            print("STILL MISSING after retry: %s" % left, flush=True)
            sys.exit(1)
        print("[round 2] retrying %d missing run(s)" % len(left), flush=True)
    print("ALL RUNS DONE for %s" % events, flush=True)
