# -*- coding: utf-8 -*-
"""Extract the per-node flooding RATE time series from a SWMM run.

Why: the 1D->2D coupling used to release each flooded node's total volume shaped
by the rainfall hyetograph, i.e. the overflow timing was assumed, not computed.
The real overflow hydrograph is needed both for a faithful release and for the
RETURN flow (surface -> sewer) of the two-way coupling.

⚠️ `pyswmm.Node.flooding` returns 0 throughout even when the report shows
flooding (checked across the 2026-07-19 peak: 843 nodes overflow in the report,
`node.flooding` = 0).  Read the written .out instead:
`Output.node_attribute(NodeAttribute.FLOODING_LOSSES, time_index=i)`.

Also derives each manhole's inlet capacity from its downstream conduit's
full-flow capacity (Manning, circular section).

Usage: python flood_series.py output/swmm_20260719_base.inp 20260719_base [run]
       (append "run" to force a fresh SWMM run; otherwise an existing .out is read)
"""
import json
import math
import os
import re
import sys
import time
ROOT = os.environ.get("EMS_SWMM_BASE",
                   os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

WORK = os.path.join(ROOT, "work")


def blocks(inp):
    txt = open(inp, encoding="utf-8", errors="replace").read()
    out = {}
    for name in ("JUNCTIONS", "CONDUITS", "XSECTIONS"):
        m = re.search(r"\[%s\](.*?)(?=\n\[)" % name, txt, re.S)
        out[name] = [l.split() for l in (m.group(1).splitlines() if m else [])
                     if l.strip() and not l.strip().startswith(";")]
    return out


def parse(inp):
    b = blocks(inp)
    elev, diam, cond = {}, {}, []
    for r in b["JUNCTIONS"]:
        if len(r) >= 2:
            try:
                elev[r[0]] = float(r[1])
            except ValueError:
                pass
    for r in b["XSECTIONS"]:
        if len(r) >= 3:
            try:
                diam[r[0]] = float(r[2])
            except ValueError:
                pass
    for r in b["CONDUITS"]:
        if len(r) >= 5:
            try:
                cond.append((r[0], r[1], r[2], float(r[3]), float(r[4])))
            except ValueError:
                pass
    return elev, diam, cond


def outlet_capacity(elev, diam, cond):
    """Manning full-flow capacity of each node's downstream conduit (m3/s)."""
    cap = {}
    for name, f, t, L, n in cond:
        d = diam.get(name)
        if not d:
            continue
        s = max((elev.get(f, 0.0) - elev.get(t, 0.0)) / max(L, 0.1), 5e-4)
        a = math.pi * d * d / 4.0
        q = (1.0 / n) * a * ((d / 4.0) ** (2.0 / 3.0)) * math.sqrt(s)
        cap[f] = max(cap.get(f, 0.0), q)
    return cap


def read_out(out, with_inflow=True):
    """-> (times, {node_id: {step: rate}}, {node_id: {step: inflow}}).

    The node inflow is needed to cap the RETURN flow by the conduit's SPARE
    capacity: a pipe already carrying its full-flow discharge cannot also
    swallow surface water, even while its node is not (yet) overflowing.
    """
    from swmm.toolkit.shared_enum import NodeAttribute
    from pyswmm import Output
    o = Output(out)
    ts = list(o.times)
    series = {}
    for i in range(len(ts)):
        d = o.node_attribute(NodeAttribute.FLOODING_LOSSES, time_index=i)
        for name, v in d.items():
            if v and v > 1e-6:
                series.setdefault(int(str(name).lstrip("N")), {})[i] = round(float(v), 6)
    inflow = {}
    if with_inflow and series:
        for nid in series:
            nm = "N%d" % nid
            try:
                s = o.node_series(nm, NodeAttribute.TOTAL_INFLOW)
            except Exception:
                break
            vals = list(s.values()) if isinstance(s, dict) else list(s)
            if not vals:
                continue
            # FULL-length series (m3/s) so the coupling can look up the spare
            # capacity at any step, not only while the node happens to flood
            inflow[nid] = [round(float(v), 5) for v in vals]
    return ts, series, inflow


def main(inp, tag, force_run=False):
    elev, diam, cond = parse(inp)
    cap = outlet_capacity(elev, diam, cond)
    rpt, out = inp.replace(".inp", ".rpt"), inp.replace(".inp", ".out")
    if force_run or not os.path.exists(out) or os.path.getmtime(out) < os.path.getmtime(inp):
        from pyswmm import Simulation
        t0 = time.time()
        with Simulation(inp, rpt, out) as sim:
            sim.step_advance(60)
            for _ in sim:
                pass
        print("  (SWMM 运行 %.0f s)" % (time.time() - t0), flush=True)

    ts, series, inflow = read_out(out)
    t0s, k = ts[0], len(ts)
    vtot = sum(sum(v.values()) for v in series.values()) * 60.0
    json.dump({"t0": t0s.isoformat(), "dt": 60, "n_steps": k,
               "flood_nodes": {str(n): series[n] for n in series},
               "node_inflow": {str(n): inflow.get(n, {}) for n in series},
               "outlet_cap": {str(n): round(cap.get("N%d" % n, 0.0), 4)
                              for n in series},
               "note": "flooding rate CMS per 60-s step (from .out "
                       "FLOODING_LOSSES); node_inflow = TOTAL_INFLOW CSS; "
                       "outlet_cap = downstream conduit full-flow capacity "
                       "(Manning, m3/s)"},
              open("%s/floodts_%s.json" % (WORK, tag), "w"))
    print("[%s] %d steps (%s ~ %s)  %d flooded nodes  total %.0f m3"
          % (tag, k, ts[0].strftime("%m-%d %H:%M"), ts[-1].strftime("%m-%d %H:%M"),
             len(series), vtot), flush=True)
    caps = sorted(c for c in cap.values() if c > 0)
    if caps:
        print("  outlet capacity: p10 %.3f  median %.3f  p90 %.3f m3/s"
              % (caps[len(caps) // 10], caps[len(caps) // 2], caps[9 * len(caps) // 10]))
    if inflow:
        rr = sorted(v / max(cap.get("N%d" % n, 1e-9), 1e-9)
                    for n, vs in inflow.items() for v in vs)
        print("  节点入流/管能力 比值 p50 %.2f  p90 %.2f（>1 = 已被占满）"
              % (rr[len(rr) // 2], rr[9 * len(rr) // 10]))
    print("-> work/floodts_%s.json" % tag)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], force_run=(len(sys.argv) > 3 and sys.argv[3] == "run"))
