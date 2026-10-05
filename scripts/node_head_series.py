# -*- coding: utf-8 -*-
"""Manhole water-surface elevation series, for the head-driven exchange (C3).

For every junction:  Z_node(t) = invert + depth(t)          from the SWMM .out
                     z_rim      = invert + maxdepth         (the manhole rim, i.e. the ground)

Why this replaces the flow-ratio criterion
------------------------------------------
The first two-way implementation decided "can this manhole accept surface water?"
by `q_in < 0.75 * Q_full` and capped the return by `Q_full - q_in`.  Both include
the return flow we are about to impose, so any I > 0 switched the inlets off and
the fixed-point iteration collapsed to the one-way solution (R = 60 397 -> 0 -> 0).

A head criterion depends only on STATE:

    dH = (z_rim + h_surface) - (Z_node)
    dH > 0  -> water enters the sewer,  Q = C*A*sqrt(2 g dH)
    dH < 0  -> the sewer is above the surface: SWMM's own flooding handles it

Nothing in it refers to the imposed flow, so the loop closes physically instead of
cancelling itself.  It also removes the arbitrary "residual capacity" cap: if the
sewer cannot take the water it surcharges, Z_node rises and dH -> 0 by itself.

Output: work/nodehead_<tag>.npz  (+ .json for the node list)
"""
import datetime as dt
import json
import os
import re
import sys
import time

import numpy as np
ROOT = os.environ.get("EMS_SWMM_BASE",
                   os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

WORK = os.path.join(ROOT, "work")


def junctions(inp):
    txt = open(inp, encoding="utf-8", errors="replace").read()
    sec = re.search(r"\[JUNCTIONS\](.*?)(?=\n\[)", txt, re.S).group(1)
    ids, inv, md = [], [], []
    for line in sec.splitlines():
        p = line.split()
        if len(p) >= 3 and re.match(r"^N\d+$", p[0]):
            try:
                ids.append(int(p[0][1:]))
                inv.append(float(p[1]))
                md.append(float(p[2]))
            except ValueError:
                continue
    return np.array(ids, np.int64), np.array(inv), np.array(md)


def main(inp, tag, force_run=False):
    ids, z_inv, z_md = junctions(inp)
    rpt, out = inp.replace(".inp", ".rpt"), inp.replace(".inp", ".out")
    if force_run or not os.path.exists(out):
        from pyswmm import Simulation
        with Simulation(inp, rpt, out) as sim:
            sim.step_advance(60)
            for _ in sim:
                pass
    from swmm.toolkit.shared_enum import NodeAttribute
    from pyswmm import Output
    # downstream conduit full-flow capacity (Manning): the SPARE capacity limits
    # how much surface water a manhole can pass on while its sewer is loaded
    import flood_series as FS
    _e, _d, _c = FS.parse(inp)
    _qf = FS.outlet_capacity(_e, _d, _c)
    q_full = np.array([float(_qf.get("N%d" % i, 0.0)) for i in ids])
    o = Output(out)
    ts = list(o.times)
    names = ["N%d" % i for i in ids]
    Z = np.zeros((ids.size, len(ts)), np.float32)
    QI = np.zeros((ids.size, len(ts)), np.float32)
    t0 = time.time()
    for i in range(len(ts)):
        d = o.node_attribute(NodeAttribute.INVERT_DEPTH, time_index=i)
        q = o.node_attribute(NodeAttribute.TOTAL_INFLOW, time_index=i)
        for k, nm in enumerate(names):
            v = d.get(nm)
            if v:
                Z[k, i] = float(v)
            v = q.get(nm)
            if v:
                QI[k, i] = float(v)
    Z += z_inv[:, None].astype(np.float32)          # depth -> water surface elevation
    np.savez_compressed("%s/nodehead_%s.npz" % (WORK, tag),
                        ids=ids, z_inv=z_inv.astype(np.float32),
                        z_rim=(z_inv + z_md).astype(np.float32), Z=Z,
                        q_in=QI, q_full=q_full.astype(np.float32))
    json.dump({"t0": ts[0].isoformat(), "dt": 60, "n_steps": len(ts),
               "n_nodes": int(ids.size), "inp": os.path.basename(inp)},
              open("%s/nodehead_%s.json" % (WORK, tag), "w"))
    print("[%s] %d 井 × %d 步  Z_node %.2f~%.2f m  井口 %.2f~%.2f m  (%.0f s)"
          % (tag, ids.size, len(ts), float(Z.min()), float(Z.max()),
             float((z_inv + z_md).min()), float((z_inv + z_md).max()),
             time.time() - t0))
    print("  零水深占比 %.1f%%；管满流能力 p10/p50/p90 = %.4f/%.4f/%.4f m3/s"
          % (100.0 * float((Z <= z_inv[:, None] + 1e-6).mean()),
             *np.percentile(q_full[q_full > 0], [10, 50, 90])))
    print("-> work/nodehead_%s.npz" % tag)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
