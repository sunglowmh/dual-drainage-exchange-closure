# -*- coding: utf-8 -*-
"""Shared metric helpers for the event-based SWMM evaluations.

Used by eval_multi_event.py (base vs B) and eval_aponded.py (Aponded sweep).
Keeping ONE implementation avoids the two scripts drifting apart.
"""
import bisect
import datetime as dt
import json
import os
import re

import numpy as np
ROOT = os.environ.get("EMS_SWMM_BASE",
                   os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

WORK = os.path.join(ROOT, "work")
OUT = os.path.join(ROOT, "output")
ME = WORK + r"/multi_event"

STN = json.load(open(WORK + r"/station_nodes.json", encoding="utf-8"))
CODE2NAME = {k: v["name"] for k, v in STN.items()}
MANHOLE = {k: v for k, v in STN.items() if not str(k).startswith("86")}


def frozen_mask(vals, run=4):
    """Samples inside a run of >= `run` identical values (telemetry hold)."""
    n = len(vals)
    keep = [True] * n
    i = 0
    while i < n:
        j = i
        while j + 1 < n and vals[j + 1] == vals[i]:
            j += 1
        if j - i + 1 >= run:
            for k in range(i, j + 1):
                keep[k] = False
        i = j + 1
    return keep


def join(mt, mv, ot, ov, tol=180):
    """Nearest-neighbour join of an observation series onto a model series."""
    md = dict(zip(mt, mv))
    out = []
    for t, v in zip(ot, ov):
        i = bisect.bisect_left(mt, t)
        best = None
        for j in (i - 1, i):
            if 0 <= j < len(mt):
                dd = abs((mt[j] - t).total_seconds())
                if dd <= tol and (best is None or dd < best[0]):
                    best = (dd, mt[j])
        if best:
            out.append((t, v, md[best[1]]))
    return out


def metrics(obs_series, sim_node, t_lo=None, t_hi=None):
    """Bias-corrected r / RMSE / peak error / peak-timing error for one station.

    The sensors have an uncalibrated zero (a pre-storm baseline of 0.65-1.5 m),
    so the offset is removed with the median before scoring -- absolute depth
    cannot be compared, only the dynamic response.
    """
    if t_lo is not None:
        obs_series = [(t, v) for t, v in obs_series if t_lo <= t <= t_hi]
    ot = [t for t, _ in obs_series]
    ov = [v for _, v in obs_series]
    pairs = join([p[0] for p in sim_node], [p[1] for p in sim_node], ot, ov)
    if len(pairs) < 15:
        return None
    o_all = np.array([p[1] for p in pairs])
    m_all = np.array([p[2] for p in pairs])
    t_all = [p[0] for p in pairs]
    keep = frozen_mask([v for _, v in obs_series])
    oidx = {t: k for k, t in enumerate(ot)}
    idx = [k for k, (t, _, _) in enumerate(pairs)
           if oidx.get(t) is None or keep[oidx[t]]]
    if len(idx) >= 10:
        o, m = o_all[idx], m_all[idx]
        n_excl = len(pairs) - len(idx)
    else:
        o, m, idx, n_excl = o_all, m_all, list(range(len(pairs))), 0
    if o.max() - o.min() < 0.10:          # flat station -> nothing to compare
        return None
    off = float(np.median(m - o))
    ma = m - off
    rmse = float(np.sqrt(np.mean((ma - o) ** 2)))
    r = float(np.corrcoef(m, o)[0, 1]) if np.std(m) > 1e-6 and np.std(o) > 1e-6 else np.nan
    nse = float(1 - np.sum((ma - o) ** 2) / np.sum((o - o.mean()) ** 2))
    pk_o = max(idx, key=lambda k: o_all[k])
    pk_m = max(idx, key=lambda k: m_all[k])
    return dict(n=len(idx), n_excl=n_excl, offset=round(off, 3), rmse=round(rmse, 3),
                r=round(r, 3) if np.isfinite(r) else None, nse=round(nse, 3),
                obs_peak=round(o_all[pk_o], 2), obs_peak_t=t_all[pk_o].strftime("%m-%d %H:%M"),
                sim_peak=round(m_all[pk_m], 2), sim_peak_t=t_all[pk_m].strftime("%m-%d %H:%M"),
                peak_err=round(m_all[pk_m] - o_all[pk_o], 2),
                dt_h=round((t_all[pk_m] - t_all[pk_o]).total_seconds() / 3600.0, 1))


def event_window(obs, sim=None):
    """Evaluation window: rain start -2 h .. rain end +8 h, capped by the
    observation coverage and truncated before a consensus second storm.

    `sim` (a sim_*.json payload) is preferred when it records `eval_window`:
    the window the simulation was built around is the authority, so the metrics
    can never end up averaging over uncovered time.
    """
    if sim and sim.get("eval_window"):
        return (dt.datetime.fromisoformat(sim["eval_window"][0]),
                dt.datetime.fromisoformat(sim["eval_window"][1]))
    import run_multi_event as R          # shared definitions, no side effects
    rain = [(dt.datetime.fromisoformat(t), v) for t, v in obs["rain"]]
    dm = int(obs["dt_min"])
    sig = [t for t, v in rain if v * 60.0 / dm >= R.THR_MM_H] or \
          [t for t, v in rain if v > 0] or [rain[0][0]]
    t_lo = sig[0] - dt.timedelta(hours=2)
    t_hi = sig[-1] + dt.timedelta(hours=R.EVAL_POST)
    o_end = R._obs_end(obs)
    if o_end and t_hi > o_end:
        t_hi = o_end
    r = R._second_storm(obs, sig[-1])
    if r and r <= t_hi:
        t_hi = max(r - dt.timedelta(minutes=30),
                   sig[-1] + dt.timedelta(hours=2))
    return t_lo, t_hi


def coverage(sim, t_lo, t_hi):
    """(ok, message) -- does the simulation actually span the metric window?"""
    w = sim.get("window")
    if not w:
        return True, ""
    a, b = dt.datetime.fromisoformat(w[0]), dt.datetime.fromisoformat(w[1])
    if a <= t_lo and b >= t_hi:
        return True, ""
    miss = max((t_lo - a).total_seconds(), (t_hi - b).total_seconds()) / 3600.0
    return False, ("模拟窗 %s~%s 未覆盖评估窗 %s~%s（缺 %.1f h）"
                   % (a.strftime("%m-%d %H:%M"), b.strftime("%m-%d %H:%M"),
                      t_lo.strftime("%m-%d %H:%M"), t_hi.strftime("%m-%d %H:%M"), miss))


def obs_manhole_series(obs, t_lo, t_hi, min_len=15):
    out = {}
    for code, rec in obs["manhole"].items():
        s = [(dt.datetime.fromisoformat(t), v) for t, v in rec["series"]
             if t_lo <= dt.datetime.fromisoformat(t) <= t_hi]
        if len(s) >= min_len:
            out[code] = s
    return out


def load_sim(event, cfg):
    """Simulation series for (event, cfg); 'A100' falls back to the legacy 'B'."""
    for c in ([cfg, "B"] if cfg == "A100" else [cfg]):
        f = "%s/sim_%s_%s.json" % (ME, event, c)
        if os.path.exists(f):
            return json.load(open(f, encoding="utf-8"))
    return None


def flood_summary(inp_or_rpt):
    """(total volume m3, node count, {node: (hours, volume m3, max rate)})."""
    txt = open(inp_or_rpt, errors="replace").read()
    i = txt.find("Node Flooding Summary")
    if i < 0:
        return None, None, {}
    sec = txt[i:]
    j = sec.find("*****", 200)
    sec = sec[:j] if j > 0 else sec
    rows, tot = {}, 0.0
    for l in sec.splitlines():
        p = l.split()
        # Node | HoursFlooded | MaxRate | day hr:min | Volume(10^6 ltr) | PondedDepth
        if len(p) >= 7 and re.match(r"^N\d+$", p[0]):
            try:
                v = float(p[5]) * 1000.0
                rows[int(p[0][1:])] = (float(p[1]), v, float(p[2]))
            except ValueError:
                continue
            tot += v
    return tot, len(rows), rows
