# -*- coding: utf-8 -*-
"""Two-way (dual-drainage) coupling: 1D SWMM <-> 2D surface, iterated.

Pass 1  SWMM (ALLOW_PONDING NO)          -> per-node flooding rate F_i(t)
Pass 1  2D with source F_i(t)            -> surface depth h(x,y)
        return inlets (manhole opening)  -> return rate R_i(t) back to the sewers
Pass 2  SWMM + [INFLOWS] R_i(t)          -> updated flooding/depths, and a
                                            convergence check against pass 1

Two physical ingredients that the one-way coupling lacked:
  * the overflow is released with its COMPUTED timing (not the rain hyetograph);
  * surface water re-enters the sewer through manhole inlets, capped by the
    downstream conduit's full-flow capacity, and only while that sewer node is
    NOT surcharged (a pressurised sewer cannot accept surface water).

Usage: python couple_2way.py 20260719            (event key)
       python couple_2way.py 20260719 1          (1 extra iteration)
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
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from solver2d import Surface2D

from pyswmm import Simulation, Nodes

WORK = os.path.join(ROOT, "work")
OUT = os.path.join(ROOT, "output")
TX, TY = 2421.762, 4753.796
CELL = float(os.environ.get("TONGJI_CELL", "2.0"))
MANNING = float(os.environ.get("TONGJI_MANNING", "0.05"))
DT = 1.0
RECORD_CB = None          # optional per-step recorder: record(t, solver) every COUPLE s
COUPLE = 60.0              # coupling interval (s)
T_END = 25200.0            # 7 h, covers the 2026-07-19 event
INLET_AREA = 0.25          # m2, manhole opening (d ~ 0.56 m)
INLET_COEF = 0.6


# --------------------------------------------------------------------- setup
def load_ctx(event):
    z = np.load(WORK + "/surf2d_z.npy")
    act = np.load(WORK + "/surf2d_active.npy")
    meta = json.load(open(WORK + "/surf2d_meta.json", encoding="utf-8"))
    gates = json.load(open(WORK + "/gates.json", encoding="utf-8"))
    sv = json.load(open(WORK + "/lake_sv.json", encoding="utf-8"))
    lake_mask = np.load(WORK + "/lake_mask.npy").astype(bool)
    net = json.load(open(WORK + "/network_v2.json", encoding="utf-8"))
    pos = {n["id"]: (n["x"] + TX, n["y"] + TY) for n in net["nodes"]}
    return z, act, meta, gates, sv, lake_mask, pos


def make_solver(z, act, meta, gates, sv, lake_mask, pos, scen=None):
    x0, y1 = float(meta["x0"]), float(meta["y1"])
    ghosts = [(int((y1 - g["y"]) / CELL), int((g["x"] - x0) / CELL),
               g["zout"], g["verdict"]) for g in gates]
    S = Surface2D(z, act, cell=CELL, manning=MANNING, gates=ghosts, dt_max=DT)
    rr, cc = np.nonzero(lake_mask & act)
    idx = S.k[rr, cc]
    scen = scen or {}
    cap = 0.7222 * float(scen.get("pump", 1.0))
    S.set_lake_outlet(idx[idx >= 0], sill_stage=sv["normal_stage"], capacity=cap)
    stage = float(scen.get("lake_m", 1.0))
    v0 = S.set_lake_init(sv["base_elev"] + stage)    # draw-down scenario
    g = scen.get("gate", "AUTO")
    if g == "GD4":
        S.lake_open_t = 4 * 3600.0                   # manual: opens after 4 h
    elif g == "GCL":
        S.lake_open_t = 1e9                          # never opened
    return S, pos, x0, y1, v0


def cell_of(pos, nid, x0, y1, shape, k):
    if nid not in pos:
        return None
    gx, gy = pos[nid]
    c, r = int((gx - x0) / CELL), int((y1 - gy) / CELL)
    if 0 <= r < shape[0] and 0 <= c < shape[1]:
        kk = k[r, c]
        return int(kk) if kk >= 0 else None
    return None


# ------------------------------------------------------------------ sources
class SourceFromSeries:
    """Per-node flooding rate -> 2D source array (buffer reused each step)."""

    def __init__(self, S, cells, rates, dt_couple=COUPLE):
        self.S = S
        self.cells = cells
        self.rates = rates
        self.dt = dt_couple
        self.buf = np.zeros(S.n)
        self.last = None

    def __call__(self, t):
        i = int(t / self.dt)
        if i >= self.rates.shape[1]:
            self.buf[:] = 0.0
            self.last = i
            return self.buf
        if i != self.last:
            self.buf[:] = 0.0
            # np.add.at: several nodes may share a cell, and `buf[cells] = ...`
            # would silently keep only the last one (under-injection)
            np.add.at(self.buf, self.cells, self.rates[:, i])
            self.last = i
        return self.buf


def prep_inlets_head(S, pos, x0, y1, npz, off, n_couple, area=INLET_AREA):
    """Inlets for the HEAD-DRIVEN exchange: every junction with a 2D cell.

    Head mode does not need the "which nodes were flooding" filter -- the
    criterion dH = (rim + h_surface) - Z_node switches each inlet on and off by
    itself.  `off` aligns the coupling steps with the 1D record.
    Returns (cells, rim_z, Z_node[j], ids).
    """
    d = np.load(npz)
    ids0, z_rim, Z = d["ids"], d["z_rim"], d["Z"]
    QI, QF = d["q_in"], d["q_full"]
    free = np.full(n_couple, 1.0e9, np.float32)      # "no data" -> never returns
    cells, rim, head, cap, ids = [], [], [], [], []
    j = np.arange(n_couple) + off
    for m, nid in enumerate(ids0):
        kk = cell_of(pos, int(nid), x0, y1, S.k.shape, S.k)
        if kk is None:
            continue
        row = free.copy()
        crow = np.zeros(n_couple, np.float32)
        good = j < Z.shape[1]
        row[good] = Z[m, j[good]]
        # SPARE capacity of the downstream conduit right now; where the sewer is
        # already full (q_in >= Q_full) the manhole cannot pass surface water on
        crow[good] = np.maximum(0.0, float(QF[m]) - QI[m, j[good]])
        cells.append(kk)
        rim.append(float(z_rim[m]))
        head.append(row)
        cap.append(crow)
        ids.append(int(nid))
    return (np.array(cells, np.int64), np.array(rim, float),
            np.array(head, np.float32), ids, np.array(cap, np.float32))


def prep_sources(S, fs, pos, x0, y1, off, n_couple, spare_frac=0.75):
    """Build (cells, rates[cell, couple_step], ok, caps, ids).

    `off` = number of 60-s steps between the 1D window start and the 2D window
    start; without it the 2D run would see the pre-rain hours only and would
    compress the whole overflow into the last minutes of the simulation.
    """
    fld = fs["flood_nodes"]
    infl_all = fs.get("node_inflow", {})
    ok_rows, cells, caps, ids, rates, cap_rows = [], [], [], [], [], []
    for nid_s, d in fld.items():
        nid = int(nid_s)
        kk = cell_of(pos, nid, x0, y1, S.k.shape, S.k)
        if kk is None:
            continue
        q_full = float(fs["outlet_cap"].get(nid_s, 0.0))
        arr = infl_all.get(nid_s)
        arr = np.asarray(arr, float) if arr else None
        r = np.zeros(n_couple)
        ok = np.ones(n_couple, bool)
        cap_row = np.full(n_couple, q_full)
        for i in range(n_couple):
            step = i + off
            q_in = float(arr[step]) if (arr is not None and step < arr.size) else 0.0
            cap_row[i] = max(0.0, q_full - q_in)
            ok[i] = q_in < spare_frac * q_full
        for step_s, q in d.items():
            i = int(step_s) - off
            if 0 <= i < n_couple:
                r[i] = q
                ok[i] = False          # overflowing node cannot accept water
        cells.append(kk)
        rates.append(r)
        ok_rows.append(ok)
        cap_rows.append(cap_row)
        caps.append(q_full)
        ids.append(nid)
    return (np.array(cells, np.int64), np.array(rates, np.float32),
            np.array(ok_rows, bool), np.array(caps), ids,
            np.array(cap_rows, np.float32))


# ------------------------------------------------------------- return inflows
def _w(path, data, is_json=True, tries=8):
    """Write text/JSON, retrying transient Windows share violations.

    A `PermissionError` here does NOT mean the file is read-only -- it means some
    other handle (search indexer, AV, a killed sibling) held it for a moment.
    Losing a completed 27-minute iteration to that is unacceptable, so retry,
    then write an "_alt" neighbour rather than dying.
    """
    for i in range(tries):
        try:
            with open(path, "w", encoding="utf-8") as f:
                if is_json:
                    json.dump(data, f, ensure_ascii=False, indent=1)
                else:
                    f.write(data)
            return path
        except PermissionError as e:
            last = e
            time.sleep(1.5 * (i + 1))
    alt = (path[:-5] + "_alt.json") if (is_json and path.endswith(".json")) \
        else (path[:-4] + "_alt" + path[-4:])
    with open(alt, "w", encoding="utf-8") as f:
        if is_json:
            json.dump(data, f, ensure_ascii=False, indent=1)
        else:
            f.write(data)
    print("  ⚠ %s 被占用(%s)，改写 %s" % (path, last, alt), flush=True)
    return alt


def collect_return(ret_hist, ids):
    """Cumulative snapshots -> per-node volume PER coupling interval (m3)."""
    series = {}
    prev = {}
    for t, acc in ret_hist:
        for nid, cum in acc.items():
            inc = cum - prev.get(nid, 0.0)
            prev[nid] = cum
            if inc > 1e-9:
                series.setdefault(nid, {})[int(t)] = \
                    series.setdefault(nid, {}).get(int(t), 0.0) + inc
    return series


def relax(imposed, R, alpha):
    """Fixed-point relaxation: I <- I + alpha * (R(I) - I).

    Needed because the two-way exchange is implicit: returning surface water
    makes the sewer overflow again, which in turn forbids return at those
    moments.  With alpha = 1 the first iterate over-shoots (1D flooding rose
    from 69 992 to 79 828 m3); alpha < 1 damps the iteration without moving the
    fixed point.
    """
    out = {}
    for nid in set(imposed) | set(R):
        a, b = imposed.get(nid, {}), R.get(nid, {})
        d = {}
        for k in set(a) | set(b):
            d[k] = a.get(k, 0.0) + alpha * (b.get(k, 0.0) - a.get(k, 0.0))
        out[nid] = d
    return out


def write_return_inp(src_inp, out_inp, imposed, t0, ids):
    """Write the imposed surface->sewer return flow into SWMM as [INFLOWS]."""
    txt = open(src_inp, encoding="utf-8", errors="replace").read()
    lines_ts, lines_in = [], []
    for nid in ids:
        pts = imposed.get(nid)
        if not pts or sum(pts.values()) < 1.0:      # ignore < 1 m3 over the event
            continue
        for t in sorted(pts):
            clock = t0 + dt.timedelta(seconds=t)
            lines_ts.append("RET_%-10d %s  %.5f"
                            % (nid, clock.strftime("%m/%d/%Y %H:%M"), pts[t] / COUPLE))
        lines_in.append("%-14s FLOW           RET_%-10d FLOW     1.0      1.0      0.0"
                        % ("N%d" % nid, nid))
    if not lines_in:
        return None
    txt = re.sub(r"(\[TIMESERIES\].*?)(?=\n\[)",
                 lambda m: m.group(1) + "\n" + "\n".join(lines_ts) + "\n",
                 txt, count=1, flags=re.S)
    if "[INFLOWS]" in txt:
        txt = re.sub(r"\[INFLOWS\].*?(?=\n\[)",
                     "[INFLOWS]\n" + "\n".join(lines_in) + "\n", txt, count=1, flags=re.S)
    else:
        txt = txt.replace("[TIMESERIES]",
                          "[INFLOWS]\n" + "\n".join(lines_in) + "\n\n[TIMESERIES]", 1)
    _w(out_inp, txt, is_json=False)
    return len(lines_in)


def run_1d(inp, tag, want_flood=True):
    """Run SWMM; return (flood series, station depths, pump flows).

    ⚠️ The flooding series comes from the written .out
    (`NodeAttribute.FLOODING_LOSSES`), NOT from `pyswmm.Node.flooding`, which
    returns 0 throughout even when the report shows flooding.
    """
    from swmm.toolkit.shared_enum import NodeAttribute
    from pyswmm import Links, Output

    stn = json.load(open(WORK + "/station_nodes.json", encoding="utf-8"))
    nodes_of = sorted({v["node"] for v in stn.values()})
    sdepth = {n: [] for n in nodes_of}
    pump = {"PS1": [], "PS2": []}
    t0, k = None, 0
    rpt, out = inp.replace(".inp", ".rpt"), inp.replace(".inp", ".out")
    with Simulation(inp, rpt, out) as sim:
        sim.step_advance(60)
        nodes, links = Nodes(sim), Links(sim)
        for _ in sim:
            if t0 is None:
                t0 = sim.current_time
            for n in nodes_of:
                sdepth[n].append(round(nodes["N%d" % n].depth, 4))
            for p in pump:
                pump[p].append(round(links[p].flow, 4))
            k += 1
    series = {}
    if want_flood:
        o = Output(out)
        for i in range(min(k, len(list(o.times)))):
            d = o.node_attribute(NodeAttribute.FLOODING_LOSSES, time_index=i)
            for name, v in d.items():
                if v and v > 1e-6:
                    # key WITHOUT the "N" prefix, so prep_sources can int() it
                    series.setdefault(int(str(name).lstrip("N")), {})[str(i)] = \
                        round(float(v), 6)
    res = dict(t0=t0.isoformat(), dt=60, n_steps=k, stn_depth=sdepth, pumps=pump,
               flood_nodes=series)
    _w("%s/tw_%s.json" % (WORK, tag), res)
    vtot = sum(sum(d.values()) for d in series.values()) * 60.0
    return res, vtot


# ------------------------------------------------------------------- driver
def main(event, extra_iter=0, use_inlets=True, max_hours=None, alpha=0.6,
         runtag=""):
    ftag = os.environ.get("TONGJI_FTSUFFIX", "base")
    fs = json.load(open("%s/floodts_%s_%s.json" % (WORK, event, ftag), encoding="utf-8"))
    obs = json.load(open("%s/multi_event/obs_%s.json" % (WORK, event), encoding="utf-8"))
    z, act, meta, gates, sv, lake_mask, pos = load_ctx(event)
    t1d = dt.datetime.fromisoformat(fs["t0"])
    # 2D window: from 30 min before the rain to 8 h after it stops.  Aligning the
    # two clocks is essential -- the 1D window starts 6 h before the rain, so
    # taking the first N steps of the series would release the overflow before
    # it rains.
    rain = [(dt.datetime.fromisoformat(t), v) for t, v in obs["rain"]]
    wet = [t for t, v in rain if v > 0]
    t_2d = wet[0] - dt.timedelta(minutes=30)
    t_end = wet[-1] + dt.timedelta(hours=8)
    t_end = min(t_end, t1d + dt.timedelta(seconds=fs["n_steps"] * 60))
    if max_hours:
        t_end = min(t_end, t_2d + dt.timedelta(hours=max_hours))
    t_span = (t_end - t_2d).total_seconds()
    off = int((t_2d - t1d).total_seconds() // 60)
    n_couple = int(t_span // 60) + 1
    print("event %s | 1D 窗口 %s (%d 步) | 2D 窗口 %s ~ %s (%.1f h, off=%d) | 淹没节点 %d"
          % (event, t1d.strftime("%m-%d %H:%M"), fs["n_steps"],
             t_2d.strftime("%m-%d %H:%M"), t_end.strftime("%m-%d %H:%M"),
             t_span / 3600, off, len(fs["flood_nodes"])))

    m_s = re.search(r"_L(\d\d)_P(\d\d\d)_(AUTO|GD4|GCL)", runtag or "")
    scen = {"lake_m": int(m_s.group(1)) / 10.0, "pump": int(m_s.group(2)) / 100.0,
            "gate": m_s.group(3)} if m_s else {"lake_m": 1.0, "pump": 1.0, "gate": "AUTO"}
    print("情景 %s  湖初位 %.2f m  泵 %.2fx  闸门 %s"
          % (runtag or "(基准)", scen["lake_m"], scen["pump"], scen["gate"]), flush=True)

    head_mode = use_inlets and ("HEAD" in (runtag or "").upper())
    head_npz = "%s/nodehead_%s_base.npz" % (WORK, event)
    if head_mode and not os.path.exists(head_npz):
        raise SystemExit("缺少 %s：先跑 scripts/node_head_series.py" % head_npz)

    inp1 = "%s/swmm_%s_base.inp" % (OUT, event)
    flood0 = sum(sum(d.values()) for d in fs["flood_nodes"].values()) * 60.0
    print("pass 1 (1D, 单向基准): 漫溢总量 %.0f m3" % flood0)

    tag = event
    imposed = {}          # relaxed surface->sewer inflow (node -> {step: m3})
    hist = []
    for it in range(extra_iter + 1):
        # ---------------- 2D with the CURRENT 1D overflow (+ return inlets)
        S, pos, x0, y1, v_lake0 = make_solver(z, act, meta, gates, sv, lake_mask,
                                             pos, scen)
        cells, rates, ok, caps, ids, capm = prep_sources(
            S, fs, pos, x0, y1, off, n_couple)
        scen_tag = (("2way-head" if "HEAD" in (runtag or "").upper() else "2way")
                    if use_inlets else "1way-aligned")
        if use_inlets and head_mode:
            # C3: head-driven exchange at every manhole (no flow-ratio mask, no
            # spare-capacity cap -- both were functions of the imposed flow and
            # made the iteration collapse onto the one-way solution).
            cells_i, rim_i, head_i, ids, cap_i = prep_inlets_head(
                S, pos, x0, y1, head_npz, off, n_couple)
            n_in = S.set_return_inlets(
                [(nid, kk, 0.0, INLET_AREA) for nid, kk in zip(ids, cells_i)],
                coef=INLET_COEF, step=COUPLE, head=head_i, rim=rim_i, cap=cap_i)
            print("  入口模式：水头驱动 + 剩余容量封顶（%d 个检查井，判据 dH>0）"
                  % n_in, flush=True)
        elif use_inlets:
            n_in = S.set_return_inlets(
                [(nid, kk, cap, INLET_AREA) for nid, kk, cap in zip(ids, cells, caps)],
                coef=INLET_COEF, ok=ok, step=COUPLE, cap=capm)
            print("  入口模式：流量比（旧判据，仅作对照）", flush=True)
        else:
            n_in = 0                       # one-way, but with the SAME flood timing
        src = SourceFromSeries(S, cells, rates)
        # columns the solver actually applies: it uses int((t+dt/2)/COUPLE) for
        # t < t_span, so the last column (index n_couple-1) is never reached.
        # Counting it inflated the "released" volume by up to 1 292 m3 (one
        # peak-time column) and showed up as a fake -1.13 % closure error.
        n_app = max(1, int(t_span // COUPLE))
        released = float(rates[:, :n_app].sum()) * COUPLE

        # A finished 2D solve is checkpointed at once: iterating costs 25-30 min
        # and the driver has twice been killed by a transient share violation
        # seconds AFTER the solve, which used to discard the whole iteration.
        ckf = "%s/tw_%s_i%d%s%s_ckpt.json" % (WORK, tag, it,
                                              "" if use_inlets else "_1way", runtag)
        ck = None
        if os.path.exists(ckf):
            try:
                c = json.load(open(ckf, encoding="utf-8"))
                if c.get("scen") == scen_tag and abs(c.get("span_s", 0) - t_span) < 1:
                    ck = c
            except (ValueError, OSError):
                ck = None
        if ck:
            R = {int(k): {int(s): v for s, v in d.items()} for k, d in ck["R"].items()}
            summary = ck["summary"]
            print("  [续跑] 复用 iter %d 的 2D 结果：回流 %.0f m3，蓄存 %.0f m3"
                  % (it, sum(sum(d.values()) for d in R.values()),
                     summary.get("stored_m3", float("nan"))), flush=True)
        else:
            t_start = time.time()
            S.run(t_span, src_of=src, dt=DT,
                  record=RECORD_CB if RECORD_CB is not None else (lambda t, s: None),
                  record_every=COUPLE)
            el = time.time() - t_start
            # the SOLVER's own injection count is authoritative: recomputing it
            # from `rates` can disagree with what was actually applied (index vs
            # interval boundaries) and would show up as a fake closure error
            IN = S.src_vol + v_lake0      # initial landscape water is an inflow too
            OUTV = S.out_vol + S.lake_out_vol + S.ret_vol + S.volume() + S.clip_vol
            print("  2D iter %d: 注入 %.0f (算出 %.0f, 1D 事件 %.0f) + 初始湖 %.0f | "
                  "校门 %.0f | 闸门-泵站 %.0f | "
                  "**回流管网 %.0f** | 蓄存 %.0f | 削截 %.0f | 闭合差 %+.4f%%  (%.0f s, %d inlets)"
                  % (it, S.src_vol, released, flood0, v_lake0, S.out_vol, S.lake_out_vol,
                     S.ret_vol, S.volume(), S.clip_vol,
                     100 * (OUTV - IN) / IN, el, n_in), flush=True)
            try:
                np.save("%s/h2d_tw_%s_i%d%s%s.npy"
                        % (WORK, tag, it, "" if use_inlets else "_1way", runtag), S.h)
            except (PermissionError, OSError) as e:
                print("  (跳过 .npy 写入：%s)" % e, flush=True)
            R = collect_return(S.ret_hist, ids)
            # land / lake split of the maximum depth: Table 9 has to state which one
            # it reports, and the campus land depth is the decision-relevant number.
            _isl = lake_mask[S.rr, S.cc].astype(bool)
            _hmax_lake = float(S.h[_isl].max()) if _isl.any() else float("nan")
            _hmax_land = float(S.h[~_isl].max()) if (~_isl).any() else float("nan")
            summary = {"scen": scen_tag, "iter": it,
                       "window": [t_2d.isoformat(), t_end.isoformat()],
                       "released_m3": released, "injected_m3": S.src_vol,
                       "lake_init_m3": v_lake0, "flood_event_m3": flood0,
                       "gate_out_m3": S.out_vol, "lake_out_m3": S.lake_out_vol,
                       "return_m3": S.ret_vol, "stored_m3": S.volume(),
                       "clipped_m3": S.clip_vol, "hmax_m": float(S.h.max()),
                       "hmax_land_m": _hmax_land, "hmax_lake_m": _hmax_lake,
                       "wet_ha": float((S.h >= 0.02).sum() * CELL ** 2 / 1e4),
                       "n_inlets": int(n_in), "runtag": runtag, "alpha": alpha}
            _w(ckf, {"R": {str(k): {str(s): v for s, v in d.items()}
                           for k, d in R.items()},
                     "summary": summary, "scen": scen_tag, "span_s": t_span,
                     "saved": dt.datetime.now().isoformat()})
        _w("%s/h2d_tw_%s_i%d%s%s_summary.json"
           % (WORK, tag, it, "" if use_inlets else "_1way", runtag), summary)
        if not use_inlets:
            print("  (单向修正对齐基准；不作 1D 重跑)")
            break
        # ---------------- relax the return inflow, then re-run the 1D model
        v_R = sum(sum(d.values()) for d in R.values())
        imposed = relax(imposed, R, alpha)
        v_I = sum(sum(d.values()) for d in imposed.values())
        inp2 = "%s/swmm_%s_return%s_i%d.inp" % (OUT, tag, runtag, it)
        n_ts = write_return_inp(inp1, inp2, imposed, t_2d, ids)
        if not n_ts:
            print("  无回流超过阈值，终止"); break
        print("  回流 R=%.0f m3 -> 施加 I=%.0f m3（松弛 α=%.2f，%d 条时序）"
              % (v_R, v_I, alpha, n_ts), flush=True)
        res, vtot = run_1d(inp2, "%s_i%d%s" % (tag, it, runtag))
        if head_mode and it < extra_iter:
            # next iteration's inlets must see the head of THIS 1D solution
            import node_head_series as NHS
            _nm = "%s%s_i%d" % (event, runtag, it)
            NHS.main(inp2, _nm)
            head_npz = "%s/nodehead_%s.npz" % (WORK, _nm)
            print("  Z_node 已更新 <- %s" % os.path.basename(head_npz), flush=True)
        print("  pass 2 (1D + 回流): 漫溢总量 %.0f m3  (基准 %.0f，%+.1f%%)"
              % (vtot, flood0, 100 * (vtot - flood0) / flood0), flush=True)
        hist.append(dict(iter=it, R=round(v_R, 1), imposed=round(v_I, 1),
                         flood=round(vtot, 1), n_ts=n_ts,
                         ret_2d=round(summary.get("return_m3", 0.0), 1)))
        _w("%s/tw_%s%s_hist.json" % (WORK, tag, runtag), hist)
        # next iteration uses the updated flood series and sewer-state mask
        fs = {"t0": res["t0"], "flood_nodes": res["flood_nodes"],
              "outlet_cap": fs["outlet_cap"], "n_steps": res["n_steps"]}
    print("done")


if __name__ == "__main__":
    def _num(x, default):
        return float(x) if x not in ("", None) else default

    ev = sys.argv[1]
    it = int(sys.argv[2]) if len(sys.argv) > 2 else 0
    inl = not (len(sys.argv) > 3 and sys.argv[3] == "oneway")
    mh = _num(sys.argv[4], None) if len(sys.argv) > 4 else None
    al = _num(sys.argv[5], 0.6) if len(sys.argv) > 5 else 0.6
    rt = sys.argv[6] if len(sys.argv) > 6 else ""
    main(ev, it, inl, mh, al, rt)
