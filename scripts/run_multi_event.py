# -*- coding: utf-8 -*-
"""Cross-event runs: storms x {base, A<area>} configs.

base : v7 settings (ALLOW_PONDING NO, PS1 2.20/0.90)
A<n> : v8 settings with ALLOW_PONDING YES, Aponded = n m2, PS1 1.60/0.90
B    : alias of A100 (kept for the runs already on disk)

Usage: python run_multi_event.py 20250627:base 20250627:A50 20260719:A400 ...
"""
import datetime as dt
import json
import re
import sys
import time
import os

from pyswmm import Simulation, Nodes, Links
ROOT = os.environ.get("EMS_SWMM_BASE",
                   os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

WORK = os.path.join(ROOT, "work")
OUT = os.path.join(ROOT, "output")
TEMPLATE = OUT + r"/swmm_v7_50a.inp"
ME = WORK + r"/multi_event"

STN = json.load(open(WORK + r"/station_nodes.json", encoding="utf-8"))
NODES = sorted({v["node"] for v in STN.values()})
NODE_NAME = {v["node"]: k for k, v in STN.items()}
PUMPS = ["PS1", "PS2"]

# Window rule v4.  Anchored on the SIGNIFICANT rainfall (intensity >= THR_MM_H),
# not on "any nonzero value" -- traces of 0.1-0.2 mm per 10 min (0.6-1.2 mm/h
# drizzle) otherwise keep the window open for many hours.  20250730 is the
# textbook case: the storm was over by 07-30 17:25 (96.4 % of the total fallen)
# while the record dribbles on to 07-31 09:25.
#
#   PRE  3 h -- the model has no dry-weather flow, so the sewers start empty;
#               3 h settles the initial transient;
#   POST 10 h -- bracketed from above by (a) the OBSERVATION coverage and (b) the
#               next UNMEASURED storm (2026-07-19's next response starts 07-20
#               20:01, 25.3 h after the rain ended -- +10 h leaves 15.3 h).
#
# Two invariants make the tail safe:
#   * EVAL_END = min(rain end + EVAL_POST, observation end)   -- never compare
#     against time that has no data;
#   * T_END   >= EVAL_END + 1 h and <= min(rain end + POST_H, observation end
#     + 2 h) -- the evaluation window is always strictly inside the run.
# A "second storm" guard may pull T_END back, but ONLY when the consensus rise
# lies beyond the evaluation window; a rise inside it truncates EVAL_END (with a
# warning) instead of silently corrupting the comparison.
PRE_H, POST_H, EVAL_POST, THR_MM_H = 3.0, 10.0, 8.0, 3.0

# guard thresholds: a rise counts only when it is large enough to be a storm
# (>= RISE_M) on a device that actually had water (post-rain min >= WET_MIN), and
# when at least N_VOTE devices agree within VOTE_WIN -- single-probe artefacts
# (a dry sensor returning 0.00 -> 0.60 m) must not move the window.
WET_MIN, RISE_M, RISE_FRAC = 0.10, 0.20, 0.30
N_VOTE, VOTE_WIN = 2, dt.timedelta(hours=2.0)


def _obs_end(obs):
    ends = [dt.datetime.fromisoformat(s[-1][0])
            for g in ("manhole", "road")
            for rec in obs.get(g, {}).values()
            for s in [rec.get("series") or []] if s]
    return max(ends) if ends else None


def _second_storm(obs, rain_end, verbose=False):
    """Time of the first CONSENSUS unexplained rise, or None.

    A device votes at the first moment its level climbs above its own post-rain
    minimum by `max(RISE_M, RISE_FRAC * amplitude)`, where `amplitude` is the
    device's own storm range.  The RELATIVE part matters: on 20260719 the two
    voting devices recover 97-98 % of their storm amplitude (a real second
    storm), whereas on 20250627 two devices move only 19-23 % while six others
    are still declining -- an operational/local effect, not rainfall.

    Devices whose post-rain minimum is below WET_MIN are skipped: a probe sitting
    at 0.00 m cannot distinguish "no water" from "baseline" (20250730's 13号
    "rose 0.60 m" simply because the well refilled).
    """
    votes, weak = [], []
    for grp in ("manhole", "road"):
        for rec in obs.get(grp, {}).values():
            s = [(dt.datetime.fromisoformat(a), b) for a, b in (rec.get("series") or [])]
            post = [(a, b) for a, b in s if a > rain_end]
            if len(post) < 10:
                continue
            tm, vm = min(post, key=lambda z: z[1])
            if vm < WET_MIN:
                continue
            amp = max(b for _, b in s) - min(b for _, b in s)
            thr = max(RISE_M, RISE_FRAC * amp)
            t_vote = None
            for a, b in post:
                if a > tm and b - vm >= thr:
                    t_vote = a
                    break
            if t_vote:
                votes.append(t_vote)
            else:
                top = max((b for a, b in post if a > tm), default=vm)
                if top - vm >= RISE_M:                 # coordinated but small
                    weak.append((max(post, key=lambda z: z[1])[0],
                                 top - vm, thr))
    votes.sort()
    for i in range(len(votes) - N_VOTE + 1):
        if votes[i + N_VOTE - 1] - votes[i] <= VOTE_WIN:
            if verbose:
                print("    [守卫] 共识抬升 %d 站（±2 h 内）：%s"
                      % (N_VOTE, votes[i + N_VOTE - 1].strftime("%m-%d %H:%M")))
            return votes[i + N_VOTE - 1]
    if verbose and weak:
        print("    [守卫] 存在抬升但未达阈值（幅度 %.2f m < 需要 %.2f m，"
              "不计为后续暴雨）：最多 %d 站" % (weak[0][1], weak[0][2], len(weak)))
    return None


def build(event, cfg, post_h=None, verbose=True):
    """post_h: override the post-rain buffer (h).  Use ~24 h only for dedicated
    recession studies."""
    obs = json.load(open("%s/obs_%s.json" % (ME, event), encoding="utf-8"))
    rain = [(dt.datetime.fromisoformat(t), v) for t, v in obs["rain"]]
    dt_min = int(obs["dt_min"])
    sig = [t for t, v in rain if v * 60.0 / dt_min >= THR_MM_H]
    wet = sig or [t for t, v in rain if v > 0] or [rain[0][0]]
    o_end = _obs_end(obs)

    t0 = wet[0] - dt.timedelta(hours=PRE_H)
    eval_end = wet[-1] + dt.timedelta(hours=EVAL_POST)
    if o_end and eval_end > o_end:
        eval_end = o_end
    t1 = wet[-1] + dt.timedelta(hours=post_h if post_h else POST_H)
    if o_end and t1 > o_end + dt.timedelta(hours=2):
        t1 = o_end + dt.timedelta(hours=2)
    if t1 < eval_end + dt.timedelta(hours=1):
        t1 = eval_end + dt.timedelta(hours=1)

    notes = []
    r = _second_storm(obs, wet[-1], verbose=verbose)
    if r:
        if r <= eval_end:                       # contamination inside the window
            eval_end = max(r - dt.timedelta(minutes=30),
                           wet[-1] + dt.timedelta(hours=2))
            notes.append("评估窗被后续降雨污染 -> 截断至 %s（抬升 %s）"
                         % (eval_end.strftime("%m-%d %H:%M"), r.strftime("%m-%d %H:%M")))
        if r < t1:
            t1 = max(r - dt.timedelta(hours=1), eval_end + dt.timedelta(hours=1))
            notes.append("窗口尾提前到 %s" % t1.strftime("%m-%d %H:%M"))
    if verbose:
        print("  [window v4] %s ~ %s (%.1f h) | 有效雨 %s~%s | 评估窗止 %s%s"
              % (t0.strftime("%m-%d %H:%M"), t1.strftime("%m-%d %H:%M"),
                 (t1 - t0).total_seconds() / 3600.0,
                 wet[0].strftime("%m-%d %H:%M"), wet[-1].strftime("%m-%d %H:%M"),
                 eval_end.strftime("%m-%d %H:%M"),
                 "  ⚠ " + "；".join(notes) if notes else ""))
    dt_min = int(obs["dt_min"])

    txt = open(TEMPLATE, encoding="utf-8").read()
    txt = txt.replace("START_DATE           09/15/2026", "START_DATE           %s" % t0.strftime("%m/%d/%Y"))
    txt = txt.replace("START_TIME           00:00:00", "START_TIME           %s" % t0.strftime("%H:%M:%S"))
    txt = txt.replace("END_DATE             09/15/2026", "END_DATE             %s" % t1.strftime("%m/%d/%Y"))
    txt = txt.replace("END_TIME             06:00:00", "END_TIME             %s" % t1.strftime("%H:%M:%S"))
    txt = txt.replace("REPORT_START_DATE    09/15/2026", "REPORT_START_DATE    %s" % t0.strftime("%m/%d/%Y"))
    txt = txt.replace("REPORT_START_TIME    00:00:00", "REPORT_START_TIME    %s" % t0.strftime("%H:%M:%S"))

    # rain gauge interval must match the data
    txt = re.sub(r"^RG1\s+INTENSITY\s+0:01",
                 "RG1             INTENSITY 0:%02d" % dt_min, txt, flags=re.M)

    lines = []
    # Aggregate onto a uniform grid: the 2026-07-19 gauge logs at ~64 s spacing,
    # so formatting to whole minutes can repeat a timestamp and SWMM then raises
    # ERROR 173 (data out of sequence).
    grid = {}
    for t, depth in rain:
        if t < t0 or t > t1:
            continue
        k = t.replace(second=0, microsecond=0,
                      minute=(t.minute // dt_min) * dt_min)
        grid[k] = grid.get(k, 0.0) + depth
    for t in sorted(grid):
        lines.append("RAIN          %s  %.2f" % (t.strftime("%m/%d/%Y %H:%M"),
                                                grid[t] * 60.0 / dt_min))
    ts = ("[TIMESERIES]\n;;observed gauge %s, %d-min, intensity mm/h\n" % (event, dt_min)
          + "\n".join(lines) + "\n")
    txt = re.sub(r"\[TIMESERIES\].*?(?=\n\[)", ts, txt, flags=re.S)

    if cfg != "base":
        # ponding configs: B == A100; A<n> = ponding with Aponded = n m2.
        # PS1 start/shutoff (1.60 / 0.90) stays the same in all ponding runs so
        # that the ONLY difference across A50/A100/A400 is the ponded area.
        m_a = re.match(r"^(?:A(\d+)|B)$", cfg)
        if not m_a:
            raise ValueError("unknown cfg %r (use base | B | A<area>)" % cfg)
        area = 100 if cfg == "B" else int(m_a.group(1))
        txt = txt.replace("ALLOW_PONDING        NO", "ALLOW_PONDING        YES")
        txt, n_pond = re.subn(
            r"^(N\d+)\s+(-?[\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+0\s*$",
            lambda m: "%s      %s      %s      %s      %s      %d" % (m.groups() + (area,)),
            txt, flags=re.M)
        txt = txt.replace(
            "PS1            N900100        N900101        PS1C           ON       2.20     0.90",
            "PS1            N900100        N900101        PS1C           ON       1.60     0.90")
        print("  ponding: Aponded=%d m2 on %d junctions, PS1 1.60/0.90" % (area, n_pond))
    path = "%s/swmm_%s_%s.inp" % (OUT, event, cfg)
    open(path, "w", encoding="utf-8").write(txt)
    return path, t0, t1, eval_end


def _dump_json(path, obj, tries=6):
    """Write JSON, retrying on a transient Windows share violation.

    A killed sibling process can still hold the handle for a few seconds, which
    surfaced as PermissionError on sim_20250730_base.json (2026-09-21) and cost
    a completed 15-minute run.
    """
    for i in range(tries):
        try:
            with open(path, "w") as f:
                json.dump(obj, f)
            return path
        except PermissionError:
            if i == tries - 1:
                alt = path.replace(".json", "_alt.json")
                with open(alt, "w") as f:
                    json.dump(obj, f)
                print("  ⚠ %s 被占用，改写 %s" % (path, alt), flush=True)
                return alt
            time.sleep(3.0 * (i + 1))


def run(event, cfg):
    inp, t0, t1, ev_end = build(event, cfg)
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
    _dump_json("%s/sim_%s_%s.json" % (ME, event, cfg),
               {"event": event, "cfg": cfg,
                "window": [t0.isoformat(), t1.isoformat()],
                "eval_window": [t0.isoformat(), ev_end.isoformat()],
                "nodes": {NODE_NAME[n]: hist[n] for n in NODES if n in NODE_NAME},
                "pumps": pumps})
    print("[%s/%s] %.0f s  window %s -> %s" % (event, cfg, el, t0, t1), flush=True)
    for n in NODES:
        if n not in NODE_NAME:
            continue
        s = hist[n]
        pk = max(s, key=lambda z: z[1])
        print("    %-28s peak %.2f @%s" % (STN[NODE_NAME[n]]["name"][:28], pk[1], pk[0][5:16]), flush=True)
    for p in PUMPS:
        on = sum(1 for _, f in pumps[p] if f > 1e-4)
        print("    %s on %d min, maxQ %.3f" % (p, on, max(f for _, f in pumps[p])), flush=True)


if __name__ == "__main__":
    for arg in sys.argv[1:]:
        ev, cfg = arg.split(":")
        run(ev, cfg)
