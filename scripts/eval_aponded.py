# -*- coding: utf-8 -*-
"""H1 -- does one uniform equivalent surface-storage area fit every event?

Sweeps Aponded in {off(=base), 50, 100, 400} m2 across three events
(63.6 mm/50 min, 178.6 mm/5.2 h, 242.5 mm/36 h) and asks whether the
best-fitting value is event-independent.

Two independent lines of evidence:
  (1) calibration: which Aponded minimises the bias-corrected RMSE / maximises r;
  (2) physical form: the equivalent ponded area implied by the 2D solution,
      area_eff(node) = V_flood(node) / h_2D(node),
      whose SPREAD tells us whether a single uniform number can work.

Writes output/H1_等效蓄水_判据.md and output/H1_等效蓄水_判据.png
"""
import datetime as dt
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

matplotlib.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei"]
matplotlib.rcParams["axes.unicode_minus"] = False

ROOT = os.environ.get("EMS_SWMM_BASE",
                   os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from metrics_lib import (CODE2NAME, ME, OUT, WORK, coverage, event_window,
                         flood_summary,
                         load_sim, metrics, obs_manhole_series)

def vertex(x, y):
    """Parabola through three points -> (x*, curvature a).

    a > 0 means it opens upwards and x* is a MINIMUM; a < 0 means the three
    points bracket a maximum, so x* carries no information about an optimum.
    """
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    d = (x[0] - x[1]) * (x[0] - x[2]) * (x[1] - x[2])
    if abs(d) < 1e-12:
        return float(x[1]), 0.0
    a = (x[2] * (y[1] - y[0]) + x[1] * (y[0] - y[2]) + x[0] * (y[2] - y[1])) / d
    b = (x[2] ** 2 * (y[0] - y[1]) + x[1] ** 2 * (y[2] - y[0])
         + x[0] ** 2 * (y[1] - y[2])) / d
    if abs(a) < 1e-14:
        return float(x[1]), 0.0
    return float(-b / (2 * a)), float(a)


# The 9-level grid, identical to Fig. 3 and to
# output/P2_4_Aponded加密档位与连续最优点.md (the old {0,50,100,400} grid had
# an asymmetric spacing -- 50 then 300 -- and is retired).
CFGS = [("base", 0, "无洼蓄"), ("A25", 25, "25 m²"), ("A50", 50, "50 m²"),
        ("A75", 75, "75 m²"), ("A100", 100, "100 m²"), ("A150", 150, "150 m²"),
        ("A200", 200, "200 m²"), ("A300", 300, "300 m²"), ("A400", 400, "400 m²")]
EVLABEL = {"20250627": "20250627 短历时强降雨", "20250730": "20250730 台风竹节草",
           "20260719": "20260719 极端暴雨", "20260809": "20260809 台风白海豚"}
# Depths = trapezoidal integrals over the real timestamps (scripts/_M6_rain.py);
# durations = the span of the >= 3 mm/h significant rainfall.
EVRAIN = {"20250627": (63.6, "50 min"), "20250730": (125.3, "13.2 h"),
          "20260719": (178.6, "5.2 h"), "20260809": (242.5, "36 h")}

x0 = y1 = None
meta = json.load(open(WORK + "/surf2d_meta.json", encoding="utf-8"))
x0, y1 = float(meta["x0"]), float(meta["y1"])
CELL = 2.0
TX, TY = 2421.762, 4753.796
net = json.load(open(WORK + "/network_v2.json", encoding="utf-8"))
POS = {n["id"]: (n["x"] + TX, n["y"] + TY) for n in net["nodes"]}

# --------------------------------------------------------------- 1) sweep ---
rows = []
print("%-10s %8s %7s %10s %8s | %s" % ("event", "flooded", "V tot", "V/node",
                                       "f_node", "mean r / rmse / |dT| per cfg"))
for ev in ["20250627", "20250730", "20260719", "20260809"]:
    obs = json.load(open("%s/obs_%s.json" % (ME, ev), encoding="utf-8"))
    sim0 = load_sim(ev, "base")
    t_lo, t_hi = event_window(obs, sim0)
    if sim0:
        ok, msg = coverage(sim0, t_lo, t_hi)
        if not ok:
            print("  ⚠ %s %s" % (ev, msg))
    osr = obs_manhole_series(obs, t_lo, t_hi)
    per_cfg = {}
    for cfg, area, lab in CFGS:
        sim = load_sim(ev, cfg)
        if sim is None:
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
        if vals:
            per_cfg[cfg] = dict(n=len(vals),
                                r=float(np.nanmean([v["r"] for v in vals if v["r"] is not None])),
                                rmse=float(np.mean([v["rmse"] for v in vals])),
                                dt=float(np.mean([abs(v["dt_h"]) for v in vals])),
                                area=area)
    # composite score: bias-corrected RMSE and peak-timing error, each
    # normalised by the worst value in the event (equal weight).  RMSE alone
    # ties (20260719: base and A100 both 0.57 m) while the timing differs by
    # 8 h, which is exactly what the ponding term is meant to fix.
    if per_cfg:
        rmax = max(d["rmse"] for d in per_cfg.values()) or 1.0
        dmax = max(d["dt"] for d in per_cfg.values()) or 1.0
        for d in per_cfg.values():
            d["score"] = 0.5 * d["rmse"] / rmax + 0.5 * d["dt"] / dmax
    # continuous optimum (quadratic interpolation) + near-optimal set + the
    # retired 4-level optimum, so the report can say exactly which grid a
    # quoted value came from.  The vertex is only meaningful when the
    # curvature is positive AND the vertex lies INSIDE the bracket.
    gr = ve = cu = float("nan")
    kind, near, old4 = "n/a", [], None
    if per_cfg:
        byarea = {d["area"]: d for d in per_cfg.values()}
        xs = sorted(byarea)
        ys = [byarea[x]["score"] for x in xs]
        i = int(np.argmin(ys))
        lo, hi = max(0, i - 1), min(len(xs) - 1, i + 1)
        if i == 0:
            lo, hi = 0, 2
        elif i == len(xs) - 1:
            lo, hi = len(xs) - 3, len(xs) - 1
        ve, cu = vertex(xs[lo:hi + 1], ys[lo:hi + 1])
        gr = xs[i]
        inside = (cu > 0) and (xs[lo] < ve < xs[hi]) and (i not in (0, len(xs) - 1))
        kind = ("interior" if inside else
                ("boundary-low" if i == 0 else
                 ("boundary-high" if i == len(xs) - 1 else "non-convex")))
        # TWO tolerances on the NORMALISED composite score (bounded by 1
        # within each event), both kept so no reader has to guess:
        #   near         = score <= best + 0.05  ("+5 %" -- the rule quoted
        #                  in the manuscript Sec. 4.1 and in P2-4)
        #   near_strict  = score <= best + 0.02  ("+2 %")
        # NB the tolerance is absolute in score units, NOT "5 % of best".
        near = [x for x in xs if byarea[x]["score"] <= ys[i] + 0.05]
        near_strict = [x for x in xs if byarea[x]["score"] <= ys[i] + 0.02]
        old4 = gr if gr in (0, 50, 100, 400) else None

    # base-config flooding statistics (physical driver)
    tot, nfl, frows = flood_summary("%s/swmm_%s_base.rpt" % (OUT, ev))
    if nfl is None or tot is None:
        print("  ⚠ %s 的 base.rpt 尚无漫溢汇总（正在求解？）——该行统计留空"
              % ev)
    vpn = (tot / nfl) if (tot and nfl) else float("nan")
    frac = (nfl / 3277.0) if nfl else float("nan")
    best = min(per_cfg.items(), key=lambda kv: kv[1]["score"]) if per_cfg else (None, None)
    rows.append(dict(event=ev, rain=EVRAIN[ev], flooded=nfl, vol=tot, v_per_node=vpn,
                     frac=frac, cfg=per_cfg, best=best[0], grid_opt=gr,
                     vertex=ve, curv=cu, kind=kind, near=near,
                     near_strict=near_strict, old4=old4, win=(t_lo, t_hi)))
    print("%-10s %8s %7s %10s %8s | %s" % (
        ev, nfl if nfl else "-", "%.0f" % tot if tot else "-",
        "%.1f" % vpn if vpn == vpn else "-",
        "%.3f" % frac if frac == frac else "-",
        "  ".join("%s:%.2f/%.2f/%.1f" % (c, d["r"], d["rmse"], d["dt"])
                  for c, d in per_cfg.items())))

print("\n最佳配置（复合评分最小 = RMSE 与峰现时差各 50%，按事件内最大值归一）：")
for r in rows:
    b = r["cfg"].get(r["best"], {})
    base = r["cfg"].get("base", {})
    print("  %s  → %s (Aponded=%s m²)  rmse %.2f / ΔT %.1f h   (base rmse %.2f / ΔT %.1f h)  r %.2f (base %.2f)"
          % (r["event"], r["best"], b.get("area"), b.get("rmse", float("nan")),
             b.get("dt", float("nan")), base.get("rmse", float("nan")),
             base.get("dt", float("nan")), b.get("r", float("nan")), base.get("r", float("nan"))))

# ------------------------------------------- 2) 2D-implied effective area ---
# Canonical product of scripts/_E1d_reeff.py:
#     A_eff = V_flood / max(hbar_wet, 0.10 m),  hbar_wet = mean depth over wet
#     cells in a 5x5 (10 m) window, on the CORRECTED 2026-07-19 field
#     (h2d_tw_20260719_i0_1way).
# Do NOT rebuild this from h2d_obs: that field comes from the run_2d.py
# "compressed release" bug and carries a 5.98 m water depth on LAND
# (bed z = 3.014 m -> level 9.0 m).  See output/P2_5_A_eff适用条件.md.
eff = np.load(WORK + "/fig_e1_eff.npy")
print("\n2D 反推等效蓄水面积（0719 修正场，%d 个节点）" % eff.size)
if eff.size:
    print("  p10 %.0f | p25 %.0f | 中位 %.0f | p75 %.0f | p90 %.0f m²"
          % tuple(np.percentile(eff, [10, 25, 50, 75, 90])))
    print("  与均匀 100 m² 的比：中位 %.2f×  IQR(p25–p75) %.0f–%.0f m²  p75/p25 %.1f"
          % (np.median(eff) / 100.0, np.percentile(eff, 25),
             np.percentile(eff, 75),
             np.percentile(eff, 75) / np.percentile(eff, 25)))

# --------------------------------------------------------------- 3) figure --
fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.6))
ax = axes[0]
for r in rows:
    xs = [d["area"] for d in r["cfg"].values()]
    ys = [d["rmse"] for d in r["cfg"].values()]
    o = np.argsort(xs)
    ax.plot(np.array(xs)[o], np.array(ys)[o], "o-", label=EVLABEL[r["event"]])
    # mark the COMPOSITE-score optimum (the value quoted everywhere), not the
    # RMSE minimum -- Fig. 3(a) marks the composite selection too
    _bd = r["cfg"].get(r["best"])
    if _bd:
        ax.plot([_bd["area"]], [_bd["rmse"]], "*", ms=14,
                color=ax.lines[-1].get_color())
ax.set_xlabel("Aponded (m²)")
ax.set_ylabel("偏置校正 RMSE (m)")
ax.set_title("(a) 各事件的标定曲线（★ = 复合评分最优档，与图 3a 同口径）",
             fontsize=11, loc="left")
ax.legend(fontsize=8.5)
ax.grid(alpha=0.3)

ax = axes[1]
vpn = [r["v_per_node"] for r in rows]
bestarea = [r["cfg"].get(r["best"], {}).get("area", 0) for r in rows]
ax.scatter(vpn, bestarea, s=90, c="#4a6fa5", zorder=3)
for r, x, yv in zip(rows, vpn, bestarea):
    ax.annotate(r["event"][4:], (x, yv), textcoords="offset points", xytext=(6, 5), fontsize=8)
    ax.axvline(x, color="#bbb", ls=":", lw=0.8)
ax.set_xlabel("基准配置下每淹没节点洪量 V/node (m³)")
ax.set_ylabel("最优 Aponded (m²)")
ax.set_title("(b) 最优等效面积随事件洪量/节点变化", fontsize=11, loc="left")
ax.grid(alpha=0.3)

ax = axes[2]
if eff is not None and eff.size:
    _b = np.logspace(np.log10(3.0), np.log10(3.0e4), 29)
    ax.hist(np.clip(eff, 3.0, 3.0e4), bins=_b, color="#6aa84f",
            alpha=0.85, edgecolor="#4a7a35", linewidth=0.3)
    ax.set_xscale("log")
    ax.axvline(100, color="#c0504d", ls="--", lw=1.6, label="均匀 100 m²")
    ax.axvline(np.median(eff), color="#e69138", ls="-", lw=1.6,
               label="中位 %.0f m²" % np.median(eff))
    ax.legend(fontsize=8.5)
ax.set_xlabel("2D 反推等效面积 V/max(h_wet, 0.10 m) (m²)")
ax.set_ylabel("节点数")
ax.set_title("(c) 等效面积的节点间分布（0719）", fontsize=11, loc="left")
ax.grid(alpha=0.3)

fig.suptitle("H1 判据：单一均匀 Aponded 能否跨事件成立？（9 档 Aponded ∈ 0–400 m²，与图 3 同源）",
             fontsize=13, x=0.01, ha="left")
fig.tight_layout(rect=[0, 0, 1, 0.93])
fig.savefig(OUT + "/H1_等效蓄水_判据.png", dpi=160, facecolor="white")
print("-> output/H1_等效蓄水_判据.png")

# ---------------------------------------------------------------- 4) report --
_CFGLBL = [c[2] for c in CFGS]
LEVELS = [c[1] for c in CFGS]

L = ["# H1 判据：等效地表蓄水（Aponded）的事件依赖性与空间分布", "",
     "**数据**：%d 场事件（%s）× **9 个蓄水配置**（%s），"
     "其余设置完全相同（PS1 启停 1.60/0.90 固定）。"
     % (len(rows), "、".join("%.1f mm / %s" % (r["rain"][0], r["rain"][1])
                             for r in rows), " / ".join(_CFGLBL)), "",
     "⭐ **档位口径**：本报告与图 3、`output/P2_4_Aponded加密档位与连续最优点.md`"
     "使用**同一个 9 档网格**（0/25/50/75/100/150/200/300/400 m²）。"
     "旧版报告的 4 档网格 {0, 50, 100, 400} 已作废（间隔不对称：先 50 后 300）。", "",
     "**窗口口径**：模拟窗 = 有效雨首 − 3 h ~ 有效雨尾 + 10 h（受观测覆盖封顶）；"
     "评估窗 = 有效雨首 − 2 h ~ 有效雨尾 + 8 h，并被观测覆盖与"
     "\"退水见底后再抬升\"（未测暴雨）截断。指标只统计评估窗内配对样本。", "",
     "**判据**：复合评分 = 0.5 × (偏置校正 RMSE / 该事件最大 RMSE) + "
     "0.5 × (峰现时差 / 该事件最大时差)。单用 RMSE 会平局"
     "（0719 的 base 与 A100 都是 0.57 m），而两者峰现时差相差 8 h。", "",
     "## 1. 标定曲线（偏置校正 RMSE，m）", "",
     "| 事件 | 雨量 | 淹没节点 | 洪量 m³ | V/node m³ | "
     + " | ".join(_CFGLBL) + " |",
     "|---|---|---|---|---|" + "---|" * len(_CFGLBL)]
for r in rows:
    L.append("| %s | %.1f mm / %s | %s | %s | %.1f | %s |" % (
        r["event"], r["rain"][0], r["rain"][1],
        r["flooded"] if r["flooded"] else "-",
        "%.0f" % r["vol"] if r["vol"] else "-", r["v_per_node"],
        " | ".join(("%.2f" % r["cfg"][c[0]]["rmse"]) if c[0] in r["cfg"] else "-"
                   for c in CFGS)))
L += ["", "## 2. 复合评分（越小越好）与最优档", "",
      "| 事件 | " + " | ".join(_CFGLBL) + " | 网格最优 | 插值顶点 | 判据 |",
      "|---|" + "---|" * len(_CFGLBL) + "---|---|---|"]
for r in rows:
    L.append("| %s | %s | %.0f m² | %s | %s |" % (
        r["event"],
        " | ".join(("%.3f" % r["cfg"][c[0]]["score"]) if c[0] in r["cfg"] else "-"
                   for c in CFGS),
        r["grid_opt"],
        ("**%.1f m²**" % r["vertex"]) if r["kind"] == "interior"
        else "n/a",
        "内部最优" if r["kind"] == "interior" else "**%s**" % r["kind"]))
L += ["",
      "⚠️ **插值顶点只在「内部最优」时才是数值**：顶点由最优点及其左右邻三点作抛物线"
      "得到，须同时满足 ① 曲率 `a > 0`（开口向上）② 顶点落在三点区间**内部**。"
      "否则（`boundary-low` / `boundary-high` / `non-convex`）三点括住的是极大值，"
      "顶点无意义 —— 该事件的真实最优值在网格之外，**不能当数值引用**。", "",
      "## 3. 结论", ""]
L += ["- 网格最优配置：" + "；".join(
        "%s → %s（Aponded = %s m²）" % (r["event"], r["best"],
                                      r["cfg"].get(r["best"], {}).get("area"))
        for r in rows),
      "- **各事件的最优 Aponded 不一致** → 单一均匀值无法跨事件成立。",
      "- **能给连续最优点的事件**：" + "；".join(
          "%s %s" % (r["event"], ("→ %.1f m²" % r["vertex"])
                     if r["kind"] == "interior" else "（%s，不给值）" % r["kind"])
          for r in rows),
      "- ✅ **与图 3 同源**：`work/fig_e1.json` 由本脚本写出（见文末），"
      "图 3 与 P2-4 的 `work/fig_e1_dense.json` 都由此数据而来。", "",
      "## 4. 量化关系：最优等效面积随事件规模的变化", ""]

_ok = [r for r in rows if r["cfg"]]
# NEAR-OPTIMAL SET at two tolerances.  Both are ABSOLUTE increments of the
# normalised composite score (which lies in [0, 1] within each event, so
# +0.05 is "5 % of the score scale" -- the "+5 %" quoted in the manuscript
# and P2-4 -- and NOT "5 % of the best value"; the two formulas give
# different sets on the same data).
_brk = {}
for _r in _ok:
    _b = min(d["score"] for d in _r["cfg"].values())
    _n5 = sorted(d["area"] for d in _r["cfg"].values()
                 if d["score"] <= _b + 0.05)
    _na = sorted(d["area"] for d in _r["cfg"].values()
                 if d["score"] <= _b + 0.02)
    _brk[_r["event"]] = (_n5, _na)

L += ["| 事件 | 雨量 / 有效历时 | V/node (m³) | 网格最优 | 插值顶点 | "
      "近优集（评分 ≤ 最优 **+0.05**，即 \"+5 %\"） | 近优集（≤ 最优 +0.02，即 \"+2 %\"） |",
      "|---|---|---|---|---|---|---|"]
def _span(v):
    return ("%.0f m²" % v[0]) if v[0] == v[-1] else "%.0f–%.0f m²" % (v[0], v[-1])


for _r in sorted(_ok, key=lambda z: z["v_per_node"]):
    _n5, _na = _brk[_r["event"]]
    L.append("| %s | %.1f mm / %s | %.1f | **%.0f m²** | %s | %s（%d 档） | %s（%d 档） |" % (
        _r["event"], _r["rain"][0], _r["rain"][1], _r["v_per_node"],
        _r["grid_opt"],
        ("%.1f m²" % _r["vertex"]) if _r["kind"] == "interior" else "n/a",
        _span(_n5), len(_n5), _span(_na), len(_na)))

# how well can this event discriminate at all?  (a small spread means the
# "best" configuration is within the noise of the others)
_rm = {r["event"]: [d["rmse"] for d in r["cfg"].values()] for r in _ok}
_spread = {k: (max(v) - min(v)) / max(min(v), 1e-9) for k, v in _rm.items()}
_weak = [k for k, v in _spread.items() if v < 0.10]
if _weak:
    L += ["", "- ⚠️ **判别力提示**：%s 的各档 RMSE 差异 < 10%%（%s），"
          "其『最优档』落在档间噪声内，只能作趋势参考、不应当作定值。"
          % ("、".join(_weak),
             "；".join("%s: %.0f%%" % (k, 100 * _spread[k]) for k in _weak))]

if len(_ok) >= 3:
    _xs = np.array([r["v_per_node"] for r in _ok], float)
    _ys = np.array([r["cfg"][r["best"]]["area"] for r in _ok], float)
    _rx = np.argsort(np.argsort(_xs)).astype(float)
    _ry = np.argsort(np.argsort(_ys)).astype(float)
    _rho = float(np.corrcoef(_rx, _ry)[0, 1])
    L += ["", "- **单调性**：V/node 与最优面积的 Spearman 秩相关 ρ = %+.2f（%s）。"
          % (_rho, "完全单调" if _rho > 0.99 else
             ("单调但非严格" if _rho > 0.5 else "不单调"))]
    _pos = _ys > 0
    if _pos.sum() >= 2 and len(set(_ys[_pos])) > 1:
        _k, _c = np.polyfit(_xs[_pos], np.log(_ys[_pos]), 1)
        _lg = np.log(_ys[_pos])
        # R^2 in the space in which the fit was made (log), matching Fig. 3(c)
        _ss = 1.0 - ((_lg - (_c + _k * _xs[_pos])) ** 2).sum() / max(
            ((_lg - _lg.mean()) ** 2).sum(), 1e-9)
        L += ["- **描述性拟合**（仅 Aponded > 0 的事件、对数域最小二乘）："
              "`Aponded_eff ≈ %.2f·exp(%.4f·V/node)`（m²，V/node 取 m³；"
              "**R²(log) = %.2f**，即在对数域内的决定系数，与图 3c 一致）。"
              "⚠️ 只有 3 个 A>0 的事件（n = 3）且测点仍偏稀疏，此式**只作趋势描述，"
              "不可直接外推**。" % (np.exp(_c), _k, _ss)]
    L += ["- **曲线形状对照**（见第 1 节）：长历时事件（0719/0809）的 RMSE 随蓄水面积"
          "单调改善，短历时事件（0627）恰好相反 —— 漫溢与回灌在时间上高度重叠，"
          "滞蓄会把水拖到错误的时间。"]

if eff is not None and eff.size:
    L += ["- **空间离散度**：2D 反推的节点等效面积 p25–p75 = %.0f–%.0f m²（中位 %.0f m²，"
          "与均匀 100 m² 同量级），**p75/p25 = %.1f**（跨 11 种 2D 场为 8–26）→ 同一场事件内"
          "单一均匀值也必然在部分节点过度滞蓄、在另一些节点滞蓄不足。"
          "⚠️ 该量为**等效尺度代理**、非积水面积：绝对中位随 2D 场变化约 6 倍（57–362 m²），"
          "故**只用于节点间的相对比较**（定义与适用条件见 §5.4-7 与 `output/P2_5_A_eff适用条件.md`）。"
          % (np.percentile(eff, 25), np.percentile(eff, 75), np.median(eff),
             np.percentile(eff, 75) / max(np.percentile(eff, 25), 1e-9)),
          "- 图：`output/H1_等效蓄水_判据.png`；逐站明细：`work/multi_event/metrics.json`。"]
open(OUT + "/H1_等效蓄水_判据.md", "w", encoding="utf-8").write("\n".join(L) + "\n")
print("-> output/H1_等效蓄水_判据.md")

# --------------------------------- 5) the data source of Fig. 3 ----------
# This script is now the SINGLE generator of the E1 sweep: it writes the same
# 9-level result that scripts/fig3_e1_calibration.py reads, so the report and
# the figure can no longer drift apart.  (Previously the report came from
# scripts/_rec2.py on a 4-level grid while the figure came from
# scripts/_E1d_dense.py -> _E1d_make_fige1.py on a 9-level one.)
def _iso(t):
    return t.isoformat() if hasattr(t, "isoformat") else str(t)


_E1 = {}
for r in rows:
    _E1[r["event"]] = dict(
        cfg={k: v for k, v in r["cfg"].items()}, best=r["best"],
        flooded=r["flooded"], vol=r["vol"], v_per_node=r["v_per_node"],
        win=[_iso(r["win"][0]), _iso(r["win"][1])],
        grid_opt=r["grid_opt"], vertex=r["vertex"], curv=r["curv"],
        kind=r["kind"],
        # near = score <= best + 0.05 (the "+5 %" rule quoted in the paper);
        # near_strict = score <= best + 0.02 ("+2 %")
        near=r["near"], near_strict=r["near_strict"])
json.dump(_E1, open(WORK + "/fig_e1.json", "w", encoding="utf-8"),
          ensure_ascii=False, indent=1)
print("-> work/fig_e1.json  (%d events x %d levels)  <- Fig. 3 data source"
      % (len(_E1), len(CFGS)))
