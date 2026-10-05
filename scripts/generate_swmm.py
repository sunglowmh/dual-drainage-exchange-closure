"""Generate a SWMM 5.2 .inp model for the storm-sewer network extracted from CAD.

All hydraulically-unknown attributes (diameter, invert, roughness) come from the
network builder as documented placeholders; land-use is assumed.
"""
import argparse
import json
import math
from datetime import datetime, timedelta
import os

import numpy as np
from scipy.spatial import cKDTree
ROOT = os.environ.get("EMS_SWMM_BASE",
                   os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

NET = os.path.join(ROOT, "work", "network_v2.json")

# --------------------------------------------------------------------------- #
# Pumping stations, located from the yellow "泵站" / "小包围泵站" text in the
# drawing and parameterised from the user's survey notes:
#   主泵站   at (575.7, 1094.9): gate from the lake into the pump chamber, then
#            a force main to the Guokang Rd municipal outfall.
#            75 kW, head 6 m, 2600 m3/h = 0.7222 m3/s.
#   小包围站 at (910.7, 962.6): small package station, holding tank
#            3.5 x 4.5 x 2.8 m (= 44.1 m3), pump on at 2.2 m, off at 0.9 m.
# --------------------------------------------------------------------------- #
PUMPS = [
    # 主泵站 (= 图中"内部泵站")：A27 留学生楼西侧，75 kW / 扬程 6 m / 2600 m3/h；
    # 南部河湖水经闸门进泵室。接纳 国康路雨水排口7 与 排口8。
    {"name": "PS1", "x": 575.7, "y": 1094.9, "area": 240.0,
     "q": 0.7222, "head": 6.0, "inlets": ["国康路雨水排口7", "国康路雨水排口8"],
     "discharge": "内部泵站接国康路排口", "fm_dia": 0.40},
    # 小包围泵站：国康路 11 号门东南侧绿地；蓄水池 3.5x4.5x2.8 m，
    # 2.2 m 启泵 / 0.9 m 停泵。接纳 国康路雨水排口6。
    {"name": "PS2", "x": 910.7, "y": 962.6, "area": 15.75,
     "q": 0.0500, "head": 4.0, "inlets": ["国康路雨水排口6"],
     "discharge": "小包围泵站接国康路排口", "fm_dia": 0.40},
]
# the small station's documented duty: 2.2 m switch-on and 0.9 m switch-off in a
# 2.8 m tank; the same water depths are applied to the main station, whose
# control levels were not surveyed
START_DEPTH, STOP_DEPTH = 2.20, 0.90


# --------------------------------------------------------------------------- #
# Shanghai design storm  (DB31/T 1043-2017《暴雨强度公式与设计雨型标准》)
#   q = 5544 (P^0.3 - 0.42) / (t + 10 + 7 lgP)^(0.82 + 0.07 lgP)   [L/(s*ha)]
#   i [mm/min] = q / 166.67 ,  valid for t = 5..180 min, P = 2..100 a
# --------------------------------------------------------------------------- #
def idf(P: float):
    """Return (A, b, n) of i(t) = A / (t + b)^n  with i in mm/min."""
    lgP = math.log10(P)
    A = 5544.0 * (P ** 0.3 - 0.42) / 166.67
    b = 10.0 + 7.0 * lgP
    n = 0.82 + 0.07 * lgP
    return A, b, n


def intensity_mmph(P: float, t_min: float) -> float:
    A, b, n = idf(P)
    return A / ((t_min + b) ** n) * 60.0


def chicago(P: float, dur: int = 120, step: int = 1, r: float = 0.4):
    """Chicago hyetograph -> [(offset_minutes, mm/h)] (1-minute steps)."""
    A, b, n = idf(P)
    tb = dur * r            # time from start to peak
    ta = dur - tb           # time from peak to end
    out = []
    for k in range(1, int(tb / step) + 1):
        x = (k * step) / r
        i = A * ((1 - n) * x + b) / ((x + b) ** (n + 1))   # mm/min
        out.append([k * step, i * 60.0])
    for k in range(1, int(ta / step) + 1):
        x = (k * step) / (1 - r)
        i = A * ((1 - n) * x + b) / ((x + b) ** (n + 1))
        out.append([int(tb) + k * step, i * 60.0])
    out.sort()
    return [(t, v) for t, v in out]


# --------------------------------------------------------------------------- #
def build_subcatchments_gis(nodes, lu, eligible=None, min_area=0.01):
    """Per-node catchment area / imperviousness taken from the GIS partition.

    `lu` is `landuse_v3.json`: an exact 2 m rasterisation of boundary.shp,
    green.shp, river.shp, buildings.shp clipped to the study area and split
    between the network nodes.

    Two merging passes, both preserving the total area exactly:
      1. nodes that may not receive runoff (DN150/200 inlet connections) hand
         their area to the nearest eligible node;
      2. nodes whose share is below `min_area` hand their area to the nearest
         node that already carries a real catchment — with 3 300 nodes over
         68.9 ha many nodes own only a few 2 m cells, and simply dropping them
         would silently discard ~10 % of the rainfall.
    """
    xy = np.array([[n["x"], n["y"]] for n in nodes])
    n = len(nodes)
    a = np.array([lu.get(x["id"], {}).get("area_ha", 0.0) for x in nodes], dtype=float)
    imp = np.array([x.get("imperv", 0.63) for x in nodes], dtype=float)

    def merge(src_mask, dst_idx):
        """Move the area of src_mask into dst_idx, area-weighting imperviousness."""
        if not src_mask.any():
            return
        num = a * imp
        np.add.at(num, dst_idx, a[src_mask] * imp[src_mask])
        np.add.at(a, dst_idx, a[src_mask])
        a[src_mask] = 0.0
        imp[:] = np.divide(num, a, out=np.full(n, 0.63), where=a > 0)

    if eligible is not None and len(eligible) < n:
        elig = np.array(sorted(eligible))
        drop = np.array([i for i in range(n) if i not in eligible])
        if len(drop):
            merge(np.isin(np.arange(n), drop),
                  elig[cKDTree(xy[elig]).query(xy[drop])[1]])
    for _ in range(4):
        small = (a > 0) & (a < min_area)
        if not small.any():
            break
        big = np.nonzero(a >= min_area)[0]
        if len(big) == 0:
            break
        src = np.nonzero(small)[0]
        merge(np.isin(np.arange(n), src),
              big[cKDTree(xy[big]).query(xy[src])[1]])

    subs = [(i, round(float(a[i]), 4)) for i in range(n) if a[i] > 0]
    return subs, float(a.sum()), a * 1e4, dict(enumerate(np.round(imp, 4)))


# --------------------------------------------------------------------------- #
# pipe sizing: rational method for the 2-year criterion (placeholder diameters)
# --------------------------------------------------------------------------- #
def q_Lps_ha(P: float, t_min: float) -> float:
    lgP = math.log10(P)
    return 5544.0 * (P ** 0.3 - 0.42) / ((t_min + 10.0 + 7.0 * lgP) ** (0.82 + 0.07 * lgP))


def size_diameter(Q, n=0.013, S=0.002, dmin=0.30, dmax=2.00):
    if Q <= 0:
        return dmin
    for D in [x / 1000.0 for x in range(300, 2001, 100)]:
        A = math.pi * D * D / 4.0
        R = D / 4.0
        Qc = (1.0 / n) * A * (R ** (2.0 / 3.0)) * math.sqrt(S)
        if Qc >= Q:
            return D
    return dmax


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--p", type=float, default=50.0, help="return period (years)")
    ap.add_argument("--dur", type=int, default=120, help="storm duration (min)")
    ap.add_argument("--out", default=os.path.join(ROOT, "output", "swmm_storm.inp"))
    args = ap.parse_args()

    net = json.load(open(NET, encoding="utf-8"))
    nodes = net["nodes"]
    links = net["links"]
    nid = {n["id"]: n for n in nodes}

    # ---- pumping stations -------------------------------------------------
    # Each station is  [in-campus outfall(s)] -> chamber(STORAGE) -PUMP->
    # throttle node -> force main -> municipal discharge outfall (YS9 / YS10).
    # The force main matters: without it SWMM sees only the small static lift and
    # the pump would run far beyond its rated duty.  With it the pipe friction
    # forms a proper system curve and the pump settles near the given duty point.
    REAL_OF = {o["name"]: o for o in json.load(
        open(os.path.join(ROOT, "work", "outfalls_real.json"), encoding="utf-8"))}
    ps_list, ps_out_ids = [], []
    off = 900100
    for ps in PUMPS:
        used = [n for n in nodes
                if n["outfall"] and n.get("name") in ps["inlets"]]
        if not used:
            print(f"WARNING: {ps['name']} found no in-campus outfall")
            continue
        for n in used:
            n["outfall"] = False           # becomes the chamber's inlet
            n["inlet_of_pump"] = ps["name"]
        # Chamber levels.  The documented duty is a 2.2 m switch-on and a 0.9 m
        # switch-off depth; the chamber invert is placed so that the STOP level
        # sits (a) 0.10 m below the lowest incoming pipe invert, so the sewer can
        # still discharge by gravity when the pump cuts out, and (b) at least
        # 0.15 m ABOVE the discharge outfall stage.  (b) matters: with the stop
        # level below the outfall the head across the pump goes negative and SWMM
        # extrapolates the three-point curve to absurd flows.
        low_inlet = min(n["invert"] for n in used)
        stage = REAL_OF[ps["discharge"]]["campus_invert"]
        stop = round(max(low_inlet - 0.10, stage + 0.15), 3)
        inv = round(stop - STOP_DEPTH, 3)
        g = max(n["ground"] for n in used)
        maxd = round(g - inv, 2)
        if maxd < START_DEPTH + 0.30:
            print(f"WARNING: {ps['name']} chamber only {maxd:.2f} m deep "
                  f"(switch-on depth is {START_DEPTH:.2f} m)")
        sid, tid, oid = off, off + 1, off + 2
        off += 3
        nodes.append({"id": sid, "x": ps["x"], "y": ps["y"], "level": -1,
                      "ground": round(g, 3), "invert": round(inv, 3),
                      "depth": maxd, "outfall": False, "anchor": False,
                      "imperv": 0.95, "deg": len(used) + 1, "diam": 0.5,
                      "storage": True, "area": ps["area"],
                      "name": ps["name"] + "集水池"})
        nodes.append({"id": tid, "x": ps["x"] + 3.0, "y": ps["y"] + 3.0,
                      "level": -1, "ground": round(g, 3),
                      "invert": round(stage, 3), "depth": 6.0,
                      "outfall": False, "anchor": False, "imperv": 0.95,
                      "deg": 2, "diam": ps["fm_dia"],
                      "name": ps["name"] + "压力出流节点"})
        do = REAL_OF[ps["discharge"]]
        d_inv = do["campus_invert"]
        nodes.append({"id": oid, "x": do["x"], "y": do["y"], "level": -2,
                      "ground": do["ground"], "invert": d_inv,
                      "depth": round(do["ground"] - d_inv, 3),
                      "outfall": True, "anchor": True, "imperv": 0.95,
                      "deg": 1, "diam": do["dia"], "name": do["name"],
                      "mode": "pump_discharge", "municipal_invert": d_inv,
                      "campus_invert": d_inv, "dia_mm": int(do["dia"] * 1000),
                      "well": do.get("well")})
        for n in used:
            d = math.hypot(n["x"] - ps["x"], n["y"] - ps["y"]) + 1.0
            links.append({"a": n["id"], "b": sid, "len": round(d, 2),
                          "diam": n["diam"], "diam_src": "pump", "rough": 0.013,
                          "pts": None})
        # force main, laid along the straight line to the municipal connection
        fm = math.hypot(do["x"] - ps["x"], do["y"] - ps["y"])
        links.append({"a": tid, "b": oid, "len": round(max(fm, 10.0), 2),
                      "diam": ps["fm_dia"], "diam_src": "forcemain",
                      "rough": 0.013, "pts": None})
        ps_list.append({"name": ps["name"], "sump": sid, "throttle": tid,
                        "out": oid, "q": ps["q"], "head": ps["head"],
                        "start": START_DEPTH, "stop": STOP_DEPTH,
                        "start_level": round(inv + START_DEPTH, 2),
                        "stop_level": stop, "area": ps["area"],
                        "maxd": maxd, "invert": round(inv, 3),
                        "fm": round(fm, 1), "fm_dia": ps["fm_dia"],
                        "n_up": len(used),
                        "inlets": [n["name"] for n in used],
                        "discharge": do["name"], "d_inv": d_inv})
        ps_out_ids.append(oid)
    for p in ps_list:
        print("  %s: %d 进水口 %s -> 集水池(inv %.2f, 深 %.2f, 启 %.2f m/停 %.2f m 深度"
              " = 标高 %.2f/%.2f) -PUMP-> 出水管 DN%d L=%.0f m -> %s (inv %.3f)" % (
                  p["name"], p["n_up"], p["inlets"], p["invert"], p["maxd"],
                  p["start"], p["stop"], p["start_level"], p["stop_level"],
                  p["fm_dia"] * 1000, p["fm"], p["discharge"], p["d_inv"]))

    lu = {r["id"]: r for r in json.load(
        open(os.path.join(ROOT, "work", "landuse_v3.json"), encoding="utf-8"))}
    subs, total_ha, area_m2, sub_imperv = build_subcatchments_gis(nodes, lu)
    print(f"subcatchments (pass 1): {len(subs)}  total area: {total_ha:.2f} ha")

    # ---- which nodes may receive surface runoff? ---------------------------
    # A node is a "collector" if it touches a conduit of DN300 or larger; DN150/200
    # are inlet-connection stubs which must not carry catchment runoff.
    ids = [n["id"] for n in nodes]
    pos = {i: k for k, i in enumerate(ids)}
    up = {i: [] for i in ids}
    maxd = {i: 0.0 for i in ids}
    for l in links:
        d = l.get("diam") or 0.0
        if d:
            maxd[l["a"]] = max(maxd[l["a"]], d)
            maxd[l["b"]] = max(maxd[l["b"]], d)
    # nodes without any drawing diameter: treat as collectors (unknown, usually mains)
    # NOTE: `eligible` holds array indices, not node ids (ids are no longer dense
    # after the topology clean-up)
    eligible = {k for k, i in enumerate(ids)
                if maxd[i] == 0.0 or maxd[i] >= 0.30}
    if not eligible:
        eligible = set(range(len(nodes)))
    subs, total_ha, area_m2, sub_imperv = build_subcatchments_gis(
        nodes, lu, eligible=eligible)
    print(f"collector nodes: {len(eligible)}; subcatchments (pass 2): {len(subs)}, "
          f"total area: {total_ha:.2f} ha")

    # upstream catchment area (used for reporting and for the pipes whose
    # diameter could not be read off the drawing)
    for l in links:
        up[l["b"]].append(l["a"])          # a -> b  (upstream -> downstream)
    order = sorted(ids, key=lambda i: -(nodes[pos[i]]["level"] or 0))
    acc_area = {}
    for i in order:
        acc_area[i] = area_m2[pos[i]] / 10000.0 + sum(acc_area.get(u, 0.0) for u in up[i])

    PSI, T_DES, P_DES = 0.65, 15.0, 2.0
    q_des = q_Lps_ha(P_DES, T_DES)
    n_dwg = 0
    for l in links:
        if l.get("diam_src") == "drawing":
            n_dwg += 1
        else:
            Q = PSI * q_des * acc_area.get(l["a"], 0.01) / 1000.0    # m3/s
            l["diam"] = size_diameter(Q)
        l["acc_ha"] = round(acc_area.get(l["a"], 0.0), 4)

    from collections import Counter
    print(f"diameters from drawing: {n_dwg} / {len(links)} "
          f"({100*n_dwg/len(links):.1f}%), computed: {len(links)-n_dwg}")
    print("diameter distribution:",
          dict(sorted(Counter(l["diam"] for l in links).items())))

    series = chicago(args.p, args.dur)
    rain_depth = sum(v / 60.0 for _, v in series)  # mm (1-min steps)
    print(f"design storm P={args.p:g}a, {args.dur}min -> {rain_depth:.1f} mm")

    start = datetime(2026, 9, 15, 0, 0)
    end = start + timedelta(hours=6)
    fmt = "%m/%d/%Y"
    fmt2 = "%H:%M:%S"

    L = []
    w = L.append
    w("[TITLE]")
    w(";; 同济大学校区雨水管网 SWMM 模型")
    w(";; 由 AutoCAD 图纸《同济排水平面图.dwg》自动提取重建")
    w(";; 地面高程/管底标高来自图纸高程点与排口标注识别；管径来自图纸标注")
    w("")
    w("[OPTIONS]")
    w("FLOW_UNITS           CMS")
    w("INFILTRATION         HORTON")
    w("FLOW_ROUTING         DYNWAVE")
    w("LINK_OFFSETS         DEPTH")
    w("MIN_SLOPE            0")
    w("ALLOW_PONDING        NO")
    w("SKIP_STEADY_STATE    NO")
    w(f"START_DATE           {start.strftime(fmt)}")
    w("START_TIME           00:00:00")
    w("REPORT_START_DATE    " + start.strftime(fmt))
    w("REPORT_START_TIME    00:00:00")
    w(f"END_DATE             {end.strftime(fmt)}")
    w(f"END_TIME             {end.strftime(fmt2)}")
    w("SWEEP_START          01/01")
    w("SWEEP_END            12/31")
    w("DRY_DAYS             0")
    w("REPORT_STEP          00:01:00")
    w("WET_STEP             00:01:00")
    w("DRY_STEP             00:05:00")
    w("ROUTING_STEP         0:00:05")
    w("VARIABLE_STEP        0.75")
    w("LENGTHENING_STEP     0")
    w("MIN_SURFAREA         1.167")
    w("MAX_TRIALS           12")
    w("HEAD_TOLERANCE       0.0015")
    w("SYS_FLOW_TOL         5")
    w("LAT_FLOW_TOL         5")
    w("MINIMUM_STEP         0.5")
    w("THREADS              4")
    w("")
    w("[EVAPORATION]")
    w("CONSTANT        0.0")
    w("DRY_ONLY        NO")
    w("")
    w("[RAINGAGES]")
    w(";;Name           Format    Interval SCF      Source")
    w(f"RG1             INTENSITY 0:01     1.0      TIMESERIES RAIN")
    w("")
    w("[SUBCATCHMENTS]")
    w(";;Name           RainGage   Outlet        Area     %Imperv  Width    %Slope   CurbLen")
    for i, ha in subs:
        n = nodes[i]
        name = f"S{n['id']}"
        width = round(1.7 * math.sqrt(ha * 10000.0), 1)
        imp = round(float(sub_imperv.get(i, n.get("imperv", 0.63))) * 100.0, 1)
        w(f"{name:<15}RG1        N{n['id']:<12}{ha:<9.4f}{imp:<9.1f}{width:<8.1f}"
          f"0.5      0")
    w("")
    w("[SUBAREAS]")
    w(";;Subcatchment   N-Imperv N-Perv   S-Imperv S-Perv   %Zero    RouteTo  %Routed")
    for i, ha in subs:
        w(f"S{nodes[i]['id']:<14}0.012    0.150    1.5      3.0      25.0     OUTLET   100.0")
    w("")
    w("[INFILTRATION]")
    w(";;Subcatchment   MaxRate  MinRate  Decay    DryTime  MaxInfil")
    for i, ha in subs:
        w(f"S{nodes[i]['id']:<14}76.2     12.7     4.0      7.0      0")
    w("")
    w("[JUNCTIONS]")
    w(";;Name           Elevation  MaxDepth InitDepth SurDepth  Aponded")
    for n in nodes:
        if n["outfall"] or n["invert"] is None or n.get("storage"):
            continue
        w(f"N{n['id']:<14}{n['invert']:<11.3f}{n['depth']:<10.2f}0        0.30      0")
    w("")
    w("[STORAGE]")
    w(";;Name           Elevation  MaxDepth InitDepth Shape       A1         A2         A3")
    for p in ps_list:
        w(f"N{p['sump']:<14}{p['invert']:<11.3f}{p['maxd']:<9.2f}0        FUNCTIONAL "
          f"{p['area']:<10.1f}0          0")
    w("")
    w("[OUTFALLS]")
    w(";;Name           Elevation  Type       Stage Data       Gated    Route To")
    w(";; The five 赤峰路 outlets discharge by gravity into the municipal storm")
    w(";; sewer (DN700): at every one the campus pipe invert is higher than the")
    w(";; municipal invert (+0.04 to +2.26 m), so the tailwater is a FIXED stage")
    w(";; set to the documented municipal pipe invert.  The three in-campus")
    w(";; 国康路 outfalls (YS6/YS7/YS8) feed the two pumping stations instead, and")
    w(";; the stations discharge through a force main to YS9 / YS10 on 国康路")
    w(";; (municipal DN1000).  Those two are written as FIXED-stage outfalls at")
    w(";; their documented invert.")
    n_free = n_pump = 0
    for n in nodes:
        if not n["outfall"] or n["invert"] is None:
            continue
        if n.get("mode") == "pumped":
            n_pump += 1
        stage = n.get("municipal_invert")
        if stage is None:
            w(f"N{n['id']:<14}{n['invert']:<11.3f}FREE                        NO")
        else:
            w(f"N{n['id']:<14}{n['invert']:<11.3f}FIXED      {stage:<8.3f}              NO")
        n_free += 1
    print(f"outfalls written: {n_free} "
          f"(FIXED at the documented municipal invert: "
          f"{sum(1 for n in nodes if n['outfall'] and n.get('municipal_invert') is not None)}; "
          f"pump discharges: {sum(1 for n in nodes if n.get('mode') == 'pump_discharge')})")
    if n_pump:
        print(f"WARNING: {n_pump} pumped outfall(s) were not picked up by a pump "
              f"station and are left as FREE discharge")
    w("")
    w("[CONDUITS]")
    w(";;Name           From Node      To Node        Length     Roughness  InOffset   OutOffset  InitFlow   MaxFlow")
    for k, l in enumerate(links):
        w(f"C{k+1:<13}N{l['a']:<14}N{l['b']:<14}{l['len']:<11.2f}"
          f"{l['rough']:<10.4f}0          0          0          0")
    w("")
    w("[XSECTIONS]")
    w(";;Link           Shape        Geom1            Geom2      Geom3      Geom4      Barrels")
    for k, l in enumerate(links):
        w(f"C{k+1:<13}CIRCULAR     {l['diam']:<16.3f}0          0          0          1")
    w("")
    w("[PUMPS]")
    w(";;Name           From Node      To Node        Pump Curve     Status   Startup  Shutoff")
    w(";; Startup/Shutoff are DEPTHS above the wet-well invert (NOT elevations).")
    w(";; The pump drives the throttle node, not the outfall directly: the force")
    w(";; main between them provides the pipe friction that makes the operating")
    w(";; point settle near the rated duty instead of running wide open.")
    for p in ps_list:
        w(f"{p['name']:<15}N{p['sump']:<14}N{p['throttle']:<14}{p['name'] + 'C':<15}"
          f"ON       {p['start']:<8.2f} {p['stop']:.2f}")
    w("")
    w("[CURVES]")
    w(";;Name           Type       X-Value    Y-Value")
    w(";;!! For pumps SWMM evaluates the curve as  flow = Y(X)  with the X column")
    w(";;   read as HEAD, not as flow (verified empirically with pump_sweep.py).")
    w(";;   Points are therefore written as (shutoff-free head, flow):")
    w(";;       0            -> maximum flow")
    w(";;       duty head    -> rated flow")
    w(";;       1.33 x head  -> zero flow")
    for p in ps_list:
        q, h = p["q"], p["head"]
        # SWMM allows the Type only on the first record of a curve; subsequent
        # records of the same curve carry just Name X Y.
        w(f"{p['name'] + 'C':<15}PUMP3      {0.0:<10.4f} {1.33 * q:.4f}")
        w(f"{p['name'] + 'C':<15}           {h:<10.4f} {q:.4f}")
        w(f"{p['name'] + 'C':<15}           {1.33 * h:<10.4f} 0.0000")
    w("")
    w("[CONTROLS]")
    w("; Pump duty is handled by the Startup/Shutoff depths in [PUMPS], which")
    w("; provides proper hysteresis (on above the start level, off below the stop")
    w("; level).  An ELSE-based RULE would re-trip the pump every time step at the")
    w("; setpoint, so no rule is emitted here.")
    w("")
    w("[COORDINATES]")
    w(";;Node           X-Coord            Y-Coord")
    for n in nodes:
        w(f"N{n['id']:<14}{n['x']:<19.3f}{n['y']:.3f}")
    w("")
    w("[VERTICES]")
    w(";;Link           X-Coord            Y-Coord")
    for k, l in enumerate(links):
        pts = l.get("pts") or []
        for p in pts[1:-1]:
            w(f"C{k+1:<14}{p[0]:<19.2f}{p[1]:.2f}")
    w("")
    w("[TIMESERIES]")
    w(";;Name           Date       Time       Value")
    base = start
    for t, v in series:
        ts = base + timedelta(minutes=t - 1)
        w(f"RAIN           {ts.strftime(fmt)} {ts.strftime(fmt2)}  {v:.2f}")
    w("")
    w("[REPORT]")
    w("INPUT      NO")
    w("CONTROLS   NO")
    w("SUBCATCHMENTS ALL")
    w("NODES ALL")
    w("LINKS ALL")
    w("")
    w("[TAGS]")
    w("")
    w("[MAP]")
    xs = [n["x"] for n in nodes]
    ys = [n["y"] for n in nodes]
    w(f"DIMENSIONS {min(xs):.3f} {min(ys):.3f} {max(xs):.3f} {max(ys):.3f}")
    w("Units      None")

    text = "\n".join(L)
    with open(args.out, "w", encoding="utf-8") as f:
        f.write(text)
    print("written", args.out)


if __name__ == "__main__":
    main()
