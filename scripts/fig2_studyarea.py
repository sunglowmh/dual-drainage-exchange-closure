# -*- coding: utf-8 -*-
"""Figure 2 -- study area: the 68.9 ha pumped campus, its drainage network, the
ten outlet chambers and the observation network, on the 2 m terrain.

(a) plan view                    (b) hydraulic boundary / outlet chain schematic

Drawing rules (revised 2026-09-29, second pass)
----------------------------------------------
* Only the conduits read from the as-built plan are drawn (n = 3 221, the links
  that carry a polyline).  The 50 synthetic connectors that wire each branch and
  each disconnected piece of the network to an outlet chamber are NOT drawn as
  pipes.  One dotted routing arrow per outlet is drawn instead, from the network
  node the plan wires that chamber to: that gap is the only place where the plan
  has no reach.
* ONE symbol per outlet class, so that (a) and (b) agree and the class does not
  rest on colour alone:
      ^  gravity outlet -> Chifeng Rd, municipal dn700         tags 1-5
      s  pumped route   -> Guokang Rd, municipal dn1000        tags 6-10
* Every outlet carries its tag, in the plan's own numbering (1-8) continued over
  the two named pumped chambers (9, 10).
* Roads are NAMED but not drawn.  Their positions still come from the survey,
  never from a basemap -- the names are offset onto the campus kerb:
      Chifeng Rd  the polyline through the five gravity outlets, extended along
                  its end tangents.  The plan names each of them as discharging
                  into Chifeng Rd, and the plan's own road-name label plots 2.1 m
                  off this line.
      Guokang Rd / Siping Rd / Miyun Rd   the northern / eastern / western campus
                  edge, lightly smoothed.  These are kerb lines, so the plan's
                  own street labels plot 14-28 m outside them, as expected.
  No centreline is plotted: any line long enough to read as a carriageway cuts
  through the buildings that line it, and the roads are only a frame of
  reference here.
* No path-effect halo anywhere: a halo would turn the text into vector outlines
  and the PDF would stop being searchable.  Labels sit on a white rounded box.
"""
import json
import sys
import os

ROOT = os.environ.get("EMS_SWMM_BASE",
                   os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection, PatchCollection
from matplotlib.lines import Line2D
from matplotlib.patches import Polygon as MPoly, FancyBboxPatch
from shapely.geometry import Point, LineString
from shapely.ops import unary_union
from shapely.prepared import prep

import figstyle as F
from gis_layers import load as gisload, TX, TY

F.use_style()
PT = 7.0                                   # the journal floor, honoured everywhere
BOX = dict(boxstyle="round,pad=0.16", facecolor="white", edgecolor="none",
           alpha=0.86)

NET = json.load(open(F.WORK + "/network_v2.json", encoding="utf-8"))
OUTF = json.load(open(F.WORK + "/outfalls_real.json", encoding="utf-8"))
STN = json.load(open(F.WORK + "/stations_shcoord.json", encoding="utf-8"))
META = json.load(open(F.WORK + "/surf2d_meta.json", encoding="utf-8"))
TXT = json.load(open(F.WORK + "/all_text.json", encoding="utf-8"))

X0, Y1, CELL = META["x0"], META["y1"], META["CELL"]
lake = np.load(F.WORK + "/lake_mask.npy")

# Design at the journal's double-column width (190 mm) so that nothing is scaled
# up at production: 7 pt in the file stays 7 pt on the page.
W, H = 189.6, 156.0
fig = plt.figure(figsize=(W * F.MM, H * F.MM))
ax = fig.add_axes([0.003, 0.201, 0.672, 0.770])      # re-set to the data aspect
axb = fig.add_axes([0.700, 0.297, 0.298, 0.672])

# ================================================================== basemap
bound = gisload("boundary")
bb = bound[0].bounds
PAD = 48
ax.set_xlim(bb[0] - PAD, bb[2] + 34)
ax.set_ylim(bb[1] - 26, bb[3] + 30)
ax.set_aspect("equal")
ax.set_xticks([])
ax.set_yticks([])
for s in ax.spines.values():
    s.set_visible(False)

# set_aspect("equal") letterboxes any axes box whose shape differs from the data
# aspect, which would push the panel title out of line with panel (b).  Size the
# box to the data instead, then axes fractions coincide with the drawn map.
AX_X, AX_W, AX_TOP = 0.003, 0.672, 0.971
_xr = ax.get_xlim()[1] - ax.get_xlim()[0]
_yr = ax.get_ylim()[1] - ax.get_ylim()[0]
ax.set_position([AX_X, AX_TOP - (AX_W * W) * (_yr / _xr) / H, AX_W,
                 (AX_W * W) * (_yr / _xr) / H])
MM_PER_M = (AX_W * W) / _xr
M_PER_PT = 0.3528 / MM_PER_M


def add_polys(layer, fc, ec="none", lw=0.0, z=2):
    patches = []
    for g in gisload(layer):
        ps = g.geoms if g.geom_type == "MultiPolygon" else [g]
        for p in ps:
            if p.geom_type == "Polygon" and p.area > 4.0:
                patches.append(MPoly(list(p.exterior.coords)))
    ax.add_collection(PatchCollection(patches, facecolor=fc, edgecolor=ec,
                                      linewidth=lw, zorder=z))


add_polys("green", "#DCEDD8", z=3)
add_polys("buildings", "#DEDEDE", z=2)

_gx = X0 + (np.arange(lake.shape[1]) + 0.5) * CELL
_gy = Y1 - (np.arange(lake.shape[0]) + 0.5) * CELL
ax.contourf(_gx, _gy, lake.astype(float), levels=[0.5, 1.5],
            colors=["#8FC8E8"], zorder=4)
ax.contour(_gx, _gy, lake.astype(float), levels=[0.5], colors=["#3F7EA8"],
           linewidths=0.6, zorder=5)

BLD_P = prep(unary_union([g.buffer(0) for g in gisload("buildings")]))

# =========================================== text placement (all in CAD metres)
BOXES = []          # placed label rectangles: (centre, u, v, half_w, half_h)


def _rect(c, s, fs, rot_deg=0.0, u=(1.0, 0.0)):
    u = np.asarray(u, float)
    v = np.array([-u[1], u[0]])
    hw = 0.5 * len(s) * fs * 0.53 * M_PER_PT + 1.2 * M_PER_PT
    hh = 0.5 * fs * 1.10 * M_PER_PT + 1.2 * M_PER_PT
    return np.asarray(c, float), u, v, hw, hh


def _overlap(a, b):
    (c1, u1, v1, w1, h1), (c2, u2, v2, w2, h2) = a, b
    for axx in (u1, v1, u2, v2):
        p = abs((c2 - c1) @ axx)
        if p > (w1 * abs(u1 @ axx) + h1 * abs(v1 @ axx)
                + w2 * abs(u2 @ axx) + h2 * abs(v2 @ axx)):
            return False
    return True


LABELS = []         # deferred text: (x, y, s, fs, colour, weight, rot, z)


def place(x, y, s, fs, colour, weight, z, cands, avoid_buildings=True):
    """Place a label at the first candidate offset that collides with nothing."""
    best = None
    for dx, dy in list(cands) + [(0, 0)]:
        c = (x + dx + TX, y + dy + TY)
        box = _rect(c, s, fs)
        clash = sum(1 for o in BOXES if _overlap(box, o))
        bad = 0
        if avoid_buildings:
            bad = sum(1 for t in np.linspace(-1, 1, 7)
                      if BLD_P.contains(Point(*(c + t * box[3] * np.array([1.0, 0.0])))))
        score = clash * 100 + bad
        if best is None or score < best[0]:
            best = (score, c, box)
        if score == 0:
            break
    score, c, box = best
    used = (c[0] - TX - x, c[1] - TY - y)
    if score:
        print("  ! label %-6r placed with score %d" % (s, score))
    BOXES.append(box)
    LABELS.append((c[0], c[1], s, fs, colour, weight, 0.0, z))
    return used


# ================================================================== network
POS = {str(n["id"]): (n["x"] + TX, n["y"] + TY) for n in NET["nodes"]}
DEG = {}
for l in NET["links"]:
    if l.get("pts"):
        DEG[str(l["a"])] = DEG.get(str(l["a"]), 0) + 1
        DEG[str(l["b"])] = DEG.get(str(l["b"]), 0) + 1
DRAWN = [k for k in POS if DEG.get(k, 0) > 0]
NODE_OF = np.array([POS[k] for k in DRAWN])
segs = [(POS[str(l["a"])], POS[str(l["b"])]) for l in NET["links"]
        if l.get("pts") and str(l["a"]) in POS and str(l["b"]) in POS]
ax.add_collection(LineCollection(segs, colors="#4A6B8A", linewidths=0.22,
                                 zorder=6))

# node class = the drawing symbol the node sits on (work/node_symbol_class.json)
NT = json.load(open(F.WORK + "/node_symbol_class.json", encoding="utf-8"))
CLS = {"circle_large": ("#2B4A68", "Node on a hatched circle symbol (1 414)"),
       "circle_small": ("#5B9BD5", "Node on a small plain circle (257)"),
       "rect": ("#8A4F9E", "Node on a grate symbol, 1.0 x 0.5 m (1 093)"),
       "none": ("#C9C9C9", "Node on no symbol (513)")}
for cls, (c, lab) in CLS.items():
    pts = [POS[i] for i in POS if NT.get(str(i)) == cls]
    ax.scatter([p[0] for p in pts], [p[1] for p in pts], s=0.85, c=c, marker="o",
               linewidths=0, zorder=8)
ax.add_patch(MPoly(list(bound[0].exterior.coords), fill=False, ec="black",
                   lw=0.9, zorder=9))

# ================================================== outlets & pumping stations
MODE = {"gravity": (F.OI["vermillion"], "^",
                    "Gravity outlet to Chifeng Rd (dn700)"),
        "to_pump": (F.OI["orange"], "s",
                    "Outlet on the pumped route to Guokang Rd (dn1000)"),
        "pump_discharge": (F.OI["orange"], "s", None)}
TAG = {"赤峰路雨水排口1": "1", "赤峰路雨水排口2": "2", "赤峰路雨水排口3": "3",
       "赤峰路雨水排口4": "4", "赤峰路雨水排口5": "5",
       "国康路雨水排口6": "6", "国康路雨水排口7": "7", "国康路雨水排口8": "8",
       "小包围泵站接国康路排口": "9", "内部泵站接国康路排口": "10"}
TAGCAND = {
    "1": [(24, -20), (24, 20), (-26, -18), (0, -32), (0, 32)],
    "2": [(24, -22), (24, 20), (-26, -20), (0, -34), (0, 32)],
    "3": [(24, -24), (-26, -22), (24, 20), (0, -36)],
    "4": [(-26, -22), (24, -24), (-26, 20), (0, -36)],
    "5": [(24, -22), (24, 20), (-26, -20), (0, 30)],
    "6": [(24, 18), (24, -20), (0, 28), (-26, 16)],
    "7": [(-28, 8), (-28, -6), (-28, 20), (-46, 0), (30, 0), (0, 26)],
    "8": [(-28, -6), (-28, 8), (-46, 0), (30, 0), (0, -26)],
    "9": [(24, 18), (0, 28), (24, -20), (-26, 18)],
    "10": [(-28, 30), (-28, 0), (30, 0), (0, 30), (0, -30)],
}
PS = {"PS1": (575.7, 1094.9), "PS2": (910.7, 962.6)}
PS_OFF = {"PS1": [(40, 0), (40, 16), (40, -16), (0, 32), (-64, 0)],
          "PS2": [(-30, -26), (30, -18), (0, -30), (0, 28)]}

# the wet wells first: they anchor the whole northern cluster
for nm, (x, y) in PS.items():
    place(x, y, nm, PT, "#1A1A1A", "normal", 14, PS_OFF[nm], avoid_buildings=False)

# ---- dotted routing links: how each chamber reaches the network ------------
ROUTE = {}
for e in OUTF:
    if e["mode"] == "gravity":
        d = np.hypot(NODE_OF[:, 0] - (e["x"] + TX), NODE_OF[:, 1] - (e["y"] + TY))
        i = int(np.argmin(d))
        ROUTE[e["name"]] = (POS[DRAWN[i]], float(d[i]), "out")
    else:
        px, py = PS[e["pump"]]
        ROUTE[e["name"]] = ((px + TX, py + TY), None,
                            "in" if e["mode"] == "to_pump" else "out_of_pump")

for e in OUTF:
    tgt, dist, kind = ROUTE[e["name"]]
    src = (e["x"] + TX, e["y"] + TY)
    # Scale the head and the end-shrinks with the link length.  These links run
    # 5-59 m, i.e. 0.6-6.5 mm on the page, so a fixed 1.6-3.6 pt shrink would
    # swallow a short link whole and leave a shapeless mark instead of an arrow
    # ending on its target.  Each shrink is capped at a quarter of the link, so
    # at least half the shaft is always drawn.
    Lpt = np.hypot(tgt[0] - src[0], tgt[1] - src[1]) * MM_PER_M / 0.3528
    kw = dict(color="#8A8A8A", lw=0.8, ls=(0, (1.9, 1.5)), zorder=6.5,
              head=max(1.6, min(4.0, 0.45 * Lpt)))
    sa = min(0.5, 0.12 * Lpt)                 # shrink at the tail
    sb = min(1.6, 0.25 * Lpt)                 # shrink at the head
    sw = min(3.6, 0.25 * Lpt)                 # shrink against a wet-well dot
    # gaussian_arrow(ax, xy, xytext) puts the HEAD at xy, so the DOWNSTREAM end
    # must be passed first.  Getting this backwards had every pumped-route arrow
    # in (a) and all five arrows in (b) pointing upstream.
    if kind == "out":                      # gravity chamber: node -> chamber
        F.gaussian_arrow(ax, src, tgt, shrinkA=sa, shrinkB=sb, **kw)
    elif kind == "in":                     # campus chamber -> wet well
        F.gaussian_arrow(ax, tgt, src, shrinkA=sa, shrinkB=sw, **kw)
    else:                                  # wet well -> pumped chamber
        F.gaussian_arrow(ax, src, tgt, shrinkA=sw, shrinkB=sa, **kw)

handles = []
for mode in ("gravity", "to_pump"):
    c, m, lab = MODE[mode]
    handles.append(Line2D([], [], marker=m, color="none", markerfacecolor=c,
                          markeredgecolor="white", markeredgewidth=0.4,
                          markersize=4.0, label=lab))
for e in OUTF:
    c, m, _ = MODE[e["mode"]]
    ax.scatter([e["x"] + TX], [e["y"] + TY], marker=m, s=20, c=c,
               edgecolors="white", linewidths=0.4, zorder=12)
    place(e["x"], e["y"], TAG[e["name"]], PT, "#7A3500", "normal", 14,
          TAGCAND[TAG[e["name"]]], avoid_buildings=False)
for nm, (x, y) in PS.items():
    ax.scatter([x + TX], [y + TY], marker="o", s=30, facecolors="white",
               edgecolors="black", linewidths=0.9, zorder=13)
handles.append(Line2D([], [], marker="o", color="none", markerfacecolor="white",
                      markeredgecolor="black", markersize=4.4,
                      label="Pumping station PS1 / PS2"))

# ================================================ observation network
MAN = [k for k in STN if k.startswith("86794807")]
ROADSENS = [k for k in STN if not k.startswith("86794807")]
ax.scatter([STN[k]["x"] for k in MAN], [STN[k]["y"] for k in MAN], marker="o",
           s=15, c=F.OI["blue"], edgecolors="white", linewidths=0.4, zorder=12.5)
ax.scatter([STN[k]["x"] for k in ROADSENS], [STN[k]["y"] for k in ROADSENS],
           marker="s", s=14, c=F.OI["sky"], edgecolors="white", linewidths=0.4,
           zorder=12.5)
handles.append(Line2D([], [], marker="o", color="none",
                      markerfacecolor=F.OI["blue"], markeredgecolor="white",
                      markersize=4, label="Ultrasonic manhole sensor (n = 12)"))
handles.append(Line2D([], [], marker="s", color="none",
                      markerfacecolor=F.OI["sky"], markeredgecolor="white",
                      markersize=4, label="Buried road-depth sensor (n = 5)"))
handles.append(Line2D([], [], color="#4A6B8A", lw=0.7,
                      label="Sewer conduit from the plan (n = 3 221)"))
handles.append(Line2D([], [], color="#8A8A8A", lw=0.8, ls=(0, (1.9, 1.5)),
                      label="Routing link to the nearest modelled node"))
for cls in ("circle_large", "circle_small", "rect", "none"):
    c, lab = CLS[cls]
    handles.append(Line2D([], [], marker="o", color="none", markerfacecolor=c,
                          markeredgecolor="none", markersize=3.4, label=lab))

# ============================================================== roads (named)
ROADFS = PT + 0.5             # roads are indicated by their names only
BOUND_CAD = np.array(bound[0].exterior.coords) - [TX, TY]


def _extend(p, d0, d1):
    p = np.asarray(p, float)
    u0 = (p[0] - p[1]) / np.linalg.norm(p[0] - p[1])
    u1 = (p[-1] - p[-2]) / np.linalg.norm(p[-1] - p[-2])
    return np.vstack([p[0] + d0 * u0, p, p[-1] + d1 * u1])


def _smooth(p, k=3, reps=2):
    p = np.asarray(p, float).copy()
    for _ in range(reps):
        q = p.copy()
        h = k // 2
        for i in range(h, len(p) - h):
            q[i] = p[i - h:i + h + 1].mean(0)
        p = q
    return p


E = {e["name"]: e for e in OUTF}
CHIFENG = _extend([(E[n]["x"], E[n]["y"]) for n in
                   ("赤峰路雨水排口5", "赤峰路雨水排口4", "赤峰路雨水排口3",
                    "赤峰路雨水排口2", "赤峰路雨水排口1")], 62.0, 62.0)
_id = lambda cx, cy: int(np.argmin(((BOUND_CAD - [cx, cy]) ** 2).sum(1)))
_iW, _iN, _iNE, _iSE = _id(-73.1, 615.0), _id(370.0, 1214.7), \
    _id(1006.2, 927.9), _id(705.1, 278.9)
GUOKANG = _extend(_smooth(BOUND_CAD[_iN:_iNE + 1]), 34.0, 44.0)
SIPING = _extend(_smooth(BOUND_CAD[_iNE:_iSE + 1][::-1]), 34.0, 34.0)
MIYUN = _extend(_smooth(np.vstack([BOUND_CAD[_iW:], BOUND_CAD[:_iN + 1]])),
                30.0, 26.0)
# 4th field = which side of the line the campus is on: the name belongs on that
# side, next to the kerb, so that the reader sees it against the study area and
# not floating in the neighbouring block.
ROADS = [("Chifeng Rd", CHIFENG, 0.30, +1), ("Guokang Rd", GUOKANG, 0.46, -1),
         ("Siping Rd", SIPING, 0.20, -1), ("Miyun Rd", MIYUN, 0.62, -1)]

for label, pts, frac, in_side in ROADS:
    p = np.asarray(pts, float) + [TX, TY]
    # No road geometry is drawn at all: any centreline long enough to read as a
    # carriageway necessarily cuts through the buildings that line it.  The
    # names alone place the four roads.
    # ---- label: nearest clear slot along the road, campus side first ----
    seg = np.diff(p, axis=0)
    L = np.hypot(seg[:, 0], seg[:, 1])
    cum = np.concatenate([[0.0], np.cumsum(L)])
    best = None
    for cand in [frac] + [frac + s * d for d in (0.03, 0.06, 0.09, 0.12, 0.15,
                                                 0.19, 0.24, 0.30)
                          for s in (+1, -1)]:
        if not 0.03 < cand < 0.97:
            continue
        t = cand * cum[-1]
        i = int(np.clip(np.searchsorted(cum, t) - 1, 0, len(seg) - 1))
        r = (t - cum[i]) / L[i]
        mid, d = p[i] + r * seg[i], seg[i] / L[i]
        ang = np.degrees(np.arctan2(d[1], d[0]))
        if ang > 90:
            ang -= 180
        if ang < -90:
            ang += 180
        for sgn in (in_side, -in_side):
            nrm = np.array([-d[1], d[0]]) * sgn
            for offs in (0.0, 5.0, 10.0, 16.0, 22.0, 30.0, 40.0, 52.0):
                lp = mid + offs * nrm
                box = _rect(lp, label, ROADFS, u=d)
                bad = sum(1 for t2 in np.linspace(-1, 1, 9)
                          if BLD_P.contains(Point(*(lp + t2 * box[3] * d))))
                bad += sum(1 for o in BOXES if _overlap(box, o))
                bad += 0 if offs <= 22.0 else (offs - 22.0) / 15.0
                if best is None or bad < best[0]:
                    best = (bad, lp, d, ang, box)
                if bad == 0:
                    break
            if best[0] == 0:
                break
        if best[0] == 0:
            break
    bad, lp, d, ang, box = best
    if bad:
        print("  road label %-11s placed with score %.1f" % (label, bad))
    BOXES.append(box)
    LABELS.append((lp[0], lp[1], label, ROADFS, "#3F3F3F", "normal", ang, 11))

# ---- draw every deferred label -------------------------------------------
for x, y, s, fs, col, wt, rot, z in LABELS:
    ax.text(x, y, s, fontsize=F.FS(fs), color=col, fontweight=wt,
            rotation=rot, rotation_mode="anchor", ha="center", va="center",
            zorder=z, bbox=BOX)

# ---- cross-check against the plan's own street-name labels ----------------
def _clean(s):
    import re
    return re.sub(r"\{\\\\[^;]*;", "", str(s)).replace("}", "").strip()


LINES = {nm: LineString(np.asarray(p, float) + [TX, TY])
         for nm, p, _, _ in ROADS}
print("  --- the plan's own street labels vs the constructed centrelines ---")
for key, nm in (("赤峰路", "Chifeng Rd"), ("国康路", "Guokang Rd"),
                ("四平路", "Siping Rd"), ("密云路", "Miyun Rd")):
    for t in TXT:
        if _clean(t.get("text", "")) == key:
            print("      %-11s label -> %5.1f m off the centreline"
                  % (nm, LINES[nm].distance(
                      Point(float(t["x"]) + TX, float(t["y"]) + TY))))
print("  --- outlet -> network routing distance (m) ---")
for e in OUTF:
    print("      %-26s %-15s %s" % (e["name"], e["mode"],
                                    "%.1f" % ROUTE[e["name"]][1]
                                    if ROUTE[e["name"]][1] else ROUTE[e["name"]][2]))

# ====================================================== decoration & labels
F.scale_bar(ax, 100, "100 m", x=0.035, y=0.028, fs=F.FS(PT))
F.north_arrow(ax, x=0.072, y=0.895, fs=F.FS(PT))
ax.text(0.0, 1.010, "a", transform=ax.transAxes, fontsize=F.FS(PT + 1),
        fontweight="bold", va="bottom", ha="left")
ax.text(0.028, 1.010, "Campus drainage network, outlet chambers and observations",
        transform=ax.transAxes, fontsize=F.FS(PT + 0.5), fontweight="bold",
        va="bottom", ha="left")

fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.004, 0.181),
           ncol=3, fontsize=F.FS(PT), frameon=False, handletextpad=0.5,
           columnspacing=1.6, labelspacing=0.30)
fig.text(0.004, 0.008,
         "Boundary, buildings, greens and river: campus GIS over a 2 m bare-earth DEM. Sewer nodes, conduits, outlet chambers and the two wet wells:\n"
         "digitised from the as-built combined drainage plan. Outlet tags follow the plan (1-8; 9 and 10 are its two named pumped chambers). Dotted lines\n"
         "are routing links, not pipes in the plan: the ten chambers sit 5-59 m from the nearest drawn conduit. Rainfall (1 min) from one campus gauge.",
         fontsize=F.FS(PT), color="#3F3F3F", ha="left", va="bottom",
         linespacing=1.45)

# ============================================== (b) outlet chain schematic
axb.axis("off")
axb.set_xlim(0, 100)
axb.set_ylim(0, 100)
axb.text(0, 101.5, "b", fontsize=F.FS(PT + 1), fontweight="bold", va="bottom",
         ha="left")
axb.text(6, 101.5, "Hydraulic boundary", fontsize=F.FS(PT + 0.5),
         fontweight="bold", va="bottom", ha="left")
axb.text(0, 97.5, "Where the campus model stops and the municipal\nsystem begins",
         fontsize=F.FS(PT), va="top", ha="left", color="#555555",
         linespacing=1.5)

axb.add_patch(FancyBboxPatch((4, 46), 92, 40,
                             boxstyle="round,pad=0.6,rounding_size=2",
                             linewidth=0.9, edgecolor="#2E6E9E",
                             facecolor="#E8F1F8"))
axb.text(50, 82.5, "CAMPUS   68.9 ha", fontsize=F.FS(PT), fontweight="bold",
         ha="center", va="center", color="#2E6E9E")
axb.text(50, 77.6, "3 277 nodes  ·  3 271 conduits (digitised)",
         fontsize=F.FS(PT), ha="center", va="center", color="#333333")
axb.text(50, 73.0, "SWMM model: 3 283 nodes  ·  3 278 links",
         fontsize=F.FS(PT), ha="center", va="center", color="#333333")
axb.plot([10, 90], [68.8, 68.8], color="#9A9A9A", lw=0.4, ls=":", zorder=1)

axb.text(6.5, 63.4, "5 gravity outlets  (1-5)", fontsize=F.FS(PT), ha="left",
         va="center")
axb.scatter([68, 74.5, 81, 87.5, 94], [63.4] * 5, marker="^", s=26,
            c=F.OI["vermillion"], edgecolors="white", linewidths=0.35, zorder=5)
axb.text(6.5, 55.0, "5 pumped-route outlets", fontsize=F.FS(PT), ha="left",
         va="center")
axb.scatter([68, 74.5, 81, 87.5, 94], [55.0] * 5, marker="s", s=22,
            c=F.OI["orange"], edgecolors="white", linewidths=0.35, zorder=5)

for x, nm, q in ((26, "PS1", "0.72 m³/s"), (74, "PS2", "0.05 m³/s")):
    axb.add_patch(FancyBboxPatch((x - 11, 29), 22, 13,
                                 boxstyle="round,pad=0.4,rounding_size=1.2",
                                 linewidth=0.9, edgecolor="#333333",
                                 facecolor="white", zorder=6))
    axb.text(x, 37.6, nm, fontsize=F.FS(PT), fontweight="bold", ha="center",
             va="center", zorder=7)
    axb.text(x, 32.6, q, fontsize=F.FS(PT), ha="center", va="center",
             color="#555555", zorder=7)
F.gaussian_arrow(axb, (26, 42.6), (26, 50.0), color="#8A8A8A", lw=0.8, head=5)
F.gaussian_arrow(axb, (74, 42.6), (74, 50.0), color="#8A8A8A", lw=0.8, head=5)
axb.text(29.0, 46.3, "7, 8", fontsize=F.FS(PT), color="#555555", ha="left",
         va="center")
axb.text(77.0, 46.3, "6", fontsize=F.FS(PT), color="#555555", ha="left",
         va="center")

axb.add_patch(FancyBboxPatch((4, 6), 92, 15,
                             boxstyle="round,pad=0.6,rounding_size=2",
                             linewidth=0.9, edgecolor="#7F7F7F",
                             facecolor="#F2F2F2"))
axb.text(50, 17.0, "MUNICIPAL STORM SEWER", fontsize=F.FS(PT),
         fontweight="bold", ha="center", va="center", color="#4D4D4D")
axb.text(26, 10.0, "Chifeng Rd   dn700", fontsize=F.FS(PT), ha="center",
         va="center")
axb.text(74, 10.0, "Guokang Rd   dn1000", fontsize=F.FS(PT), ha="center",
         va="center")

F.gaussian_arrow(axb, (8, 23.0), (8, 45.8), color=F.OI["vermillion"],
                 lw=0.8, head=5)
axb.text(5.2, 34.4, "gravity", fontsize=F.FS(PT), color=F.OI["vermillion"],
         ha="center", va="center", rotation=90)
F.gaussian_arrow(axb, (26, 22.4), (26, 28.6), color="#8A3B00", lw=0.8, head=5)
F.gaussian_arrow(axb, (74, 22.4), (74, 28.6), color="#8A3B00", lw=0.8, head=5)
axb.text(29.0, 25.5, "10", fontsize=F.FS(PT), color="#8A3B00", ha="left",
         va="center")
axb.text(77.0, 25.5, "9", fontsize=F.FS(PT), color="#8A3B00", ha="left",
         va="center")

F.save(fig, "Fig2_study_area", tight=False)
