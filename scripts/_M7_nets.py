# -*- coding: utf-8 -*-
"""M7 -- conceptual (synthetic) sewer networks for the cross-system test of kappa(Lambda).

Three rule-generated networks, NONE of them digitised from a drawing:

  N1  "comb/grid"    10 laterals x 10 nodes + trunk, 80 m spacing, lateral slope 2 promille
  N2  "random tree"  seeded random dendritic network, edge lengths 60-150 m
  N3  "sparse flat"  6 laterals x 6 nodes + trunk, 150 m spacing, lateral slope 0.5 promille

Shared conventions (identical to the campus model):
  CMS / DYNWAVE / ALLOW_PONDING NO, Horton infiltration, n = 0.013,
  MaxDepth 2.5 m, SurDepth 0.30 m, Aponded 0, one INTENSITY rain gage,
  one FREE (FIXED-stage at invert) outfall, dendritic (every junction has
  exactly one outgoing conduit, no cycles).

Diameters: D = clip(scale * (0.15 + 0.16*sqrt(A_c)), 0.20, 2.00) m, where A_c is
the contributing area (ha) of the conduit's upstream node.  The global scale is
calibrated per network by the caller so that the flooding onset falls INSIDE the
16-member load family (the same condition that makes the campus family
informative); the claim under test is where the transition sits in Lambda, which
does not depend on where the onset sits in return-period space.
"""
import math
import random
import os
ROOT = os.environ.get("EMS_SWMM_BASE",
                   os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

OPTIONS = """FLOW_UNITS           CMS
INFILTRATION         HORTON
FLOW_ROUTING         DYNWAVE
LINK_OFFSETS         DEPTH
MIN_SLOPE            0
ALLOW_PONDING        NO
SKIP_STEADY_STATE    NO
START_DATE           09/15/2026
START_TIME           00:00:00
REPORT_START_DATE    09/15/2026
REPORT_START_TIME    00:00:00
END_DATE             09/15/2026
END_TIME             06:00:00
SWEEP_START          01/01
SWEEP_END            12/31
DRY_DAYS             0
REPORT_STEP          00:01:00
WET_STEP             00:01:00
DRY_STEP             00:05:00
ROUTING_STEP         0:00:05
VARIABLE_STEP        0.75
LENGTHENING_STEP     0
MIN_SURFAREA         1.167
MAX_TRIALS           12
HEAD_TOLERANCE       0.0015
SYS_FLOW_TOL         5
LAT_FLOW_TOL         5
MINIMUM_STEP         0.5
THREADS              4
"""


class Net:
    def __init__(self, name, label):
        self.name, self.label = name, label
        self.junc = {}      # name -> invert elevation (m)
        self.xy = {}
        self.edges = []     # (cname, upstream, downstream, length)
        self.subs = {}      # node -> (area_ha, imperv, width_m)
        self.outfall = None  # (name, elevation)
        self.out_conduit = None

    # ---------- construction helpers ----------
    def add_junction(self, j, elev, x, y, area, imperv, width):
        self.junc[j] = elev
        self.xy[j] = (x, y)
        self.subs[j] = (area, imperv, width)

    def add_edge(self, a, b, length):
        self.edges.append(("C%d" % (len(self.edges) + 1), a, b, length))

    # ---------- checks ----------
    def check_forest(self):
        out = {}
        for _, a, b, _l in self.edges:
            assert a not in out, "node %s has two outgoing conduits" % a
            out[a] = b
        assert set(out) == set(self.junc), "some junction has no outgoing conduit"
        for j in self.junc:                       # no cycles
            seen, v = set(), j
            while v in out:
                assert v not in seen, "cycle at %s" % v
                seen.add(v)
                v = out[v]
            assert v == self.outfall[0], "%s does not reach the outfall" % j
        # every edge downhill (or flat)
        for _, a, b, _l in self.edges:
            if b != self.outfall[0]:
                assert self.junc[a] >= self.junc[b] - 1e-9, "adverse slope %s->%s" % (a, b)
        assert len(self.edges) == len(self.junc), "E != V_junctions"

    # ---------- contributing areas ----------
    def downstream(self):
        d = {}
        for _, a, b, _l in self.edges:
            d[a] = b
        return d

    def contrib_ha(self):
        """Contributing subcatchment area (ha) of every junction, by leaf-to-root accumulation."""
        d = self.downstream()
        own = {j: self.subs[j][0] for j in self.junc}
        indeg = {j: 0 for j in self.junc}
        for _, a, b, _l in self.edges:
            if b in indeg:
                indeg[b] += 1
        acc = dict(own)
        queue = [j for j in self.junc if indeg[j] == 0]
        seen = 0
        while queue:
            j = queue.pop()
            seen += 1
            b = d[j]
            if b in acc:                       # trunk/outfall-side junction
                acc[b] += acc[j]
                indeg[b] -= 1
                if indeg[b] == 0:
                    queue.append(b)
        assert seen == len(self.junc), "topological walk missed nodes"
        return acc, d

    # ---------- INP text ----------
    def inp(self, rain, scale, title):
        self.check_forest()
        acc, d = self.contrib_ha()

        def diam(node):
            D = scale * (0.15 + 0.16 * math.sqrt(acc[node]))
            return min(2.0, max(0.20, round(D, 3)))

        L = []
        L.append("[TITLE]")
        L.append(";; M7 conceptual network %s -- %s (scale %.2f)" % (self.name, title, scale))
        L.append(";; rule-generated; NOT digitised from any drawing")
        L.append("[OPTIONS]")
        L.append(OPTIONS.rstrip("\n"))
        L.append("[EVAPORATION]")
        L.append("  CONSTANT        0.0")
        L.append("  DRY_ONLY        NO")
        L.append("[RAINGAGES]")
        L.append(";;Name           Format    Interval SCF      Source")
        L.append("RG1             INTENSITY 0:01     1.0      TIMESERIES RAIN")
        L.append("[SUBCATCHMENTS]")
        L.append(";;Name           RainGage   Outlet        Area     %Imperv  Width    %Slope   CurbLen")
        for j, (a, im, w) in sorted(self.subs.items()):
            L.append("%-14s RG1        %-12s %-8.4g %-8.4g %-8.4g %-8.4g 0" % ("S" + j[1:], j, a, im, w, 0.5))
        L.append("[SUBAREAS]")
        L.append(";;Subcatchment   N-Imperv N-Perv   S-Imperv S-Perv   %Zero    RouteTo  %Routed")
        for j in sorted(self.subs):
            L.append("%-14s 0.012    0.150    1.5      3.0      25.0     OUTLET   100.0" % ("S" + j[1:]))
        L.append("[INFILTRATION]")
        L.append(";;Subcatchment   MaxRate  MinRate  Decay    DryTime  MaxInfil")
        for j in sorted(self.subs):
            L.append("%-14s 76.2     12.7     4.0      7.0      0" % ("S" + j[1:]))
        L.append("[JUNCTIONS]")
        L.append(";;Name           Elevation  MaxDepth InitDepth SurDepth  Aponded")
        for j, e in sorted(self.junc.items()):
            L.append("%-14s %-9.4f 2.50     0        0.30      0" % (j, e))
        L.append("[OUTFALLS]")
        L.append(";;Name           Elevation  Type       Stage Data       Gated    Route To")
        nm, el = self.outfall
        L.append("%-14s %-9.4f FIXED      %-9.4f NO" % (nm, el, el))
        L.append("[CONDUITS]")
        L.append(";;Name           From Node      To Node        Length     Roughness  InOffset   OutOffset  InitFlow   MaxFlow")
        for cn, a, b, ln in self.edges:
            L.append("%-14s %-14s %-14s %-9.2f 0.0130    0          0          0          0" % (cn, a, b, ln))
        L.append("[XSECTIONS]")
        L.append(";;Link           Shape        Geom1            Geom2      Geom3      Geom4      Barrels")
        for cn, a, b, ln in self.edges:
            L.append("%-14s CIRCULAR     %-9.4f       0          0          0          1" % (cn, diam(a)))
        L.append("[COORDINATES]")
        for j, (x, y) in sorted(self.xy.items()):
            L.append("%-14s %-9.1f %-9.1f" % (j, x, y))
        nm, el = self.outfall
        L.append("%-14s %-9.1f %-9.1f" % (nm, 0.0, 0.0))
        L.append("[TIMESERIES]")
        L.append(";;Name           Date       Time       Value")
        for t, v in rain:
            mm, ss = divmod(int(t) * 60, 60)
            L.append("RAIN           09/15/2026 %02d:%02d:00  %.2f" % (mm // 60, mm % 60, v))
        L.append("[REPORT]")
        L.append("INPUT      NO")
        L.append("CONTROLS   NO")
        L.append("SUBCATCHMENTS ALL")
        L.append("NODES ALL")
        L.append("LINKS ALL")
        return "\n".join(L) + "\n"


# ------------------------------------------------------------------ N1: comb
def build_N1():
    """10 laterals x 10 nodes + trunk, 80 m spacing, lateral slope 2 permille."""
    net = Net("N1", "grid")
    SP = 80.0
    cover = 1.5

    def inv(x, y):
        return 0.0005 * x + 0.002 * y - cover

    for i in range(1, 11):                        # laterals, x = 80..800
        for k in range(1, 11):                    # along y = 80..800
            j = "J%d_%d" % (i, k)
            net.add_junction(j, inv(SP * i, SP * k), SP * i, SP * k, 0.64, 60.0, 80.0)
        net.add_junction("T%d" % i, inv(SP * i, 0.0), SP * i, 0.0, 0.64, 60.0, 80.0)
    for i in range(1, 11):
        for k in range(10, 1, -1):                # k=10 -> 9 -> ... -> 2 -> 1
            net.add_edge("J%d_%d" % (i, k), "J%d_%d" % (i, k - 1), SP)
        net.add_edge("J%d_1" % i, "T%d" % i, SP)   # lateral -> trunk
    for i in range(10, 1, -1):                    # trunk: high x -> low x
        net.add_edge("T%d" % i, "T%d" % (i - 1), SP)
    inv_out = inv(SP, 0.0) - 0.001 * SP           # last drop to the outfall
    net.outfall = ("OUT", inv_out)
    net.edges.append(("C%d" % (len(net.edges) + 1), "T1", "OUT", SP))
    return net


# ------------------------------------------------------------------ N2: random tree
def build_N2(seed=20261005):
    """Seeded random dendritic network: ~120 junctions, edges 60-150 m,
    invert RISING by 2 promille of the path length away from the outfall
    (upstream is higher), so every edge is downhill."""
    net = Net("N2", "random tree")
    rng = random.Random(seed)
    RISE = 0.002
    area, imperv, width = 0.50, 70.0, 70.0
    path = {0: 0.0}                               # internal id -> path length from outfall
    xy = {0: (0.0, 0.0)}
    parent = {}
    nid = 0
    while nid < 120:
        if nid == 0:
            p = 0                                 # the first node alone connects to the outfall
        else:
            p = rng.randrange(1, nid + 1)         # all others attach to existing junctions
        nid += 1
        L = rng.uniform(60.0, 150.0)
        ang = rng.uniform(0.0, 2 * math.pi)
        x, y = xy[p][0] + L * math.cos(ang), xy[p][1] + L * math.sin(ang)
        path[nid] = path[p] + L
        xy[nid] = (x, y)
        parent[nid] = p
    for v in range(1, 121):
        name = "J%d" % v
        net.add_junction(name, RISE * path[v], xy[v][0], xy[v][1], area, imperv, width)
        down = "J%d" % parent[v] if parent[v] else "OUT"
        net.add_edge(name, down, path[v] - path[parent[v]])
    inv_out = RISE * path[1] - 0.05               # below the single trunk node
    net.outfall = ("OUT", inv_out)
    return net


# ------------------------------------------------------------------ N3: sparse flat
def build_N3():
    """6 laterals x 6 nodes + trunk, 150 m spacing, lateral slope 0.5 promille."""
    net = Net("N3", "sparse flat")
    SP = 150.0
    cover = 1.5

    def inv(x, y):
        return 0.00025 * x + 0.0005 * y - cover

    for i in range(1, 7):
        for k in range(1, 7):
            j = "J%d_%d" % (i, k)
            net.add_junction(j, inv(SP * i, SP * k), SP * i, SP * k, 2.25, 50.0, 150.0)
        net.add_junction("T%d" % i, inv(SP * i, 0.0), SP * i, 0.0, 2.25, 50.0, 150.0)
    for i in range(1, 7):
        for k in range(6, 1, -1):
            net.add_edge("J%d_%d" % (i, k), "J%d_%d" % (i, k - 1), SP)
        net.add_edge("J%d_1" % i, "T%d" % i, SP)
    for i in range(6, 1, -1):
        net.add_edge("T%d" % i, "T%d" % (i - 1), SP)
    inv_out = inv(SP, 0.0) - 0.0005 * SP
    net.outfall = ("OUT", inv_out)
    net.edges.append(("C%d" % (len(net.edges) + 1), "T1", "OUT", SP))
    return net


def summary(net):
    a = sum(v[0] for v in net.subs.values())
    return dict(name=net.name, label=net.label, junctions=len(net.junc),
                conduits=len(net.edges), area_ha=round(a, 2))


if __name__ == "__main__":
    import io
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "gen", os.path.join(ROOT, "scripts", "generate_swmm.py"))
    gen = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(gen)
    rain = gen.chicago(2.0, 120, 1, 0.4)
    for build in (build_N1, build_N2, build_N3):
        net = build()
        print(summary(net))
        txt = net.inp(rain, 1.0, "test")
        print("  INP lines:", len(txt.split("\n")))
        io.open(os.path.join(ROOT, "work", "_M7_test_%s.inp") % net.name,
                "w", encoding="utf-8").write(txt)
    print("self-check done")
