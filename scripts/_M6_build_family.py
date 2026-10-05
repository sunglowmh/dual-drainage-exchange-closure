# -*- coding: utf-8 -*-
"""M6 -- build the controlled design-storm family (corrected generator).

Uses the SAME generator as the model itself (`generate_swmm.chicago`, the
Shanghai standard DB31/T 1043-2017), so a rebuilt member is identical to an
existing one up to the rain block.

Template = output/swmm_v7_50a.inp.  Its JUNCTIONS / CONDUITS / XSECTIONS /
SUBCATCHMENTS / SUBAREAS / PUMPS / STORAGE / OUTFALLS / OPTIONS / COORDINATES
sections are byte-identical to v7_2a and v7_5a (verified), i.e. the only thing
that changes across the family is the hyetograph.  The family therefore varies
the LOAD alone.

Self-check: P = 2, 5 and 50 must reproduce the rain series of swmm_v7_{2,5,50}a.
"""
import importlib.util
import json
import re
import sys
import os
ROOT = os.environ.get("EMS_SWMM_BASE",
                   os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

OUT = os.path.join(ROOT, "output")
WORK = os.path.join(ROOT, "work")
TEMPLATE = OUT + "/swmm_v7_50a.inp"

spec = importlib.util.spec_from_file_location(
    "gen", os.path.join(ROOT, "scripts", "generate_swmm.py"))
gen = importlib.util.module_from_spec(spec)
spec.loader.exec_module(gen)
print("generator loaded; i(2a, 10 min) = %.1f mm/h" % gen.intensity_mmph(2, 10))


def series(txt):
    m = re.search(r"\[TIMESERIES\](.*?)(?=\n\[|\Z)", txt, re.S)
    v = []
    for l in m.group(1).splitlines():
        q = l.split()
        if len(q) >= 4 and q[0] == "RAIN":
            v.append(float(q[3]))
    return v


tpl = open(TEMPLATE, encoding="utf-8").read()
m = re.search(r"\[TIMESERIES\](.*?)(?=\n\[|\Z)", tpl, re.S)
lines = m.group(1).splitlines()
idx = [i for i, l in enumerate(lines)
       if l.strip() and not l.strip().startswith(";") and l.split()[0] == "RAIN"]
tpl_s = [float(lines[i].split()[3]) for i in idx]
prefix = lines[idx[0]].split()[:3]          # RAIN  09/15/2026  00:00:00
print("template: %d points, peak %.2f mm/h, depth %.2f mm"
      % (len(tpl_s), max(tpl_s), sum(tpl_s) / 60.0))

print()
print("self-check against the members already on disk:")
for P, tag in ((2.0, "v7_2a"), (5.0, "v7_5a"), (50.0, "v7_50a")):
    ref = [v for _, v in gen.chicago(P, 120, 1, 0.4)]
    act = series(open("%s/swmm_%s.inp" % (OUT, tag), encoding="utf-8").read())
    if len(ref) == len(act):
        print("  P=%-4g vs %-8s max|diff| = %.4f mm/h   (%d pts)"
              % (P, tag, max(abs(a - b) for a, b in zip(ref, act)), len(ref)))
    else:
        print("  P=%-4g vs %-8s LENGTH MISMATCH %d vs %d" % (P, tag, len(ref), len(act)))

P_LIST = ([float(x) for x in sys.argv[1:]] if len(sys.argv) > 1
          else [0.25, 0.5, 1.0, 3.0, 10.0, 20.0, 100.0])
made = []
print()
for P in P_LIST:
    vals = [v for _, v in gen.chicago(P, 120, 1, 0.4)]
    assert len(vals) == len(idx), (P, len(vals), len(idx))
    new = list(lines)
    for i, v in zip(idx, vals):
        p = new[i].split()
        new[i] = "%-15s %-10s %-9s %.2f" % (p[0], p[1], p[2], v)
    body = "\n".join(new)
    txt = tpl[:m.start(1)] + body + tpl[m.end(1):]
    txt = re.sub(r"P=50a", "P=%ga" % P, txt, count=1)
    tag = "%g" % P
    path = "%s/swmm_M6_P%sa.inp" % (OUT, tag)
    open(path, "w", encoding="utf-8").write(txt)
    s = series(txt)
    made.append(dict(P=P, tag="P" + tag, file=path.split("/")[-1],
                     peak=round(max(s), 2), depth_mm=round(sum(s) / 60.0, 2)))
    print("  wrote %-26s P=%-6g peak %7.2f mm/h  depth %7.2f mm"
          % (path.split("/")[-1], P, max(s), sum(s) / 60.0))

json.dump(made, open(WORK + "/_M6/family.json", "w", encoding="utf-8"), indent=1)
print("-> work/_M6/family.json")
