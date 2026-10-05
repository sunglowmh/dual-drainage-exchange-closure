# -*- coding: utf-8 -*-
"""E3 load probes (replacing the 200/500/2000-a return-period probes).

Rationale
---------
The old high-return-period members (P = 200/500/2000 a) required extrapolating
the Shanghai design IDF (DB31/T 1043-2017) far outside its calibrated range,
which invites the objection that a 2000-year storm is being presented as a
design scenario.  The purpose of those members was never a design scenario: it
was to locate the SHAPE of the export-load curve, i.e. to test whether the
system has a hard export ceiling.

The same question can be asked without any reference to return period, by
scaling the 100-year hyetograph by a factor f.  That keeps the load axis
meaningful ("the design storm, multiplied") instead of making a climate-
statistics claim.

Design: f = 1.2 / 1.4 / 1.6.  These reproduce the load range previously covered
by P = 200 / 500 / 2000 a (compare work/_M6/analysis.json: ww was
75 920 / 89 174 / 111 656 m3).

Template = output/swmm_v7_50a.inp, identical to the design-storm family members
outside [TIMESERIES]; only the RAIN values change.

Self-check: f = 1.0 must reproduce output/swmm_M6_P100a.inp.
"""
import importlib.util
import json
import re
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
base = [v for _, v in gen.chicago(100.0, 120, 1, 0.4)]
assert len(base) == len(idx), (len(base), len(idx))
print("100 a hyetograph: %d points, peak %.2f mm/h, depth %.2f mm"
      % (len(base), max(base), sum(base) / 60.0))

print()
print("self-check: f = 1.0 vs output/swmm_M6_P100a.inp")
ref = series(open(OUT + "/swmm_M6_P100a.inp", encoding="utf-8").read())
print("   max|diff| = %.4f mm/h  (%d pts)"
      % (max(abs(a - b) for a, b in zip(base, ref)), len(ref)))

FAC = [1.2, 1.4, 1.6]
made = []
print()
for f in FAC:
    vals = [v * f for v in base]
    new = list(lines)
    for i, v in zip(idx, vals):
        p = new[i].split()
        new[i] = "%-15s %-10s %-9s %.2f" % (p[0], p[1], p[2], v)
    txt = tpl[:m.start(1)] + "\n".join(new) + tpl[m.end(1):]
    txt = re.sub(r"P=50a", "100a x%g load probe" % f, txt, count=1)
    tag = "L%d" % round(f * 100)
    path = "%s/swmm_M6_%s.inp" % (OUT, tag)
    open(path, "w", encoding="utf-8").write(txt)
    s = series(txt)
    made.append(dict(fac=f, tag=tag, file=path.split("/")[-1],
                     peak=round(max(s), 2), depth_mm=round(sum(s) / 60.0, 2)))
    print("  wrote %-24s x%-4g  peak %7.2f mm/h  depth %7.2f mm"
          % (path.split("/")[-1], f, max(s), sum(s) / 60.0))

json.dump(made, open(WORK + "/_M6/loadprobe.json", "w", encoding="utf-8"),
          indent=1, ensure_ascii=False)
print("-> work/_M6/loadprobe.json")
