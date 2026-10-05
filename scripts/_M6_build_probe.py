# -*- coding: utf-8 -*-
"""M6 -- capacity probes: push the load far past the standard's validity range
just to locate the PLATEAU of the export-load curve, i.e. the effective export
capacity C_eff.  These runs are NOT used for any design statement.

Usage: python _M6_build_probe.py 200 500 2000
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

P_LIST = [float(a) for a in sys.argv[1:]] or [200.0, 500.0, 2000.0]

tpl = open(TEMPLATE, encoding="utf-8").read()
m = re.search(r"\[TIMESERIES\](.*?)(?=\n\[|\Z)", tpl, re.S)
lines = m.group(1).splitlines()
idx = [i for i, l in enumerate(lines)
       if l.strip() and not l.strip().startswith(";") and l.split()[0] == "RAIN"]


def series(txt):
    mm = re.search(r"\[TIMESERIES\](.*?)(?=\n\[|\Z)", txt, re.S)
    return [float(l.split()[3]) for l in mm.group(1).splitlines()
            if len(l.split()) >= 4 and l.split()[0] == "RAIN"]


for P in P_LIST:
    vals = [v for _, v in gen.chicago(P, 120, 1, 0.4)]
    new = list(lines)
    for i, v in zip(idx, vals):
        p = new[i].split()
        new[i] = "%-15s %-10s %-9s %.2f" % (p[0], p[1], p[2], v)
    txt = tpl[:m.start(1)] + "\n".join(new) + tpl[m.end(1):]
    txt = re.sub(r"P=50a", "probe P=%ga" % P, txt, count=1)
    tag = "P%ga" % P
    path = "%s/swmm_M6_%s.inp" % (OUT, tag)
    open(path, "w", encoding="utf-8").write(txt)
    s = series(txt)
    print("  wrote %-28s peak %8.2f mm/h  depth %8.2f mm"
          % (path.split("/")[-1], max(s), sum(s) / 60.0))
print("done")
