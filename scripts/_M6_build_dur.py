# -*- coding: utf-8 -*-
"""M6 -- the DURATION family: same return period, different storm duration.

The intensity family holds the duration fixed, so it cannot tell whether the
controlling quantity is the event VOLUME (V_in) or the event mean RATE
(V_in / T).  Varying the duration at a fixed return period separates the two:

  * if the response collapses on  V_in / (C T)  (a rate ratio), the duration
    family will move the state along the same axis as the intensity family;
  * if it collapses on  V_in / V_ref  (a volume ratio), it will not.

Usage: python _M6_build_dur.py 5 30 60 240 360
       (return period, then durations in minutes)
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

P = float(sys.argv[1]) if len(sys.argv) > 1 else 5.0
DURS = [int(x) for x in sys.argv[2:]] or [30, 60, 240, 360]

tpl = open(TEMPLATE, encoding="utf-8").read()
m = re.search(r"\[TIMESERIES\](.*?)(?=\n\[|\Z)", tpl, re.S)
lines = m.group(1).splitlines()
idx = [i for i, l in enumerate(lines)
       if l.strip() and not l.strip().startswith(";") and l.split()[0] == "RAIN"]

made = []
for dur in DURS:
    vals = [v for _, v in gen.chicago(P, dur, 1, 0.4)]
    new = list(lines)
    # the template has 120 rain rows; replace them with `dur` rows
    first = idx[0]
    hdr = [l for l in lines[:first] if True]
    rain = ["%-15s %-10s %-9s %.2f" % ("RAIN", "09/15/2026", "%02d:%02d:00" % (t // 60, t % 60), v)
            for t, v in enumerate(vals)]
    rest = lines[idx[-1] + 1:]
    newblock = lines[:first] + rain + rest
    body = "\n".join(newblock)
    txt = tpl[:m.start(1)] + body + tpl[m.end(1):]
    # extend the simulation so the event can drain (start + dur + 12 h)
    end_min = dur + 12 * 60
    txt = re.sub(r"(END_DATE\s+)\S+", r"\g<1>09/15/2026", txt, count=1)
    txt = re.sub(r"(END_TIME\s+)\S+", r"\g<1>%02d:%02d:00" % (end_min // 60, end_min % 60),
                 txt, count=1)
    txt = re.sub(r"P=50a", "P=%ga dur=%dmin" % (P, dur), txt, count=1)
    tag = "D%d" % dur
    path = "%s/swmm_M6_%s.inp" % (OUT, tag)
    open(path, "w", encoding="utf-8").write(txt)
    dep = sum(vals) / 60.0
    made.append(dict(tag=tag, P=P, dur_min=dur, depth_mm=round(dep, 2), peak=round(max(vals), 2)))
    print("  wrote %-22s P=%ga dur=%3d min  peak %7.2f mm/h  depth %7.2f mm"
          % (path.split("/")[-1], P, dur, max(vals), dep))

json.dump(made, open(WORK + "/_M6/duration_family.json", "w", encoding="utf-8"), indent=1)
print("-> work/_M6/duration_family.json")
