"""Run a SWMM model and summarise surcharge / flooding (waterlogging) results.

Usage:  python run_model.py <model.inp> [csv_out]
"""
import re
import sys
import time
import os

from pyswmm import Simulation
ROOT = os.environ.get("EMS_SWMM_BASE",
                   os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

inp = sys.argv[1]
csv_out = sys.argv[2] if len(sys.argv) > 2 else None
rpt = inp.replace(".inp", ".rpt")
out = inp.replace(".inp", ".out")

t0 = time.time()
with Simulation(inp, rpt, out) as sim:
    sim.step_advance(3600)          # only sample hourly; routing step is internal
    for _ in sim:
        pass
print(f"[run] {inp}  finished in {time.time()-t0:.1f} s")

txt = open(rpt, errors="replace").read()


_STARS = re.compile(r"\*{5,}")


def section(title):
    """Return the text of a report section (from its title up to the next one)."""
    i = txt.find(title)
    if i < 0:
        return ""
    m1 = _STARS.search(txt, i)
    if not m1:
        return txt[i:]
    m2 = _STARS.search(txt, m1.end())
    return txt[i:m2.start() if m2 else len(txt)]


# continuity
m = re.search(r"Continuity Error \(%\) \.*\s*(-?[\d.]+)", txt)
print(f"[continuity] routing error = {m.group(1) if m else '?'} %")

# ---- node surcharge ----
body = section("Node Surcharge Summary")
lines = [l for l in body.splitlines() if re.match(r"\s*N\d+\s+JUNCTION", l)]
if lines:
    hrs = [float(l.split()[2]) for l in lines]
    above = [float(l.split()[3]) for l in lines]
    print(f"[surcharge] {len(lines)} nodes surcharged; "
          f"max height above crown = {max(above):.3f} m; "
          f"max hours = {max(hrs):.2f} h; total node-hours = {sum(hrs):.1f}")
else:
    print("[surcharge] none")

# ---- node flooding / ponding ----
body = section("Node Flooding Summary")
rows = []
for l in body.splitlines():
    p = l.split()
    if len(p) >= 7 and re.match(r"^N\d+$", p[0]):
        rows.append(p)

out_lines = []
if rows:
    tot_vol = sum(float(p[-2]) for p in rows)
    max_pond = max(float(p[-1]) for p in rows)
    max_hrs = max(float(p[1]) for p in rows)
    print(f"[flooding] {len(rows)} nodes flooded; total volume = {tot_vol:.3f} x10^6 l "
          f"({tot_vol*1000:.0f} m3); max ponded depth = {max_pond:.3f} m; "
          f"max duration = {max_hrs:.2f} h")
    print("\n  top 25 flooding nodes (node, hours, max rate CMS, total 10^6 l, ponded m)")
    for p in sorted(rows, key=lambda x: -float(x[-1]))[:25]:
        print(f"   {p[0]:<8} {p[1]:>5} {p[2]:>9} {p[-2]:>10} {p[-1]:>9}")
    if csv_out:
        with open(csv_out, "w", encoding="utf-8") as f:
            f.write("node,hours_flooded,max_rate_CMS,total_volume_1e6l,ponded_depth_m\n")
            for p in rows:
                f.write(",".join(p) + "\n")
        print(f"  -> written {csv_out}")
else:
    print("[flooding] no node flooding reported")

# ---- outfall loading ----
body = section("Outfall Loading Summary")
ofl = []
for l in body.splitlines():
    p = l.split()
    if len(p) >= 5 and re.match(r"^N\d+$", p[0]):
        ofl.append(p)
if ofl:
    tot = sum(float(p[-1]) for p in ofl)
    print(f"\n[outfalls] {len(ofl)} outfalls; total discharged = {tot:.3f} x10^6 l "
          f"({tot*1000:.0f} m3); peak = {max(float(p[3]) for p in ofl):.3f} CMS")
