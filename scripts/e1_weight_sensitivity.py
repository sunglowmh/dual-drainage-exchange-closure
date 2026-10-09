# -*- coding: utf-8 -*-
"""v3 review F4b -- composite-score weight sensitivity (pure re-analysis).

The composite score S = w * rmse/rmse_max + (1-w) * |dT|/|dT|_max uses w = 0.5
("arbitrary weighting", conceded in the manuscript).  Re-derive, from
work/fig_e1_dense.json (the same per-level metrics behind Table 3):

  (1) the per-event optimum under w = 0.3 / 0.5 / 0.7 (w on the RMSE term);
  (2) the +0.05 near-optimal set width under each w;
  (3) the two end-member transfer costs quoted in the paper
      (0809 -> 0627 and 0627 -> 0809, composite score) under each w.

If optima and costs barely move, one sentence in the manuscript closes the
issue; no new figure, no new table in the main text.
"""
import json

import numpy as np

WORK = r"E:/work/tongji/swmm_model/work"
D = json.load(open(WORK + "/fig_e1_dense.json", encoding="utf-8"))
EVENTS = ["20250627", "20250730", "20260719", "20260809"]
SHORT = {"20250627": "0627", "20250730": "0730", "20260719": "0719", "20260809": "0809"}

res = {}
for w in (0.3, 0.5, 0.7):
    opt, near = {}, {}
    for ev in EVENTS:
        r = D[ev]
        xs, per = r["xs"], r["per"]
        rm = np.array([per[str(x)]["rmse"] for x in xs])
        dtt = np.array([per[str(x)]["dt"] for x in xs])
        s = w * rm / rm.max() + (1 - w) * dtt / dtt.max()
        i = int(np.argmin(s))
        opt[ev] = xs[i]
        near[ev] = [x for x, v in zip(xs, s) if v <= s[i] + 0.05]
    # end-member transfer costs: apply donor's optimum, score on recipient
    def cost(donor, recv):
        sd = opt[donor]
        rm = np.array([D[recv]["per"][str(x)]["rmse"] for x in D[recv]["xs"]])
        dtt = np.array([D[recv]["per"][str(x)]["dt"] for x in D[recv]["xs"]])
        s = w * rm / rm.max() + (1 - w) * dtt / dtt.max()
        j = D[recv]["xs"].index(sd)
        own = s[int(np.argmin(s))]
        return 100.0 * (s[j] / own - 1.0)
    res[w] = dict(opt=opt, near=near,
                  cost_0809_to_0627=cost("20260809", "20250627"),
                  cost_0627_to_0809=cost("20250627", "20260809"))

for w in (0.3, 0.5, 0.7):
    r = res[w]
    print("w(RMSE)=%.1f  optima: %s" % (w, {SHORT[k]: v for k, v in r["opt"].items()}))
    print("           near-opt widths: %s"
          % {SHORT[k]: (max(v) - min(v)) for k, v in r["near"].items()})
    print("           transfer 0809->0627 %+.0f %%  |  0627->0809 %+.0f %%"
          % (r["cost_0809_to_0627"], r["cost_0627_to_0809"]))

json.dump(res, open(WORK + "/fig_e1_weight_sens.json", "w", encoding="utf-8"))
print("saved work/fig_e1_weight_sens.json")
