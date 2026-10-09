# dual-drainage exchange-closure test suite

Code and model inputs for: *Sewer–surface exchange closures in dual-drainage
flood modelling: event dependence, identifiability and decision consequences*
(submitted to **Environmental Modelling & Software**).

## Software (EMS template fields)

- **Name**: dual-drainage exchange-closure test suite
- **Developer and contact information**: Maohui Zheng, Tongji University
  (zmh@tongji.edu.cn; authors: M. Zheng, C. Liu, K. Zhang)
- **First available**: 2026
- **Hardware required**: any x86-64 workstation; all SWMM runs are single-node;
  the two-dimensional solver is single-process (each run of the 13.7 h 0719
  window takes ≈ 20–45 min at the 2 m grid)
- **Software required**: US EPA SWMM engine via pyswmm 2.x / swmm-toolkit 0.17;
  Python ≥ 3.10 with NumPy, SciPy, matplotlib (see `requirements.txt`)
- **Program language**: Python 3; SWMM 5.2 input files
- **Program size**: ≈ 6.1 kLOC across 29 scripts
- **Availability and cost**: free, MIT licence; archived at
  **https://github.com/sunglowmh/dual-drainage-exchange-closure** (version v1.0.1).
  The repository is the archive of record.

## What maps to what (paper → code)

| Paper element | Script(s) |
|---|---|
| Design-storm / load-family generation (Chicago r = 0.4) | `scripts/generate_swmm.py`, `scripts/run_multi_event.py`, `scripts/run_model.py` |
| 1D → 2D release series | `scripts/flood_series.py` |
| Two-dimensional diffusive-wave solver | `scripts/solver2d.py` |
| One-way / two-way coupling (closures, E3/E4 scenarios) | `scripts/couple_2way.py`, `scripts/node_head_series.py` |
| Metrics, observation audit helpers | `scripts/metrics_lib.py` |
| E1 robustness (Table S4): four-criterion transfer costs | `scripts/e1_four_criteria_matrix.py` |
| E1 robustness (Table S4): composite-score weight sensitivity | `scripts/e1_weight_sensitivity.py` |
| E1 robustness (Table S4): manhole-only ponding-area sensitivity | `scripts/e1_manhole_subset_sensitivity.py`, `data/node_symbol_class.json` |
| E1 equivalent-area inversion (A_eff) | `scripts/eval_aponded.py` |
| E2 identical-twin identifiability audit | `scripts/identical_twin_sampling.py` |
| E3 controlled load family (13 design storms + load probes) | `scripts/_M6_build_family.py`, `_M6_build_loadprobe.py`, `_M6_build_dur.py`, `_M6_build_probe.py`, `_M6_final.py` |
| Conceptual networks N1–N3 (cross-system replication) | `scripts/_M7_nets.py`, `scripts/_M7_family.py` |
| Figures 1–8 | `scripts/figstyle.py`, `scripts/fig1_framework.py` … `scripts/fig8_M7_networks.py` |

`_M6*` = controlled load family (experiment E3); `_M7*` = conceptual networks;
the names are kept for traceability with the supplementary reproduction index.

## Layout expected at run time

Scripts read and write a workspace with this layout (the archive ships two
example model inputs under `data/`):

```
<EMS_SWMM_BASE>/            # set the EMS_SWMM_BASE environment variable
├── work/                   # floodts_*.json, obs_*.json (restricted), products
├── output/                 # SWMM .inp/.rpt/.out, figures
└── scripts/                # this archive
```

Without `EMS_SWMM_BASE` the scripts fall back to the repository's parent
directory. Example: regenerate the two-dimensional baseline of the paper —

```
python scripts/couple_2way.py 20260719 0 oneway        # one-way, full window
```

## Data policy

- `data/swmm_v7_50a.inp` — the 3 277-node model template (design storm 50 a);
  all event and family members are derived from it by the scripts above.
- `data/swmm_20260719_base.inp` — the 2026-07-19 observed-event model.
- `data/station_nodes.json` — sensor-to-node metadata (station identifiers,
  coordinates, matching distance 1.8–2.3 m); derived quantities only, no raw
  measurements.
- The 1-min rain-gauge and manhole-level sensor records are **not included**:
  they are owned by the campus facilities operator and were provided under a
  data-sharing agreement that permits analysis and publication of derived
  quantities but not redistribution of the raw series (Option C, declared in
  the manuscript's Software and data availability section). Aggregated event
  statistics sufficient to reproduce the calibration and the audit are shipped
  with the repository under `data/aggregates/` (see below).

## Archive contents (`data/aggregates/`)

Aggregated event statistics — no raw sensor series. See
`data/aggregates/MANIFEST.md` for the full list.

- `e1_grid_per_event.json` — per-event, per-level composite score, RMSE, r and
  peak-timing offset over the 9-area grid (Fig. 3; Table S4 row 1)
- `e1_four_criteria.json` — per-event optimum under four criteria + transfer-cost
  matrix (Section 4.1)
- `e1_four_levels.json`, `e1_weight_sensitivity.json`, `e1_manhole_subset.json`
  — E1 robustness variants (Table S4)
- `observation_audit_per_station.json` — per-station observation-audit metrics
- `alarm_log_audit.json` — alarm-log audit (counts, gaps, truncated peaks)
- `water_balance_per_run.json` — per-run 1D network water balance

The script and run identifier behind every reported number are listed in
Supplementary Material S1.

## Known environment note

With pyswmm 2.1.0 + swmm-toolkit 0.17.0 we observed `File Error 435` when
reading back SWMM binary `.out` files written by the same environment
(both for freshly written and for older, structurally valid files). Scripts
that only *write* `.out` via `pyswmm.Simulation` are unaffected; parsing
`.out` may require a swmm-toolkit/pyswmm version pair verified against each
other, or the documented plain-text `.rpt` summaries.

## Licence

MIT — see `LICENCE`.
