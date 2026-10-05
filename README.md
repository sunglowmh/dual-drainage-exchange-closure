# dual-drainage exchange-closure test suite

Code and model inputs for: *Sewer–surface exchange closures in dual-drainage
flood modelling: event dependence, identifiability and decision consequences*
(submitted to **Environmental Modelling & Software**).

## Software (EMS template fields)

- **Name**: dual-drainage exchange-closure test suite
- **Developer and contact information**: *(to be completed: author names, corresponding author, e-mail)*
- **First available**: 2026
- **Hardware required**: any x86-64 workstation; all SWMM runs are single-node;
  the two-dimensional solver is single-process (each run of the 13.7 h 0719
  window takes ≈ 20–45 min at the 2 m grid)
- **Software required**: US EPA SWMM engine via pyswmm 2.x / swmm-toolkit 0.17;
  Python ≥ 3.10 with NumPy, SciPy, matplotlib (see `requirements.txt`)
- **Program language**: Python 3; SWMM 5.2 input files
- **Program size**: ≈ 5.7 kLOC across 26 scripts
- **Availability and cost**: free, MIT licence; archived at
  **https://github.com/…** *(to be completed)*, DOI **10.5281/zenodo.XXXXXX** *(to be completed)*

## What maps to what (paper → code)

| Paper element | Script(s) |
|---|---|
| Design-storm / load-family generation (Chicago r = 0.4) | `scripts/generate_swmm.py`, `scripts/run_multi_event.py`, `scripts/run_model.py` |
| 1D → 2D release series | `scripts/flood_series.py` |
| Two-dimensional diffusive-wave solver | `scripts/solver2d.py` |
| One-way / two-way coupling (closures, E3/E4 scenarios) | `scripts/couple_2way.py`, `scripts/node_head_series.py` |
| Metrics, observation audit helpers | `scripts/metrics_lib.py` |
| E1 equivalent-area inversion (A_eff) | `scripts/eval_aponded.py` |
| E2 identical-twin identifiability audit | `scripts/oracle_sampling.py` |
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
- `data/station_nodes.json` — sensor-to-node metadata (names, coordinates,
  matching distance 1.8–2.3 m). **Author check before publication**: contains
  station location descriptions; no raw measurements.
- The 1-min rain-gauge and manhole-level sensor records are **not included**:
  they are owned by the campus facilities operator and were provided under a
  data-sharing agreement that permits analysis and publication of derived
  quantities but not redistribution of the raw series (Option C, declared in
  the manuscript's Software and data availability section). Aggregated event
  statistics sufficient to reproduce the calibration and the audit are shipped
  with the Zenodo archive.

## Known environment note

With pyswmm 2.1.0 + swmm-toolkit 0.17.0 we observed `File Error 435` when
reading back SWMM binary `.out` files written by the same environment
(both for freshly written and for older, structurally valid files). Scripts
that only *write* `.out` via `pyswmm.Simulation` are unaffected; parsing
`.out` may require a swmm-toolkit/pyswmm version pair verified against each
other, or the documented plain-text `.rpt` summaries.

## Licence

MIT — see `LICENCE`.
