# Aggregated statistics shipped with the archive

These files are the aggregated event statistics behind the analyses of the
paper.  They contain no raw sensor series (which remain with the data owner)
and no raw model output, only per-event / per-run summary quantities.

| File | What it is | Reproduces |
|---|---|---|
| `e1_grid_per_event.json` | Per-event, per-level composite score, bias-corrected RMSE, correlation r and peak-timing offset over the 9-area grid | Fig. 3; Table S4 row 1 |
| `e1_four_criteria.json` | Per-event optimum under each of four criteria + the transfer-cost matrix | Table S4 row 1; Section 4.1 |
| `e1_four_levels.json` | The coarse 4-level grid variant | Section 4.1 sensitivity |
| `e1_weight_sensitivity.json` | Composite-score optimum vs weight w(RMSE) = 0.3 / 0.5 / 0.7 | Table S4 row 2 |
| `e1_manhole_subset.json` | E1 repeated with the ponded area applied only to the 1 414 inspection-manhole nodes | Table S4 row 3 |
| `observation_audit_per_station.json` | Per-station observation-audit metrics (n, RMSE, r, NSE, observed/simulated peaks) per event and configuration | Section 3.4; Section 4.2 |
| `alarm_log_audit.json` | Alarm-log audit: alarm counts by type, logging gaps, truncated peaks | Section 2.3 |
| `e1_robustness.json` | Censored-station re-analysis (station 01 excluded), composite-weight sweep and the 4- vs 9-level normalisation check | Table S4 rows 4-5; Section 5.1 |
| `e1_robustness.json` | Censored-station re-analysis (station 01 excluded), composite-weight sweep and the 4- vs 9-level normalisation check | Table S4 rows 4-5; Section 5.1 |
| `water_balance_per_run.json` | Per-run 1D network water balance (overflow volume, number of flooding nodes) | Section 4.2 |

The script and run identifier behind every reported number are listed in
Supplementary Material S1.
