# Evaluator (maintainer-only — not agent-facing)

Scores a submitted `run_qp.py` against held-out ARPES.

```
hidden/L{1,2,3}/{K,Mg}/   held-out metals per level: element_config.json, grid.csv,
                          arpes_reference.csv (scoring target), + level data files
gold/<El>_gold.csv        gold occupied bands (QP + KS) from reference/gold_runner.jl;
                          used for band-count cap and the KS-baseline audit
source_data/<El>/         raw grids + ARPES (the build source for hidden/ and packet)
validate_submission.py    the scorer
```

## Run

```bash
python validate_submission.py --submission-dir <dir> --level 2 --json result.json
```

`--level {1,2,3}` selects `hidden/L<level>`. The validator runs
`python <dir>/run_qp.py --element-config ... --grid ... --out ...` per hidden
element, then scores.

## Scoring guards (why low RMSE means "did the physics")

1. **Flooding guard** — predictions are capped to the gold occupied-band count
   per k-point. A submission that floods the energy window (to make nearest-band
   matching trivial) is `REJECTED_FLOODING`. This was the fatal exploit in the
   original scorer.
2. **Nearest-band matching** is then safe (n_pred ≈ n_occ = 1–3 bands).
3. **KS-baseline gate** — thresholds (PASS < 0.30 eV) sit below what bare KS
   scores (~0.41–0.61 eV), so a no-physics submission FAILs. `ks_baseline_rmse_eV`
   is reported for audit.

## Maintainer audit checklist

- `run_qp.py` must not read any `arpes_reference.csv`.
- Same code path for public and hidden elements; no per-element branches/hardcodes.
- `method.md` (if provided) describes a parameter-free correction, not a fit.
- Low RMSE + a hidden-element hardcode = invalid.

## Tests

```bash
python -m pytest tests/ -q
```
