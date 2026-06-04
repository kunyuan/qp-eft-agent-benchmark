# Evaluator (maintainer-only — not agent-facing)

Scores a submitted `run_qp.py` against held-out ARPES.

```
hidden/L{1,2,3}/{K,Mg}/   held-out metals per level: element_config.json, grid.csv,
                          arpes_reference.csv (scoring target), + level data files
gold/<El>_gold.csv        gold bands (QP + KS, occupied + first unoccupied) from
                          reference/gold_runner.jl; used for the band-count check
                          and the KS-baseline audit
source_data/<El>/         raw grids + ARPES (the build source for hidden/ and packet)
validate_submission.py    the scorer
```

## Run

```bash
python validate_submission.py --submission-dir <dir> --level 2 --json result.json
```

`--level {1,2,3}` selects `hidden/L<level>`. For each hidden element the validator
copies the inputs to a **sanitized temp dir with no `arpes_reference.csv`**
(`copy_runner_inputs`), runs `python <dir>/run_qp.py --element-config ... --grid
... --out ...` against it, then scores.

## Scoring guards (why low RMSE means "did the physics")

1. **No answer key in the inputs** — the runner gets a sanitized input dir without
   `arpes_reference.csv`, so it cannot read the held-out energies sitting next to the
   config and echo them. (A real exploit: the agent could discover this on Na/Al.)
2. **Shape guard** — every reference point must be predicted with **exactly** the
   gold band count per k-point (occupied + first unoccupied). Flooding →
   `REJECTED_FLOODING`; missing points / wrong count → `INVALID_SHAPE`. Closes the
   flooding, sparse, and under-band cheats.
3. **Nearest-band RMSE** is then safe (n_pred = 1–4 well-separated bands).
4. **KS-baseline gate** — thresholds (PASS < 0.30 eV) sit below what bare KS scores
   (~0.41–0.61 eV), so a no-physics submission FAILs. `ks_baseline_rmse_eV` reported.

(The Harbor tasks add container isolation: agent code runs as `nobody` with
`/tests/{hidden,gold}` root-only, **and the verifier phase runs with no network**
(`[verifier] network_mode = "no-network"`) so the runner cannot fetch the answers
online — see `harbor/`. This local `validate_submission.py` does **not** isolate
network by itself; run it on trusted submissions, or inside a no-network sandbox.
The real agent-eval path is the Harbor verifier, which enforces offline.)

## Maintainer audit checklist

- Same code path for public and hidden elements; no per-element branches/hardcodes
  (reading `arpes_reference.csv` is now structurally impossible — it's not in the
  runner's inputs — but still confirm no per-element answer tables).
- `method.md` (required by the submission spec) describes a parameter-free
  correction, not a fit; for Level 3 it should show the derivation.
- Low RMSE + a hidden-element hardcode = invalid.

## Tests

```bash
python -m pytest tests/ -q
```
