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
4. **Tight per-element thresholds (maintainer-only)** — `PASS_THRESHOLD_EV = {"K":
   0.17, "Mg": 0.21}`, calibrated to the gold (0.139 / 0.187) + a small margin;
   overall PASS requires *every* hidden element to clear its own bar. These sit well
   below the bare-KS baselines (~0.41–0.61 eV, so a no-physics submission FAILs;
   `ks_baseline_rmse_eV` reported) and just above the gold, so a cruder approximation
   that only roughly matches the magnitude does not pass. The agent never sees these
   numbers — the task states a physical ~0.1 eV target instead.

(The Harbor tasks add container isolation: agent code runs as `nobody` with
`/tests/{hidden,gold}` root-only — see `harbor/`.)

## Maintainer audit checklist

- Same code path for public and hidden elements; no per-element branches/hardcodes
  (reading `arpes_reference.csv` is now structurally impossible — it's not in the
  runner's inputs — but still confirm no per-element answer tables).
- `method.md` (required by the submission spec) describes a parameter-free
  correction, not a fit; for Level 3 it should show the derivation.
- Low RMSE + a hidden-element hardcode = invalid.
- **Vertex diagonal check (Level 3).** From the submission's derived per-channel
  coupling potential `V_vertex(r)`, compute `∫ u_c(r)² V_vertex(r) dr`. The
  consistent value is **0**: the dynamical vertex must be a zero-mean fluctuation
  against the core orbital — its static mean is already counted in `DeltaE_c`
  and in the PSP, so a large nonzero diagonal means the coupling double-counts
  the static sector (and was typically locked in by numerical coincidence on the
  public elements rather than derived). One-line classifier over observed
  failure modes: gold `(V_H_c − J_c)` → 0; `J_c`-as-coupling → `J_c ≠ 0`;
  bare-nuclear `Z/r` → `Z⟨1/r⟩_c`, huge. See
  `experiments/L3_fable5/REPORT.md` and
  `maintainer_sources/NOTE-vertex-derivation-Li.md` (§7).

## Tests

```bash
python -m pytest tests/ -q
```
