# QP-EFT Agent Benchmark

A benchmark for **agentic theoretical-physics reasoning + scientific coding**,
built on *Kohn–Sham Hamiltonian from Effective Field Theory: Quasiparticle Band
Narrowing from Frozen Core Dynamics* (arXiv:2604.25199).

The task: predict the **occupied quasiparticle band energies** of simple metals.
A bare Kohn–Sham band (from DFTK) overestimates the ARPES bandwidth by 20–35%
for alkali metals; the agent must reconstruct and implement the parameter-free
frozen-core correction `E_QP - E_F = z_core(n,k) (E_KS - E_F)` and beat bare KS
against held-out ARPES.

Theory is **not** graded directly. It is tested *implicitly*: only code that
encodes the correct physics produces correct numbers — so we withhold the
formula (at the harder levels) and grade the predictions.

## Difficulty ladder (same physics, withhold more)

| Level | Given | Withheld | Tests |
|-------|-------|----------|-------|
| **L1** | full `z_core` formula + precomputed `f_c(K)` tables | — | DFTK wiring, applying the correction, generalization |
| **L2** | the formula; atomic core data (`u_c`, `V_H_c`) | the form factors | implementing the form-factor quadrature |
| **L3** | only the physical setup + atomic core data | the formula | deriving `z_core` from the EFT |

All levels are scored identically against held-out ARPES.

## What the agent gets vs what is hidden

- **Agent-facing:** `agent_packet/` — README + `levels/L{1,2,3}/` with the public
  development elements **Na, Al** (config, grid, ARPES, level data) and the task
  + theory/setup docs. Hand the agent ONE level's directory.
- **Hidden (maintainer-only):** `evaluator/` — held-out metals **K, Mg** per
  level, the gold band references, and the scorer. The agent develops against
  Na/Al only; its generic `run_qp.py` is run on the concealed metals.

## Submission & scoring

The agent submits `run_qp.py`:

```bash
python run_qp.py --element-config element_config.json --grid grid.csv --out qp_bands.csv
```

```bash
python evaluator/validate_submission.py --submission-dir <dir> --level 2 --json result.json
```

Scoring is nearest-band RMSE vs ARPES, with predictions **capped to the gold
occupied-band count** per k-point (flooding is rejected). **PASS < 0.30 eV**,
PARTIAL 0.30–0.40, FAIL otherwise. Bare KS scores ~0.4–0.6 eV → FAIL, so the
correction is necessary to pass (the report states the KS baseline for audit).

## Validation (gold reference vs real ARPES)

| El | bare KS RMSE | QP RMSE | verdict |
|----|--------------|---------|---------|
| Na | 0.413 | 0.078 | PASS |
| K  | 0.614 | 0.139 | PASS |
| Mg | 0.434 | 0.187 | PASS |
| Al | 0.414 | 0.248 | PASS |

End-to-end (reference submission → evaluator, hidden K+Mg, level 2): overall
PASS, mean RMSE 0.163 eV. The reference solution (`reference/`) reproduces the
paper's Table I exactly and is the gold used to calibrate thresholds and band
counts. See `DESIGN.md` for the full design and `reference/README.md` for
provenance.

## Layout

```text
agent_packet/levels/L{1,2,3}/{Na,Al}/   agent-facing packets (one level per run)
evaluator/
  hidden/L{1,2,3}/{K,Mg}/               held-out metals (config+grid+data+ARPES)
  gold/                                 gold occupied-band references
  source_data/                          raw grids+ARPES (maintainer source)
  validate_submission.py                scorer (--level)
reference/                              gold solution, atomic solver, generators
environment/                            pinned Julia project (DFTK 0.7.25)
```

## Reproducing the benchmark data

```bash
julia --project=environment reference/gen_packet_data.jl   # atomic data + f_c tables
python3 reference/build_packet.py                          # assemble levels + hidden
```
