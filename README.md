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
| **L4** | only the problem + public ARPES + Li as a closed-form anchor element (mandatory Li-first analytic stage) | the formula AND all atomic data; no prescribed approximations | open frontier: own atomic solver, complete contraction enumeration on Li, declared/controlled approximations, match or **beat** the published leading order |

All levels are scored identically against held-out ARPES.

## What the agent gets vs what is hidden

- **Agent-facing:** `agent_packet/` — README + `levels/L{1,2,3,4}/` with the public
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

Scoring is nearest-band RMSE vs held-out ARPES, with **tight per-element thresholds**
calibrated to the gold (hidden `K < 0.17`, `Mg < 0.21` eV; gold scores 0.139 / 0.187
and faithful implementations reproduce it to ~0.001 eV). Overall PASS requires
*every* hidden element to clear its own bar. Bare KS scores ~0.4–0.6 eV → FAIL, so
the correction is necessary (the report states the KS baseline for audit). The
agent is **not** told the grading bar — the instruction states a *physical* target
(reach ~0.1 eV agreement with experiment, the level of the many-body reference
eDMFT) with public Li/Na/Ca LDA/eDMFT/experiment anchors to calibrate against.

Anti-cheat guards (so a low RMSE means "did the physics"):
- The runner is handed a **sanitized input directory with no `arpes_reference.csv`**
  — it cannot read the answers sitting next to the config. In the Harbor tasks it
  also runs as an unprivileged user with the hidden/gold data root-only.
- Predictions must cover **every** reference point with **exactly** the gold band
  count (occupied bands + the first unoccupied band) per k-point — flooding,
  sparse, and under-band submissions are rejected.
- Develops against Na/Al only; the held-out metals are concealed, so memorized
  per-element answers can't be hardcoded into generic code.

## Run as a Harbor task (containerized agent harness)

Each level is also packaged as a [Harbor](https://github.com/harbor-framework/harbor)
task under `harbor/qp-eft-L{1,2,3,4}-*/` — a pinned container (Julia + DFTK), the
oracle, and the verifier (see `harbor/README.md`):

```bash
harbor run --agent oracle --path harbor/qp-eft-L1-apply          # check the task is solvable
harbor run --agent claude-code --path harbor/qp-eft-L2-formfactor # run a real agent
```

Fresh solver agents (given only one level's packet) reproduce the gold and pass
the concealed K/Mg: L1 (verified across 3 independent agents) and L2 (which must
implement the form-factor quadrature itself) both reach overall PASS ≈ 0.163 eV.

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
agent_packet/levels/L{1,2,3,4}/{Na,Al}/ agent-facing packets (one level per run)
evaluator/
  hidden/L{1,2,3,4}/{K,Mg}/             held-out metals (config+grid+data+ARPES)
  gold/                                 gold band references (occupied + first unoccupied)
  source_data/                          raw grids+ARPES (maintainer source)
  validate_submission.py                scorer (--level; sanitizes inputs, anti-cheat guards)
harbor/qp-eft-L{1,2,3,4}-*/             the four levels as Harbor tasks (env/oracle/verifier)
reference/                              gold solution, atomic solver, generators
environment/                            pinned Julia project (DFTK 0.7.25)
```

## Reproducing the benchmark data

```bash
julia --project=environment reference/gen_packet_data.jl   # atomic data + f_c tables
python3 reference/build_packet.py                          # assemble agent_packet + hidden
python3 reference/build_harbor.py                          # assemble the four Harbor tasks
```
