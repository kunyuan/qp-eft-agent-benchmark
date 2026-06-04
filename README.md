# EFT-QP Agent Benchmark Harness

This repository is an evaluator harness for testing whether a problem-solving
agent can reconstruct and implement the frozen-core quasiparticle correction
for simple-metal band narrowing.

Do not give this full repository to a problem-solving agent. Give only
`agent_packet/`. The `evaluator/` directory contains hidden validation data.

## Benchmark Idea

The public development set is Na and Al. The agent sees their grids,
experimental ARPES references, and material configs. It must use these to:

1. reproduce the Kohn-Sham baseline with DFTK/LDA/GTH;
2. diagnose why the Na Kohn-Sham bandwidth is too large while Al is nearly OK;
3. derive a parameter-free frozen-core quasiparticle correction;
4. implement a generic runner:

```bash
python run_qp.py --element-config element_config.json --grid grid.csv --out qp_bands.csv
```

The hidden validation set is K and Mg. The evaluator calls the same runner on
hidden configs and grids, then compares the output to held-out ARPES data.

## Directory Layout

```text
agent_packet/          agent-facing task packet; safe to hand to solvers
evaluator/             hidden validation data and scorer
environment/           Julia environment declaration for DFTK work
maintainer_sources/    optional paper/source materials for benchmark maintainers
```

## Running The Hidden Validator

After an agent submits a directory containing `run_qp.py`, run:

```bash
python evaluator/validate_submission.py \
  --submission-dir /path/to/submission \
  --json result.json
```

The validator passes hidden `element_config.json` and `grid.csv` files to the
runner. The runner must write a CSV with columns:

```text
element,point_id,t,E_pred_eV
```

The scorer matches each held-out ARPES point to the nearest submitted occupied
band at the same `point_id`.

## Scoring

Per hidden element:

| RMSE (eV) | Verdict |
| --- | --- |
| `< 0.20` | PASS |
| `0.20 - 0.40` | PARTIAL |
| `>= 0.40` | FAIL |

The overall verdict uses the mean hidden-set RMSE.

## Leakage Rules

- Agent-facing solvers must receive only `agent_packet/`.
- Hidden `K/` and `Mg/` folders must not be included in prompts, context, or
  mounted workspaces for the problem-solving agent.
- The full paper and old benchmark README should not be shown to solvers,
  because they reveal validation trends and intended-method diagnostics.
- During validation, the agent code may receive hidden element configs and
  grids as ordinary runtime inputs. It must not receive hidden ARPES references.

