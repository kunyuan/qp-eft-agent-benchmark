# Evaluator Notes

This directory is not agent-facing.

Hidden validation elements live under `hidden/`. Each element folder contains:

- `grid.csv`: runtime grid passed to the submitted solver;
- `element_config.json`: runtime material and core-model input;
- `arpes_reference.csv`: held-out scoring reference, never passed to solver.

Run a submitted solver with:

```bash
python evaluator/validate_submission.py --submission-dir /path/to/submission --json result.json
```

The validator calls:

```bash
python /path/to/submission/run_qp.py \
  --element-config evaluator/hidden/<El>/element_config.json \
  --grid evaluator/hidden/<El>/grid.csv \
  --out /tmp/<El>_qp_bands.csv
```

Then it scores nearest-band RMSE against `arpes_reference.csv`.

## Maintainer Audit Checklist

- Confirm the submitted `run_qp.py` does not read files from this repository's
  `evaluator/hidden/*/arpes_reference.csv`.
- Confirm the code path for hidden elements is the same as for Na/Al.
- Confirm `method.md` describes a parameter-free correction, not a fitted
  bandwidth scaling.
- Treat a low hidden RMSE plus a hardcoded hidden-element branch as invalid.

