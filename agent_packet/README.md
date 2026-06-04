# Agent Task Packet: Frozen-Core Quasiparticle Bands

You are given an offline scientific benchmark. Your goal is to build a generic
pipeline that predicts occupied quasiparticle band energies for simple metals.

You may use only the information in this packet and the software environment
provided by the benchmark runner. Do not use external experimental band data.

## Public Development Systems

You are given public data for:

- `data/public/Na/`
- `data/public/Al/`

Each element directory contains:

- `grid.csv`: k-points as `point_id,t`, where `t` is fractional position along
  the named high-symmetry path.
- `arpes_reference.csv`: public ARPES quasiparticle targets for development.
- `element_config.json`: crystal, DFTK, and frozen-core model inputs.

The hidden evaluator will later test the same code on additional simple metals.
Do not special-case Na or Al.

## Required Submission

Submit a directory containing:

- `run_qp.py`: command-line runner with this interface:

```bash
python run_qp.py --element-config element_config.json --grid grid.csv --out qp_bands.csv
```

- `method.md`: short derivation and implementation notes.
- Any source files needed by `run_qp.py`.

The output CSV must contain:

```text
element,point_id,t,E_pred_eV
```

Use one row per occupied band at each k-point. Repeat `point_id` when multiple
occupied bands are present.
Do not emit dense candidate-energy grids or other extra bands: the evaluator
checks that the number of rows per k-point is physically consistent with the
element.

## Suggested Workflow

1. Read `sources/theory_background.md`.
2. Run the DFTK baseline described in `tasks/01_dftk_baseline.md`.
3. Compare Na/Al Kohn-Sham predictions with public ARPES data using
   `tasks/02_na_al_failure_diagnosis.md`.
4. Derive the frozen-core quasiparticle correction using
   `tasks/03_theory_derivation.md` and `sources/formula_sheet.md`.
5. Implement the generic runner described in `tasks/04_qp_correction_code.md`.
6. Package the final submission as described in `tasks/05_submission.md`.

## Non-Negotiable Rules

- No empirical fitting to ARPES.
- No per-element tuning.
- No hardcoded public ARPES values in `run_qp.py`.
- No hidden-data assumptions.
- A converged LDA/PBE Kohn-Sham band alone is not a quasiparticle prediction
  for this benchmark.
