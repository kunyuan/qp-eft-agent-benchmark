# Task 05: Final Submission

Submit a directory containing:

```text
run_qp.py
method.md
<supporting source files>
```

Before submission, test locally on the public elements:

```bash
python run_qp.py \
  --element-config agent_packet/data/public/Na/element_config.json \
  --grid agent_packet/data/public/Na/grid.csv \
  --out Na_qp_bands.csv

python run_qp.py \
  --element-config agent_packet/data/public/Al/element_config.json \
  --grid agent_packet/data/public/Al/grid.csv \
  --out Al_qp_bands.csv
```

Then compare public outputs against public ARPES data using your own scorer.

Your hidden score is based on additional simple metals. The evaluator will call
the same `run_qp.py` interface with hidden configs and hidden grids.

