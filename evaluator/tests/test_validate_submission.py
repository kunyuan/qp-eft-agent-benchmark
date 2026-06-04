import csv
import json
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
VALIDATOR = ROOT / "evaluator" / "validate_submission.py"


def write_csv(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerows(rows)


def test_validator_scores_nearest_band_predictions(tmp_path):
    hidden = tmp_path / "hidden"
    element_dir = hidden / "X"
    write_csv(element_dir / "grid.csv", [["point_id", "t"], [1, 0.0], [2, 1.0]])
    write_csv(
        element_dir / "arpes_reference.csv",
        [["point_id", "t", "E_expt_eV"], [1, 0.0, -2.0], [2, 1.0, -1.0]],
    )
    (element_dir / "element_config.json").write_text('{"element": "X"}\n')

    submission = tmp_path / "submission"
    submission.mkdir()
    runner = submission / "run_qp.py"
    runner.write_text(
        "\n".join(
            [
                "import argparse, csv",
                "p = argparse.ArgumentParser()",
                "p.add_argument('--element-config')",
                "p.add_argument('--grid')",
                "p.add_argument('--out')",
                "args = p.parse_args()",
                "with open(args.out, 'w', newline='') as f:",
                "    w = csv.writer(f)",
                "    w.writerow(['element', 'point_id', 't', 'E_pred_eV'])",
                "    w.writerow(['X', 1, 0.0, -2.5])",
                "    w.writerow(['X', 1, 0.0, -2.02])",
                "    w.writerow(['X', 2, 1.0, -1.01])",
            ]
        )
        + "\n"
    )

    result_json = tmp_path / "result.json"
    proc = subprocess.run(
        [
            sys.executable,
            str(VALIDATOR),
            "--submission-dir",
            str(submission),
            "--hidden-dir",
            str(hidden),
            "--json",
            str(result_json),
        ],
        text=True,
        capture_output=True,
    )

    assert proc.returncode == 0, proc.stderr
    result = json.loads(result_json.read_text())
    assert result["overall"]["verdict"] == "PASS"
    assert result["per_element"]["X"]["n_scored"] == 2
    assert result["per_element"]["X"]["rmse_eV"] < 0.03


def test_validator_rejects_missing_runner(tmp_path):
    hidden = tmp_path / "hidden"
    element_dir = hidden / "X"
    write_csv(element_dir / "grid.csv", [["point_id", "t"], [1, 0.0]])
    write_csv(element_dir / "arpes_reference.csv", [["point_id", "t", "E_expt_eV"], [1, 0.0, -1.0]])
    (element_dir / "element_config.json").write_text('{"element": "X"}\n')

    submission = tmp_path / "submission"
    submission.mkdir()

    proc = subprocess.run(
        [
            sys.executable,
            str(VALIDATOR),
            "--submission-dir",
            str(submission),
            "--hidden-dir",
            str(hidden),
        ],
        text=True,
        capture_output=True,
    )

    assert proc.returncode != 0
    assert "run_qp.py" in proc.stderr
