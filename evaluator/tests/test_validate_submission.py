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


def write_gold(gold_dir, element, rows):
    # rows: [element, point_id, t, E_pred_eV, E_KS_eV, Delta]
    gold_dir.mkdir(parents=True, exist_ok=True)
    write_csv(gold_dir / f"{element}_gold.csv",
              [["element", "point_id", "t", "E_pred_eV", "E_KS_eV", "Delta"]] + rows)


def run_validator(submission, hidden, gold, result_json=None):
    cmd = [
        sys.executable,
        str(VALIDATOR),
        "--submission-dir",
        str(submission),
        "--hidden-dir",
        str(hidden),
        "--gold-dir",
        str(gold),
    ]
    if result_json is not None:
        cmd += ["--json", str(result_json)]
    return subprocess.run(cmd, text=True, capture_output=True)


def test_validator_scores_nearest_band_predictions(tmp_path):
    hidden = tmp_path / "hidden"
    element_dir = hidden / "X"
    write_csv(element_dir / "grid.csv", [["point_id", "t"], [1, 0.0], [2, 1.0]])
    write_csv(
        element_dir / "arpes_reference.csv",
        [["point_id", "t", "E_expt_eV"], [1, 0.0, -2.0], [2, 1.0, -1.0]],
    )
    (element_dir / "element_config.json").write_text('{"element": "X"}\n')
    gold = tmp_path / "gold"
    write_gold(gold, "X", [["X", 1, 0.0, -2.5, -2.7, 0.08], ["X", 1, 0.0, -2.02, -2.2, 0.09],
                           ["X", 2, 1.0, -1.01, -1.1, 0.09]])

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
    proc = run_validator(submission, hidden, gold, result_json)

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


def test_validator_rejects_flooding_submission(tmp_path):
    # An agent that floods each k-point with a dense energy grid (to make
    # nearest-band matching trivial) is rejected, not scored.
    hidden = tmp_path / "hidden"
    element_dir = hidden / "X"
    write_csv(element_dir / "grid.csv", [["point_id", "t"], [1, 0.0], [2, 0.5]])
    write_csv(element_dir / "arpes_reference.csv",
              [["point_id", "t", "E_expt_eV"], [1, 0.0, -2.0], [2, 0.5, -1.0]])
    (element_dir / "element_config.json").write_text('{}\n')
    gold = tmp_path / "gold"
    write_gold(gold, "X", [["X", 1, 0.0, -2.0, -2.3, 0.1], ["X", 2, 0.5, -1.0, -1.2, 0.1]])

    submission = tmp_path / "submission"
    submission.mkdir()
    (submission / "run_qp.py").write_text(
        "\n".join([
            "import argparse, csv",
            "p = argparse.ArgumentParser()",
            "p.add_argument('--element-config'); p.add_argument('--grid'); p.add_argument('--out')",
            "a = p.parse_args()",
            "w = csv.writer(open(a.out, 'w', newline=''))",
            "w.writerow(['element','point_id','t','E_pred_eV'])",
            "[w.writerow(['X', pid, 0.0, round(-5.0+0.05*i, 3)]) for pid in (1,2) for i in range(120)]",
        ]) + "\n"
    )

    proc = run_validator(submission, hidden, gold)
    assert proc.returncode != 0, proc.stdout
    result = json.loads(proc.stdout)
    assert result["per_element"]["X"]["verdict"] == "REJECTED_FLOODING", result


def test_validator_does_not_expose_sibling_arpes_to_runner(tmp_path):
    hidden = tmp_path / "hidden"
    element_dir = hidden / "X"
    write_csv(element_dir / "grid.csv", [["point_id", "t"], [1, 0.0]])
    write_csv(element_dir / "arpes_reference.csv", [["point_id", "t", "E_expt_eV"], [1, 0.0, -2.0]])
    (element_dir / "element_config.json").write_text('{"element": "X"}\n')
    gold = tmp_path / "gold"
    write_gold(gold, "X", [["X", 1, 0.0, -2.0, -2.2, 0.1]])

    submission = tmp_path / "submission"
    submission.mkdir()
    (submission / "run_qp.py").write_text(
        "\n".join([
            "import argparse, csv",
            "from pathlib import Path",
            "p = argparse.ArgumentParser()",
            "p.add_argument('--element-config'); p.add_argument('--grid'); p.add_argument('--out')",
            "a = p.parse_args()",
            "arpes = Path(a.grid).parent / 'arpes_reference.csv'",
            "with open(arpes, newline='') as inp, open(a.out, 'w', newline='') as out:",
            "    r = csv.DictReader(inp)",
            "    w = csv.writer(out)",
            "    w.writerow(['element','point_id','t','E_pred_eV'])",
            "    [w.writerow(['X', row['point_id'], row.get('t', ''), row['E_expt_eV']]) for row in r]",
        ]) + "\n"
    )

    proc = run_validator(submission, hidden, gold)

    assert proc.returncode != 0
    assert "runner failed" in proc.stderr


def test_validator_does_not_expose_element_name_in_runner_paths(tmp_path):
    hidden = tmp_path / "hidden"
    element_dir = hidden / "X"
    write_csv(element_dir / "grid.csv", [["point_id", "t"], [1, 0.0]])
    write_csv(element_dir / "arpes_reference.csv", [["point_id", "t", "E_expt_eV"], [1, 0.0, -2.0]])
    (element_dir / "element_config.json").write_text('{"element": "X"}\n')
    gold = tmp_path / "gold"
    write_gold(gold, "X", [["X", 1, 0.0, -2.0, -2.2, 0.1]])

    submission = tmp_path / "submission"
    submission.mkdir()
    (submission / "run_qp.py").write_text(
        "\n".join([
            "import argparse, csv, sys",
            "p = argparse.ArgumentParser()",
            "p.add_argument('--element-config'); p.add_argument('--grid'); p.add_argument('--out')",
            "a = p.parse_args()",
            "if 'X' in a.element_config or 'X' in a.grid or 'X' in a.out:",
            "    sys.exit(3)",
            "w = csv.writer(open(a.out, 'w', newline=''))",
            "w.writerow(['element','point_id','t','E_pred_eV'])",
            "w.writerow(['case', 1, 0.0, -2.0])",
        ]) + "\n"
    )

    proc = run_validator(submission, hidden, gold)

    assert proc.returncode == 0, proc.stderr


def test_validator_rejects_missing_reference_points(tmp_path):
    hidden = tmp_path / "hidden"
    element_dir = hidden / "X"
    write_csv(element_dir / "grid.csv", [["point_id", "t"], [1, 0.0], [2, 0.5]])
    write_csv(element_dir / "arpes_reference.csv",
              [["point_id", "t", "E_expt_eV"], [1, 0.0, -2.0], [2, 0.5, -1.0]])
    (element_dir / "element_config.json").write_text('{"element": "X"}\n')
    gold = tmp_path / "gold"
    write_gold(gold, "X", [["X", 1, 0.0, -2.0, -2.2, 0.1], ["X", 2, 0.5, -1.0, -1.2, 0.1]])

    submission = tmp_path / "submission"
    submission.mkdir()
    (submission / "run_qp.py").write_text(
        "\n".join([
            "import argparse, csv",
            "p = argparse.ArgumentParser()",
            "p.add_argument('--element-config'); p.add_argument('--grid'); p.add_argument('--out')",
            "a = p.parse_args()",
            "w = csv.writer(open(a.out, 'w', newline=''))",
            "w.writerow(['element','point_id','t','E_pred_eV'])",
            "w.writerow(['X', 1, 0.0, -2.0])",
        ]) + "\n"
    )

    proc = run_validator(submission, hidden, gold)

    assert proc.returncode != 0, proc.stdout
    result = json.loads(proc.stdout)
    assert result["per_element"]["X"]["verdict"] == "INVALID_SHAPE", result
    assert result["per_element"]["X"]["n_missing"] == 1


def test_validator_rejects_missing_occupied_bands(tmp_path):
    hidden = tmp_path / "hidden"
    element_dir = hidden / "X"
    write_csv(element_dir / "grid.csv", [["point_id", "t"], [1, 0.0]])
    write_csv(element_dir / "arpes_reference.csv", [["point_id", "t", "E_expt_eV"], [1, 0.0, -2.0]])
    (element_dir / "element_config.json").write_text('{"element": "X"}\n')
    gold = tmp_path / "gold"
    write_gold(gold, "X", [["X", 1, 0.0, -2.5, -2.7, 0.08], ["X", 1, 0.0, -2.0, -2.2, 0.09]])

    submission = tmp_path / "submission"
    submission.mkdir()
    (submission / "run_qp.py").write_text(
        "\n".join([
            "import argparse, csv",
            "p = argparse.ArgumentParser()",
            "p.add_argument('--element-config'); p.add_argument('--grid'); p.add_argument('--out')",
            "a = p.parse_args()",
            "w = csv.writer(open(a.out, 'w', newline=''))",
            "w.writerow(['element','point_id','t','E_pred_eV'])",
            "w.writerow(['X', 1, 0.0, -2.0])",
        ]) + "\n"
    )

    proc = run_validator(submission, hidden, gold)

    assert proc.returncode != 0, proc.stdout
    result = json.loads(proc.stdout)
    assert result["per_element"]["X"]["verdict"] == "INVALID_SHAPE", result
    assert result["per_element"]["X"]["band_count_mismatches"] == {"1": {"expected": 2, "got": 1}}
