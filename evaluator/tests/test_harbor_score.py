import csv
import json
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
TASK = ROOT / "harbor" / "qp-eft-L2-formfactor"
SCORER = TASK / "tests" / "score.py"
HIDDEN = TASK / "tests" / "hidden"
GOLD = TASK / "tests" / "gold"


def write_pred(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["element", "point_id", "t", "E_pred_eV"])
        w.writerows(rows)


def run_score(pred_dir):
    proc = subprocess.run(
        [
            sys.executable,
            str(SCORER),
            "--pred-dir",
            str(pred_dir),
            "--hidden",
            str(HIDDEN),
            "--gold",
            str(GOLD),
        ],
        text=True,
        capture_output=True,
    )
    return proc, json.loads(proc.stdout)


def test_harbor_score_rejects_missing_reference_points(tmp_path):
    pred = tmp_path / "pred"
    for el in ("K", "Mg"):
        with (HIDDEN / el / "arpes_reference.csv").open(newline="") as f:
            row = next(csv.DictReader(f))
        write_pred(pred / f"{el}_out.csv", [[el, row["point_id"], row["t"], row["E_expt_eV"]]])

    proc, result = run_score(pred)

    assert proc.returncode != 0, proc.stdout
    assert result["overall"]["verdict"] == "INVALID_SHAPE"
    assert result["per_element"]["K"]["n_missing"] == 9
    assert result["per_element"]["Mg"]["n_missing"] == 87


def test_harbor_score_rejects_missing_occupied_bands(tmp_path):
    pred = tmp_path / "pred"
    for el in ("K", "Mg"):
        ref = {
            int(r["point_id"]): float(r["E_expt_eV"])
            for r in csv.DictReader((HIDDEN / el / "arpes_reference.csv").open(newline=""))
        }
        gold_bands = {}
        with (GOLD / f"{el}_gold.csv").open(newline="") as f:
            for r in csv.DictReader(f):
                gold_bands.setdefault(int(r["point_id"]), []).append(float(r["E_pred_eV"]))
        rows = []
        for pid, bands in gold_bands.items():
            if pid in ref:
                rows.append([el, pid, "", min(bands, key=lambda v: abs(v - ref[pid]))])
        write_pred(pred / f"{el}_out.csv", rows)

    proc, result = run_score(pred)

    assert proc.returncode != 0, proc.stdout
    assert result["overall"]["verdict"] == "INVALID_SHAPE"
    assert result["per_element"]["Mg"]["verdict"] == "INVALID_SHAPE"
    assert result["per_element"]["Mg"]["band_count_mismatches"]
