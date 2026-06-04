#!/usr/bin/env python3
"""Self-contained verifier scorer (mirrors evaluator/validate_submission.py).

Reads predictions already produced by test.sh, scores them against the hidden
ARPES with the flooding guard + KS-baseline gate, writes result.json, and exits
0 iff overall verdict is PASS.

  score.py --pred-dir <dir> --hidden <dir> --gold <dir> --json <out.json>
"""
from __future__ import annotations
import argparse, csv, json, math
from pathlib import Path

PASS_RMSE_EV = 0.30
PARTIAL_RMSE_EV = 0.40
FLOOD_POINT_FRAC = 0.20
FLOOD_TOTAL_RATIO = 1.5


def read_reference(path: Path):
    with path.open(newline="") as f:
        return {int(r["point_id"]): float(r["E_expt_eV"]) for r in csv.DictReader(f)}


def read_predictions(path: Path):
    out = {}
    if not path.exists():
        return out
    with path.open(newline="") as f:
        rows = csv.DictReader(f)
        if rows.fieldnames is None or "point_id" not in rows.fieldnames \
                or "E_pred_eV" not in rows.fieldnames:
            return out
        for r in rows:
            v = (r.get("E_pred_eV") or "").strip()
            if v:
                out.setdefault(int(r["point_id"]), []).append(float(v))
    return out


def read_gold(path: Path):
    qp, ks = {}, {}
    with path.open(newline="") as f:
        for r in csv.DictReader(f):
            pid = int(r["point_id"])
            qp.setdefault(pid, []).append(float(r["E_pred_eV"]))
            ks.setdefault(pid, []).append(float(r["E_KS_eV"]))
    return qp, ks


def verdict(rmse):
    if rmse is None:
        return "NO_PREDICTION"
    if rmse < PASS_RMSE_EV:
        return "PASS"
    if rmse < PARTIAL_RMSE_EV:
        return "PARTIAL"
    return "FAIL"


def _rmse(a):
    return math.sqrt(sum(x * x for x in a) / len(a))


def _nearest_rmse(ref, bands):
    res = [min(bands[p], key=lambda v: abs(v - ref[p])) - ref[p]
           for p in ref if p in bands and bands[p]]
    return _rmse(res) if res else None


def score_element(ref, pred, gold_qp, gold_ks):
    n_occ = {p: len(b) for p, b in gold_qp.items()}
    over = [p for p, b in pred.items() if len(b) > n_occ.get(p, 1) + 1]
    total_pred = sum(len(b) for b in pred.values())
    total_gold = sum(n_occ.values()) or 1
    if len(over) > FLOOD_POINT_FRAC * max(len(pred), 1) or total_pred > FLOOD_TOTAL_RATIO * total_gold:
        return {"verdict": "REJECTED_FLOODING", "rmse_eV": None,
                "n_points_over_band_cap": len(over), "total_predicted_bands": total_pred,
                "total_gold_bands": total_gold}
    scored = sorted(set(ref) & set(pred))
    missing = sorted(set(ref) - set(pred))
    unexpected = sorted(set(pred) - set(ref))
    mismatches = {
        str(p): {"expected": n_occ.get(p, 0), "got": len(pred.get(p, []))}
        for p in sorted(ref)
        if p in pred and len(pred.get(p, [])) != n_occ.get(p, 0)
    }
    if missing or unexpected or mismatches:
        return {"verdict": "INVALID_SHAPE", "rmse_eV": None,
                "n_scored": len(scored), "n_missing": len(missing),
                "n_unexpected": len(unexpected), "n_points_over_band_cap": len(over),
                "band_count_mismatches": mismatches,
                "note": "predictions must cover every reference point with exactly the gold occupied-band count"}
    if not scored:
        return {"verdict": "NO_PREDICTION", "rmse_eV": None, "n_scored": 0}
    res = [min(pred[p], key=lambda v: abs(v - ref[p])) - ref[p] for p in scored]
    rmse = _rmse(res)
    return {"verdict": verdict(rmse), "rmse_eV": round(rmse, 6),
            "mean_signed_eV": round(sum(res) / len(res), 6),
            "n_scored": len(scored), "n_missing": len(set(ref) - set(pred)),
            "n_points_over_band_cap": len(over),
            "ks_baseline_rmse_eV": round(r, 6) if (r := _nearest_rmse(ref, gold_ks)) else None}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pred-dir", required=True, type=Path)
    ap.add_argument("--hidden", required=True, type=Path)
    ap.add_argument("--gold", required=True, type=Path)
    ap.add_argument("--json", type=Path)
    a = ap.parse_args()

    per = {}
    for eldir in sorted(p for p in a.hidden.iterdir() if p.is_dir()):
        el = eldir.name
        ref = read_reference(eldir / "arpes_reference.csv")
        pred = read_predictions(a.pred_dir / f"{el}_out.csv")
        gqp, gks = read_gold(a.gold / f"{el}_gold.csv")
        per[el] = score_element(ref, pred, gqp, gks)

    rmses = [r["rmse_eV"] for r in per.values() if r.get("rmse_eV") is not None]
    mean = sum(rmses) / len(rmses) if rmses else None
    rejected = any(r.get("verdict") == "REJECTED_FLOODING" for r in per.values())
    invalid = any(r.get("verdict") == "INVALID_SHAPE" for r in per.values())
    empty = any(r.get("verdict") == "NO_PREDICTION" for r in per.values())
    if rejected:
        overall = "REJECTED_FLOODING"
    elif invalid:
        overall = "INVALID_SHAPE"
    elif empty:
        overall = "NO_PREDICTION"
    else:
        overall = verdict(mean)
    result = {"per_element": per,
              "overall": {"mean_rmse_eV": round(mean, 6) if mean is not None else None,
                          "verdict": overall,
                          "thresholds_eV": {"pass_below": PASS_RMSE_EV, "partial_below": PARTIAL_RMSE_EV}}}
    text = json.dumps(result, indent=2, sort_keys=True)
    print(text)
    if a.json:
        a.json.write_text(text + "\n")
    return 0 if overall == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
