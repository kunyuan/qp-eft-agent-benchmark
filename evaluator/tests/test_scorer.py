"""Behavioral tests for the rewritten scorer: gold passes, KS-only fails,
flooding is rejected. Uses the real K/Mg gold + ARPES data."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "evaluator"))
import validate_submission as V  # noqa: E402

GOLD = ROOT / "evaluator" / "gold"
HIDDEN = ROOT / "evaluator" / "hidden" / "L2"  # ARPES is level-independent


def _load(el):
    ref = V.read_reference(HIDDEN / el / "arpes_reference.csv")
    gold_qp, gold_ks = V.read_gold(GOLD / f"{el}_gold.csv")
    return ref, gold_qp, gold_ks


def test_gold_qp_passes_for_each_hidden_element():
    for el in ("K", "Mg"):
        ref, gold_qp, gold_ks = _load(el)
        r = V.score_element(ref, gold_qp, gold_qp, gold_ks, element=el)
        assert r["verdict"] == "PASS", (el, r)
        assert r["rmse_eV"] < V.pass_threshold(el)


def test_bare_ks_fails():
    # submitting the uncorrected KS bands must NOT pass (the whole point)
    for el in ("K", "Mg"):
        ref, gold_qp, gold_ks = _load(el)
        r = V.score_element(ref, gold_ks, gold_qp, gold_ks, element=el)
        assert r["verdict"] in ("FAIL", "PARTIAL"), (el, r)
        assert r["rmse_eV"] >= V.pass_threshold(el)


def test_tight_thresholds_separate_gold_from_crude_model():
    # the per-element bars (K 0.17, Mg 0.21) separate faithful frozen-core impls
    # (gold ~0.139/0.187 -> PASS) from a cruder state-independent single-Z model
    # that mis-fits the k-dependence (~0.188/0.220 -> NOT pass, lands in PARTIAL).
    assert V.verdict(0.139, "K") == "PASS"
    assert V.verdict(0.187, "Mg") == "PASS"
    assert V.verdict(0.188, "K") != "PASS"    # just over the 0.17 bar -> PARTIAL
    assert V.verdict(0.220, "Mg") != "PASS"   # just over the 0.21 bar -> PARTIAL
    # bare KS is clearly wrong -> FAIL
    assert V.verdict(0.614, "K") == "FAIL" and V.verdict(0.434, "Mg") == "FAIL"
    # an element without a calibrated bar falls back to the generic 0.30
    assert V.verdict(0.25, "Xx") == "PASS" and V.verdict(0.41, "Xx") == "FAIL"


def test_flooding_is_rejected():
    # a submission that floods each point with a dense energy grid is rejected,
    # even though nearest-band matching would otherwise give ~0 RMSE.
    el = "K"
    ref, gold_qp, gold_ks = _load(el)
    flood = {pid: [round(-5.0 + 0.05 * i, 3) for i in range(120)] for pid in gold_qp}
    r = V.score_element(ref, flood, gold_qp, gold_ks)
    assert r["verdict"] == "REJECTED_FLOODING", r


def test_flooding_would_have_cheated_without_guard():
    # sanity: confirm the flood really would score ~0 under naive nearest-match
    el = "K"
    ref, gold_qp, gold_ks = _load(el)
    flood = {pid: [round(-5.0 + 0.05 * i, 3) for i in range(120)] for pid in gold_qp}
    naive = V._nearest_rmse(ref, flood)
    assert naive < 0.05, naive  # exploit is real -> guard is necessary


def test_sparse_predictions_rejected():
    # dropping reference points (predict only the easy/deep ones) -> INVALID_SHAPE
    el = "K"
    ref, gold_qp, gold_ks = _load(el)
    keep = sorted(gold_qp)[:3]
    sparse = {pid: gold_qp[pid] for pid in keep}
    r = V.score_element(ref, sparse, gold_qp, gold_ks)
    assert r["verdict"] == "INVALID_SHAPE" and r["n_missing"] > 0, r


def test_wrong_band_count_rejected():
    # emitting fewer bands than the gold (drop the first-unoccupied) -> INVALID_SHAPE
    el = "Mg"  # multi-band: gold has 2-4 bands per point
    ref, gold_qp, gold_ks = _load(el)
    under = {pid: bands[:1] for pid, bands in gold_qp.items()}  # one band per point
    r = V.score_element(ref, under, gold_qp, gold_ks)
    assert r["verdict"] == "INVALID_SHAPE" and r["n_band_count_mismatches"] > 0, r


if __name__ == "__main__":
    test_gold_qp_passes_for_each_hidden_element()
    test_bare_ks_fails()
    test_flooding_is_rejected()
    test_flooding_would_have_cheated_without_guard()
    print("all scorer behavior tests passed")
