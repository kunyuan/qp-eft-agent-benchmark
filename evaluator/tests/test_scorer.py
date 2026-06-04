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
        r = V.score_element(ref, gold_qp, gold_qp, gold_ks)
        assert r["verdict"] == "PASS", (el, r)
        assert r["rmse_eV"] < V.PASS_RMSE_EV


def test_bare_ks_fails():
    # submitting the uncorrected KS bands must NOT pass (the whole point)
    for el in ("K", "Mg"):
        ref, gold_qp, gold_ks = _load(el)
        r = V.score_element(ref, gold_ks, gold_qp, gold_ks)
        assert r["verdict"] in ("FAIL", "PARTIAL"), (el, r)
        assert r["rmse_eV"] >= V.PASS_RMSE_EV


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


if __name__ == "__main__":
    test_gold_qp_passes_for_each_hidden_element()
    test_bare_ks_fails()
    test_flooding_is_rejected()
    test_flooding_would_have_cheated_without_guard()
    print("all scorer behavior tests passed")
