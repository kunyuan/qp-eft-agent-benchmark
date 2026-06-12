# Guided session addendum — the vertex derived exactly, and the convention that hid it

After probe 3's evaluation (K FAIL, calibrated E0 = π), the maintainer ran one
**guided** follow-up session in the same workspace: a maintainer note
([`MAINTAINER_NOTE.md`](MAINTAINER_NOTE.md), archived verbatim) diagnosed the
operator-class error from the agent's own logged symptoms and named the
amplitude class (the spectator-mediated transfer objects G/U the agent had
itself written in its channel-C algebra), with an inline comment at the error
site in `derivation_notes.md`. **This session is guided, not a clean probe**;
the guided-vs-derived split is stated explicitly in the agent's final report
and below.

## What the agent then derived (15 minutes)

The per-shell hybridization vertex with the zero-diagonal subtraction:

```
M_c(p) = G_c(p) − U_c·χ̃_c(p),   G_c(p) = ⟨p|V_H[|χ_c|²]|χ_c⟩,  U_c = ⟨χ_cχ_c|v|χ_cχ_c⟩
```

plus, unguided: a Brillouin/collapse lemma proving no one-body vertex exists,
the reference-mismatch argument selecting the same-shell pair term, Li closed
forms, multi-shell generalization, and every sensitivity variant.

**This is the gold vertex.** Maintainer verification:

- **Li, p → 0 (closed form):** agent M(0) = (33/81 − 5/8)ζ · χ̃(0) =
  **−1.88207** Ha·Bohr^{3/2}; gold f(0) (maintainer derivation note, PR #15) =
  **−1.88207**. Identical to five digits.
- **Na 2s, full K range (its own atomic solver vs the published gold table):**

  | K (Bohr⁻¹) | agent M_2s(K) | gold f_c(K) | ratio |
  |---|---|---|---|
  | 0.02 | −5.068 | 1.4503 | −3.495 |
  | 0.50 | −4.170 | 1.1900 | −3.504 |
  | 1.00 | −2.261 | 0.6398 | −3.535 |
  | 2.00 | +0.353 | −0.1040 | −3.389 |

  Constant ratio ≈ **−√4π = −3.5449** (±1.5% over the range where f is
  significant): the same function, an overall sign (irrelevant in |·|²) and a
  spherical-harmonic normalization convention apart.

## The agent's verdict — and why it was honestly wrong

The agent assembled its derived vertex in **its own λ convention**
(Q ∝ 1/(4πΩ); orbital-eigenvalue denominators) and found per-shell effective
scales 0.62/0.31/0.37 Ha against the gate-demanded 3.35/3.18/2.65 — a 5–10×
shortfall with a closed-form ceiling proof on Li — and, per the note's
disagree branch, **reported the disagreement honestly and refused to update
production**. It even ran a convention audit: the locked E0 = π model
reproduces its gates through the identical code path — "no hidden contraction
factor."

**Why the audit was structurally blind:** E0 = π had been calibrated *inside*
the same convention, so the convention factor was absorbed into the constant —
making the wrong normalization self-consistent. This mechanism retroactively
explains the whole series: every "universal calibrated constant" (L3's
Z²/z_val, probe 2's √(9π), probe 3's π) was the correct vertex × an
unrecognized embedding-normalization factor.

## The campaign's final accounting

With the operator class pointed out, this model generation derives the
complete gold physics — exactly, to the published-table level. What separates
its assembly from PASS is the **absolute normalization of the atomic-to-Bloch
contraction** (where √4π, 1/Ω, and the ΔSCF-vs-orbital-eigenvalue denominator
choices land) — bookkeeping the gold fixes operationally, that the maintainer
derivation note (PR #15 §6) itself left as "operationally calibrated," and
that no packet artifact pins at L3/L4 (the L1/L2 f_c tables exist precisely to
pin it there). The frontier the benchmark now measures, stated exactly:

1. *(clean probes)* finding the spectator-mediated vertex class without
   guidance — unsolved in three attempts; and
2. *(this addendum)* anchoring the embedding normalization without
   calibration — unsolved by the agent **and not yet derived from first
   principles by anyone in this project**.

Point 2 is a completeness question for the written materials (the paper's
Sec IV/V treatment of the contraction normalization) as much as a model
capability gap. Recommended protocol change (reopening the withdrawn P4 with
four-fold evidence, in corrected form): not a dimensional-analysis check (the
agent's convention was dimensionally self-consistent) but a rule —
**"calibrated constants absorb normalization factors and make wrong
conventions self-consistent; the absolute normalization must be anchored by
one calibration-free physical quantity"** — plus, optionally, exposing one
worked contraction example at a single k-point as packet documentation
(equivalent to one row of an f_c table; design decision pending).

Artifacts: [`vertex_derived.py`](vertex_derived.py),
[`lamderived_{Li,Na,Al}.csv`](lamderived_Li.csv), updated
[`NOTES.md`](NOTES.md) / [`derivation_notes.md`](derivation_notes.md) (§5),
[`MAINTAINER_NOTE.md`](MAINTAINER_NOTE.md). Production `solution/` unchanged
(the disagree branch); the evaluated submission remains the one scored in
[`result_L4.json`](result_L4.json).
