# L4 probe 3 — fable5 under the full revised protocol: complete enumeration in one hour, best-ever public scores, and the calibration debt paid on K

Third L4 probe, Claude **Fable 5**, 2026-06-11/12 — first run under the complete
revised protocol (PR #18 Li-first + PR #20: single-field enumeration algebra
with the leg-assignment clause, anchor-calibration prohibition, vertex-diagonal
bookkeeping discipline, working protocol in the task spec) AND the first with
the lengthened harness watchdog (30-min stall limit: zero thinking-stretch
kills; the two session deaths were platform socket drops during a service
outage). 3 sessions, ~4.5 h net.

Artifacts: [`NOTES.md`](NOTES.md) (the agent's log — the whole story),
[`derivation_notes.md`](derivation_notes.md), [`solution/`](solution/),
[`result_L4.json`](result_L4.json).

## Verdict

| element | KS baseline | **submission** | bar | gold | mean signed | verdict |
|---------|------------:|---------------:|-----|-----:|------------:|---------|
| **K** (hidden) | 0.614 | **0.395** | 0.17 | 0.139 | **+0.367** | FAIL |
| **Mg** (hidden) | 0.434 | **0.258** | 0.21 | 0.187 | +0.148 | PARTIAL |

Overall **FAIL** (mean 0.326). Public: Li z_Γ = 0.771 (anchor 0.747), Na
**0.094**, Al **0.219 — the first probe to beat the gold's 0.248 on a public
element**. Shape checks clean; e2e from config+grid-only inputs verified for
all three dev elements.

## What the revised protocol demonstrably fixed

1. **The leg-assignment clause worked on first contact**: within ~50 minutes
   the agent produced the 16-leg-assignment table of the quartic Coulomb term
   with a ΔN_c column, kept the number-non-conserving classes (its own
   annotation on K6: *"This is the assignment an enumeration that conserves
   core number separately would miss"*), and enumerated channels A/B/C — the
   core-hole sector entered a probe's enumeration for the first time in the
   series (probe 1 never reached it in 20 h; probe 2's completeness proof
   excluded it by construction).
2. **The three-anchor gate killed every wrong vertex within minutes**: bare
   nuclear −Z/r (the L3 killer) died on first cross-element test (Li 6× small,
   Na 3× big, Al 7× big); the pure-overlap/PK family died on the Li anchor;
   the atomic-potential vertex died on the 7:1 element span. Hypothesis-test
   cycles that took the L3-era probes whole runs now cost minutes each.
3. **The disclosure regime worked**: the final ledger states E0 = π E_h as
   **CALIBRATED** ("derivation incomplete at that point per SETUP wording"),
   the single non-derived constant, with the per-shell required magnitudes
   (3.35/3.18/2.65 Ha) shown rather than hidden, and the caveat
   *"universality across diffuse n=3 concealed cores is an assumed
   extrapolation"* — the K failure was self-predicted in the submission.

## How it still failed — the trilogy's invariant

Shipped model: λ_nk = Σ_c E0²·Q_c,nk/(ε_nk−ε_c)², E0 = π Ha. Structurally
the same family as probe 2's 9π·S²/Δ² and, at one remove, the L3 probe's
Z²/z_val: **a universal calibrated constant × per-shell overlap kinematics /
energy-difference²**. The constant absorbs what should be the derived vertex
content — the orthogonalized core-potential strength (V_H,c − J_c), whose
per-shell variation is exactly the 25% spread (3.35→2.65 Ha) the single π
flattens. On compact 2s cores the flattening is benign (hence the best-ever
public scores); on K's diffuse 3s it is the whole error: mean signed +0.367 eV
of overcorrection.

**The overfitting signature across probes 2→3 is the cleanest evidence**: probe
3 fits the three public anchors better than probe 2 (Na 0.094 vs 0.085 ≈ tie;
Al 0.219 vs 0.248-matching; Li gate inside tolerance) yet does *worse* on both
hidden elements (K 0.395 vs 0.304; Mg 0.258 vs 0.240). Tightening a
calibrated constant against compact-core anchors deepens the commitment to the
wrong diffuse-core scaling.

It also (again) parked channel B's q→0 singularity as "regularized by valence
screening, expect small" — the unfinished contraction whose proper treatment
(the divergent long-wavelength piece is the statically-counted orbital
average; subtracting it is the −J_c physics) is where the derived vertex
lives.

## Series standings (hidden set; gold 0.139 / 0.187)

| probe | protocol | public best | K | Mg | verdict |
|-------|----------|------------|--:|---:|---------|
| L3 fable5 | L3 (formula-adjacent) | Na 0.078† | 0.464 | 0.199 | FAIL |
| L4 probe 1 | pre-revision | Na 0.235 | — | — | interrupted |
| L4 probe 2 | Li-first | Na 0.085 / Al 0.248 | 0.304 | 0.240 | FAIL |
| **L4 probe 3** | **full revision** | **Na 0.094 / Al 0.219** | **0.395** | **0.258** | **FAIL** |

† L3 provided the propagator structure; not comparable.

## The measured capability boundary (three-fold replicated)

With every scaffolding, epistemic, and infrastructure obstacle removed, the
model now reliably: performs the complete single-field contraction enumeration,
self-computes all atomic inputs, runs honest multi-anchor falsification at
minutes-per-hypothesis, maintains audit-grade ledgers, and predicts its own
failure domains. What it has not done in three independent attempts under
three escalating protocols: **carry the antisymmetrized contraction of its own
channel B through the q→0 bookkeeping to the closed-form vertex** — it
substitutes a calibrated universal constant at precisely that point, every
time. That step — the birth of (V_H,c − J_c) — is the sharpest measured edge
of this model generation, and the standing definition of what passing L4 by
derivation requires.

## Infrastructure note

First probe under the 30-minute stall watchdog
(`CLAUDE_ASYNC_AGENT_STALL_TIMEOUT_MS=1800000`): multiple 15–27-minute
derivation stretches completed without kills; both session deaths were
platform socket drops (status page: "Minor Service Outage"), resumed
losslessly via the DESIGN §5d template with zero physics paraphrase in the
resume prompts (contamination-free relay achieved for the first time in the
series).
