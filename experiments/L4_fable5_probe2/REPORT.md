# L4 probe 2 — fable5 under Li-first: complete enumeration, a conjectured prefactor, and a self-predicted failure

Second probe of **Level 4**, Claude **Fable 5**, 2026-06-11 — the first run of
the **Li-first protocol** (PR #18) on a fresh agent: new sandbox, zero
inherited artifacts, contamination-free resume templates (sessions resumed only
by pointing at the agent's own logs). 3 sessions, ~2 h 40 min net working time
(session 1 cut by an account usage limit; session 2 by an idle-wait watchdog
kill — its background jobs completed anyway; session 3 finalized).

Artifacts: [`NOTES.md`](NOTES.md) (the agent's own log — read this first),
[`derivation_notes.md`](derivation_notes.md), [`solution/`](solution/)
(`run_qp.py`, `method.md`, `ks_band.jl`), [`result_L4.json`](result_L4.json)
(evaluator output).

## Verdict

| element | KS baseline | **submission** | bar | gold | mean signed | verdict |
|---------|------------:|---------------:|-----|-----:|------------:|---------|
| **K** (hidden) | 0.614 | **0.304** | 0.17 | 0.139 | **+0.246** | FAIL |
| **Mg** (hidden) | 0.434 | **0.240** | 0.21 | 0.187 | +0.155 | PARTIAL |

Overall **FAIL** (mean 0.272 eV); `beats_gold_all` false. Shape checks clean.
Public elements: Na **0.0848** (baseline 0.078), Al **0.2481** (baseline
0.248, matched to three digits); Li anchor closed by construction
(calibration).

## What the Li-first protocol changed (the A/B result)

| | probe 1 (pre-revision) | **probe 2 (Li-first)** |
|---|---|---|
| "fluctuation channels are too small" established | ~20 h | **31 min**, with a Li⁺ polarizability sum-rule cross-check |
| contraction enumeration | asserted (D²+X²), never derived | **complete from the start** ({A,B,C} × {DD,DX,XX}, completeness proof for the 1s² core, λ per class) |
| bare −Z/r vertex (killed L3 probe) | n/a | **rejected at birth** by the three-anchor Z-pattern |
| Li closed forms | n/a | exact match to the published variational treatment (Z_s = 43/16, E = −7.2227 Ha, J = 1.6797 Ha) — derived independently |

## The remaining gap — and how it failed

After honestly killing every natural one-body vertex it could derive (GTH-NL,
ΔV = V_AE−V_GTH, kinetic, −Z/r, contact, and two more), the agent found the
empirical fingerprint of the correct coupling — **vertex ≈ (universal ~5.3 Ha)
× core overlap, with 1/Δ² denominators** — and shipped

```
λ_kn = Ξ² Σ_a |⟨χ_a|ψ̃_kn⟩|² / Δ_a²,   Ξ² = 9π Ha²  (CONJECTURAL — disclosed)
```

calibrated on the Li anchor (28.250 vs 9π = 28.274, 0.08%), s-channels only
(zero-range argument + Na band-shape confirmation). This is structurally
isomorphic to the gold (per-s-channel overlap-type amplitude squared over an
excitation energy squared) but with the **vertex content guessed, not
derived**: gold's form factor weights the overlap by the orthogonalized core
potential (V_H,c − J_c), whose effective strength in the overlap region is
exactly the "unexplained universal ~5.3 Ha" the agent measured but did not
derive. The two parameterizations agree on compact 2s cores (hence Na/Al
match) and diverge on diffuse cores — S²/Δ² grows much faster than the gold's
F²/ΔE² as the core shell becomes shallow and extended.

**The agent predicted its own failure mode before submission**: its
self-initiated Ca cross-check (config self-built, eDMFT 3.24 eV anchor from
the README table) overcorrected ~5× on Ca's diffuse 3s semicore; it diagnosed
the point-core limit breaking for large-radius s shells, disclosed it as a
validity-domain caveat (ledger L1.10), and **refused to patch it** ("any
Ca-fitted form factor would violate no-tuning"). Hidden K has the same
diffuse 3s semicore: mean signed error +0.246 eV — overcorrection, exactly as
predicted. A submission that ships with an honest falsification of itself.

## Ledger / audit notes

- Vertex diagonal: **nonzero** (Ξ = 5.317 Ha, channel-independent), argued
  consistent via Phillips–Kleinman energy bookkeeping (⟨χ_a|V_PK(ε_v)|χ_a⟩ =
  Δ_a; NC-PSP absorbs the first-order energy dependence). Structurally
  different from the gold's zero-mean fluctuation vertex — the divergence
  point to adjudicate in any follow-up.
- No fitting headroom used beyond the single disclosed constant; e2e verified
  byte-identical from config+grid-only inputs; element-blind; own AE-LDA
  atomic solver (eigenvalues validated against LDA tables).
- Contamination ledger: resume prompts for sessions 2–3 contained no physics
  paraphrase (pure pointers to its own logs). The Li-first SETUP itself states
  the anchor and the "undershoot ⇒ missing contraction" logic — intentional,
  user-approved protocol design (PR #18), uniform for all future probes.

## Cross-probe standings (hidden set)

| probe | route | K | Mg | verdict |
|-------|-------|--:|---:|---------|
| L3 fable5 | posited bare −Z/r vertex | 0.464 | 0.199 | FAIL |
| L4 probe 1 | direct channel, fully audited | — | — | interrupted at the exchange doorstep |
| **L4 probe 2** | complete enumeration + conjectured prefactor | **0.304** | **0.240** | **FAIL (self-predicted domain)** |
| gold | derived (V_H,c − J_c) vertex | 0.139 | 0.187 | PASS |

## The one-line lesson

The Li-first stage moved the frontier from *"which channel?"* (now solved
structurally in minutes) to *"derive the vertex content"* — the ~5.3 Ha that
probe 2 measured, conjectured, and disclosed is precisely the orthogonalized
core-potential strength that must come out of the antisymmetrized
integrate-out. That is now the sharpest remaining test L4 poses.
