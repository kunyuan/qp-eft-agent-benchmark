# L3 probe — agent "fable5": an 80%-real integrate-out, felled by a posited vertex

A single fresh, no-web agent (Claude **Fable 5**, 2026-06-10) run on **Level 3**
(derive the frozen-core quasiparticle correction from the physical setup — no
formula given), host-run protocol (no Docker): sandboxed packet copy, the
maintainer repo renamed away during the run, evaluator applied afterwards.

Solution preserved here: [`solution/run_qp.py`](solution/run_qp.py),
[`solution/qp_correction.py`](solution/qp_correction.py),
[`solution/ks_dump.jl`](solution/ks_dump.jl),
[`solution/method.md`](solution/method.md) (fable5's own derivation, verbatim).
Evaluator output: [`result_L3.json`](result_L3.json).

---

## TL;DR

- **FAIL overall**: hidden **Mg 0.199 eV PASS** (bar 0.21), hidden **K 0.464 eV
  FAIL** (bar 0.17, KS baseline 0.614). Shape checks clean (no missing points,
  no over-cap) — a genuine physics failure, not an engineering one.
- It produced a **real step-by-step integrate-out** for the pole structure
  (one pole per core s-channel; static part identified with the PSP and dropped
  as double counting; the `O(iω)` remainder giving a state-resolved
  `z_core(ν,k) = 1/(1+w)`), and **correctly killed the direct
  (core-polarization) channel quantitatively** (`w_direct ~ 3×10⁻⁴`,
  variance-bounded) — identifying the exchange (Fock) channel as the carrier.
  Both of these are the paper's actual structure.
- **The fatal step:** at the vertex it stopped contracting and **posited** a
  coupling instead — "the energy transfer ΔE_c is sudden relative to all
  screening scales, so the vertex is the unscreened nuclear `−Z/r`", normalized
  per valence electron: `w = (Z²/z_val)·|⟨φ_c|1/r|ψ̃_νk⟩|²/ΔE_c²`. Its
  `method.md` claims this matching step "is not fixed by a closed-form
  integral". **That claim is false** — the Gaussian integrate-out closes it:
  the exchange contraction's vertex is `u_c(r)·[V_H,c(r) − J_c]` (the core
  orbital's own Hartree potential, orbital-average `J_c` subtracted as the
  double-counting removal), with no `Z²` and no `1/z_val`.
- **Why Na/Al could not save it:** on Z=11/13 the wrong `Z²/z_val·(1/r)` vertex
  and the correct `(V_H,c − J_c)` vertex are nearly degenerate — its required
  Na/Al coupling ratio (~3.8) is reproduced to ~5% by the wrong form. Public
  Na **0.078** (= gold), Al **0.222** (better than gold's 0.248). K (Z=19)
  breaks the degeneracy: `Z²` over-couples the diffuse 3s/3p core →
  systematically over-narrowed band (`mean_signed = +0.415 eV`).
- Its own scan machinery (`work/scan*.py`) was already sweeping candidate
  potential weights inside the same sine-transform form factor — the correct
  vertex was within reach of the apparatus it built; the *selection* failed,
  not the machinery.

## Results

### Hidden metals (the benchmark verdict)

| element | KS baseline | **QP RMSE** | bar | mean signed | verdict |
|---------|------------:|------------:|-----|------------:|---------|
| **K** (hidden)  | 0.614 | **0.464** | 0.17 | +0.415 | FAIL |
| **Mg** (hidden) | 0.434 | **0.199** | 0.21 | +0.083 | PASS |

Overall mean RMSE 0.331 eV → **FAIL** (every element must clear its bar).

### Public dev elements (its self-check, end-to-end through `run_qp.py`)

| element | bare KS | corrected | gold |
|---------|--------:|----------:|-----:|
| Na (31 pts) | 0.413 | **0.078** | 0.078 |
| Al (61 pts) | 0.414 | **0.222** | 0.248 |

## What it derived (its physics)

Full text in [`solution/method.md`](solution/method.md). The chain:

1. Split `ψ = ψ̃ + Σ φ_c b_c` (oblique, non-orthogonal — correctly identified
   as the pseudopotential split); integrate the quadratic core modes out →
   one-pole induced potential per channel, `Σ_c(iω) = |Λ_c|²/(iω − ε_c)`.
2. Static part = level repulsion already inside `V_PSP` → dropped
   (double-counting subtraction). `O(iω)` part → `z_core⁻¹ = 1 + Σ_c Λ_c²/ΔE_c²`.
   *(Both steps match the paper's Term 1–3 / Term 4 split.)*
3. Channel selection: computed the direct (density-response) channel in full
   with closure bounds → negligible; antisymmetry essential; the effect is the
   exchange channel ("the valence hole virtually hops into the core shell").
   *(Also correct.)*
4. **Vertex (the divergence):** posited `Λ_c ∝ Z·⟨φ_c|1/r|ψ̃⟩` (sudden/unscreened
   heuristic) with an EFT-flavored `1/z_val` spectral-weight normalization —
   instead of finishing the Wick contraction, which yields
   `Λ_c = ⟨φ_c (V_H,c − J_c)|ψ̃⟩` with O(1 Ha) scale and mild Z-growth.

## Comparison across L3 probes

| agent | coupling form | state-resolved? | hidden K | hidden Mg | verdict |
|-------|---------------|-----------------|---------:|----------:|---------|
| 16 of 17 early agents | scalar single-`z` (various) | no | ~0.19–0.32 | mixed | FAIL |
| r5 | real-space variance overlap | partially | best-of-6 | INVALID_SHAPE | FAIL |
| clean2 | `z = 1−λ`, on-site closure | no (single z) | 0.190 | 0.211 | PARTIAL |
| **ex2** | `Σ_G c_nk(G)`·(sine-transform of full `V_H^core`) | **yes** | **0.135** | **0.193** | **PASS** |
| **fable5 (this run)** | `Σ_G c_nk(G)`·(sine-transform of `Z/r`), `Z²/z_val` | **yes** | 0.464 | **0.199** | **FAIL** |
| gold | `Σ_G c_nk(G)`·(sine-transform of `V_H,c − J_c`) | yes | 0.139 | 0.187 | PASS |

Notable: fable5 is the **second** agent ever to reach a state-resolved
momentum-space form factor, and its *static-sector* work (channel exclusion,
double-counting logic, Al beating gold on the public split) is the strongest
seen at L3 — but ex2 passed with a *less* principled potential because full
`V_H^core` happens to scale benignly with Z, while fable5's `Z²` does not.
The benchmark's K/Mg split did exactly its job: K is where coupling-strength
errors on diffuse cores explode.

## Lesson (for SETUP_L3 / prompt design)

The `SETUP_L3` requirement "*do a real step-by-step integrate-out — don't
posit/guess the coupling*" remains the live discriminator. fable5 obeyed it
for the pole structure and channel selection but violated it **once**, at the
vertex, substituting a physically-narrated guess for the final contraction —
and that single posit is the entire failure. A sharper phrasing ("the vertex
itself must be the closed-form result of the contraction; any 'matching'
constant you cannot derive is a red flag") would target exactly this failure
mode. Conversely: an agent's own admission in `method.md` that a step is "not
fixed by a closed-form integral" is a reliable audit signal that the
derivation left the controlled path.

## Run facts

- Model: Claude Fable 5 (`claude-fable-5`), fresh subagent, no web, packet-only
  sandbox (`qp_eft_L3_fable5_run/`), maintainer repo hidden during the run.
- Wall clock ≈ 4 h 52 min; ≈ 315k output tokens; 58 tool calls. Work pattern:
  long (20–30 min) reasoning stretches punctuated by bursts of scripted
  numerics (3 rounds of coupling-form scans + GTH pseudopotential analytic
  reconstruction for double-counting accounting).
- Evaluator: `validate_submission.py --level 3` on host (DFTK 0.7.25 pinned,
  same protocol as the README's host-run validation). Note
  `solution/run_qp.py` resolves the Julia project via an absolute path,
  overridable with `QP_JULIA_PROJECT`.
