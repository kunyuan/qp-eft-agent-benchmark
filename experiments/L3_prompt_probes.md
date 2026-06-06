# L3 prompt probes — 14 agent attempts across three prompt generations

This note records every fresh, no-web agent run on Level 3 (derive the frozen-core
quasiparticle correction) and what they reveal about the prompt and the benchmark.
All runs are honest agents given the packet + local DFTK only, audited for no web
use, scored through the hardened evaluator on the hidden metals K and Mg
(per-element bars `K < 0.17`, `Mg < 0.21`; PASS requires every element under its bar).

## TL;DR

- **0 / 14 pass.** No agent has cleared the hidden K/Mg bars.
- The single thing they all miss is the **coupling form factor
  `f_c(K) = √(4π)/K ∫ u_c (V_H_c − J_c) sin(Kr) dr`** — the *orthogonalised
  fluctuation* `(V_H_c − J_c)`, momentum-resolved — and the resulting
  **state-dependence** of `z_core`. Every agent substitutes a scalar
  (`J_c`, `⟨V_H_c²⟩`, `J_c²`, a field variance, …) and collapses to a single `z`.
- The closest is **clean2** (reframed prompt): single-`z` with coupling `J_c`,
  K 0.190 / Mg 0.211 (Mg only 0.001 over) — and that is a *numerical coincidence*
  (see "structure vs coupling" below), not better physics.
- The newest prompt (gives the paper's Eq (1)/(2)) makes agents get the **structure
  and bookkeeping right** (clean `z_core⁻¹ = 1 + Σ|coupling|²/ΔE²`, parameter-free,
  no fitted constants) but the **coupling is still wrong** — so the failure is now
  cleanly localised to the one irreducible step.

## The three prompt generations

| gen | what the prompt gave | data |
|-----|----------------------|------|
| **old** (r1–r6) | "valence QP is `G≈z/(iω−H_KS)`; the missing physics is dynamical and lives in the core; derive it" — with a "uniform/Galilean" framing that steered toward a single `z` | `u_c`, `V_H_c`, `ΔE_c` |
| **reframed** (clean1–6) | EFT "integrate out the core" framing; "given" = coherent `G` for an inhomogeneous gas; mechanism left open ("do not assume whether it renormalises `z`, shifts, or is state-dependent"); eDMFT narrative | `u_c`, `V_H_c`, `ΔE_c` |
| **new** (new1–2) | the paper's **Eq (1)** tree-level propagator + **Eq (2)** QP energy + explicit `H_KS`; a "how to derive the effective valence action" roadmap (all-electron → integrate core out → `g_v⁻¹=g_0⁻¹+δV_pp` → static=`V_PSP` / dynamic=`z_core`); closure + s-channel stated as controlled approximations | **`u_c` only** (`V_H_c` dropped — derivable from `u_c`, and handing it leaks the coupling) |

## The 14 attempts

| run | prompt | structure | coupling kernel | Na | Al | **K** | **Mg** | verdict |
|-----|--------|-----------|-----------------|----|----|-------|--------|---------|
| r1–r4, r6 | old | single-`z`(-ish) | assorted scalars | — | — | — | — | FAIL (0/5) |
| r5 | old | variance-ish | state-ish variance | 0.077 | — | 0.260 | INVALID_SHAPE | FAIL (see `L3_r5/`) |
| clean1 | reframed | **state-dep** | variance + AE-restored | 0.194 | — | 0.306 | 0.361 | FAIL |
| **clean2** | reframed | single-`z` | `J_c` (= ⟨V_H_c⟩) | 0.089 | 0.209 | **0.190** | **0.211** | **PARTIAL (closest)** |
| clean3 | reframed | single-`z` | `⟨V_H_c²⟩` | 0.076 | 0.324 | 0.308 | 0.266 | FAIL |
| clean4 | reframed | single-`z` | `J_c²` | 0.079 | 0.300 | 0.318 | 0.251 | FAIL |
| clean5 | reframed | **state-dep** | `Σ_G c_nk(G) φ̂_c` + **fitted K≈48** | 0.098 | 0.228 | 0.239 | 0.596 | FAIL |
| clean6 | reframed | single-`z` | field variance (atomic) | 0.077 | 0.334 | 0.296 | 0.277 | FAIL |
| new1 | new | single-`z` | `J_c²` | 0.080 | 0.296 | 0.320 | 0.249 | FAIL |
| new2 | new | single-`z` | `J_c²` | 0.080 | 0.296 | 0.320 | 0.249 | FAIL |

(RMSE in eV, nearest-band vs ARPES. `—` = not retained.)

## The gold, for reference

```
z_core(n,k) = 1 / (1 + Σ_c |F_c(n,k)|² / ΔE_c²)        (state-dependent)
F_c(n,k)    = Σ_G c_nk(G) f_c(|k+G|)                    (Bloch overlap of the form factor)
f_c(K)      = √(4π)/K ∫ u_c(r) [V_H_c(r) − J_c] sin(Kr) dr
J_c         = ∫ u_c(r)² V_H_c(r) dr
ε_QP − ε_F  = z_core (ε_KS − ε_F)
```

## Why everyone misses it

The form factor needs **three** things at once, and each agent gets a different one
or two — never all three:

1. **a spatial/momentum form factor** (the K-resolved sine transform), not a scalar
   moment;
2. **weighting by the core Hartree potential `V_H_c`** (the core's mean-field
   response), not the bare Coulomb field — agents who do a naive 2nd-order closure
   get the *field variance* `⟨1/r²⟩ − V_H²` (clean1/clean6) instead;
3. **the `−J_c` orthogonalisation** that makes the coupling a zero-mean fluctuation
   (`∫u_c²(V_H_c − J_c) = 0`) — agents subtract the static part *in frequency*
   (`Σ(ω) − Σ(E_F)`, which they all do) but miss this *spatial* subtraction.
   Several use `J_c` — the *subtracted constant itself* — as the coupling.

On top of that, the **state-dependence collapses** because the dev set (Na/Al, both
tight 2s cores) genuinely has < 1 % `z`-variation across the band; agents that test
Bloch-weighting on Na/Al correctly find it negligible and drop it — the variation
only matters on the hidden K (diffuse 3s). The unrepresentative dev set is a trap
no prompt wording fixes.

## Structure vs coupling — what the new prompt changed

Giving the paper's Eq (1)/(2) made the new-prompt agents get the **second-order
structure right** — both new1 and new2 derived the *correct* shape
`z_core⁻¹ = 1 + Σ_c |U_c|²/ΔE_c²`, parameter-free, prefactor exactly 1, with no
fitted constants and no `1−λ`-vs-`1/(1+λ)` convention slips. They also correctly
**rebuilt `V_H_c` from `u_c`** (the dropped column did not trip them). But the
coupling is still `U_c = J_c` (the self-Coulomb), not the form factor `f_c`:

| | coupling structure | coupling | K | Mg |
|---|---|---|---|---|
| gold | `\|F_c\|²/ΔE²` | form factor `f_c` (state-dep) | pass | pass |
| new1/new2 | `\|U_c\|²/ΔE²` ✓ (correct 2nd-order shape) | `J_c` (self-Coulomb) ✗ | 0.320 | 0.249 |
| clean2 | `J_c¹/ΔE²` ✗ (not even the right power) | `J_c` | 0.190 | 0.211 |

So the new prompt produced a **more principled but numerically worse** result than
clean2. clean2 was "best" only by a **coincidental cancellation** (it used `J_c¹`
*and* linearised `z = 1−λ` instead of `1/(1+λ)`; the weak coupling and the
over-correcting linearisation roughly cancel — see `L3_clean2/REPORT.md §2`). The
new prompt removes that coincidence and exposes the real error: *the coupling is
`J_c`, and it should be the `(V_H_c − J_c)` form factor.*

## Conclusion for the benchmark

- **0/14 is the right frontier signal.** The barrier is a single, well-defined
  physics step (the orthogonalised, momentum-resolved `(V_H_c − J_c)` form factor +
  its state-dependence), not prompt ambiguity or missing data.
- **Scaffolding helps everything *except* that step.** The new prompt fixes the
  framework (Eq 1/2), the bookkeeping (slope-not-value, parameter-free), and removes
  the leak (`V_H_c`) — and the failure cleanly concentrates on the one irreducible
  step. This confirms that step cannot be hinted without giving it away (it would
  turn L3 into L2).
- If a guided-derivation rung is wanted, it belongs in a separate **L2.5** that
  hands the state-dependent `F_c = Σ_G c_nk(G) f_c(|k+G|)` *structure* (not the
  `(V_H_c − J_c)` weight), leaving the true L3 untouched.
