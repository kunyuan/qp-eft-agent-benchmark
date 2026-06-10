# Maintainer note — first-principles derivation of the (V_H,c − J_c) vertex (Li worked example)

**Maintainer-only.** This note derives the gold dynamical vertex from the bare
electron–ion action via the dual-fermion integrate-out, with every checkpoint
verified numerically against the paper (arXiv 2604.25199) and against the
benchmark's own published artifacts (the L1 `fc_table` files). It exists
because the vertex is the discriminating step of Level 3: of the 18 fresh
agents probed so far, every failure mode is a deviation from one of the four
checkpoints below. Companion post-mortem:
[`../experiments/L3_fable5/REPORT.md`](../experiments/L3_fable5/REPORT.md).

Do not place any of this in a solving agent's context.

---

## 0. Starting point: the electron–ion action

One Li ion (Z = 3) plus the electron field, imaginary time:

```
S[ψ̄,ψ] = ∫dτ { Σ_σ ∫d³r ψ̄_σ (∂_τ − ∇²/2 − Z/r − μ) ψ_σ
              + ½ Σ_σσ′ ∫∫ ψ̄_σ ψ̄′_σ′ |r−r′|⁻¹ ψ′_σ′ ψ_σ }
```

Orbital basis {φ_1s, φ_s}: variational 1s, φ = √(α³/π)e^{−αr},
α = Z − 5/16 = 2.6875.

## 1. Reference system: push the valence up

Dual fermion needs an exactly solvable reference. Define
S_ref = S + ψ̄Rψ with R = Λ·P_val (valence projector × large scale Λ):
an interacting *atom* whose valence levels sit at E₀ + Λ, so the low-energy
Fock space is just the core configurations (Li⁺: 1s²) — finite-dimensional,
exactly solvable. Physical system and reference differ by the ONE-BODY term
−ψ̄Rψ only.

## 2. Dual-fermion transform: integrate the core out

Hubbard–Stratonovich on the one-body difference (dual field φ_d), then
integrate ψ exactly inside the reference:

```
e^{ψ̄Rψ} ∝ ∫Dφ_d e^{−φ̄_d R⁻¹ φ_d + φ̄_d ψ + ψ̄ φ_d},
⟨e^{φ̄_d ψ + ψ̄ φ_d}⟩_ref = exp{ φ̄_d G_ref φ_d + (higher cumulants) }
```

Truncating at the Gaussian level (higher reference cumulants enter the valence
sector at higher order in 1/Λ) and rescaling φ_d → g_c⁻¹ψ gives the paper's
Eqs (14)–(17):

```
g_v⁻¹ = g₀⁻¹ + δV_pp,
δV_pp = g₀⁻¹R⁻¹g₀⁻¹  +  Σ_c  +  Σ_c G_ref Σ_c
        └─ projector ─┘
        (diverges in the core sector: projects core states out)
```

No physics has been approximated yet; everything sits in the reference's exact
Green's function G_ref.

## 3. Solve the Li reference (1st order in v; sufficient to O(1/Λ²))

| state | energy | value (Ha) | paper |
|---|---|---:|---:|
| ξ₀ = 1s² | E₀ = 2ε + J | −7.223 | −7.219 |
| core hole a_c\|ξ₀⟩ | ε | −4.451 | −4.449 |
| excitation | ΔE = \|ε\| − J | **2.771** | 2.770 |

with ε = α²/2 − Zα, J = ⟨φφ|v|φφ⟩ = 5α/8 = 1.680. **First appearance of J:
in the energy denominator.** The ground state carries an O(1/Λ) admixture

```
|Φ₀⟩ = |ξ₀⟩ − (1/Λ) Σ_{s,c} M_sc |ξ₁^{sc}⟩ + O(Λ⁻²),
M_sc = ⟨ξ₁^{sc}|V̂_ee|ξ₀⟩.
```

## 4. 1/Λ expansion of G_ref, valence sector

**Particle part** (add an electron to a pushed-up valence level) → the static
core mean field Σ^HF_st = Σ_c n_c[⟨cs|ct⟩ − ⟨cs|tc⟩] = **Term 1** (Hartree,
spin factor 2: 2u·δ(r−r′)) **+ Term 2** (same-spin exchange kernel −φφ′v).

**Hole part** (the crucial one): a_t|Φ₀⟩ lands on the admixture → core-hole
intermediate states:

```
(G_h)_st(iω) = (1/Λ²) Σ_c M_sc M*_tc / (iω + ΔE_c)
```

Slater–Condon for M_sc (excite 1s↑ → s↑):

- j = 1s↓ (opposite-spin partner): exchange dies by spin; the direct term
  survives → ⟨s|u|c⟩ with u(r) = (1−e^{−2αr})/r − αe^{−2αr}, the Hartree
  potential of ONE core electron.
- j = c (self-orbital): ⟨sc‖cc⟩ = ⟨sc|cc⟩ − ⟨sc|cc⟩ ≡ 0 —
  **antisymmetry kills self-scattering exactly** (direct = exchange for the
  self-orbital; this is the "exchange origin" of the subtraction).

So the bare vertex function is a(r) ≡ u(r)φ(r), with M_sc = ⟨s|a⟩.

## 5. Where −J is born: from matrix elements to the operator kernel

Decompose along φ:

```
a = ã + Jφ,   ã ≡ (u − J)φ,   ⟨φ|ã⟩ = J − J = 0   (numerically 1e−13)
```

On the orthogonal valence subspace |a⟩⟨a| ≡ |ã⟩⟨ã| — the contraction itself
does not distinguish the two. But the core sector of δV_pp is owned by the
diverging projector (core states are projected out; their statics are already
in ΔE_c and in the PSP), so the pole kernel must have **zero action in the
core sector**. The unique representative satisfying that is ã:

```
Term 4 = φ(r)φ(r′) [u(r) − J][u(r′) − J] / (iω + ΔE_1s)        (★)
```

The displaced pieces do not vanish; the identity
|a⟩⟨φ| + |φ⟩⟨a| − J|φ⟩⟨φ| = |ã⟩⟨φ| + |φ⟩⟨ã| + J|φ⟩⟨φ| shows they land in the
static projection sector — **Term 3**: −φφ′[u(r) + u(r′) − J]. Terms 1–4
together are the paper's Eq (37).

**J appears twice from a single partition** — in the denominator
(ΔE = |ε| − J) and in the vertex subtraction (u − J). Omitting either is the
same bookkeeping error (double counting) seen from two sides; the mechanism
that mandates zero self-coupling is antisymmetry (§4, j = c cancellation),
and the value J is pinned by the zero-diagonal condition ⟨φ|vertex|φ⟩ = 0.

## 6. To z_core and the crystal (numerical closure)

Subtract Term 4's static value at ε_F (absorbed by the PSP fit), linearize
the Dyson equation:

```
z_core(νk) = 1 / (1 + Σ_c |F^c_νk|²/ΔE_c²),
F^c_νk = Σ_G c_νk(G) f^c_{|k+G|},
f^c_K = (√4π / K) ∫ u_c [V_H,c − J_c] sin(Kr) dr
```

Verified in this note's preparation (scripts run against the repo's own data):

| check | derived | reference |
|---|---|---|
| Li ΔE, J, E₀ | 2.771 / 1.680 / −7.223 Ha | paper: 2.770 / 1.680 / −7.219 |
| ⟨φ\|(u−J)\|φ⟩ | ~1e−13 | 0 (zero-diagonal condition) |
| Li z_Γ (BCC, Ω = 145.9 Bohr³) | ≈ 0.73–0.75 | paper: ≈ 0.75 |
| Na 2s f(0)/f(0.5)/f(1.0) from L3 raw u_c (own Poisson V_H,c, J_2s = 1.146 Ha) | 1.444 / 1.184 / 0.636 | L1 `fc_table_2s.csv`: 1.451 / 1.190 / 0.640 (≤0.5%) |

**Normalization convention (a real trap):** the gold contraction uses the
DFTK cell-normalized c_νk(G) directly against f — NO 1/√Ω. Inserting the
"normalized plane-wave matrix element" 1/√Ω makes Li z_Γ ≈ 0.997 (correction
vanishes) — numerically the same trap that made agent fable5 reject the
correct vertex family as "30× too weak" on Na. The L1/L2 f_c tables exist
precisely to pin this convention; at L3 it is a fourth silent failure point.

## 7. The four L3 checkpoints (audit taxonomy)

Every step of the derivation is forced — there is no place to "posit". The
four points where probed agents have actually deviated:

1. **Same-spin Slater–Condon carried to the end** (§4): the j = c
   cancellation is antisymmetry's veto on self-scattering. Skipping it →
   wrong vertex family (fable5: bare −Z/r).
2. **Non-orthogonality awareness** (§5): smooth/pseudized states overlap φ_c
   (Na 2s: ⟨ψ̃|φ_c⟩ ≈ 0.17), so the kernel's φ_c-diagonal matters. The audit
   question "what is ⟨φ_c|your vertex|φ_c⟩, and what should it be?" probes
   exactly the one number the contraction leaves free. Gold: 0. ex2's
   V_H^core: ~few Ha. fable5's Z/r: Z⟨1/r⟩_c, huge. One question separates
   all three.
3. **ΔE_c decomposition** (§3): the packet's DeltaE_c_Ha already contains
   −J_c (Li: ΔE = |ε| − J, verifiable in closed form). Treating ΔE_c as a
   black box while using a bare vertex double-counts J.
4. **Contraction convention** (§6): where the volume factor lands in F is
   calibration, not physics; at L3 nothing in the packet pins it.

Suggested SETUP_L3 sharpening (not yet applied): require the submission's
method.md to state the value of ⟨φ_c|vertex|φ_c⟩ and justify it. This is a
one-line machine-checkable invariant (∫u_c²(V_H,c − J_c)dr = 0) that cannot
be satisfied by narrative.

---

*Provenance: derivation follows the structure of arXiv 2604.25199 Secs I–II
and Eq (37); all numerical checks were re-run independently from the repo's
packet data and the paper's closed forms (2026-06-10). Grew out of the
fable5 L3 post-mortem discussion.*
