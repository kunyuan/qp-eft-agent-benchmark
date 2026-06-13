# A first-principles effective field theory for the valence electron of Li
### (what "integrating out the atomic 1s and 2s" can and cannot mean)

Companions: AUDIT-EFT_for_pseudopotential.md (1s machinery),
NOTE-2s-channel-derivation.md (choices 1/2/3, unbiased J's, ED-DMFT).

## 1. Obstruction theorem (why there is no "remove-the-2s" local EFT)

An EFT with local-in-time Wilson coefficients requires every integrated-out
pole to lie above the cutoff. The atomic 2s sector's charge-fluctuation poles
sit at ±U/2 = ±2.39 eV around μ — *inside* the band (E_F = 4.74 eV,
E_c = 7.5 eV). The gradient expansion of any "2s-removed" effective action
has convergence radius U/2 and fails at band energies; its zeroth-order
resummation is Hubbard-I, which falsely makes Li a Mott insulator. Removing
the 2s wholesale also removes the electron the EFT is supposed to describe.
By contrast the 1s sector's poles sit at ΔE_c ≈ 70 eV: ω/ΔE_c ≤ 0.1 across
the band, and a local expansion is exponentially good. **The two shells are
EFT-inequivalent; the construction must be two-stage.**

## 2. Two-stage Wilsonian construction

**Stage 1 (cutoff Λ₁ ≈ 30 eV): integrate the 1s multiplet sector.**
Operator tower and matched coefficients (all audited):
- static one-body δV_pp(r,r′)/v — the four-term kernel; the relevant operator;
- iω-linear one-body term, coefficient λ_1s = f²/(vΔE_c²) = 0.003 — the
  leading irrelevant correction, suppression (ω/ΔE_c);
- induced four-fermion (core-polarization mediated) — O((V/ΔE_c)⁴), dropped.

**Stage 2 (cutoff Λ₂ ≲ 1–2 eV < U/2): integrate (i) band states outside a
shell around the Fermi surface, (ii) the high-energy *part of the local
spectral weight* (the 2p1h/1p2h continuum and atomic-multiplet satellites).**
What survives is the quasiparticle field ψ̃ on the Fermi surface — a
Shankar–Polchinski EFT:

    S = Σ_{k,ω} ψ̃†[ -iω + v_F^*(k−k_F) ]ψ̃
        + (marginal) Landau function F(θ,θ′) + BCS channel
        + irrelevant towers O(ω/Λ₂, (k−k_F)²/...)

Power counting: standard FS scaling; quadratic terms and forward/BCS
four-fermion marginal, everything else irrelevant.

## 3. Matching (the Wilson coefficients, computed)

| coefficient | matching object | value (Li, KS free-e band) | sensitivity/real bands |
|---|---|---|---|
| field normalization z (c = √z ψ̃ + incoh.) | 1/(1+λ_ω), λ_ω = (U/E_F)²J_ω, J_ω = 0.139 | **0.88** | m_KS ≈ 1.35–1.5 → z = 0.76–0.80; choice-3 bracket 0.79–0.87 |
| Fermi velocity v_F^* (bandwidth) | (1+λ_k)/(1+λ_ω)·v_F^KS | **m*/m = 1.14** | real bands → 1.26–1.32 |
| **λ_k ≡ ∂Σ/∂ε_k at k_F (NEW)** | CRN-paired MC, ±4%/±10% windows | **−0.0004 ± 0.0004** | Σ(k,0) flat through FS; curvature only near band edges |
| Landau F (tree) | U·ρ(E_F) | 0.75 | 2nd-order matching = next work item |
| λ_1s (stage 1) | f²/(vΔE_c²) | 0.003 | — |

Cross-element (KS free-e bands): λ_ω = 0.141/0.278/0.455 → z = 0.88/0.78/0.69,
m*/m = 1.14/1.28/1.46 for Li/Na/K.

## 4. The two locality statements (the campaign's physics, in one place)

The single most clarifying outcome:
- **Locality of FORM is emergent and excellent:** λ_k ≈ 0 at the FS, so
  m*/m = 1/z to better than 1% — a k-independent self-energy ansatz (the
  paper's ε_K → ε_K/(1+Δ); eDMFT) has the *right structure* for this model.
- **Locality of VALUE is wrong by 3.3×:** computing λ_ω with local/DOS
  phase space (impurity rung) gives 0.46–0.60 where the momentum-conserving
  value is 0.14 (free-e) – 0.26 (real bands).

So the paper's Stage-2 operator was the right operator with the wrong
matching: the coefficient must be the momentum-conserving **valence**
(2s-fluctuation) λ_ω — not a core quantity (λ_1s = 0.003 is negligible) and
not an impurity-level value.

## 5. Predictions and falsifiable contrasts

With real-band matching for Li: z ≈ 0.78–0.82, m*/m = γ*/γ ≈ 1.25–1.30.
- Specific heat / dHvA mass enhancement: ≈ 1.26 (cf. experimental
  γ-enhancement for Li ~ 1.4–1.6 incl. phonons; electronic part consistent).
- ARPES: occupied-bandwidth narrowing factor m/m* ≈ 0.78; QP weight z ≈ 0.8.
- **Sharp contrast with the core mechanism:** our missing weight 1−z ≈ 0.2
  is a *smooth 2p1h continuum within a few eV* of the QP (no sharp line);
  the published z_core = 0.75 would instead put ~26% weight in a satellite
  ~104 eV below E_F. Valence-band XPS/EELS lineshapes distinguish these
  unambiguously.
- Kohn–Luttinger: the marginal BCS coupling induced at 2nd order is
  attractive in ℓ ≥ 1 — exponentially small Tc, irrelevant in practice.

Validity window: |ω|, v_F|k−k_F| ≲ 1–2 eV; the EFT breaks at the
charge-fluctuation scale U/2 ≈ 2.4 eV where the integrated-out continuum
re-enters. Remaining matching work items: Landau F at 2nd order, real KS
bands (top priority — dominates current error), full γ⁽⁴⁾ frequency
structure, U×V_inter (screening/collective) class.
