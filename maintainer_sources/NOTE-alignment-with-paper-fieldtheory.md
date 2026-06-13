# Aligning the valence EFT with the paper's field theory

The previous note (NOTE-valence-EFT.md) wrote stage 2 in the
Shankar–Polchinski **Fermi-surface** normal form (field = shell modes, z =
field normalization). The paper arXiv 2604.25199 uses a different — and for
this problem more natural — normal form:

> a **full-space electron field**, LO = the KS Hamiltonian obtained by
> resumming the LDA class of self-energy, and corrections entering as a
> **multiplicative spectral-weight (dilation) factor** z(ν,k) = 1/(1+Δ(ν,k))
> on the propagator.

The two are the same EFT in different normal forms, related by a field
redefinition. This note makes the map explicit and re-expresses every
campaign result in the paper's variables. Companions: NOTE-valence-EFT.md,
NOTE-2s-channel-derivation.md, AUDIT-EFT_for_pseudopotential.md.

## 1. The field redefinition connecting the two normal forms

Near the QP pole the dressed propagator is

    G(k,ω)^{-1} = ω(1+Δ(k)) − (ε_k^{KS} − μ) − δΣ_static(k)            (paper)

           ⇔   G(k,ω) ≈ z(k) / (ω − E_k^*) + G_incoh,
                z(k) = 1/(1+Δ(k)),  E_k^* = [ε_k^{KS}−μ+δΣ(k)]/(1+Δ(k)).

The FS-EFT field of the previous note is c(k) = √z(k)·ψ̃(k) + (incoherent).
So **the paper's dilation z(ν,k) is identically the field normalization of
the valence EFT** — kept as a propagator factor instead of absorbed into the
field. Nothing physical differs; the paper's form is preferred here because it
states the result band-wide (all ν,k), not just on the FS.

## 2. Term-by-term identification

| paper object | EFT meaning | campaign value (Li) |
|---|---|---|
| LO KS Hamiltonian ε_k^{KS} | stage-2 LO (resummed LDA class of V_inter) | the approved approximation |
| dilation denominator Δ(ν,k) = Σ_c \|F_c\|²/ΔE_c² | **Σ over integrated-out channels of** \|coupling\|²/(energy)² | see §3 |
| z(ν,k) = 1/(1+Δ) | QP weight = field normalization | 0.74–0.90 (dispersing, §4) |
| frozen-core projection (Eqs. 88–94) | the **core** term in the Δ-sum | Δ_core ≈ 0.003 (after the metric fix) |
| k-independent rescaling ε_K→ε_K/(1+Δ) | the ansatz δΣ_static ≈ 0, single Δ | **validated**: δΣ(k)≈0 (§3), but Δ disperses |

## 3. The Δ-sum: which channels, and the two errors

Paper's denominator, read as a sum over everything integrated out:

    Δ(k) = Δ_core(k)  +  Δ_2s(k)  +  …
         = [Σ_c |F_c|²/(v ΔE_c²)]  +  [(U/E_F)² J_ω(k)]  +  (inter-site)

- **Δ_core**: the paper's only retained term. Correctly normalized (metric of
  Eqs. 88–89) it is 0.003 — negligible. The published 0.75 is this term with
  the missing 1/v (the documented bug), which numerically — and
  coincidentally — landed near the valence answer below.
- **Δ_2s**: the **valence charge-fluctuation channel**, U = I−A (ΔSCF),
  momentum-conserving phase space J_ω(k). This is the physical origin of the
  dilation. It is exactly the kind of term the paper's z(ν,k) is built to
  hold; it was simply never summed.
- δΣ_static(k): the paper sets it to zero. Confirmed: CRN-MC gives
  ∂Σ/∂ε_k|_{k_F} = −0.0004 ± 0.0004, and Σ(k,0) is flat to <1% across the
  occupied band. **The paper's k-independent-rescaling ansatz is correct in
  form.** (Locality of FORM ✓.)

The single non-form error is the matching VALUE: Δ must be the
momentum-conserving valence number, not the (3.3× larger) local/impurity one,
and not the core one.

## 4. Band-wide dilation in the paper's variables (the concrete new result)

The paper writes z(ν,k) as a function but evaluates a single dilation. The
campaign now supplies the full dispersion. CRN-MC, one-band KS model
(J_ω(k) ≡ E_F²·⟨2p1h⟩, errors ~3e-4):

```
 k/k_F   (ε_k−μ)/E_F    J_ω(k)    z(k) free-e    z(k) real-band (m*=1.4)
 0.00      −1.00        0.173       0.85              0.74
 0.36      −0.87        0.171       0.85              0.75
 0.54      −0.71        0.168       0.85              0.75
 0.72      −0.48        0.157       0.86              0.76
 0.90      −0.19        0.146       0.87              0.77
 1.00       0.00        0.139       0.88              0.78
 1.13       0.28        0.124       0.89              0.80
 1.26       0.59        0.107       0.90              0.82
```

z increases from band bottom to top — deeper hole states have more 2p1h decay
phase space. The **occupied-bandwidth narrowing** (the ARPES observable),
W*/W_KS = ∫ over occupied (1+δΣ)/(1+Δ):

    free-electron bands:  W*/W_KS = 0.85   (single-Δ would give 0.88)
    real Li KS (m*≈1.4):  W*/W_KS = 0.74   (single-Δ would give 0.78)

The z(k) dispersion adds ~3% narrowing beyond the single-number dilation — a
genuine, computable refinement of the paper's normal form (its z(ν,k) holds
this for free; one just evaluates it at each k).

## 5. Net reconciliation

The paper's field theory is, structurally, **the correct EFT for this
problem**:
- LO KS from resummed LDA ✓
- corrections as a multiplicative dilation z(ν,k) on the KS propagator ✓
- k-independent static part (form of the rescaling) ✓ (δΣ≈0, derived)

Its two defects are both in the **matching of the dilation**, not the field
theory:
1. **Channel**: the dilation was matched to the core term only; the dominant
   contribution is the valence 2s charge-fluctuation channel (never summed).
2. **Bug**: the core term itself was computed without the Eqs. 88–89 metric
   (missing 1/v), inflating Δ_core from 0.003 to ~0.34 and accidentally
   producing a number near the (uncomputed) valence value.

Corrected and matched to the valence channel, the paper's own z(ν,k)
formalism gives z_Li ≈ 0.78–0.82 with band dispersion 0.74→0.82 and
occupied-bandwidth narrowing W*/W ≈ 0.74 — in the measured ARPES range and
with a smooth few-eV incoherent tail (not a 104-eV satellite). Same field
theory; right channel; no bug.
