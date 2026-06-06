# Method (Level 3): the leading frozen-core dynamical quasiparticle correction

**No web / no external lookup was used.** The derivation below is reconstructed
from the physical setup in `instruction.md` and the provided atomic data only;
the DFTK API was discovered by introspection of the pinned local install.

---

## 1. What the static pseudopotential keeps, and what it throws away

The rigorous starting point is the all-electron Hamiltonian. A frozen-core
pseudopotential integrates the core electrons out **statically**: the LDA / GTH
Kohn–Sham band that DFTK produces *is* that static reduction `H_KS`. The piece it
discards is the core's **frequency-dependent** response: a valence electron, via
the Coulomb interaction, virtually excites a core particle–hole pair (excitation
energy `ΔE_c`) and the core relaxes back. This is the leading (second-order),
core-induced contribution to the valence self-energy, and it is exactly the
*dynamical* physics that eDMFT (which keeps the full frequency dependence) needs
and that static G0W0 / hybrids / meta-GGA miss. Because the static reduction is
already in `H_KS`, the dynamical part must be derived from the **same** reduction
(not bolted on) to avoid double-counting the core.

On the valence side we are given that the Green's function is a sharp coherent
quasiparticle, `G(iω) ≈ z / (iω − H_KS)`, with weight `z`. So the only question
is how the dynamical core term modifies `z` and the band energies.

## 2. Second-order core self-energy and closure

For a valence Bloch state, second-order perturbation theory in the residual
valence–core Coulomb coupling gives a self-energy with poles at the core
excitation energies. Concentrating (as the task licenses) all of channel `c`'s
spectral weight at the single high-energy scale `ω = ΔE_c` (closure: `ΔE_c` is
several Hartree, far above the valence Fermi energy), the channel contributes one
effective bosonic mode:

```
    Σ_c(ω) = W_c / (ω − E_F + ΔE_c).
```

The mode's coupling spectral function `S_c(ω)` is the core density response of
orbital `c` weighted by the valence–core Coulomb interaction. Its **integrated
weight** (zeroth moment / closure sum rule) is the Coulomb self-energy of the core
orbital — i.e. the total strength with which orbital `c` can scatter a valence
electron through the Coulomb interaction:

```
    ∫ (dω/π) S_c(ω) = ⟨u_c² V_H_c⟩ ≡ ∫₀^∞ u_c(r)² V_H_c(r) dr     (Hartree),
```

where `u_c(r)=r R_c(r)` is the all-electron core orbital (`∫u_c² dr = 1`, one
electron) and `V_H_c(r)` is precisely that orbital's single-electron Hartree
potential — both supplied as `atomic_core_<c>.csv`. Putting all of this weight at
the closure pole, `S_c(ω) = π W_c δ(ω−ΔE_c)` with `W_c = ⟨u_c² V_H_c⟩`.

## 3. Consequence for the quasiparticle: a weight, not a shift

Expanding `Σ_c(ω)` about the Fermi level (`E_F` lies far below `ΔE_c`):

* the **constant** part `W_c/ΔE_c` is a static, state-independent shift — it is
  part of the static core–valence physics already folded into the pseudopotential,
  and it cancels in energies measured relative to `E_F`;
* the **slope** is the dynamical content,
  `Σ_c'(E_F) = − W_c/ΔE_c²` (equivalently the Kramers–Kronig weight integral
  `∫(dω/π) S_c(ω)/ω² = W_c/ΔE_c²`).

The quasiparticle weight is `z = 1/(1 − Σ'(E_F)) ≈ 1 − Σ_c (W_c/ΔE_c²)`. Hence

```
    z = 1 − λ,     λ = Σ_c ⟨u_c² V_H_c⟩ / ΔE_c² ,
```

and the quasiparticle band energy **relative to the Fermi level** is the KS energy
rescaled by the weight (the standard coherent-quasiparticle result, linearised
because the occupied bandwidth ≪ `ΔE_c`):

```
    E_QP − E_F = z · (E_KS − E_F).
```

This is a **band-narrowing**: `z < 1` pulls every occupied state toward `E_F`,
shrinking the occupied bandwidth while leaving `E_F` fixed — exactly the ARPES
puzzle (20–35 % narrowing for the alkalis, only a few % for Al).

## 4. Is the correction state-dependent?

**It is a single weight `z` per element — k-independent (one factor per element),
not per-Bloch-state.** This follows from the derivation: the core fluctuation is
**on-site and high-energy** (`ΔE_c` ≫ valence scale), so after closure the
dynamical self-energy is a **local** operator. A local self-energy is the same for
every Bloch state — it renormalises the whole occupied manifold by one factor and
does not disperse. This is precisely the eDMFT picture of a local (momentum-
independent) dynamical self-energy on the correlated, core-adjacent orbital. The
empirical KS-vs-ARPES ratio confirms it: for Na, `E_expt/E_KS ≈ 0.81` is flat
across the Γ→N path (a single weight), not an energy-dependent shift.

## 5. Why this is parameter-free and generalises

Every ingredient of `λ = Σ_c ⟨u_c² V_H_c⟩/ΔE_c²` is fixed:

* `⟨u_c² V_H_c⟩` is read directly from `atomic_core_<c>.csv` (a pure atomic
  number — the core orbital's Hartree self-energy);
* `ΔE_c` is read from `core_model.json`;
* the KS band, `E_F`, and the coefficients come from the pinned DFTK run.

There is no fitted constant, no fit to ARPES, and no per-element branch — the same
code path runs on every element. The `1/ΔE_c²` denominator makes the deep inner
core (1s) automatically negligible (its `ΔE_c` is tens of Hartree) and lets the
outer (2s) core dominate; the element dependence then rides mainly on `ΔE_2s`,
which is why metals with a softer outer core (small `ΔE_2s`, e.g. Na) narrow much
more than those with a stiffer one (Al). The hidden metals have a different core
structure, but the same atomic quantities `{u_c, V_H_c, ΔE_c}` enter identically.

## 6. Self-test (public elements)

Using the pinned LDA/GTH setup (`Ecut`, `kgrid`, `smearing` from the config):

| element | λ (derived) | z | KS-vs-ARPES RMSE | QP-vs-ARPES RMSE |
|---------|------------:|------:|-----------------:|-----------------:|
| Na | 0.162 | 0.838 | 0.413 eV | **0.089 eV** |
| Al | 0.047 | 0.953 | 0.414 eV | **0.209 eV** |

Γ-point occupied-band depth cross-check (independent of any ARPES fit):

* Na: KS 3.27 eV → QP **2.74 eV** (experiment 2.65–2.78; LDA 3.30, eDMFT 2.84).
* Al: KS 11.18 eV → QP **10.66 eV** (experiment ≈ 10.58).

The correction removes essentially all of the LDA bandwidth overestimate and lands
at the eDMFT/experiment level (~0.1 eV), with no element-specific tuning.

## Output format

`run_qp.py` runs one SCF, computes bands at `k_frac = t · endpoint_frac` for every
grid row, and emits, per point, every occupied band (`E_KS < E_F`) plus the single
lowest unoccupied band, as `E_pred_eV = z·(E_KS − E_F)` (eV, relative to `E_F`,
ascending), in columns `element,point_id,t,E_pred_eV`.
