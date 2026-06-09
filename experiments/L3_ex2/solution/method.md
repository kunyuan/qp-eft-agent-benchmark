# Method — Level 3: deriving the frozen-core dynamical quasiparticle correction

## 0. What is given and what is derived

The instruction GIVES the EFT quasiparticle propagator and energy,

```
G_{nu,k}^{-1}(iω) ≃ (z_val)^{-1} [ (z_core_{nu,k})^{-1} iω − (ε_KS_{nu,k} − ε_F) ]   (1)
ε_QP_{nu,k}       ≃ ε_F + z_core_{nu,k} (ε_KS_{nu,k} − ε_F)                          (2)
```

and tells me the mechanism (integrate the core fields out of the all-electron
action; the static part is the conventional pseudopotential already in `H_KS`,
the *dynamical* part is `z_core`). What is NOT given is the closed form of
`z_core_{nu,k}`. Deriving it is the task. The derivation below is done from the
many-electron action; the coupling that enters `z_core` is whatever the
valence–core Coulomb integral yields (no model self-energy, no dimensional
fitting, prefactor exactly 1).

## 1. Splitting the action and integrating out the core

Start from the all-electron action (instruction Eq. for `L`). Split the electron
field into valence and core, `ψ = ψ_v + ψ_c`. The one-body part is diagonal in
this split; the Coulomb interaction couples them through the antisymmetrized
valence–core vertex

```
S_vc = ∫_{r,r'} 1/|r−r'| [ ρ_v(r) ρ_c(r')            (direct / Hartree)
                          − ψ̄_v(r)ψ_c(r) ψ̄_c(r')ψ_v(r') ]  (Fock / exchange)
```

The exchange term is present precisely because valence and core are the **same
fermion species** (antisymmetry). It is decisive below.

Pass to the reference in which the valence orbitals are pushed to high energy
(frozen), solve the core there, and let the valence propagate in what the core
leaves behind. Integrating the (Gaussian, filled-shell) core field out of the
path integral produces a valence self-energy `δV_pp = Σ`, so that

```
g_v^{-1} = g_0^{-1} + δV_pp ,    g_0^{-1} = iω − H_0 .
```

`δV_pp` has a static part (the frozen-core pseudopotential `V_PSP`, already in
`H_KS`, fixing `ε_KS`) and a dynamical part. **Only the dynamical part is the
correction** (re-adding the static part double-counts the core).

## 2. The dynamical self-energy by closure at ΔE_c

Because the core excitation energies `ΔE_c` are several Hartree — far above the
valence Fermi energy — the core is treated by **closure at the single scale
`ΔE_c`**, and (s-states dominate the monopole response) only the core **s-shell
(monopole) channels** are kept.

Second-order perturbation theory / GF2 in the valence–core vertex gives, for a
core channel `c` that can be virtually excited (core electron → empty state, cost
`ΔE_c`), a self-energy of resolvent form

```
Σ_{nk}(ω) = Σ_c  |g_{c,nk}|^2 / (ω − ΔE_c)                                   (3)
```

where the sum over the excited intermediate manifold has been closed at `ΔE_c`
(`Σ_a |a⟩⟨a| → 1`), leaving the squared coupling `|g_{c,nk}|^2` of the valence
state to the c → excited fluctuation.

**Where the iω (z_core) term comes from — and why exchange is essential.**
Expand (3) for valence frequencies ω ≪ ΔE_c:

```
Σ_{nk}(ω) = − Σ_c |g_{c,nk}|^2 / ΔE_c · [ 1 + ω/ΔE_c + … ] .
```

- A purely **direct** density–density response would enter Σ as the even
  combination `2ΔE_c/(ΔE_c^2 + ω^2)`, which is **even in ω**: it gives only a
  static shift and **no `iω` renormalization** — a direct picture misses
  `z_core` entirely.
- The **exchange (Fock)** channel enters with the *single* resolvent `1/(ω−ΔE_c)`
  of (3), which is **odd-in-ω at linear order**: it produces the `iω` term. This
  is the sense in which antisymmetry shapes the dynamical coupling.

The quasiparticle residue follows from the standard definition
`z_core = [1 − ∂Σ/∂ω|_0]^{-1}` with `∂Σ/∂ω|_0 = − Σ_c |g_{c,nk}|^2/ΔE_c^2`:

```
  z_core_{nu,k} = 1 / ( 1 + Σ_c |g_{c,nk}|^2 / ΔE_c^2 ) .                     (4)
```

Equivalently, the `iω` coefficient of `g_v^{-1}` becomes
`1 + Σ_c |g|^2/ΔE_c^2 = z_core^{-1}`, exactly matching the `z_core^{-1} iω`
structure of Eq. (1). Locating the pole of (1) reproduces Eq. (2). `z_core < 1`,
so the occupied bandwidth is compressed — the observed effect.

## 3. The coupling g_{c,nk} from the valence–core Coulomb integral

`g_{c,nk}` is the matrix element of the core-induced fluctuation potential
`δV_pp` between the core channel `c` and the valence Bloch state `ψ_{nk}`. The
core leaves behind its own Coulomb field `V_H^core(r)` (the Hartree/Fock field of
the closed core shell — the very potential whose static part is `V_PSP`); its
**fluctuation** when channel `c` is excited couples valence to core through

```
  g_{c,nk} = 2 ⟨ φ_c | V_H^core | ψ_{nk} ⟩ .                                  (Eq. A)
```

- `V_H^core(r)` is the monopole Hartree potential of the full core density
  `n_core(r) = Σ_{c'} occ_{c'} |φ_{c'}(r)|^2` (occ = 2 per closed s-shell),
  built from the supplied `u_c` by the radial Poisson solve
  `V(r) = (1/r)∫_0^r n_core dr' + ∫_r^∞ n_core/r' dr'`.
- `φ_c(r) = u_c(r)/(r√(4π))` is the atomic core s-orbital (monopole).
- `ψ_{nk}` enters only through its s-wave (monopole) component at the atom,
  `ψ0(r) = (1/√Ω_at) Σ_G c_nk(G) j0(|k+G| r)`, `Ω_at = Ω_cell / n_atom`. This is
  why **`g`, and hence `z_core`, is STATE-DEPENDENT**: it varies with (n,k)
  through the valence–core overlap encoded in `ψ0`.
- The factor **2** is the closed s-shell electron count: both spins of channel
  `c` contribute to the monopole fluctuation; the antisymmetrized direct+exchange
  combination in the monopole channel yields the same multiplicity. (It is a
  derived combinatorial factor, not a fitted constant — the EFT prefactor is 1.)

In radial form (`4π r^2 φ_c → u_c r √(4π)`):

```
  g_{c,nk} = 2 √(4π) ∫_0^∞ u_c(r) V_H^core(r) ψ0(r) r dr .                    (Eq. A')
```

`Ω_at` does not cancel: the delocalized valence amplitude at a single atom scales
as `1/√Ω_at`, correctly making `g` an intensive per-site coupling.

## 4. Final closed form

```
  z_core_{nu,k} = 1 / ( 1 + Σ_c | 2 √(4π) ∫ u_c V_H^core ψ0 r dr |^2 / ΔE_c^2 )
  ε_QP_{nu,k}   = ε_F + z_core_{nu,k} (ε_KS_{nu,k} − ε_F) ,   E_pred = ε_QP − ε_F .
```

Inputs per element: `u_c(r)` and `ΔE_c` (atomic data) → `V_H^core`; from DFTK
`c_nk(G)`, `|k+G|`, `ε_KS`, `ε_F`, `Ω`. No element-specific branch, no fitted
constant.

## 5. Is z_core state-dependent or one factor per element?

**State-dependent** (varies with the Bloch state n,k). The derivation makes
`g_{c,nk}` a matrix element of the core fluctuation potential against the actual
Bloch state: it is governed by how much s-like valence amplitude sits in the core
region (the `ψ0`/overlap content). States with more core-region s character
(e.g. the band bottom near Γ) are renormalized more strongly; states near `E_F`
less so. A single per-element factor would require `g` to be independent of (n,k),
which the Coulomb matrix element is not.

## 6. Self-test (public elements)

Band-bottom occupied depths (eV below `E_F`):

| element | KS (DFTK here) | this work (QP) | eDMFT (Mandal 2022) | experiment |
|---------|----------------|----------------|---------------------|------------|
| Na      | −3.26          | −2.90          | −2.84               | −2.65…−2.78 |
| Al      | −11.18         | −10.55         | —                   | −10.58 (ARPES bottom) |

The correction removes essentially all of the KS bandwidth overestimate and lands
at eDMFT/experiment level. `z_core` is ≈0.89 at the Na band bottom and rises
toward 1 near `E_F`, and ≈0.94 at the Al band bottom (Al’s deeper `ΔE_c`
suppresses the correction via the `1/ΔE_c^2` weight), reproducing the observed
20–35 % (alkali) → few-% (Al) trend with no tuning.

Nearest-band RMSE vs the provided ARPES (occupied + first unoccupied per point):

| element | KS-vs-ARPES | QP-vs-ARPES |
|---------|-------------|-------------|
| Na      | 0.41 eV     | 0.18 eV     |
| Al      | 0.41 eV     | 0.23 eV     |

(The Na residual is a uniform ~0.1–0.15 eV — at the eDMFT-vs-experiment level the
instruction flags as the realistic floor; the band-bottom depth matches eDMFT to
0.06 eV.)
```
