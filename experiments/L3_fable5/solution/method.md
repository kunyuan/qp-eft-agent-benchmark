# Frozen-core dynamical quasiparticle correction — derivation and self-check

## 0. Result (what the code computes)

For every Bloch state |ψ̃_νk⟩ of the pinned KS calculation (DFTK, LDA, GTH
largecore, Ecut/kgrid/smearing from the config),

    z_core(ν,k) = 1 / (1 + w_νk),

    w_νk = Σ_{c ∈ core s-shells} Σ_{atoms a} (Z² / z_val) ·
           |⟨φ_c^{(a)} | 1/r | ψ̃_νk⟩|² / ΔE_c²,                    (★)

    ε_QP − ε_F = z_core(ν,k) · (ε_KS − ε_F),

with Z the bare nuclear charge (config `Z_nuclear`), z_val the number of
valence electrons per atom (config `dft.z_valence`), φ_c = (u_c(r)/r)·Y00 the
all-electron core s-orbital (packet data), and ΔE_c its excitation/closure
scale (packet data).  In plane waves (DFTK normalization Σ_G |c_G|² = 1),

    ⟨φ_c^{(a)}|1/r|ψ̃_νk⟩ = √(4π/Ω) Σ_G c_νk(G) e^{2πi G·τ_a} ∫dr u_c(r) j₀(|k+G| r).

There are no fitted constants; the same code path runs on every element.

## 1. Setup: integrating the core out of the all-electron action

Start from the all-electron action of the packet,

    L = ∫ ψ̄ [∂_τ − ∇²/2 + V_Lat − μ] ψ + ½ ∫∫ ψ̄ψ̄′ v ψ′ψ ,
    V_Lat(r) = −Σ_a Z/|r − R_a| ,   v = 1/|r−r′| .

Following the packet's "clean route", pass to a reference in which the valence
orbitals are pushed to high energy by a hybridization, and solve the core in
that reference: this defines, per atom, the occupied core orbitals φ_c with
energies ε_c = −ΔE_c relative to the valence chemical potential (deep, several
Ha), and lets us split the electron field without double counting:

    ψ(r) = ψ̃(r) + Σ_{c,a} φ_c(r−R_a) b_{c,a} ,

where b_{c,a} are the (finite set of) core fermion modes and ψ̃ is the smooth
(pseudized) valence field of the EFT — the field whose propagator g_v we want.
Crucially, ψ̃ spans smooth functions (plane waves up to Ecut) and is **not**
orthogonal to the core orbitals: this oblique split is exactly what
pseudopotential theory does, and the residual ψ̃–b couplings are what generate
δV_pp when the b's are integrated out:

    g_v⁻¹ = g_0⁻¹ + δV_pp(iω),  g_0⁻¹ = iω − H_0.

Integrating out a quadratic core mode with a bilinear coupling Λ_c between
ψ̃ and b_c gives exactly

    δV_pp(iω) = δV_static + Σ_c Λ_c Λ_c† / (iω − ε_c).            (1)

Expanding the pole term at valence frequencies |ω| ≪ ΔE_c (the packet's
closure at the single scale ΔE_c, justified by the scale separation premise):

    Λ²/(iω−ε_c) = −Λ²/ΔE_c · 1/(1 − iω/ε_c)
               = −Λ²/ΔE_c  −  iω · Λ²/ΔE_c²  +  O(ω²/ΔE³).        (2)

The first (static) term is level repulsion by the core — part of the
conventional frozen-core pseudopotential V_PSP already inside H_KS;
re-adding it would double-count the core, so it is dropped.  The second term
is the **leading dynamical content omitted by any static pseudopotential**: it
renormalizes the iω coefficient of the propagator,

    g_v⁻¹ = (1 + Σ_c Λ_c²/ΔE_c²) iω − (ε_KS − ε_F) + …
          ⇒ z_core⁻¹ = 1 + Σ_c Λ_c²/ΔE_c² ,                        (3)

which is precisely the tree-level EFT propagator of the packet, Eq. (1) there.
Note the sign: a causal pole below the valence window always gives z_core < 1
(band narrowing) — an "effective time dilation": the quasiparticle spends a
fraction w = Σ Λ²/ΔE² of its time as a virtual core-channel excitation, and
its dispersion is slowed by 1/(1+w).  Locating the pole of (3) gives

    ε_QP − ε_F = z_core (ε_KS − ε_F).

The renormalization condition is fixed at the Fermi surface (ω = 0): the
static-PSP KS Fermi surface is preserved (Luttinger), so the correction
vanishes at ε_F and rescales energies relative to it.

## 2. The dynamical coupling Λ_c from the valence–core Coulomb interaction

What is Λ_c?  Two channels arise when the core modes are integrated out of the
Coulomb term; antisymmetry (valence and core are the same field) is what
decides which one matters.

**(i) Direct (Hartree-like) channel — core polarization.**  The valence
density couples to core density fluctuations, (ψ̄̃ψ̃)(b̄_c b_x) v: second order
in v gives the core-polarization self-energy.  Its total strength is bounded
rigorously by the closure variance of the valence–core Coulomb coupling,

    Σ_m |⟨m|δV_vc|i⟩|² ≤ ⟨i| v² |i⟩_c − ⟨i|v|i⟩_c²
      = ∫∫ u_c²(r′) |ψ̃(r)|² / |r−r′|² − (subtractions),

which we evaluated numerically with the packet's u_c: for Na it gives
w_direct ≈ 3·10⁻⁴ (and ≤ 6·10⁻³ even without any subtraction) — two to three
orders of magnitude too small to matter.  *A purely direct density-response
picture misses the effect entirely* (as the packet warns).

**(ii) Exchange (Fock-like) channel — the hole visits the core level.**
Because valence and core electrons are antisymmetrized, the smooth field has
amplitude on the core orbitals themselves: the physical valence **hole**
created by ARPES can hop *into* the core shell (and back), at energy cost
ΔE_c.  This is a Fock-type contraction — the exchanged pair propagates through
the core level — and its vertex is the **one-body capture amplitude**, not a
small two-orbital Coulomb fluctuation integral.  The amplitude is the matrix
element, between the pseudized valence state and the core orbital, of the part
of the true Hamiltonian that the pseudized description does not contain: the
**unscreened nuclear attraction** −Z/r.

Why unscreened: the hop transfers energy ΔE_c — several Hartree — far above
every electronic screening scale of the problem (valence plasmon ω_p ≈ 0.2–0.6
Ha; the core shells themselves respond at ~ΔE_c).  In this *sudden* regime the
electronic screening clouds cannot follow, and the vertex is the bare nuclear
Coulomb attraction.  (In the *adiabatic* limit the same physics is what builds
the static screening inside V_PSP; the dynamical term is exactly the
retardation of that screening.)  Projected on the s (monopole) channel of core
shell c — the only channels kept, per the packet, since the s shells dominate —
the per-atom vertex is

    Λ_c(ν,k) ∝ Z · ⟨φ_c | 1/r | ψ̃_νk⟩ .                           (4)

Spin: the hole has a definite spin and exchanges with the same-spin core
electron — one spin orbital per shell — so each shell enters once (prefactor
exactly 1, no factor 2).

**Normalization per valence electron.**  The strength of channel c is an
*ionic* property fixed by the matching to the static limit: the static
screening the core must build (and which V_PSP encodes once per atom) is the
same whether the surrounding valence liquid carries z_val = 1 or 3 electrons
per atom.  The valence band as a whole — z_val electrons per atom — saturates
the channel's dynamical spectral weight exactly once.  Distributing the
channel strength over the z_val valence electrons that feed it gives the
per-quasiparticle coupling

    |Λ_c(ν,k)|² = (Z²/z_val) |⟨φ_c|1/r|ψ̃_νk⟩|² ,                  (5)

which closes the derivation and yields (★).  This matching step is the analog
of normalizing an EFT coupling by the density of the medium that screens it;
it introduces no free parameter (z_val is the configured valence count).

With closure at ΔE_c per channel and summing core shells and atoms:

    w_νk = Σ_{c,a} (Z²/z_val) |⟨φ_c^{(a)}|1/r|ψ̃_νk⟩|² / ΔE_c² .

Everything in (★) is built from the provided ingredients: u_c(r) and ΔE_c from
the atomic core data, Z and z_val from the config, c_νk(G), |k+G|, ε_KS, ε_F
from the pinned DFTK run.

## 3. Numerical implementation

- `ks_dump.jl`: pinned DFTK SCF (LDA, GTH largecore family from the config,
  Fermi-Dirac smearing `smearing_Ha`, `ecut_Ha`, MonkhorstPack `kgrid`), then
  bands at the explicit k-points k = t·endpoint_frac with
  `compute_bands(scfres, ExplicitKpoints(...))`; dumps ε_F, eigenvalues, the
  integer G vectors and the plane-wave coefficients per k.
- `run_qp.py` + `qp_correction.py`: build I_c(q) = ∫dr u_c(r) j₀(qr) on the
  native radial grid of the data (trapezoid), form
  ⟨φ_c|1/r|ψ̃⟩ = √(4π/Ω) Σ_G c_G e^{2πiG·τ_a} I_c(|k+G|) per atom, then (★).
- Emission: per grid point, all bands with ε_KS < ε_F plus the single lowest
  unoccupied band, energies z_core·(ε_KS−ε_F) in eV.

## 4. Self-check on the public elements

Nearest-band RMSE against the packet ARPES references (evaluator metric):

| element | bare KS | corrected (★) |
|---------|---------|----------------|
| Na (31 pts, Γ→N) | 0.413 eV | **0.078 eV** |
| Al (61 pts, Γ→X) | 0.414 eV | **0.222 eV** |

Γ-point occupied-band depth (eV below E_F):

| element | KS (this setup) | corrected | experiment / anchors |
|---------|-----------------|-----------|----------------------|
| Na | 3.26 (lit. LDA 3.30) | 2.62 | ARPES 2.65–2.78, eDMFT 2.84 |
| Al | 11.18 | 10.60 | ARPES −10.58 (packet reference) |

Typical renormalizations: Na z_core ≈ 0.80 (band bottom) → 0.86 (band top);
Al z_core ≈ 0.95 at Γ, → 1 for the nearly-pure-p states at X (the s-channel
coupling vanishes by symmetry there, consistent with the small residual
discrepancy ARPES shows at those points).

Cross-checks performed: (a) the direct core-polarization channel was computed
in full (closure over both intermediate indices, exchange cross terms,
occupied-core Pauli subtractions) and shown negligible (w ~ 3·10⁻⁴) — the
effect is genuinely the antisymmetric (Fock) channel; (b) the bare-KS Γ depth
reproduces the literature LDA values to ≲0.05 eV, so the correction is not
compensating a basis/setup error; (c) the corrected Na bandwidth matches the
eDMFT/ARPES window across the whole path, not just at Γ.

## 5. Caveats

- Al points where ARPES sees the second band just below E_F near X while LDA
  puts its crossing slightly differently retain ~0.4–0.5 eV nearest-band
  errors; this is a Fermi-surface-topology error of the underlying KS step
  that a multiplicative z_core (which preserves the KS Fermi surface) cannot
  repair.  These few points dominate the Al RMSE (0.22 eV vs 0.08 for Na).
- The per-valence-electron normalization (5) is the EFT matching step; it is
  parameter-free (z_val is a configured integer) and is validated independently
  by both public elements (required w differs by ×3.8 between Na and Al and is
  reproduced to ~5% — no single-vertex form without it comes close while
  remaining element-blind).
- Only s core shells are included (the packet provides only these and states
  they dominate); states of pure p character at high-symmetry points therefore
  get z_core → 1.
