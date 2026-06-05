# Physical setup (Level 3)

## The puzzle

For decades, angle-resolved photoemission (ARPES) on the simple alkali metals
posed a puzzle: the measured occupied bandwidth is 20-35% NARROWER than the
free-electron / local-density (LDA) prediction (only a few % for Al). The
discrepancy resisted the usual fixes — a comprehensive study (Mandal et al.) found
that G0W0, hybrid functionals, and meta-GGA all fail to reproduce the ARPES
bandwidths; only embedded dynamical mean-field theory (LDA+eDMFT), which keeps the
full frequency dependence of the self-energy, achieves quantitative agreement. The
missing physics is therefore DYNAMICAL — and, as the rest of this setup makes
precise, it lives in the core.

Your task is to DERIVE, from first principles, the leading parameter-free
correction that produces this narrowing, then implement it and predict the band
energies.

## What is already understood: the valence quasiparticle (given)

The valence side is not in doubt. For an inhomogeneous electron gas — electrons
interacting through their Coulomb repulsion in a fixed external one-body potential
— accurate many-body calculation confirms, as a high-precision result you may take
as given, that in the Fermi region (occupied states, at and below the Fermi level)
the single-particle Green's function is a sharp coherent quasiparticle,

    G(iw_n)  ~=  z / (i w_n  -  H_KS) ,

with the dispersion set by an effective single-particle (Kohn-Sham) Hamiltonian
`H_KS` and a quasiparticle weight `z`. (This is why ARPES sees sharp peaks and why
KS eigenvalues are nearly the quasiparticle energies.)

## Where the puzzle lives: the core (the catch)

A real metal is NOT an inhomogeneous electron gas — besides the valence electrons
it has CORE electrons, and the rigorous starting point is the ALL-ELECTRON
Hamiltonian. A conventional (frozen-core) pseudopotential integrates the core out
STATICALLY: it keeps the static core-valence physics but discards the frequency
dependence — exactly the dynamical effect eDMFT showed is needed. (Bolting a
separate many-body correction onto a pseudopotential instead double-counts the
core, which was already integrated out implicitly — so derive it from one
consistent reduction.)

So treat the core correctly. Starting from the all-electron Hamiltonian, integrate
the core electrons out. The reduction has

- a STATIC part = the conventional pseudopotential your DFTK run already uses (so
  `H_KS` and the bands you compute ARE that static reduction); and
- a DYNAMICAL part = the core's frequency-dependent response — a valence electron
  virtually excites a core particle-hole pair (energy `DeltaE_c`) and the core
  relaxes back — absent from any static pseudopotential.

## The task

Derive that dynamical part from first principles — it is the leading
(second-order) core-induced contribution, with the core response controlled by
`DeltaE_c` — then derive its consequence for the quasiparticle band energies
(relative to the Fermi level). Do NOT assume the mechanism: whether the dynamical
term renormalizes the weight `z`, shifts the energies, or something else; how it
combines with the coherent `G` above; how large the effect is; and whether it acts
identically on every Bloch state — are all for you to work out from the
derivation, not to posit. One simplification is controlled: `DeltaE_c` is several
Hartree, far above the valence Fermi energy, so the core excitations are
high-energy modes you may integrate out by closure at the single scale `DeltaE_c`.

## What you have

Per core s-channel `c`: the all-electron radial core orbital `u_c(r)` and its
single-orbital Hartree potential `V_H_c(r)` (`atomic_core_<c>.csv`), and the core
excitation energy `DeltaE_c` (`core_model.json`). From your DFTK calculation: the
KS eigenvalues, the Fermi level, and the plane-wave coefficients `c_nk(G)` of
each Bloch state (`|k+G|` is the Cartesian length in Bohr^-1).

Your correction must be parameter-free (no fitted constant, no fit to ARPES) and
use the same code path for every element. It is graded on held-out metals whose
core structure differs from the public Na/Al — matching Na/Al is necessary but
NOT sufficient, so a form tuned to those two will not generalize.
