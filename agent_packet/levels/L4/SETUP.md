# Physical setup (Level 4 — open frontier)

## The puzzle and the goal

For the simple metals, Kohn-Sham LDA overestimates the occupied bandwidth seen in
ARPES by 20-35% (alkali) down to a few % (Al). Ordinary DFT functionals (LDA, GGA,
hybrids, G0W0) do not fix it; only a full many-body treatment (eDMFT) reproduces the
measured bandwidths. The missing physics is a DYNAMICAL effect of the core
electrons.

This level is the OPEN rung: nothing about the answer's structure is given — no
propagator ansatz, no decomposition roadmap, no prescribed approximations, and NO
atomic data. You choose the reference, the channels, and every truncation; each
must be declared, controlled, and justified. A published leading-order frozen-core
treatment of this problem exists and is the calibration baseline: through this
exact pinned pipeline it reaches nearest-band RMSE 0.078 eV (Na) and 0.248 eV (Al)
on the public development elements. Matching that level on the concealed metals
means you reconstructed the leading-order physics; beating it there means physics
beyond the published treatment. Both are in scope — the second is the point of
this level.

## Starting point

The rigorous starting point is the ALL-ELECTRON action — valence AND core electrons
in the lattice potential, with the full Coulomb interaction,

    L = int_{r,sigma} bar_psi_{r,sigma} [ d_tau - nabla^2/2 + V_Lat - mu ] psi_{r,sigma}
        + (1/2) int_{r,r'} ( bar_psi_{r,sigma} bar_psi_{r',sigma'}
                              psi_{r',sigma'} psi_{r,sigma} ) / |r - r'| .

Integrate the core fields out of the path integral to obtain an effective action
for the valence electrons alone, and DERIVE — do not posit — the correction it
implies for the quasiparticle band energies. One bookkeeping fact you must
respect: the pinned KS step already contains a static frozen-core pseudopotential
(`V_PSP`), so the static part of whatever your integrate-out yields is already
counted — your correction must consist only of what static pseudopotentials omit,
and re-adding any static piece double-counts the core.

## What you have — and what you must build yourself

`element_config.json` (lattice, atom positions, `Z_nuclear`, `dft.z_valence`, the
pinned DFT settings), `grid.csv`, and the public ARPES for Na and Al. There is NO
atomic core data in this packet: whatever atomic inputs your derivation needs —
core orbitals, core potentials, excitation energies — you compute yourself (e.g.
with your own radial atomic solver; the core occupation follows from `Z_nuclear`
and `z_valence`). Choose the level of atomic theory and justify it. Any core
excitation energy you use must be defined by an explicit expression in your own
reference — state it, do not treat it as a black box.

## Required consistency ledger (in method.md)

1. An approximation ledger: every truncation (reference choice, channel
   selection, closure/spectral treatment, frozen vs relaxed core orbitals,
   screening) declared, with a controlled error estimate.
2. The explicit expression for every core excitation energy you use.
3. For each core channel: the diagonal matrix element of your derived coupling
   against that channel's orbital, and why its value is consistent with your
   energy bookkeeping.
4. If you claim accuracy beyond the published leading order, attribute the gain
   to the specific physical term that produced it.

## Rules

Parameter-free: no fitted constant. If you need to tune a constant to match the
public elements, your derivation is incomplete — re-derive, do not fit. The same
code path runs on every element; it is graded on held-out metals whose core
structure differs from the public Na/Al, so matching Na/Al is necessary but not
sufficient.

Your derivation must remain valid across core-shell types — shells both tighter
and more diffuse than the public elements' (e.g. n = 3 s-channels of much larger
radial extent): the shell- and Z-dependence of the coupling must come out of the
derivation, not be calibrated on the tight cores of the development set.

Derive first, evaluate second. Do not use public-element magnitudes to SELECT
between candidate couplings: entire families of incorrect vertices are nearly
degenerate on Na/Al, so numerical agreement there can never substitute for the
derivation. Conversely, a correctly derived coupling can appear wrong by a large
factor if your Bloch-state contraction conventions are inconsistent — if the
magnitude looks off, audit the conventions end-to-end against the public ARPES
before abandoning the derivation.
