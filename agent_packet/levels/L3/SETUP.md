# Physical setup (Level 3)

## The puzzle and the goal

For the simple metals, Kohn-Sham LDA overestimates the occupied bandwidth seen in
ARPES by 20-35% (alkali) down to a few % (Al). Ordinary DFT functionals (LDA, GGA,
hybrids, G0W0) do not fix it; only a full many-body treatment (eDMFT) reproduces the
measured bandwidths. The missing physics is therefore a DYNAMICAL effect of the core
electrons. Your goal: DERIVE, from first principles, the parameter-free correction
that captures it — reaching eDMFT / experiment-level accuracy — and predict the
band energies.

## The quasiparticle propagator (the result you are reproducing)

This effective theory rests on (i) the scale separation between the core excitation
energies and the valence Fermi energy (which lets the core be integrated out), and
(ii) the locality of the interacting uniform electron gas at metallic densities
(established by high-order diagrammatic Monte Carlo), which dictates how many-body
effects are absorbed into renormalized parameters. The tree-level quasiparticle
propagator of the EFT is

    G_{nu,k}^{-1}(i omega)  ~=  (z_val)^{-1} [ (z_core_{nu,k})^{-1} i omega
                                              - ( eps_KS_{nu,k} - eps_F ) ]      (1)

where the dispersion is governed by the eigenvalues `eps_KS_{nu,k}` of the
Kohn-Sham Hamiltonian

    H_KS  =  -1/2 nabla^2  +  V_PSP  +  V_H[n]  +  V_xc[n]      (V_PSP: GTH large-core),

`z_val` is the valence quasiparticle residue (an overall prefactor that does NOT
shift the pole), and `z_core_{nu,k}` is the frozen-core renormalization, which
rescales the frequency relative to the dispersion — an effective time dilation from
virtual core excitations. Locating the pole gives the energy ARPES measures,

    eps_QP_{nu,k}  ~=  eps_F  +  z_core_{nu,k} ( eps_KS_{nu,k} - eps_F ).      (2)

Your DFTK run gives `eps_KS_{nu,k}` and `eps_F`. What is NOT given is
`z_core_{nu,k}` — deriving its closed form is the task.

## How to derive the effective valence action

The rigorous starting point is the ALL-ELECTRON action — valence AND core electrons
in the lattice potential, with the full Coulomb interaction,

    L = int_{r,sigma} bar_psi_{r,sigma} [ d_tau - nabla^2/2 + V_Lat - mu ] psi_{r,sigma}
        + (1/2) int_{r,r'} ( bar_psi_{r,sigma} bar_psi_{r',sigma'}
                              psi_{r',sigma'} psi_{r,sigma} ) / |r - r'| .

INTEGRATE THE CORE FIELDS OUT of the path integral to obtain an effective action for
the valence electrons alone. (A clean route: pass to a reference in which the valence
orbitals are frozen — pushed to high energy by a hybridization — so the core is
solved in that reference and the valence then propagates in the potential the core
leaves behind.) The result is a valence propagator

    g_v^{-1}  =  g_0^{-1}  +  delta_V_pp ,        g_0^{-1} = i omega - H_0,
                                                  H_0 = -nabla^2/2 + V_Lat - mu,

where `delta_V_pp` is the pseudopotential induced by the core. It has
- a STATIC (energy-independent) part: the conventional frozen-core pseudopotential
  `V_PSP` already inside `H_KS` — this sets the dispersion `eps_KS` in Eq. (1);
- a DYNAMICAL (frequency-dependent) part: this is the `z_core` renormalization in
  Eq. (1), and it is the leading thing static pseudopotentials omit.

Re-adding the static part would double-count the core; the correction is the
dynamical part only. Derive it — and hence `z_core_{nu,k}` and, through Eq. (2), the
band energies (relative to E_F) and the bandwidth.

Do this as a real derivation, not a model: start from the many-electron action and
integrate the core fields out STEP BY STEP to obtain `δV_pp` — and in particular its
dynamical (frequency-dependent) term — from the valence–core Coulomb interaction
itself. The coupling that enters `z_core` is whatever this integral yields; do not
posit a model self-energy or fix the coupling by dimensional analysis.

Be careful with EXCHANGE: valence and core electrons are antisymmetric, so `δV_pp`
carries a Fock/exchange contribution alongside the direct (Hartree) one, and the
antisymmetry shapes the dynamical coupling — a purely direct density-response
picture misses it.

One simplification is controlled: the core excitation energies `DeltaE_c` are several
Hartree, far above the valence Fermi energy, so treat the core excitations by closure
at the single scale `DeltaE_c`, and — since they dominate — keep only the core
s-shell (monopole) channels.

Do not posit the form of `z_core_{nu,k}`: what it is built from out of the core data,
and how it depends on the Bloch state, are for you to derive.

## What you have

Per core s-channel `c`: the all-electron radial core orbital `u_c(r)`
(`atomic_core_<c>.csv`, normalized int u_c^2 dr = 1) and its excitation energy
`DeltaE_c` (`core_model.json`). `DeltaE_c` is the full many-electron excitation
energy of the interacting core — the total-energy difference
E_core(hole in c) − E_core(ground state) — NOT a bare orbital-eigenvalue
difference; keep your coupling and your energy denominators consistent with this
definition. From your DFTK run: the KS eigenvalues, E_F, and the
plane-wave coefficients `c_nk(G)` of each Bloch state (`|k+G|` in Bohr^-1). Anything
else your derivation needs — potentials, integrals — you build from these.

Parameter-free: no fitted constant (the correct prefactor is exactly 1). If you need
to tune a constant to match the public elements, your derivation is incomplete —
re-derive, do not fit. The same code path runs on every element; it is graded on
held-out metals whose core structure differs from the public Na/Al, so matching
Na/Al is necessary but not sufficient.

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
