# Physical setup (Level 3)

For these simple metals the Kohn-Sham occupied bandwidth overestimates ARPES by
20-35% (alkali) down to a few % (Al). DFT gives no reason KS eigenvalues should
be quasiparticle energies — yet they nearly are. Your task is to find, derive,
and implement the leading parameter-free correction, then predict band energies.

What is missing is dynamical. A frozen-core (large-core) pseudopotential removes
the core electrons; its STATIC core-valence physics is already included. The
omitted piece is the dynamical response of virtual core excitations. Two facts
make the leading correction controlled:

1. Scale separation: core excitation energies `DeltaE_c` (provided per core
   s-channel in `core_model.json`) are several Hartree, far above the valence
   Fermi energy — they can be integrated out as high-energy modes.
2. The interacting uniform electron gas at metallic density is approximately
   Galilean invariant over the occupied Fermi ball, so the static tree-level
   Hamiltonian coincides with the Kohn-Sham Hamiltonian; only the quasiparticle
   pole is rescaled by the frozen-core dynamics.

You are given, per core s-channel c: the all-electron radial core orbital
`u_c(r)` and its single-orbital Hartree potential `V_H_c(r)`
(`atomic_core_<c>.csv`), and the core excitation energy `DeltaE_c`. From your
DFTK calculation you have the KS eigenvalues and the plane-wave coefficients of
each Bloch state. Derive the leading post-SCF quasiparticle correction and
implement it. It must be parameter-free and use the same code path for every
element.
