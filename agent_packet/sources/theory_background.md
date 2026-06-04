# Theory Background

## The Problem

Kohn-Sham density functional theory gives a ground-state density and a static
auxiliary Hamiltonian. Its eigenvalues are often compared with ARPES, but DFT
does not guarantee that bare Kohn-Sham eigenvalues are quasiparticle energies.

For simple metals, this becomes an observable problem: the occupied bandwidth
from a static frozen-core Kohn-Sham calculation can be too large. A correction
is needed that is not another fitted exchange-correlation functional.

## Physical Mechanism

A frozen-core pseudopotential removes core electrons from the explicit
valence calculation. The static part of the core-valence physics is already
encoded in the pseudopotential. The missing piece is dynamic: virtual core
excitations introduce a frequency-dependent self-energy.

Two approximations make the leading correction controlled:

1. Core excitation energies are much larger than the valence Fermi energy.
   This lets core excitations be integrated out as high-energy modes.
2. For metallic valence electrons, the uniform electron gas has an approximate
   Galilean-invariant quasiparticle structure over the occupied Fermi ball.
   This makes the static tree-level Hamiltonian coincide with the Kohn-Sham
   Hamiltonian, while the frozen-core dynamics rescales the pole.

The resulting leading-order prediction is post-SCF: first compute a standard
Kohn-Sham band structure, then compress occupied Kohn-Sham energies toward the
Fermi level by a frozen-core factor.

## What Counts As A Valid Theory Here

A valid method must be parameter-free. It can use:

- crystal structure and lattice constants from `element_config.json`;
- standard DFTK/LDA/GTH settings from `element_config.json`;
- the frozen-core model inputs from `element_config.json`;
- formulas in this packet.

It must not use:

- fitted bandwidth scaling factors;
- element-specific knobs chosen to match ARPES;
- hidden experimental data;
- manual nudging of output numbers.

## Implementation Target

Implement a generic runner that accepts an element config and a grid, computes
or approximates the Kohn-Sham occupied bands, applies the frozen-core
quasiparticle correction, and writes all occupied quasiparticle energies at
each requested `point_id`.

The hidden evaluator will call the same runner on other simple metals. The
code should therefore depend on input files and formulas, not on element names.

