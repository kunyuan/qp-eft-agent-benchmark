# Formula Sheet

This sheet summarizes the formulas needed for the benchmark. Use them as a
theory-to-code guide, not as adjustable fitting rules.

## Kohn-Sham To Quasiparticle Pole

Let `E_KS(n,k)` be a Kohn-Sham eigenvalue and `E_F` the Fermi level. The
leading frozen-core quasiparticle prediction is:

```text
E_QP(n,k) = E_F + z_core(n,k) * (E_KS(n,k) - E_F)
```

Occupied states have `E_QP < E_F`. The benchmark CSV expects energies relative
to the Fermi level:

```text
E_pred_eV = E_QP(n,k) - E_F
          = z_core(n,k) * (E_KS(n,k) - E_F)
```

## Frozen-Core Factor

The state-dependent expression is:

```text
z_core(n,k) = 1 / (1 + sum_c |F_c(n,k)|^2 / DeltaE_c^2)
```

where `c` labels core excitation channels and `DeltaE_c` is the core excitation
energy in Hartree.

The coherent form factor is:

```text
F_c(n,k) = sum_G c_nk(G) * f_c(|k + G|)
```

where `c_nk(G)` are plane-wave coefficients of the Kohn-Sham Bloch state and
`f_c(K)` is the radial core form factor.

## Warmup Scalar Model

The public `element_config.json` files include `core_model.z_core_gamma`. You
may use this for the scalar warmup task:

```text
E_pred_eV = z_core_gamma * E_KS_relative_eV
```

This warmup is not the full theory. A stronger implementation should compute
state-dependent `z_core(n,k)` from the available core model inputs and
Kohn-Sham wavefunction information.

## Double Counting

Do not add the static value of the core self-energy to the Kohn-Sham
Hamiltonian again. The static core effects are already represented by the
frozen-core pseudopotential and exchange-correlation potential. The correction
used here is the frequency-dependent change around the Fermi level, which
appears as the `z_core` factor.

## Multi-Band Output

For each requested k-point, output every occupied quasiparticle band. The
evaluator matches each measured ARPES point to the nearest predicted band at
the same `point_id`.

