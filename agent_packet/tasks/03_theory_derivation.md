# Task 03: Derive The Correction

Goal: write the derivation that justifies the correction you will implement.

Use:

- `sources/theory_background.md`
- `sources/formula_sheet.md`
- public Na/Al diagnostics in `element_config.json`

Your `method.md` must explain:

1. why the missing effect is dynamic frozen-core physics rather than a static
   exchange-correlation or pseudopotential refit;
2. why the static core self-energy should not be added again;
3. why the leading quasiparticle energy has the form
   `E_QP - E_F = z_core * (E_KS - E_F)`;
4. how your code obtains `z_core` from the input files;
5. why the code should generalize to additional simple metals without
   per-element fitting.

The derivation can be concise, but it must be specific enough that another
engineer can map every equation to a code path.

