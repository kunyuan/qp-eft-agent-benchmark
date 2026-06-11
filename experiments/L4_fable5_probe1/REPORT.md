# L4 probe 1 — fable5: the direct channel, fully audited, to the doorstep of exchange

First probe of **Level 4** (open frontier), Claude **Fable 5**, 2026-06-10/11.
Host-run relay protocol: 11 sessions over ~18 h wall clock, ~7 productive
(~9 h of working time), the rest zero-progress harness kills during
platform-degradation windows. **User-paused before deliverables were complete**
— archived as an *interrupted* attempt; its scientific trajectory and the
protocol lessons it forced are the value.

Artifacts: [`NOTES.md`](NOTES.md) (the agent's own anti-interruption relay
log — the authoritative timeline), [`derivation_notes.md`](derivation_notes.md),
[`solution/`](solution/) (pipeline as of pause). No evaluator verdict (not
submitted).

## TL;DR

- The agent derived and **exhaustively audited the direct (core-polarization)
  channel**: retarded kernel from the core's density–density fluctuation
  propagator, Fan–Migdal on-shell self-energy, multipole 3j decomposition,
  Fermi-sea blocking. Audits: plane-wave-limit identity vs its own v1 (3%),
  TRK f-sums (8.06/10 Na, 8.23/10 Al), self-computed core polarizabilities,
  Lmax/box/emax convergence scans, an exact-completeness rearrangement of the
  closure counterterm, and a quantified cross-exchange bound (~−20% of D²,
  <0.01 eV narrowing). It found and fixed its own Wigner–Eckart prefactor bug
  via an audit identity.
- **Honest ceiling, self-measured:** best derived narrowing ~0.28 eV of the
  ~0.65 eV Na deficit (v5 Na RMSE 0.2444 vs baseline 0.078). It re-derived the
  double-counting boundary for the GTH reference (plain-LDA-atom fit contains
  no adiabatic CPP ⇒ drop the +1/Δ counterterm) — a bookkeeping claim worth
  adjudicating against the paper.
- **The Li-first revision (PR #18, landed mid-campaign) detonated on contact:**
  one session later the agent measured λ_Γ(Li) = 0.028 vs the 0.338 anchor —
  *"factor ~12 missing. CONFIRMS SETUP warning: contraction missing"* — then
  bounded its entire channel class (λ_max ≈ 0.01–0.03 on Li), ruled out
  correlation-type mechanisms with an atomic ground truth (total CV correlation
  −0.05 eV), and identified *"CV exchange K(1s,2s) = 0.6 eV — the only
  eV-scale state-dependent CV quantity. Investigating."* **Paused at the
  doorstep of the exchange channel**, reached by its own scale accounting.
- Audit watchlist (unresolved): its variance bound claims to close
  "density/exchange contractions of the H₁c² class" — likely carries the 1/Ω
  contraction-convention suppression (the recurring trap); would have needed
  adjudication had it proceeded.

## Physics question raised for the paper (arXiv 2604.25199)

Its audited direct-channel slope λ ≈ 0.10 (Na) vs the paper's Γ₃ suppression
estimate ε_F/ΔE_c ~ 0.02–0.05: overlap with the hole channel (its unproven
"D²+X² already counted" claim), overestimate, or a real entry in the paper's
residual budget (Al 0.248)? The agent's own framework later argued the GTH
reference contains no adiabatic CPP — which sharpens, not settles, the
question. Worth a maintainer-side calculation in the paper's framework.

## Protocol lessons (all encoded in `.claude/memory/` and PR #18)

1. **600 s stream watchdog vs long thinking**: 6+ kills, always at the deepest
   derivation moments; "short bursts + externalize every step" (sessions 7+)
   was the fix and doubled as audit-quality fuel.
2. **NOTES.md relay**: added at session 4, made kills minute-level losses;
   now part of the standing protocol (pattern-001).
3. **Resume-note contamination**: sessions 3–4 resume prompts contained
   steering phrases ("what other contraction", "identical fermions /
   antisymmetry"); disclosed here; the fixed-template resume (point only at
   the agent's own logs) was adopted from probe 2 onward.
4. **Ledger blindness to missing channels** (pitfall-003): every ledger item
   passed while the dominant channel was absent — the Li-first stage is the
   structural fix and was validated within one session of landing.

## Relation to probe 2

Probe 2 ([`../L4_fable5_probe2/`](../L4_fable5_probe2/)) is the clean A/B
counterpart: fresh agent, Li-first protocol from minute one, no inherited
artifacts, pure resume templates. It reached this probe's 20-hour conclusions
in 31 minutes and went further. Read the two NOTES.md side by side.
