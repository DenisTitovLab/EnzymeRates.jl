# SUPERSEDED — EqualAI × NonequalAI Coupling Briefing

> **DO NOT IMPLEMENT THIS DOCUMENT.** It is retained only as a breadcrumb.

## Why it was superseded

This briefing was written on a **false premise**. It claimed that mixed
`:EqualAI` × `:NonequalAI` catalytic configurations are physically inconsistent
and should be **rejected** at `AllostericMechanism` construction (via a
rank/nullspace validator), with the enumeration rewritten to never produce them.

Reading the actual hand-verified mechanisms (`PK`, `m_mixed`) showed that is
wrong: those configurations are **valid**. A single `:NonequalAI` binding group
with `:EqualAI` catalysis is fine — the catalytic reverse rate is a *dependent*
parameter that legitimately differs between the A and I states (PK expects
`n_haldane = 2`; `m_mixed` expects rate 0 at equilibrium). Rejecting them would
break correct mechanisms.

The **real** bug is much smaller and lives entirely in the synth-dep machinery:
an `:EqualAI` dependent param whose Haldane/Wegscheider RHS references a
`:NonequalAI` symbol self-maps in `rename_T` and overwrites its active-state
assignment.

## What to read / implement instead

- **Contained fix (implement this now):**
  - Spec: `docs/superpowers/specs/2026-05-29-equalai-nonequalai-coupling-design.md`
  - Plan: `docs/superpowers/plans/2026-05-29-equalai-coupling-impl.md`
- **Follow-up PRs (after the parent refactor):**
  - `docs/superpowers/specs/2026-05-29-direction-symmetry-constraint-resolution.md`
    — direction-invariant dependent-parameter removal for all mechanisms.
  - `docs/superpowers/specs/2026-05-29-nonequalai-rank-validity.md`
    — the rank/nullspace algorithm for NonequalAI degeneracy (the *correct* home
    for the rank idea this briefing mis-scoped as a rejection-everything validator;
    it rejects only genuinely-degenerate configs, e.g. a lone NonequalAI in a
    pure-RE Wegscheider loop).
