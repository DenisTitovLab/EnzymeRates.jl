# Phase 2 — Enumerator to decomposed Species end-to-end (design)

**Branch:** `refactor-to-concrete-types-instead-of-symbols`
**Status:** implemented — `init_mechanisms` / `expand_mechanisms` emit decomposed `Mechanism` / `AllostericMechanism` end-to-end; the `_Raw*` working-rep family and all form-name string helpers are deleted; full suite green.
**Date:** 2026-05-26
**Prerequisite reading:**
[`2026-05-20-concrete-types-refactor-design.md`](2026-05-20-concrete-types-refactor-design.md),
[`2026-05-22-concrete-types-refactor-continuation-design.md`](2026-05-22-concrete-types-refactor-continuation-design.md),
[`../plans/2026-05-25-finish-refactor-roadmap.md`](../plans/2026-05-25-finish-refactor-roadmap.md),
and `.claude/CLAUDE.md`.

This is the detailed Phase 2 spec the roadmap deferred ("full task detail
deferred to its own brainstorm once Phase 1 lands"). Phase 1's green gate is
confirmed; this spec is grounded in the empirical diagnostic of §2, not the
roadmap's pre-Phase-1 sketch.

## 1. Goal

Make the mechanism-enumeration pipeline produce **decomposed** `Species`
(`Species([Substrate(:A)], :E)`) instead of the **opaque** form it produces
today (`Species([], :E_A)` — empty bound list, opaque conformation Symbol).
This removes the last parallel representation of enzyme forms and closes the
two outstanding primary success criteria of the continuation spec:

- **§3.1 #1** — one concrete struct family; no opaque-vs-decomposed parallel
  representation remains in `src/`.
- **§3.1 #2** — no Symbol-string dispatch; form names are not parsed back into
  structure (the parse-back helpers are deleted).

Everything else (test integrity, perf gates, chokepoint architecture)
carries forward unchanged from the prior specs.

## 2. Empirical baseline (diagnostic, 2026-05-26)

The roadmap flagged a contested question (one `_expand_*` move vs. the
`_catalytic_topologies` production boundary as the spike target) and a feared
failure (multi-substrate enumerator output crashing with "Cycle N not
proportional", continuation spec §11.1). Both were resolved empirically before
writing this spec.

**Diagnostic 1 — bi-bi derivation.** Built the bi-bi reaction
(`A[C], B[N] = P[C], Q[N]`), ran `init_mechanisms` (77 mechanisms), and
compiled + derived the 30 simplest (by step count, including the 8-parameter
multi-substrate central-complex forms):

```
OK: 30   not-proportional: 0   other-fail: 0
```

The §11.1 lossiness is **structural and size-independent** — if direction
reconstruction from `bound_metabolite` were broken, it would corrupt the small
central complexes too. It does not. The lossiness only bites *lumped,
multi-metabolite-per-step* opaque forms (the hand-written `EABEPQ <--> EQ + P`
fixtures), which Phase 1 already decomposed/deleted. **The enumerator never
generates those** — every enumerated step binds exactly one metabolite,
recorded in `bound_metabolite`. The remaining 47 (larger) init mechanisms hit
the documented compile-cost limit (CLAUDE.md "Known Issues"), which is a
compile-time concern, not a stoichiometry-correctness one.

**Consequence:** the roadmap's "Phase 2 exit gate — bi-bi derives without the
not-proportional error" is **already green**. It becomes a *regression guard*,
not a goal. Phase 2 is a behavior-preserving internal refactor, not a bug fix.

**Diagnostic 2 — where the opaque rep lives.** Structural dump of an init
mechanism's steps:

```
grp1: from=E  to=E_A   bound=[]  binds=A
grp2: from=E  to=E_Q   bound=[]  binds=Q
grp3: from=Estar to=Estar_B bound=[] binds=B
...
```

`init_mechanisms` emits **opaque** Species (empty bound lists, names `:E_A`,
`:Estar_B`, `:E_Q`). But `_make_species`
(`src/mechanism_enumeration.jl:155`) — the constructor the topology
backtracker calls — *already builds decomposed Species*
(`Species([Substrate(:A), …], :E)`). The decomposition is **lost on a
round-trip through `_RawSpec`**: `_stepspec_from_step`
(`src/mechanism_enumeration.jl:854`) renders Step → `_RawStep` Symbol form
names, the dead-end/grouping helpers run on the opaque `_RawSpec`, and
`EnzymeMechanism(spec::_RawSpec)` rebuilds Species opaquely.

**Diagnostic 3 — the keystone fact.** `name(::Species)` renders a decomposed
Species to the *same* Symbol as the opaque form:

```
name(Species([], :E_A))               == name(Species([Substrate(:A)], :E))   == :E_A
name(Species([], :Estar_B))           == name(Species([Substrate(:B)], :Estar)) == :Estar_B
name(Species([], :E_Q))               == name(Species([Product(:Q)], :E))      == :E_Q
```

Therefore swapping the enumerator from opaque to decomposed Species produces
**byte-identical form names**, hence an identical King-Altman graph, hence
identical rate equations. The refactor is fully gated by the existing
enumeration-count tests and derivation regression.

## 3. Target architecture

```
_catalytic_topologies / backtrack!   →  Vector{Vector{Step}}, decomposed   (unchanged)
   → dead-end enumeration: add CompetitiveInhibitor-bound decomposed Steps
   → equivalence grouping: structural merge of Vector{Vector{Step}}
   → Mechanism(reaction, steps)        — directly, no _RawSpec round-trip
```

The hardest part (topology generation with ping-pong / residual-atom
bookkeeping) **already emits decomposed structs**. Phase 2 deletes the
`_RawSpec` round-trip and lets the decomposed Steps flow straight through
dead-end enumeration and equivalence grouping into `Mechanism`.

### 3.1 Kill list

- `_RawSpec`, `_RawStep`, `_RawAllostericSpec` structs.
- `_mechanism_spec_from_steps`, `_stepspec_from_step` (the Step→_RawStep
  downgrade).
- `_mechanism_from_raw`, `_raw_from_mechanism` (the round-trip conversions).
- `EnzymeMechanism(spec::_RawSpec)`, `compile_mechanism(::_RawSpec)`,
  `compile_mechanism(::_RawAllostericSpec)`.
- Form-name parse/synthesis helpers: `_dead_end_form_name`,
  `_bound_mets_from_form_name`, `_parse_bound`, `_is_estar_form`, and
  `_form_name` (already dead — 0 callers).

### 3.2 Explicitly kept

- Atom helpers `_atoms_dict`, `_subtract_atoms`, `_can_pingpong`, `_add_atoms`
  — they perform real residual-atom math for the ping-pong backtracker, keyed
  by metabolite name, not form-name parsing. (Recon mislabeled these as
  opaque-Symbol helpers; the diagnostic shows they operate on atom dicts.)
- The 2-arg `EnzymeMechanism(metabolites, reactions)` constructor and the
  `EnzymeMechanism(m::Mechanism)` constructor — both used outside the
  enumerator (DSL, `compile_mechanism`). Consolidating the three construction
  paths is out of scope for this spec (see §7).

## 4. Execution approach (decision: A, inside-out)

Chosen over boundary-first (B — long red window, the §11.4 failure mode) and
parallel-then-swap (C — temporarily reintroduces a parallel representation).
Approach A keeps the suite green at every commit and front-loads the only real
risk.

## 5. The spike (gate before full commitment)

Rewrite **dead-end enumeration + `_apply_equivalence_grouping`** to operate on
decomposed `Vector{Vector{Step}}`, *still converting to `_RawSpec` at the
boundary* via the existing `_mechanism_spec_from_steps`, so the rest of the
pipeline is untouched and the change is isolated.

**Spike gate — prove on `test_mechanism_enumeration.jl`:**

1. Enumeration counts unchanged: **bi-bi = 11, ter-ter = 283** (CLAUDE.md
   verified counts).
2. `dedup!` canonicalizes identically (struct `==`/`hash`).
3. No `init_mechanisms` perf or compile-budget regression (compile-budget
   gate 750).
4. Derivation identical (rate-equation string + numerical) for a
   bi-bi-with-competitive-inhibitor mechanism.

**Open verification point (roadmap §⑥).** Today `_dead_end_form_name`
synthesizes dead-end form names *without* the Phase-1 `inh` suffix that the
`::Inh` rename gives competitive-inhibitor-bound derivation forms. Decomposed
`name(Species([…, CompetitiveInhibitor(:I)], :E))` may render differently from
`_dead_end_form_name`'s output. The spike must confirm either (a) name-equality
between the decomposed dead-end form and what the derivation expects, or (b)
that the rate equation is invariant to dead-end form-node spelling (node names
are graph labels; parameter names derive from step rep-index, and the canonical
hash is structural). Whichever holds, the spike's derivation-identity check
(gate #4) is the empirical arbiter.

**Escalation clause (roadmap + continuation §11.4).** If the opaque working-rep
proves load-bearing in a way decomposed Steps cannot express, **STOP and
escalate to Denis** — do not work around it by reverting architecture. This is
the exact failure mode a prior agent hit on this branch.

## 6. Sequencing (commits, green at each)

1. **Spike** — decomposed dead-end enumeration + equivalence grouping behind
   the `_RawSpec` boundary; pass the §5 gate. **← Denis reviews here.**
2. **Collapse the round-trip** — `init_mechanisms` builds `Mechanism` directly
   from the decomposed Steps; delete `_RawSpec`, `_mechanism_from_raw`,
   `_mechanism_spec_from_steps`, `_stepspec_from_step`,
   `EnzymeMechanism(::_RawSpec)`, `compile_mechanism(::_RawSpec)`.
3. **Sweep `_expand_*` moves** — replace the four `_dead_end_form_name` call
   sites with decomposed dead-end Step construction; delete `_RawAllostericSpec`,
   `_raw_from_mechanism`, `compile_mechanism(::_RawAllostericSpec)`.
4. **Delete form-name helpers** — `_dead_end_form_name`,
   `_bound_mets_from_form_name`, `_parse_bound`, `_is_estar_form`, `_form_name`;
   grep-confirm zero `src/` callers each. **← Denis reviews here.**
5. **Final** — dead-code sweep, LOC re-baseline + gate renegotiation (§8), §3
   success-criteria verification, README + CLAUDE.md update, PR.

## 7. Out of scope

- **Consolidating the three `EnzymeMechanism` construction paths** (2-arg
  tuple, `::_RawSpec`, `::Mechanism`) beyond deleting the `::_RawSpec` path.
  The 2-arg DSL path and the `::Mechanism` path both have non-enumerator
  callers; merging them is a separate change.
- **Parameter-naming refactor** (`:K1` → `:K_ATP`) — deferred per continuation
  spec §7.
- **Enumerator compile-cost** for large mechanisms (the 47 bi-bi forms that
  exceed the budget) — a pre-existing documented limit, not introduced or
  fixed here.

## 8. LOC treatment (decision: renegotiate to a measured target)

Per roadmap §④ and continuation spec §3.2 (LOC is the **secondary,
non-gating** criterion). The branch is at 9,421 src LOC; the §3.1 kill list
realistically removes ~400–700 LOC. The Final stage (§6.5) re-baselines the
true deletion against `wc -l src/*.jl` and documents the honest achieved number
in the PR, rather than treating ≤8,200 as a hard gate that could force a hack.
The gating criteria remain §3.1 #1–#4.

## 9. Non-negotiables (per commit)

Carried from the prior specs and the roadmap "Conventions" section:

- **TDD** — failing test → implement → green. No `--amend`; stay on the branch.
- **Test integrity** — no deletion/weakening/`@test_skip`/`@test_broken`; any
  helper-test deletion needs a `docs/superpowers/refactor-deleted-tests.md`
  §2.1 entry *in the same commit*. `bash scripts/check_test_integrity.sh main`
  EXIT = 0 (no pipe — it masks the exit).
- **Perf + compile + chokepoint gates per commit** — `test_rate_equation_performance`
  (0-alloc/<100ns), the 3 compile-budget gates, `test_chokepoint.jl`. They build
  their own small mechanisms, so run them every commit to catch a regression at
  its origin.
- **Full suite check** — `julia --project=. -e 'using Pkg; Pkg.test()' >
  /tmp/out.log 2>&1; echo "EXIT=$?"` then read the `Test Summary` line (the
  trailing echo swallows `Pkg.test`'s exit).
- **Evergreen comments** — no "Stage N"/"previously"/"legacy"/"will be" in
  code; spec/plan docs may reference stages.
- **Commit footer** — `src delta: -X / +Y net Z, cumulative: ±W` (`wc -l
  src/*.jl`, cumulative vs main 7,136); end messages with
  `Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>`.

## 10. Open questions

None at spec time. The one genuine uncertainty (dead-end form-name convention,
§5) is isolated by the spike and arbitrated empirically by the spike's
derivation-identity gate; if it cannot be satisfied, the §5 escalation clause
applies.
