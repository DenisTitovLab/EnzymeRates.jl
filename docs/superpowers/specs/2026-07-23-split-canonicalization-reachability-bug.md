# Bug: canonicalization merges Wegscheider-tied splits during enumeration, making distinct models unreachable

Date: 2026-07-23
Status: diagnosed, not yet fixed

## Summary

The split move (`_expand_split_kinetic_group`) canonicalizes each child with
`_canonical_mechanism`, which merges kinetic groups whose binding constants are
tied by a Wegscheider cycle. A single split that is individually a
parameter no-op is therefore reversed and dropped. This is correct
*deduplication* — the split fits the parent's rate equation — but it also prunes
the enumeration graph: because the intermediate single-split node never enters
the frontier, its children are never generated. Some of those children are
genuinely distinct models. They become unreachable, are never fit, and can never
be selected.

The bug conflates two separate jobs in one machine: **deduplication** (do not fit
the same rate equation twice) and **reachability** (do not prune paths to
distinct rate equations). Canonicalization does the first correctly, but doing it
during enumeration also does the second, incorrectly.

The deeper cause is that seeds carry **redundant steps** — surplus binding edges a
rapid-equilibrium segment does not need — and those edges are the cycle-closers
that create the Wegscheider ties in the first place. The primary fix (Direction A)
is therefore a canonical form that removes redundant *steps*, not one that merges
redundant *groups*: with the surplus edges gone, the existing split move has no
reverted splits and reaches the distinct models directly. Verified: 0/30 reverted
splits on step-minimized seeds versus 19/39 on the current seeds.

## Reproduction and evidence

Reaction: random-order bi-bi, `A + B ⇌ P + Q`, no regulators.

**1. A single split is Wegscheider-tied and reversed.** Take a seed's `bind A`
kinetic group — the parallel routes `E→EA`, `EB→EAB`, `EP→EAP`, which share one
dissociation constant. Split `E→EA` into its own group:

- before canonicalization: 6 kinetic groups (the split is performed)
- after `_canonical_mechanism`: 5 groups, identical to the parent
- the Wegscheider rename map on the raw split: `K_A_EB → K_A_E`

The split nominally creates `K_A_EB` (binding A to the EB form), but the
detailed-balance cycle around the A/B binding square forces `K_A_EB = K_A_E`. The
split adds no free parameter, so `_merge_tied_kinetic_groups` collapses it and
`_expand_split_kinetic_group` drops it (`child == mc`).

**2. Every single split of the random-order seed is reversed.**
`_expand_split_kinetic_group(seed)` returns zero children: no split-move edge
leaves the seed.

**3. Distinct richer models provably exist.** A fully-connected binding square
(E, EA, EB, EAB) has four rate constants — `K_A_E`, `K_B_E`, `K_B_EA`, `K_A_EB` —
constrained by one detailed-balance relation
`K_A_E · K_B_EA = K_B_E · K_A_EB`. That leaves **three** independent parameters,
versus **two** in the shared-constant collapse (`K_A_E = K_A_EB = K_A`,
`K_B_E = K_B_EA = K_B`). The third parameter is the substrate **interaction
factor** — how binding one substrate shifts the enzyme's affinity for the other —
a standard random-order rate-equation feature. It is a legitimate model with a
distinct rate equation.

**4. The distinct model is unreachable.** Reaching the three-parameter model
requires a mechanism whose binding routes sit in separate groups. Every path
there passes through an intermediate single-split, and every intermediate single
split is reversed (evidence 1 and 2). The collapsed seed cannot take even the
first step, so the interaction-factor model never enters the search.

**5. The derivation already reduces the tie; canonicalization is redundant for
correctness.** Compiling the raw 6-group split directly (`EnzymeMechanism` does
not re-merge — see the call-site audit below) derives the correct equation on its
own:

- `rate_equation_string` emits `# Wegscheider constraints: K_A_EB = K_A_E
  (substituted into v)` and the reduced `v` is byte-identical to the parent's,
  apart from that comment line.
- The comment-stripped dedup key (`_rate_eq_dedup_key`, i.e. `eq_hash`) is
  identical for the raw split and its parent, and both reduce to 5 fitted
  parameters.

So the group-merge does no work the derivation does not already do. Its only role
is to let the split move's `child == mc` check detect the duplicate at the
mechanism-structure level — a check that equation identity (`eq_hash`) performs
correctly on the raw form.

**6. The distinct richer models are valid and the derivation reduces them
correctly.** Fully atomizing the same seed (every step in its own group) and
compiling it yields a genuinely richer model, confirmed through the compiler
rather than the rename map:

- collapsed parent: 5 kinetic groups, **5 fitted parameters**
- fully atomized: 13 kinetic groups, **9 fitted parameters**, distinct `eq_hash`
- the four extra free parameters are `K_A_EB`, `K_A_EP`, `K_B_EQ`, `K_P_EQ` — the
  substrate/product interaction factors

The derivation reduces the 13-group form to 9 free parameters (not 13), so it
handles the uncollapsed input correctly and the richer model is a legitimate,
distinct rate equation. Reaching it requires a *coordinated* multi-group split;
every single-group split along the way is Wegscheider-redundant and is collapsed
by the current move, which is why the model is unreachable. (An earlier
`_build_wegscheider_rename_map` probe on the atomized form reported 14
"independent" parameters; that number is unreliable — the map-builder assumes
canonical input — and the compiler's 9 supersedes it.)

## Root cause

`_expand_split_kinetic_group` calls `_canonical_mechanism` on every child, stores
the **remade** (merged) child, and drops any child that canonicalizes back to the
parent:

```julia
child = _canonical_mechanism(_with_steps(m, _split_one_step(...)))
child == mc || push!(results, child)
```

Two things happen here, and only the first is wanted:

- **Dedup (correct).** `child == mc` recognizes a split that fits the parent's
  equation and drops it, so it is not fit twice.
- **Remake + prune (the defect).** The value stored is `_canonical_mechanism(...)`
  — the merged mechanism, not the raw split — and the tied split is deleted
  outright. Because the raw split never enters the frontier, its children are
  never generated.

`_canonical_mechanism` repeatedly applies `_merge_tied_kinetic_groups`, which
merges two kinetic groups when their binding constants map to the same
Wegscheider representative (`_build_wegscheider_rename_map`).

### Call-site audit

`_canonical_mechanism` (and `_merge_tied_kinetic_groups` beneath it) is called in
exactly **one** place: the split move, `_expand_split_kinetic_group`, for both
`Mechanism` and `AllostericMechanism`. `init_mechanisms`, `expand_mechanisms`, and
`EnzymeMechanism` (compilation) do **not** call it. So the group-merge remake is
localized to the split move — the fix surface is small — and every other
mechanism in the enumeration is already stored in its as-built form. The
"collapsed" group form is an invariant maintained incrementally (seeds are built
collapsed; the split move re-merges), not a global normalization pass.

### Two kinds of redundancy — and the one that matters

There are two distinct redundancies, and they lead to different canonical forms:

- **Redundant groups (mergeable).** Two kinetic groups whose binding constants
  the current `_canonical_mechanism` ties to the same Wegscheider representative,
  and merges. By this measure the seeds are already minimal — 0/55 bi-bi seeds
  have a mergeable group, and RE→SS never creates one; only the split move does
  (~79% of raw splits produce a mergeable group). This is the redundancy the
  current code acts on.
- **Redundant steps (removable edges).** A single step (edge) whose removal
  leaves the *equation* unchanged. In a rapid-equilibrium segment the form
  abundances are path-independent, so a parallel binding route is thermodynamic
  surplus: it can be dropped without changing the rate equation. `_canonical_-
  mechanism`'s group-merge does **not** touch these — it only merges groups, never
  removes steps.

By the *step* measure the seeds are **not** minimal. Removing each step from a
seed and checking `eq_hash`:

| Reaction | Seeds with ≥1 removable step |
|----------|------------------------------|
| uni-uni | 0 / 1 (a simple chain has no parallel routes) |
| bi-bi | 10 / 20 sampled; fully-connected seeds carry more |

### Step-minimization dissolves the reverted-split problem

Removing redundant *steps* is the canonical form that matters, because it removes
the cycle-closing edges that create the Wegscheider constraints in the first
place. Minimizing each bi-bi seed to a step-fixpoint, then comparing how the split
move behaves:

| | Reverted splits |
|--|-----------------|
| Original seeds | 19 / 39 |
| Step-minimized seeds | **0 / 30** |

After minimization every split is genuine. With the cycle-closer gone, the
remaining groups sit on a spanning structure with no cycle to force a split's new
parameter back to dependent — so the split move reaches richer models *directly*,
using the existing move, with no reverted intermediates and no frontier bloat.

This reframes the fix (see Direction A below): the primary lever is a canonical
form that removes redundant **steps**, not one that merges redundant **groups**.

**Where step-minimization must run.** It is needed wherever a move *adds* edges:
`init_mechanisms` (seeds carry surplus routes) and the edge-adding moves
(dead-end inhibitor and the allosteric moves, which add mirror steps — not yet
measured). It is **not** needed after the split or RE→SS moves: split only
regroups existing edges and RE→SS only flips a group's status, so neither adds a
removable edge. This answers the "seed-only vs after-every-move" question: seeds
plus edge-adding moves, but not split or RE→SS.

The existing code comment already recognizes the no-op splits but treats
dropping them as purely beneficial:

> "Splitting a group adds a parameter, but a Wegscheider cycle often forces that
> new parameter straight back to dependent, making the split a model-space no-op
> that fits the parent's equation. `_canonical_mechanism` merges such splits
> back, so each candidate is canonicalized and dropped when it returns to the
> parent."

Dropping them is right for *fitting*. It is wrong for *reachability*, because a
no-op split can be the only bridge to a non-no-op split.

## Impact

- **Missed models.** Interacting-site variants of random-order mechanisms — any
  model where a binding constant depends on what else is bound — are unreachable
  whenever the collapsed form ties the constants. These are real rate equations
  the package should be able to select.
- **Silent.** The search reports nothing missing; it simply never constructs
  these candidates.
- **Scope.** The gap is largest for random-order and symmetric topologies, where
  many binding routes share a constant. Ordered and dead-end-broken topologies
  are less affected, because their splits are not tied and already survive.

## Fix directions

All directions share one verified foundation — the derivation performs its own
Wegscheider reduction (findings 5–6), so a raw or step-minimized mechanism fits
correctly and its reduced equation matches the fully-connected form's. Direction A
is the primary fix: it removes the root cause. Direction B is a fallback that
tolerates the symptom instead. The parameter-level move is separate, more
extensive work.

### Direction A (primary): a step-minimal canonical form

Redefine the canonical form to remove redundant **steps** (edges), not to merge
redundant groups. A mechanism is canonical when no single step can be removed
without changing its rate equation. Physics dominates the hypothesis: a
rapid-equilibrium binding route the thermodynamics makes surplus is not part of
the model, so the canonical mechanism drops it.

The verified payoff (see *Step-minimization dissolves the reverted-split
problem*): once the surplus edges are gone, the split move has **no reverted
splits** — 0/30 on step-minimized bi-bi seeds, versus 19/39 on the current seeds.
The redundant edge was the cycle-closer that created the Wegscheider tie; without
it, every split frees a genuine parameter. So the *existing* split move reaches
the richer models directly, with no redundant intermediates to keep and no
frontier bloat. This is strictly better than Direction B: it removes the
duplicates at the source rather than carrying and deduplicating them.

The canonical form doubles as the dedup key. Two formulations that reduce to the
same physics minimize to the same steps, so keying fit-dedup on the step-minimal
form also collapses the parameter-renaming duplicates that `eq_hash` misses (the
~66% redundant-fit problem a prior investigation found).

Where it must run (measured): at `init_mechanisms` and after any move that *adds*
edges (dead-end inhibitor and the allosteric moves, which add mirror steps — the
latter not yet audited); **not** after split or RE→SS, which add no edges.

Open design points:

- **Determinism.** When two parallel routes are interchangeable, minimization must
  pick which to keep by a deterministic tie-break, or the canonical form is not
  unique. The `Step`/`Mechanism` canonical ordering is the natural basis.
- **Validity.** Removing an edge must leave a connected, atom-conserving,
  derivable mechanism. The measured runs only counted removals whose result
  compiled, so the minimizer must enforce this rather than assume it.
- **Completeness of reachability.** Minimization plus the existing split reaches
  models that split the *kept* edges. Models that need a *removed* edge back, with
  an independent constant, still require the separate add-edge move below. How
  large that residual is has not been measured.

### Direction B (fallback): keep raw forms, deduplicate fits by equation

If a correct step-minimizer proves hard to define (determinism, validity), the
symptom can be tolerated instead: keep the current fully-connected forms, stop the
split move from remaking and dropping tied splits, and separate the two jobs of
the `child == mc` check —

- **Expansion keeps raw forms.** Store the raw split, not
  `_canonical_mechanism(...)`; dedup the frontier by *structural* identity so a
  split is a distinct node from its parent and its children are reachable.
- **Fitting deduplicates by equation.** Collapse candidates by `eq_hash` before
  fitting, so the tied split and its parent are fit once.

This is obviously complete but reintroduces the frontier growth the current drop
avoids (~2/3 of bi-bi mechanisms), which must then be bounded — keep-all-fit-unique,
one-step lookahead, or a finest-representative per `eq_hash` class. Direction A
avoids this cost entirely, which is why it is preferred.

### A separate, more extensive fix: a parameter-level split move

Some interaction-factor models (finding 6: up to 9 parameters, distinct
`eq_hash`) require a *removed* edge back with an independent constant — the
residual Direction A leaves open. Reaching those is a larger change, logged here
as separate work. A move defined at the **parameter** level — "introduce one
independent parameter" by adding an edge and freeing its constant — reaches them
directly. That move is harder to define (it must reason about which edge additions
break a Wegscheider tie) and carries a completeness obligation (proving it reaches
every physics-distinct equation), so it belongs in its own spec once Direction A
lands. How large this residual is — how many distinct models need a re-added edge
rather than a split of a kept edge — has not been measured and should be, to size
this work.

### Open questions and risks

- **Cost of computing the minimal form (Direction A).** Minimization tests
  edge removals for equation preservation. Done naively (compile-and-compare per
  candidate edge) it is expensive; whether the Wegscheider rename map alone can
  identify removable edges cheaply, without compiling, needs working out.
- **`_build_wegscheider_rename_map` on non-canonical forms.** On a fully-atomized
  mechanism the rename map returned no dependents, contradicting the
  detailed-balance constraints that must exist. The map-builder may assume
  canonical input. The derivation's own reduction handles raw forms correctly
  (finding 5), but if minimization or dedup relies on the rename map over
  finer-grained forms, this path needs auditing.
- **`eq_hash` fallback robustness (Direction B).** If Direction B is taken,
  fit-dedup rests on `eq_hash`, which misses the renaming duplicates noted above.
  Direction A's step-minimal key does not have this problem, which is a further
  reason to prefer it.
- **Performance.** None of this touches `rate_equation` runtime — the change is
  entirely enumeration-time. Direction A should not grow the frontier (it removes
  duplicates); Direction B does and must be benchmarked against a real run
  (LDH/PFKP) before merge.

## Not in scope

- The `rate_equation` derivation and its 0-allocation / sub-120 ns contract.
- The RE/SS and combinatorics documentation, which describes the *current*
  behavior. Once this is fixed, the RE/SS subsection's "2⁴–2⁶" ceiling changes
  and the doc must be updated to match.
