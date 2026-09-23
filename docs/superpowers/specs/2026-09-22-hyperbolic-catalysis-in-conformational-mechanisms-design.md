# Hyperbolic catalysis inside conformational mechanisms

Date: 2026-09-22. Status: implemented.

## Goal

The enumerator combines two independent sources of concentration powers in one
mechanism: an MWC conformational equilibrium and a catalytic scheme with enough
steady-state (SS) steps to put `[A]^2` in its own rate equation (random SS
bi-bi, random SS product release in uni-bi). Sigmoidality from either source
looks the same to noisy data, so the stacked mechanism fits as well as either
alone, costs fits, widens the beam, and is sometimes selected. This design makes
the two sources exclusive: a conformational mechanism's catalytic scheme must be
hyperbolic in every substrate and product.

## The rule

A conformational mechanism (today `AllostericMechanism`; later KNF, mnemonic,
or slow-isomerization types) may only carry a catalytic scheme whose rate
equation has degree at most 1 in every substrate and product concentration.
The rule holds at every catalytic multiplicity, including `n = 1`.

Every binding of a substrate or product at its catalytic site counts, abortive
complexes included. Binding of a declared competitive inhibitor is a separate
source of non-hyperbolic terms and stays allowed inside conformational
mechanisms: an inhibitor binds a site of its own by definition, including a
substrate declared as a dead-end inhibitor. The predicate therefore leaves
out every step that touches a form carrying a declared inhibitor.

The rule is an enumeration prior, not a validity property of the equation. It
lives in the enumeration moves, never in a constructor: a hand-written
`@allosteric_mechanism` with random SS binding still derives and fits.

## The predicate

`_hyperbolic_catalysis(m)` returns `true` when the King–Altman denominator of
`m`'s catalytic scheme has degree at most 1 in every substrate and product.

Construction:

1. Drop every step whose either form carries a declared inhibitor (a
   `Regulator`). Inhibitor bindings and the catalytic steps mirrored onto
   inhibitor-bound forms leave together, so what remains is the catalytic
   scheme. Every substrate and product binding in it is scored, dead-end and
   abortive bindings included.
2. Build the rapid-equilibrium (RE) segment quotient graph over the kept
   steps: nodes are RE segments (union-find over RE steps, as in
   `_compute_re_groups`); every SS step is an edge in both directions.
3. Within each segment, the per-metabolite exponent clearing of
   `_bottomless_re_segment` gives each form the metabolites it carries beyond
   the segment's bottom form ("extras").
4. For metabolite `X`, score a directed edge `u → v` as
   `[the step binds X in that direction] + extras(source form of u → v)[X]`,
   and score a segment `S` as `max over forms in S of extras[X]`. A segment
   can carry `X` twice without any edge: a ping-pong whose second chemistry
   step is at rapid equilibrium puts free `E` one `B` above the covalent
   intermediate, and an abortive `E(B, Q)` two above it.
5. The degree of `X` is the maximum over root segments `S` of
   `score(S) + (max-weight spanning arborescence toward S)`. This is the
   `X`-exponent of a Cha-style King–Altman denominator term: the root
   segment's weight times the tree's edge weights.
6. `_hyperbolic_catalysis(m)` is `true` when every substrate and product has
   degree at most 1.

Only "at most 1" matters, so the maximum is never solved in general. The
degree exceeds 1 exactly when one segment or one edge scores 2 or more, or a
root scoring 1 has an arborescence toward it holding a scoring edge, or some
arborescence holds two scoring edges. An arborescence toward `S` containing
given edges exists iff every segment still reaches `S` once each given edge's
source keeps that edge as its only way out (a digraph has a spanning
arborescence toward `S` iff every vertex reaches `S`). That is a handful of
breadth-first searches per metabolite on a graph of at most a few dozen nodes,
so a call costs about 0.1 ms on a ter-ter parent.

An allosteric mechanism is scored on its full catalytic step set, which is its
A-state. The I-state prunes `:OnlyA` groups and the forms they disconnect.
Pruning removes edges and forms and can only shrink segments; a smaller segment
has a bottom at least as high, so every form's extras shrink or stay, and every
I-state arborescence is an A-state arborescence. The I-state degree therefore
never exceeds the A-state degree, and one evaluation per mechanism suffices.
The exactness test below checks the I-state anyway.

## Placement

One trait names the types the rule applies to:

```julia
_requires_hyperbolic_catalysis(::Mechanism) = false
_requires_hyperbolic_catalysis(::AllostericMechanism) = true
```

A future conformational type adds one method. Two moves call the predicate:

- `_expand_to_allosteric(m::Mechanism, rxn)` returns no children when
  `_hyperbolic_catalysis(m)` is `false`. One evaluation per parent, not per
  variant.
- `_expand_re_to_ss(m)` filters its emitted minimal flip sets through the
  predicate when `_requires_hyperbolic_catalysis(m)`. The check stays out of
  the `gains` predicate of `_minimal_gaining_sets`: a set that fails there is
  extended with more flips, and more SS steps never restore hyperbolicity, so
  the search would walk every superset for nothing.

No other move touches the flux-carrying SS structure. Split keeps every edge
and only changes which steps share parameters. The dead-end move adds
inhibitor bindings, which the predicate leaves out. Allo-state and
regulatory-site moves change tags and sites only. Seeds carry one SS step and
one RE segment, so they pass unless that segment's weight already carries a
power: 2 of the 7 bi-bi ping-pong seeds carry the abortive complexes E·B·Q on
both enzyme forms, have B² and Q², and are not promoted.

Reachability is preserved. Every conformational mechanism with a hyperbolic
catalytic scheme is still reached: from the RE seed by flips that stay
hyperbolic, or by promoting a hyperbolic `Mechanism`. A rejected flip child's
descendants are non-hyperbolic too, so dropping it loses nothing.

## Testing

Predicate unit tests, each fixture written inline with the macros:

- ordered SS bi-bi: `true`
- random SS bi-bi: `false`
- random RE bi-bi with one SS chemistry step: `true`
- ordered SS bi-bi with substrate A as an abortive complex on E(Q): `false`
  (A binds its catalytic site)
- ordered SS bi-bi with A declared as a dead-end inhibitor, all four
  placements: `true` (the powers come from the inhibitor site)
- uni-bi with random SS product release: `false`
- ping-pong with one SS step per half-reaction: `true`

Move tests in `test/test_mechanism_enumeration.jl`, following its three rules:

- `_expand_to_allosteric` on a uni-bi with random SS product release emits
  zero children; on an ordered SS bi-bi with an abortive complex of A on
  E(Q) it emits none, and with A declared as a dead-end inhibitor on E and
  E(Q) it emits 31 children.
- `_expand_re_to_ss` on a uni-bi with random RE product release emits 3
  children as a plain `Mechanism`, every expected child written out; the
  same parent as an allosteric mechanism emits only the hyperbolic one.

Exactness test: for every `Mechanism` in the uni-bi and ping-pong enumerations
to two expansion levels and the bi-bi enumeration to one level plus the
allosteric children of its level-1 allosteric mechanisms, derive it and check
that the structural degree of every substrate and product equals the largest
exponent of that symbol in the derived denominator. For every
`AllostericMechanism` at the same depth, check the same on both the A-state
and I-state projections. The reactions declare no inhibitors, so nothing is
left out and each mechanism is compared with its own derivation.

Performance: the predicate must not move the per-parent expansion time
measurably; the ter-ter seed is the check.

## Docs and version

- `docs/src/identify/enumeration_engine.md`: state the modeling assumption and
  that dead-end inhibition is exempt.
- `docs/src/developer.md`: the predicate, its scoring, and the trait.
- `Project.toml`: 0.6.0 → 0.7.0 (enumeration output changes).

## Out of scope

- KNF, mnemonic, and slow-isomerization mechanism types. The trait and the
  predicate are the hook they will use.
- A reaction-level opt-out. The rule is hard-wired by decision.
