# How many mechanisms? The combinatorics of enumeration

Even a simple enzyme has a surprising number of plausible mechanisms. The reason
is that a mechanism is built from several independent choices, and the counts
multiply. This page walks through those choices for an enzyme with two
substrates and two products, a bi-bi reaction, with concrete numbers. The point
is not the exact totals but how fast they grow, which is why
[`identify_rate_equation`](@ref) filters the search instead of fitting
everything.

## Binding order

Substrates can bind one after another in a fixed order, or in any order. Two
substrates give three arrangements: A first, B first, or either. Products leave
the same way, independently, so a bi-bi has 3 × 3 = 9 binding-and-release
skeletons. A third substrate raises the substrate side from 3 to 13, and a
ping-pong enzyme, which releases a product before the next substrate binds, adds
more skeletons on top.

## Dead-end complexes

An enzyme can bind a substrate and a product at the same time and get stuck: a
dead-end complex that sits off the catalytic cycle. Each pairing of a substrate
with a product may or may not be allowed, and every reactant must take part in
at least one. For a bi-bi that gives 7 patterns of dead ends. Together with
binding order, that is why `init_mechanisms` returns 55 bi-bi starting
mechanisms before any refinement, and 35,665 for a ter-ter reaction.

The `shared_catalytic_site` keyword removes patterns a chemist already knows are
impossible, such as ATP and ADP occupying one site together.

## Which forms share an affinity

A starting mechanism gives each metabolite one binding constant, whatever the
enzyme already holds. The search then asks whether the affinity for A changes
once B is bound, or once a product is bound. Each such question divides A's
binding steps into two groups with separate constants, and the questions
compose: A's affinity may depend on B alone, on P alone, or on both. For the
random-order bi-bi with dead ends, the search reaches 16 different ways of
sharing constants, from one constant per metabolite up to one constant per
binding step.

## Equilibrium or steady state

Each binding can be fast enough to sit at equilibrium, or slow enough that its
on and off rates both matter. The search flips groups of shared steps to steady
state, one group or several at a time, and each flip divides the enzyme's forms
into more regions that are in equilibrium among themselves. The number of ways
to divide the forms is the real count here: the random-order bi-bi with dead
ends has 9 enzyme forms and 1,434 distinct ways of partitioning them into
equilibrated regions. The number of steady-state variants grows with the number
of enzyme forms, which grows with substrates, products, and dead ends together.

## Competitive inhibitors

A declared competitive inhibitor can bind any set of the enzyme forms that hold
the metabolite it competes with. One inhibitor adds a median of about four
variants per mechanism; a second inhibitor binds independently and roughly
doubles that. Each further inhibitor multiplies the count again.

## Allosteric regulation

Making an enzyme allosteric adds a second conformation. Every binding group may
then bind in both conformations equally, in both with different affinities, or
in the active one only, and every regulator takes one of four states. A bi-bi
promotes to about 15 allosteric variants with no regulator declared, and each
regulator adds variants of its own; two regulators may also share a site or sit
at separate sites.

## Putting it together

The choices multiply. Binding order times dead ends gives the 55 starting
mechanisms; sharing structure, steady-state detail, inhibitors, and allosteric
states each multiply that by tens to thousands. Fitting one equation takes
seconds to minutes, so an exhaustive search of a bi-bi with one regulator would
take years on one core. The search therefore keeps only the promising
candidates at each parameter count, drops equations too dense to fit, starts
from fully regulated mechanisms when the regulators are known, and honors
`shared_catalytic_site`. See [Best mechanism selection](@ref) for those filters
and [The enumeration engine](@ref) for the moves that build the space.
