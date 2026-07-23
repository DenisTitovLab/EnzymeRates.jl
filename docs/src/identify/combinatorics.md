# How many mechanisms? The combinatorics of enumeration

The set of mechanisms EnzymeRates can build is a *product* of
independent structural choices, so it grows explosively with the number of
substrates, products, and regulators. Understanding that explosion is the key to
understanding why the package aggressively filters the search instead of fitting
everything.

Write ``n_s``, ``n_p``, ``n_i``, and ``n_a`` for the number of substrates,
products, competitive inhibitors, and allosteric regulators. The total count is,
to an order of magnitude, a product over six categories:

```math
\text{mechanisms} \;\approx\;
(\text{catalytic topologies}) \times
(\text{kinetic-group splits}) \times
(\text{dead-end forms}) \times
(\text{RE/SS options}) \times
(\text{competitive-inhibitor options}) \times
(\text{allosteric-regulator options}).
```

For an enzyme with two substrates, two products, and one regulator the count
already starts around ``10^6``, and it grows rapidly with every added substrate,
product, or regulator. Even at one minute per fit, ``10^6`` equations would take
nearly two years on a single core, or about 17 hours on a 1000-core cluster.
These are calculated estimates, not measured counts: deduplication and validity
pruning trim the raw product, but the growth *rate* is the point. Even the bare
starting set grows steeply — `init_mechanisms` returns 1, 3, 55, and 42,220
mechanisms for uni-uni, uni-bi, bi-bi, and ter-ter, before a single expansion
move runs. The subsections below give each factor a formula and a concrete bi-bi
number.

## Catalytic topologies

A catalytic topology is one binding-and-release skeleton: the order in which
substrates bind and products leave. Ligands may bind one at a time in a fixed
order, or several at once (random order), so the number of arrangements of `n`
ligands is the **ordered Bell number** (Fubini number)

```math
a(n) = \sum_{k=1}^{n} k!\, S(n, k),
```

where ``S(n,k)`` is a Stirling number of the second kind: ``a(1) = 1``,
``a(2) = 3``, ``a(3) = 13``, ``a(4) = 75``. Substrate binding and product release
are independent, so a reaction has

```math
a(n_s)\, a(n_p)
```

sequential catalytic topologies — ``a(2)\,a(2) = 9`` for a bi-bi. When the atom
inventory lets a covalent fragment persist between the two half-reactions, the
engine also builds ping-pong topologies, which push the count higher still: a
ter-ter reaction reaches 250 topologies against the sequential ``13 \times 13``.

## Kinetic-group splits

The split move (move 2 in [The enumeration engine](@ref)) carves a multi-step
kinetic group into separate groups, giving each split-off step its own rate
constant. A mechanism with ``S`` steps in ``G`` groups has ``S - G`` steps it can
split off, on the order of ``2^{S-G}`` groupings. A bi-bi seed has ``S \approx 8``
steps in ``G = 5`` groups — a handful of split variants — and splitting also drives
the RE/SS factor below: once every step sits in its own group, each can be
rapid-equilibrium or steady state on its own.

## Dead-end forms

A dead-end form binds a mix of substrates and products at the catalytic site — an
abortive complex that leaves the cycle. Which reactants may co-occupy is a
bipartite graph on ``n_s`` substrates and ``n_p`` products in which every reactant
touches at least one edge, and the number of such graphs is

```math
N(n_s, n_p) = \sum_{i, j} (-1)^{i+j} \binom{n_s}{i} \binom{n_p}{j}\,
              2^{(n_s - i)(n_p - j)},
```

giving ``N(1,1) = 1``, ``N(2,2) = 7``, and ``N(3,3) = 265``. The
`shared_catalytic_site` reaction keyword prunes exactly this factor, by forbidding
a physiological substrate/product pair (`ATP`/`ADP`, `NAD`/`NADH`) from ever
sharing the site.

## RE/SS options

Every kinetic group is treated as either rapid-equilibrium (RE) or steady state
(SS), and the RE→SS move (move 1 in [The enumeration engine](@ref)) flips each RE
group independently. A mechanism with ``g`` RE groups therefore reaches

```math
2^{g}
```

distinct RE/SS assignments, and ``g`` grows as the split move (move 2) carves the
mechanism's ``S`` steps into finer groups. Splitting can, in principle, give every
step its own group, so a fully split mechanism reaches up to

```math
2^{S}
```

assignments — on the order of ``2^9`` for the random-order bi-bi's nine steps.
Thermodynamic ties collapse some of these onto the same rate equation, so the
distinct-model count sits somewhat below the raw ``2^S``, as with every factor on
this page.

!!! warning "Known limitation"
    The search does not currently reach the full ``2^S``. A correctness bug in the
    split move's canonicalization reverts the intermediate splits needed to fully
    separate a mechanism's steps, so the reachable RE/SS assignments top out around
    ``2^4 = 16`` to ``2^6 = 64`` rather than the intended ``2^9``. The random-order
    topology is the worst case: its binding steps cannot be separated at all, so it
    stays at ``2^4 = 16``. This is tracked in the project design specs
    (`docs/superpowers/specs/2026-07-23-split-canonicalization-reachability-bug.md`)
    and does not affect the correctness of the equations that *are* enumerated —
    only the completeness of the set.

## Competitive-inhibitor options

Each declared competitive inhibitor (move 3 in [The enumeration engine](@ref))
binds some nonempty set of the substrate-competing and product-competing enzyme
forms — a second competition count, now over the ``n_i`` inhibitors and the ``K``
forms they can reach. A bi-bi has up to 9 such competition patterns. In practice
one declared inhibitor adds a median of about four dead-end binding variants per
mechanism; a second inhibitor binds independently and roughly doubles that, to
about eight. Each further inhibitor multiplies the count again, so this factor
grows geometrically in ``n_i``.

## Allosteric-regulator options

Four moves make a mechanism allosteric, and their options multiply together.

Promotion to a two-conformation MWC model (move 4) fans a single mechanism into
many allosteric variants at once — one per non-empty subset of binding groups that
can go `:OnlyA`, plus one per `(regulator, tag)` pairing. A bi-bi promotes to about
15 variants with no allosteric regulator declared, and each declared regulator adds
two more V-type pairings: 17 variants with one regulator, 19 with two.

Each of the ``n_a`` regulators then takes an allosteric state — `:OnlyA`, `:OnlyI`,
`:EqualAI`, or `:NonequalAI` — as does each catalytic group, so the states grow as
roughly ``4^{n_a}`` times the per-group choices (adding an allosteric ligand,
move 5, and relaxing a shared state to independent A/I, move 6). Adding a ligand
contributes about three variants per regulator: three with one regulator, six with
two.

Finally, regulators may share a site (move 7). With a single regulator there is
nothing to merge, so the move is inert; it activates only once a mechanism carries
two or more regulators at separate sites, where the number of ways to partition
``n_a`` regulators across sites is the Bell number ``B(n_a)`` (``B(2) = 2``,
``B(3) = 5``).

## Why the search is filtered

Fitting even a fraction of that space is infeasible, so the search is filtered,
not exhaustive. The beam keeps only the promising candidates at each parameter
count; `eq_complexity_filter` drops equations whose denominator is too dense to
fit or derive practically; required-regulator seeding skips the
partially-regulated lower shelf; and `shared_catalytic_site` removes mechanisms a
chemist already knows are wrong. See [Best mechanism selection](@ref) for the
filters that bound the search, and [The enumeration engine](@ref) for the moves
that build the space in the first place.
