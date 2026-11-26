# The enumeration engine

The enumeration engine is the core of EnzymeRates.jl: it defines the *universe*
of enzyme mechanisms that [`identify_rate_equation`](@ref) can sample. The
engine starts from the simplest mechanisms and applies a fixed set of moves that
each turn one mechanism into a slightly more complex one, so the search visits
candidates in order of increasing complexity. A rate equation the engine cannot
construct is never fit, and so can never be selected.

At a high level the engine has two halves. The starting set holds the simplest
catalytic mechanisms for a reaction — each binding order, with optional dead-end
substrate and product inhibition — and their merged and Theorell–Chance variants. A
fixed set of single moves then grows the set, each taking a mechanism and returning
slightly more complex children: splitting a shared rate constant, flipping a
rapid-equilibrium group to steady state, adding a regulator, making the enzyme
allosteric, and so on. After every step the candidates are deduplicated, since
different move sequences can sometimes reach the same mechanism. Understanding the two
halves explains exactly which equations the search can reach.

`EnzymeRates.Enumeration.candidate_mechanisms` lists what the search would propose. It is
not exported. Called with a reaction, it returns the starting set; called with the
reaction and one of its candidates, it returns that candidate's children by every move.
Each candidate holds its mechanism in the field `mechanism`:

```@example engine
using EnzymeRates
rxn = @enzyme_reaction begin
    substrates: A[C], B[N]
    products:   P[C], Q[N]
end
start = EnzymeRates.Enumeration.candidate_mechanisms(rxn)
length(start)
```

```@example engine
start[1].mechanism
```

```@example engine
children = EnzymeRates.Enumeration.candidate_mechanisms(rxn, start[1])
length(children)
```

Both calls take the size limits of [`identify_rate_equation`](@ref),
`max_param_count` and `eq_complexity_filter`, and the first also its two lists of
optional regulators; by default nothing is dropped for size and every declared regulator
is required. Every candidate returned has a well-defined rate law and is within the
limits, so it is one the search may fit.

## The starting set

The starting set is first the seeds, then their merged and Theorell–Chance variants. For
each catalytic topology the engine enumerates all substrate/product dead-end subsets,
assigns one steady-state catalytic step, and collapses binding steps that share the same
`(metabolite, RE/SS)` class into a single kinetic group; each result is a seed.
When a reaction's atom inventory allows a covalent fragment to persist between
half-reactions, it also builds ping-pong mechanisms: the modified enzyme carries a
residual, named in its form's name (`E_res_plusA_minusP`), and a step that would return
the enzyme to free `E` with an empty residual mid-cycle is rejected, since it would
split the reaction into two disconnected half-cycles. See [Ping-pong mechanisms](@ref)
for detail on these mechanisms.

A seed writes chemistry as an isomerization: the enzyme isomerizes to a
product-bound form, and each release is its own step. A binding step may change
the enzyme's conformation; a conformation change is not a chemical reaction.

### Merged and Theorell–Chance variants

The moves never reach some textbook families from a seed. Among them are the
steady-state ordered and ping-pong laws, whose rates at zero products are
V·A·B/(Kia·Kb + Kb·A + Ka·B + A·B) and V·A·B/(Kb·A + Ka·B + A·B). Each seed
therefore also yields variants that fold its chemistry isomerization into the
steps around it:

1. **Merged base.** Every chemistry isomerization merges onto its product side
   (`_merge_isomerization`): the isomerization and its group go, and every other
   step at its substrate-side form moves to its product-side form, the *merged
   complex*. The last substrate then binds straight into the merged complex
   (`EA + B ⇌ EPQ`, a fused binding), and the products leave by plain
   releases. A sequential seed merges one isomerization, a ping-pong seed two.
2. **Theorell–Chance bases.** A merged complex with exactly two steps, the binding
   into it and one release, whose other ends are two different forms, is eliminated
   (`_eliminate_form`). Its two steps become one Theorell–Chance step
   (`EA + B ⇌ EQ + P`) in a group of its own. Each such complex gives one base. The
   uni-uni merged complex `EP` gives none: its binding and its release both lead to `E`.
3. **Search.** Every step of a base is at rapid equilibrium, which makes its rate
   infinite. Each smallest set of the base's groups whose flip to steady state
   gives a valid, flux-carrying, non-degenerate mechanism is a variant. As in the
   flip move, only groups that carry flux with every step at steady state take
   part.
4. **Lumping twins.** A merged variant whose every merged complex has both of its
   steps at steady state, each alone in its group, is skipped: it has the family,
   at the same count, of the unmerged mechanism with rapid-equilibrium flanks
   around a steady-state isomerization.
5. **Duplicates.** A variant equal to a seed or to an earlier variant is dropped.

Every test but the guard reads a candidate's flags over the base's graph; the guard runs
on the candidate built as a mechanism, once it passes the others:

- **Valid.** No cycle of rapid-equilibrium steps runs the net reaction
  (`has_turnover_cycle`), and the guard accepts the candidate
  (`check_valid_rate_equation`). Wherever the engine applies the guard, it also applies the
  *binding test*, which drops an allosteric candidate in which a conformation keeps a form
  that only bindings it excludes enter, though the guard may accept it.
- **Flux-carrying.** Every steady-state group carries net flux.
- **Non-degenerate.** A degenerate mechanism's rate at zero products ignores some
  substrate or grows without bound as the substrates grow, or the same holds for
  products. Two conditions rule it out:
  - **V, a maximal rate both ways** (`_has_vmax`). Turning rapid equilibrium every
    steady-state step that takes up or gives off a substrate leaves no
    rapid-equilibrium turnover cycle, and the same holds for products.
  - **C, chemistry in equilibrium with one side at most**
    (`_chemistry_equilibrates_both_sides`). No chemistry node is left both by a
    rapid-equilibrium step that releases a substrate and by one that releases a
    product, and no Theorell–Chance step is at rapid equilibrium. A chemistry
    node is a merged complex with every form that rapid-equilibrium
    isomerizations join to it, or a set of forms joined by rapid-equilibrium
    isomerizations. The merged bi-bi variant with its A and Q groups at steady
    state fails C: its merged complex releases B through its equilibrated entry
    and P through its equilibrated exit. The variant with A and P at steady state
    passes.

Every seed with two chemistry isomerizations, such as a ping-pong seed, fails C
itself. A seed holds one steady-state step, so one isomerization is at rapid
equilibrium, and so are the bindings and releases on both its sides: the forms it
joins release a substrate and a product by rapid-equilibrium steps. In a bi-bi
ping-pong seed they release B and Q, so at zero products its rate does not depend
on B. The starting set holds, in place of such a seed, its flip children that are not
degenerate, and the search never fits the seed itself (see
[Best mechanism selection](@ref)).

The starting set therefore holds mixed parameter counts. A bi-bi reaction has 55 seeds
at 5 parameters and 184 variants. Its 164 merged variants fit 5 parameters when two of
their kinetic groups are at steady state (108 of them) and 6 when three are (56). Its 20
Theorell–Chance variants fit 5 to 7: the 4 without dead-end forms fit 5, and the 16 with
them 5, 6 or 7. When its atoms allow a ping-pong cycle it has 62 seeds and 202 variants;
its 7 ping-pong seeds fit 6 parameters and their 18 variants 5, one fewer. A uni-uni seed
has no variant.

## Starting from the required regulators

For a reaction that declares regulators, the seeds start far below the
useful region: every seed binds zero regulators, so the search must climb through
every partially-regulated mechanism before it reaches one that binds them all. When
the effectors are already known — phosphofructokinase and its five regulators, say —
that lower shelf is pure waste.

So the starting set begins where the search belongs. The engine takes the unregulated
starting set, in which each degenerate seed is already replaced by its flip children
(`_base_tier`), and grows it through the regulator-binding moves alone — go allosteric, add
an allosteric regulator, add a competitive inhibitor — under three constraints that keep
the set small: each required regulator binds at its own single-ligand site, no regulator
is `:NonequalAI`, and no partially-regulated mechanism survives. A kinetic group can be
`:NonequalAI`: move 4 tags a set of groups so where the same set tagged `:OnlyA` is
refused, and move 6, the only other way to such a version, does not run here. The result is
every fully-regulated mechanism those moves build, and nothing that binds fewer
regulators; like the seeds and variants, it holds mixed parameter counts. The beam then
refines these with the detail moves — steady-state flips, splits, and `:NonequalAI`
relaxations — as usual.

The growth starts from a degenerate seed's flip children, not from the seed, because the
engine refuses some allosteric states on the seed that it accepts on its children: on the
six-step ping-pong seed, each version with one half-reaction `:OnlyA`. Started from the
seed, the search would never build those versions, since with a regulator required it fits
no plain mechanism whose children could.

By default every declared regulator is required. `identify_rate_equation`'s
`optional_allosteric_regulators` and `optional_competitive_inhibitors` keywords move
named regulators back to optional, so the beam adds them as refinements rather than
forcing them into every seed; listing every regulator as optional recovers the
unregulated starting set of seeds and variants. Requiring a competitive inhibitor and an
allosteric regulator together multiplies the set: for the bi-bi reaction
`A[CX] + B[N] → P[C] + Q[NX]` at oligomeric states 1 and 2 it holds 236,444 mechanisms,
which take about three minutes to build, against 14,078 with only the regulator required
and 1,548 with only the inhibitor (the repository's `ROADMAP.md` lists this). When no
mechanism binds every required regulator — a uni-uni reaction whose substrate is also
declared a competitive inhibitor, say — `candidate_mechanisms` stops with an error that
names the regulators and these two keywords. Declaring a regulator's type in the reaction —
`::Activator` or `::Inhibitor` — pins it to one allosteric state and shrinks the set
further. The [Identify tutorial](@ref) works a concrete example.

## The seven expansion moves

`candidate_mechanisms(reaction, parent)` applies all seven moves to `parent` and returns the
resulting children in one list. Each move is applied in every applicable way, so one parent
yields many children — every set of rapid-equilibrium groups that can flip to steady state,
every way a group can be split by binding context, every way an inhibitor can bind, and so
on. Every child is checked to conserve atoms, and a child that the well-posedness guard
rejects (see [What the package rejects](@ref)) or that fails the binding test is dropped.
Above one catalytic subunit a child is also dropped when either conformation's law carries
a concentration to a power of its own, whichever move built it (see **Conformational
mechanisms carry hyperbolic catalytic schemes** under [Modeling choices](@ref)).

The seven moves:

### 1. Flip rapid-equilibrium groups to steady state

Turns whole kinetic groups from rapid equilibrium (RE) to steady state (SS): every
step in a group converts together, keeping one pair of rate constants per group.
The move emits one child per smallest set of groups whose joint flip divides a
rapid-equilibrium segment in two. On a seed every group cuts on its own, so each
RE group the move offers gives one child. Once splits have separated a metabolite's
binding steps, or a catalytic step from its inhibitor-bound mirror, a single group may be
bridged by an RE route through the others; the move then flips the bridging groups
together, because a steady-state step whose two ends stay in one equilibrated segment
never reaches the rate equation.

Competitive-inhibitor binding never flips, and a steady-state group must carry net flux
(both are [Modeling choices](@ref)). A group none of whose steps would carry flux with
every step at steady state is never offered for flipping. A set of flips that leaves one
of its groups without flux in the child, because an equilibrated route around it carries
the turnover, counts as failed and is extended like a set that divides no segment.

The move never flips a flank of a steady-state chain (`_chain_flank_groups`). A
*steady-state chain* is a steady-state isomerization between two forms that each have one
other step, a binding into that form, with each of the three steps alone in its group. The
two bindings are the chain's *flanks*. By
[Apparent constants on a steady-state chain](@ref), flipping a flank adds one or two
parameters and no rate law. The uni-uni seed is such a chain and has no flip child. In an
allosteric mechanism the chain is read on the catalytic scheme, whatever the states of its
groups, and its flanks are skipped too. There the rule is a modeling choice, not a proof:
a flank's flip can show through the two conformations, mostly over more than one catalytic
subunit, and the search does not build those mechanisms (the repository's `ROADMAP.md`
lists this).

**Parameter delta:** +1 per flipped group in most cases. A flip can give an
equation that is the parent's up to a reparametrization (one more constant) even
though it divides a segment, as a flank of a steady-state chain would. The move skips
those flanks; any other such child is fit once and loses to its parent on parsimony.

### 2. Split a kinetic group by binding context

Divides one kinetic group into two, so the two parts have separate rate
constants. The division is always by context: the steps whose enzyme form
already carries some other ligand Y, sits in a given conformation, or carries
a given covalent residual, against the rest. It encodes the hypothesis that
the affinity for a metabolite depends on what is bound next to it, on the
enzyme's conformation, or on its covalent state; with Y a competitive
inhibitor it separates a catalytic step from its inhibitor-bound mirror.

A single split often frees no parameter. In random-order binding, thermodynamics
ties the new constant straight back to the old one, and the rate function is
unchanged. The move therefore emits the smallest sets of simultaneous splits, at
most one per group, whose combined effect raises the number of independent
parameters; the thermodynamic constraint solver decides. Random-order bi-bi needs
the A group and the B group split together before either constant is free.

A part of a steady-state group that carries no flux, such as an abortive binding
separated from the catalytic binding whose constants it shared, is emitted at
rapid equilibrium: its two rates would enter the equation only as their ratio, and
the equilibrium form is the same family with one constant fewer. The count test
runs on that form, so a constant the thermodynamic ties pull back is still
rejected. A set whose child fails the guard or the binding test counts as failed too,
so the move extends it rather than emitting it. No split child holds a redundant
copy (see move 3). The move judges each child whole: splitting one group can
complete a copy's dwell gauge, and splitting a group that forms the productive
complexes the copy duplicates can break it, so whether a part of a copy group is
redundant depends on every group the child splits.

**Parameter delta:** at least +1 by construction. A child that would add nothing
is never emitted.

### 3. Add a competitive inhibitor binding site

Adds binding steps for a competitive inhibitor declared in the reaction, as a
dead-end complex with the enzyme. The move enumerates the combinations of enzyme
species the inhibitor can bind, emitting one child mechanism per combination,
subject to three rules:

- **Binding capacity.** The inhibitor is not added to a form that already holds
  every substrate or every product, counting productive bindings only. A
  substrate or product declared as its own competitive inhibitor binds as a copy
  at a dead-end site of its own, apart from the reactant's catalytic site, so a
  copy fills no catalytic site.
- **Mirror steps.** If the inhibitor binds two enzyme forms that a catalytic
  step already connects, a mirror step is added between the two inhibitor-bound
  forms, so the inhibitor-bound branch stays connected to the cycle. Each mirror
  inherits its counterpart's kinetic group and adds no parameter. Productive
  steps decide competition: an inhibitor that competes with a substrate or
  product binds at forms where that metabolite is free on a productive step — a
  binding's source form, and for a Theorell–Chance step `EA + B ⇌ EQ + P`
  the form `EA` for `B` and `EQ` for `P` — and never at a form that holds it
  productively, while a form that carries only its copy is a site. A form that
  already carries a competing ligand never receives the inhibitor,
  so the binding step of a competing substrate is never mirrored and the
  inhibitor-bound branch can never complete the net reaction. In a ping-pong
  mechanism it can carry out the one half-reaction whose ligands the inhibitor
  does not compete with.
- **No redundant copy.** A dead-end complex duplicates a productive complex, one
  that carries no competitive inhibitor, when a rapid-equilibrium route joins the
  two and binds exactly the copied ligand along it, so that their weights are
  proportional. A complex of the same composition formed by a steady-state binding
  is no duplicate. A substrate or product declared as its own competitive inhibitor
  then binds in a second orientation of a complex the mechanism already has. The
  move drops a copy when every complex it forms duplicates a productive complex and
  a dwell gauge absorbs its constant into the existing binding's, so that the data
  cannot separate the two. It keeps a copy whose gauge fails, as when a shared
  kinetic group pins the existing binding, since its constant may be identifiable.
  A complex that matches only an inhibitor-bound complex counts as new: the sites
  that pin each inhibitor may keep both constants apart. In an allosteric mechanism
  the test runs in every conformation the copy binds, with a kinetic group shared by
  both conformations rescaled alike in each, and every conformation holds its free
  enzyme. A copy that duplicates a complex in the active conformation but is the
  only such complex in the inactive one keeps a visible constant, and the move keeps
  it. On a uni-uni seed with one conformation every placement is redundant, and
  `candidate_mechanisms` reports the unsatisfiable requirement.

The inhibitor's own binding steps form one fresh kinetic group (one new
dissociation constant `K_R`).

On an allosteric parent each placement gives three children, the new group `:EqualAI`,
`:OnlyA` and `:OnlyI`: both conformations bind the inhibitor with one constant, only the
active one binds it, or only the inactive one does. The move drops a child in which a group
that binds a competitive inhibitor runs in neither conformation, or binds it redundantly as
above: either way the data cannot determine its constant. Take the three-step ping-pong
`E + Q <--> EQ`, `E + A <--> E_res_plusA_minusP + P`, `E_res_plusA_minusP + B <--> EQ`
with both half-reactions `:OnlyA`, so that the inactive conformation never reaches the
covalent form `E_res_plusA_minusP`. An inhibitor placed at free E gives all three children.
Placed at the covalent form alone it gives the `:EqualAI` and the `:OnlyA` child, and no
`:OnlyI` one, which would bind nowhere. Where the inactive conformation has none of the
forms the inhibitor binds, as with the covalent form here, the `:EqualAI` and the `:OnlyA`
child have one rate equation and share a fit. Both stay, since they part once move 6 lets
the inactive conformation reach such a form.

**Parameter delta:** +1 per competitive-inhibitor group added.

### 4. Promote a non-allosteric to allosteric mechanism

Converts a plain mechanism to allosteric ones, each with the plain mechanism's steps and
an allosteric state for every kinetic group ([Mechanisms with allosteric regulators](@ref)
defines the states). A chemistry group holds an isomerization, a fused binding or a
Theorell–Chance step (`_chemistry_groups`); every other group is a binding group. For each
oligomeric state of the reaction's `allowed_oligomeric_states`, the move emits three kinds
of child.

**Every set of `:OnlyA` groups.** One child per nonempty set of kinetic groups: those
groups `:OnlyA`, so that the inactive conformation does not run them, and every other group
`:EqualAI`. The uni-uni mechanism `E + S ⇌ ES`, `ES <--> EP`, `E + P ⇌ EP` has three
groups and seven such children: the inactive conformation lacks S's binding, the chemistry
or P's binding, any two of them, or all three. A ping-pong child whose `:OnlyA` groups are
one half-reaction's chemistry stops that half in the inactive conformation and leaves the
other running in both, as phosphorylation stops PDH's E1 and leaves E2 and E3 running. A
group that binds a competitive inhibitor has a third state, `:OnlyI`, in which only the
inactive conformation binds the inhibitor, so the uni-uni mechanism with an inhibitor bound
at free E has four groups and 2³ · 3 − 1 = 23 children. The engine tags no other group
`:OnlyI`. The move does not test these children. A set of `:OnlyA` groups is not always a
well-posed mechanism (see [How an excluded step is switched off](@ref)), and
`candidate_mechanisms` drops the children that fail, as it does for every move.

**A declared regulator beside a chemistry-only set.** A child whose `:OnlyA` groups are
all chemistry groups, and which tags no group `:OnlyI`, comes twice more for each declared
allosteric regulator: with the regulator at a new regulatory site, once `:OnlyA` and once
`:OnlyI`. Such a child shows nothing in a substrate or product titration that the plain
mechanism lacks, so without a regulator it could leave the search before move 5 adds one.
With a regulator `R` declared, the uni-uni child whose one `:OnlyA` group is the chemistry
comes three times: alone, with `R` `:OnlyA` and with `R` `:OnlyI`. A child with an `:OnlyA`
binding group changes the titration and comes alone.

**`:NonequalAI` sets that no `:OnlyA` set leads to.** An inactive conformation can run
every step, some groups with constants of its own: one that binds A and B more weakly and
runs the same chemistry, say. Move 6 reaches most such versions from a child above by
tagging its `:OnlyA` groups `:NonequalAI`. Where the same groups tagged `:OnlyA` are
refused, no child leads there, so this move emits the `:NonequalAI` version itself. It finds
each *minimal gaining set* (`_minimal_gaining_sets`): a set of groups that, tagged
`:NonequalAI` with every other group `:EqualAI`, adds a fitted constant beyond `L` and gives
a child the engine returns, while no smaller set does. It emits that child when the same
set tagged `:OnlyA` would not be returned. One rapid-equilibrium binding group alone never
gains here: an inactive conformation that runs a whole catalytic cycle must give the same
`Keq` as the active one, so thermodynamics gives the binding back the active conformation's
constant ([Thermodynamic constraints of MWC equations](@ref) works the uni-uni case). Two
rapid-equilibrium bindings gain one constant, the second tied to the first. A steady-state
group, a binding or chemistry, brings two rate constants of which the cycle ties back one,
so it gains one constant alone.

The six-step ping-pong with both half-reactions at steady state shows the three kinds:

```
E + A ⇌ EA
E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP
E + Q ⇌ EQ
E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP
EA <--> EP_res_plusA_minusP
EB_res_plusA_minusP <--> EQ
```

Its six groups give 63 sets of `:OnlyA` groups. Eight minimal sets gain: each steady-state
half-reaction alone, and each of the six pairs of rapid-equilibrium bindings. Tagged
`:OnlyA`, three of the pairs are refused, A with B, Q with P, and B with P, so the move
emits those three tagged `:NonequalAI`: 66 children at one subunit with no regulator. The
other five sets are among the 63, and move 6 reaches their `:NonequalAI` versions from them.
With a regulator declared, each of the three children whose `:OnlyA` groups are
half-reactions alone (the first, the second, both) comes twice more.

The move drops a child in which a group that binds a competitive inhibitor runs in neither
conformation, or binds it redundantly (see move 3). An `:OnlyI` inhibitor binds nowhere when
the `:OnlyA` groups cut every form it binds off from free E in the inactive conformation.
Nothing stands in for a dropped child. The all-`:EqualAI` version is never emitted (the
conformations would be identical and `L` would cancel), and children whose laws agree up to
renaming the constants share one key, so identify fits the law once. A parent gets no
allosteric child only when every child is dropped. None of the 5,158 starting mechanisms of
the ter-bi reaction `A[C] + B[N] + D[X] → P[CN] + Q[X]` is such a parent: with no regulator
declared they get a median of 26 allosteric children. The search applies this move only to
a plain mechanism.

A child does not always reveal `L`. Over more than one catalytic subunit the conformational
equilibrium enters the rate with the oligomeric state as its power. Over one subunit, when
the inactive conformation turns nothing over, it is a dead-end branch of the free enzyme,
and `L` can be a *phantom*, absorbed by an exact rescaling of the other constants. Whether
it is depends on what the inactive conformation keeps. Take the children that tag the
chemistry and at least one binding `:OnlyA`. In every such child of the rapid-equilibrium
ordered bi-bi seed a rescaling of the active state's constants absorbs `L`, whatever the
inactive conformation binds. In the merged ordered bi-bi whose B group mixes a fused and a
plain binding it is a phantom only where the inactive conformation binds nothing: no
rescaling absorbs an A or Q binding kept there. Binding nothing is not enough: in a
random-order merged bi-bi, whose A group takes A up at E and at EB under one rate constant,
`L` shows in every such child. The engine emits the children where `L` is a phantom all the
same, since their own children can carry laws no other mechanism has. When the search has
fitted another mechanism with the law such a child has at `L = 0`, identify reads that fit
instead of fitting the child; the `L_phantom` column in [Best mechanism selection](@ref)
gives the exact rule. A child whose inactive conformation turns the reaction over is always
fitted, though over one subunit some such children, and some with an `:OnlyI` inhibitor, fit
a constant the data cannot determine (the repository's `ROADMAP.md` lists these).

A parent whose catalytic scheme already carries concentration powers, such as a
random-order scheme with steady-state binding or a substrate that traps a steady-state
intermediate in an abortive complex, is promoted at oligomeric state 1 only; see
**Conformational mechanisms carry hyperbolic catalytic schemes** under
[Modeling choices](@ref).

In MWC terminology the A-state corresponds to the R-state and the I-state
to the T-state of the original Monod–Wyman–Changeux nomenclature; this
package uses A/I throughout.

**Parameter delta:** +1 for a child of `:OnlyA`, `:OnlyI` and `:EqualAI` groups — the
conformational equilibrium constant `L` is the sole new parameter, a phantom in some
children over one subunit. `:OnlyA` zeroes out the inactive state and `:OnlyI` the active
one, and neither adds a parameter of its own. +2 for a child with a regulator — `L` plus
the regulator's binding constant. At least +2 for a `:NonequalAI` child — `L` plus what its
set adds.

### 5. Add an allosteric ligand

Adds one of the reaction's allosteric regulators as a ligand at a new or an existing
regulatory site, with an allosteric state drawn from `{:OnlyA, :OnlyI, :NonequalAI}` (plus
`:EqualAI` only at a mixed existing site). The search applies this move only to an
allosteric mechanism.

**Parameter delta:**
- `:OnlyA` or `:OnlyI`: **+1** (one binding constant `K_A` or `K_I`).
- `:NonequalAI`: **+2** (independent `K_A` and `K_I`).

### 6. Relax kinetic groups to `:NonequalAI`

Tags kinetic groups `:NonequalAI`, from `:EqualAI`, `:OnlyA` or `:OnlyI`: the inactive
conformation then runs them with constants of its own, and the active conformation runs them
too. The search applies this move only to an allosteric mechanism. The move emits one child
per minimal gaining set (`_minimal_gaining_sets`), binding and chemistry groups alike: a set
of groups not yet `:NonequalAI` whose tagging raises the fitted count and gives a child the
engine returns, with no inhibitor bound redundantly, while no smaller set does. A set that
fails is extended by any other group. Then it tags each regulatory ligand `:NonequalAI`
alone, keeping the child when the fitted count rises. No group leaves `:NonequalAI`.
`:OnlyA` is the limit of `:NonequalAI` as the inactive conformation's affinity or rate goes
to zero, and `:OnlyI` the limit as the active conformation's affinity does, so a child that
fits the data well with a group `:OnlyA` can grow into the version that runs it weakly in
the inactive conformation.

Sets matter because thermodynamics ties constants together. Where the inactive conformation
runs a whole catalytic cycle, one rapid-equilibrium binding tagged alone gets the active
conformation's constant back and adds nothing, while two such bindings add one and a
steady-state group adds one alone (see move 4). An `:OnlyA` binding beside a tagged one can
have its affinity pinned in the inactive conformation, which the guard refuses; tagging both
frees it. Two parents of the six-step ping-pong of move 4, at one subunit, show both cases:

- With its first half-reaction `:OnlyA` and every other group `:EqualAI`, the inactive
  conformation turns nothing over, no cycle ties a new constant to an old one, and each
  of the six groups tagged alone gains: six children.
- With the bindings of A and P `:OnlyA`, A tagged beside P `:OnlyA` is refused, and so is
  P beside A, while the two tagged together gain. The bindings of B and Q and the second
  half-reaction gain alone. The first half-reaction adds nothing, since the inactive
  conformation, which binds neither A nor P, never reaches its forms. Four children.

**Parameter delta:** at least +1 by construction. A rapid-equilibrium binding adds at most
one constant and a steady-state group at most two; a thermodynamic cycle of the inactive
conformation can take one back. A regulatory ligand adds one.

### 7. Merge two regulatory sites

Combines two regulatory sites into one shared site, so their ligands compete for it rather
than binding independently. For each merged pair the move also enumerates the antagonist
forms: one ligand may switch to `:EqualAI` — binding both conformations equally, with no
allosteric effect of its own — so that it acts purely by displacing the other ligand from
the shared site. An activator that displaces an inhibitor still reads as activation, and an
inhibitor that displaces an activator still reads as inhibition, so the observable effect is
preserved. The search applies this move only to an allosteric mechanism, and a mechanism
with one regulatory site gets no child from it.

**Parameter delta:** **+0**. Each ligand keeps its one binding constant; only the
binding topology changes. The merged mechanism is a distinct rate equation, which the
beam fits alongside the independent-site form so that cross-validation can choose
between them.

## Modeling choices

The moves encode a few decisions about which mechanisms are worth fitting. They
are choices, not theorems, and this section states them so they can be revisited.

**Steady-state detail enters by whole groups.** Steps that share a rate constant
flip to steady state together, in the smallest set of groups that separates a
new equilibrated segment. Sharing structure is refined first, by splits, and
steady-state detail follows it. A mechanism in which one enzyme form binds a
metabolite at equilibrium while another binds it at steady state is reachable,
but only after a split has given the two bindings separate constants.

**A steady-state group carries flux.** A group whose every step lies off every
cycle that runs the reaction carries no net flux at steady state. Its forward and
reverse rates enter the equation only as their ratio, which the equilibrium form
already has, so the moves never leave such a group at steady state: the flip does
not make one, and the split reverts the one it would create. The test is
structural, on the graph of rapid-equilibrium segments (`_flux_carrying_groups`),
and reads each step's metabolite lists, so it holds for fused and
Theorell–Chance steps as well. It reads the whole catalytic scheme, so it misses a
steady-state step between two complexes of an `:OnlyI` inhibitor, which runs only in the
inactive conformation, where that conformation turns nothing over (the repository's
`ROADMAP.md` lists this).

**Conformational mechanisms carry hyperbolic catalytic schemes.** An MWC
conformational equilibrium and a random-order steady-state catalytic scheme
each put a metabolite concentration to a power in the rate equation, and
sigmoidal data cannot tell the two sources apart. The engine therefore never
combines them: above one catalytic subunit it returns no mechanism, whichever move built
it, in which either conformation's law carries such a power (`_hyperbolic_where_required`,
applied by both methods of `candidate_mechanisms`), and promotion to allosteric gives a
parent whose catalytic scheme already carries one no variant above one subunit. The
searches of moves 4 and 6 read the rule too: tagging one more group `:NonequalAI` can make
a conformation's law hyperbolic again, so a set whose child breaks the rule is extended.
With one catalytic subunit the conformational equilibrium
only reweights each enzyme form and adds no power of its own, so such a
mechanism may carry any catalytic scheme. The test is structural
(`_hyperbolic_catalysis`): it reads the equation's degree in a metabolite from the
rapid-equilibrium segment graph without deriving it. On an allosteric mechanism it reads
each conformation's law (`steps_running_in_allosteric_conformations`), since a
conformation that runs part of the catalytic graph can carry a power the whole graph
lacks: once the inactive conformation stops one half-reaction of a ping-pong, its abortive
complexes can put A² into its law. Every binding of a substrate or product at its
catalytic site counts, abortive complexes included, so an ordered scheme whose substrate
traps a steady-state intermediate, as pyruvate does with E·NAD⁺ in lactate dehydrogenase,
stays out of allosteric mechanisms. Binding of a declared competitive inhibitor does not
count, because an inhibitor binds a site of its own; declaring a substrate as a dead-end
inhibitor is how substrate inhibition enters an allosteric mechanism, provided the copy is
not redundant (move 3). Hand-written mechanisms are not subject to the rule: an allosteric
`@enzyme_mechanism` with random-order steady-state binding still derives and fits.

**Every step but a plain binding is chemistry.** A binding takes up exactly one
metabolite and gives off none. It is plain when its product form is its source
form with that metabolite added, the conformation free to change, and fused
otherwise. The engine keeps its own record of what each form it builds holds, and the
moves tell catalysis from binding by those records: every step but a plain binding is
chemistry, whether an isomerization, a fused binding or a Theorell–Chance step. Kinetic
groups go by stoichiometry instead: a group's steps take up and give off the same
metabolites, so a fused and a plain binding of one metabolite may share its constant,
as the B group of a merged seed with a dead-end complex does (`EA + B ⇌ EPQ` and
`EQ + B ⇌ EBQ`). A step that only gives off one metabolite is stored as the binding it
reverses, so a hand-written release that carries chemistry (`EAB <--> EQ + P`) is a
fused binding of `P` into `EAB`: a merged complex on the substrate side, which the
enumerator never builds. The moves expand only the candidates the engine built, so such
a mechanism takes no part in the search.

**Merging happens in the starting set only.** No move merges a complex or
eliminates a form, and no move rewrites a child into another mechanism of its
family. Merged and Theorell–Chance mechanisms descend from the starting
variants. The flip rule judges each flip when it is made: when a later split makes a
steady-state chain, a flank flipped before stays steady state, and the mechanism stays
as a tolerated phantom.

**Competitive-inhibitor binding stays at equilibrium.** An inhibitor bound to
two forms that a catalytic step connects carries flux through its mirror step,
so keeping its binding at rapid equilibrium is a modeling choice ("inhibitor
binding is fast"). What competition does decide is turnover: the inhibitor
competes with at least one substrate and one product, and the binding step of
a competing ligand is never mirrored onto the inhibitor-bound branch, so that
branch never completes the net reaction. Ping-pong is the one enumerated family
whose catalytic cycle has more than one chemistry step; there an inhibitor that
competes with the first half-reaction's ligands can still let the modified enzyme run
the second half with the inhibitor bound, which is the two-site picture, and the
inhibitor must leave before the next cycle.

**Splits follow binding context.** A group is divided only by whether another
ligand is already bound, by conformation, or by covalent residual. Arbitrary
partitions are not tried, a group is never split three ways in one move, and
groups never merge. A partition that no sequence of context splits produces is
unreachable.

**A child is never a provable copy of its parent.** The moves reject a child
whose rate function can be shown to equal the parent's: a split the constraint
solver ties back, a flip that leaves the segment count unchanged, a flip that
leaves a steady-state group without flux, a flip of a steady-state chain's flank in a
plain mechanism, a redundant dead-end copy, an inhibitor that binds in neither
conformation, a `:NonequalAI` tag that adds no fitted constant. A few children whose law
is the parent's up to renaming the constants survive; they share the parent's key, so
identify fits their law once, and they never win selection.
