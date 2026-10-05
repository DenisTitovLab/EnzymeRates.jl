# The enumeration engine

The enumeration engine is the core of EnzymeRates.jl: it defines the *universe*
of enzyme mechanisms that [`identify_rate_equation`](@ref) can sample. The
engine starts from the simplest mechanisms and applies a fixed set of moves that
each turn one mechanism into a slightly more complex one, so the search visits
candidates in order of increasing complexity. A rate equation the engine cannot
construct is never fit, and so can never be selected.

At a high level the engine has two halves. `EnzymeRates.init_mechanisms` builds
the starting set: the simplest catalytic mechanisms for a reaction — each binding
order, with optional dead-end substrate and product inhibition — and their merged
and Theorell–Chance variants. `EnzymeRates.expand_mechanisms` then grows the set
through a fixed set of single moves, each taking a mechanism and returning
slightly more complex children: splitting a shared rate constant, flipping a
rapid-equilibrium group to steady state, adding a regulator, making the enzyme
allosteric, and so on. After every step the candidates are deduplicated, since
different move sequences can sometimes reach the same mechanism. Both functions
are internal, but understanding them explains exactly which equations the search
can reach.

## `EnzymeRates.init_mechanisms`

`EnzymeRates.init_mechanisms(reaction)` returns the starting mechanisms of the
search: first the seeds, then their merged and Theorell–Chance variants. For each
catalytic topology it enumerates all substrate/product dead-end subsets, assigns
one steady-state catalytic step, and collapses binding steps that share the same
`(metabolite, RE/SS)` class into a single kinetic group; each result is a seed.
When a reaction's atom inventory allows a covalent fragment to persist between
half-reactions, it also builds ping-pong mechanisms: the modified enzyme stays
on conformation `:E` carrying a residual rather than a separate conformation
label, and a step that would return the enzyme to free `:E` with an empty
residual mid-cycle is rejected, since it would split the reaction into two
disconnected half-cycles. See [Ping-pong mechanisms](@ref) for detail on these
mechanisms.

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
   (`E(A) + B ⇌ E(P, Q)`, a fused binding), and the products leave by plain
   releases. A sequential seed merges one isomerization, a ping-pong seed two.
2. **Theorell–Chance bases.** A merged complex with exactly two steps, the binding
   into it and one release, is eliminated (`_eliminate_form`): the two steps
   become one Theorell–Chance step (`E(A) + B ⇌ E(Q) + P`) in a group of its
   own. Each such complex gives one base.
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

The tests run on a candidate's groups; a mechanism is built only for a variant:

- **Valid.** No cycle of rapid-equilibrium steps runs the net reaction
  (`_re_turnover_cycle`), and no rapid-equilibrium segment is bottomless
  (`_bottomless_re_segment`, the rule the `Mechanism` constructor enforces).
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
on B. The beam expands such a seed without fitting it (see
[Best mechanism selection](@ref)).

The starting set therefore holds mixed parameter counts. A bi-bi reaction has 55
seeds at 5 parameters and 184 variants: 164 merged ones, 108 at 5 parameters and
56 at 6, and 20 Theorell–Chance ones at 5 to 7. When its atoms allow a ping-pong
cycle it has 62 seeds and 202 variants; its 7 ping-pong seeds fit 6 parameters
and their 18 variants 5. A uni-uni seed has no variant.

## `EnzymeRates.seed_mechanisms`

For a reaction that declares regulators, `init_mechanisms` starts far below the
useful region: every seed binds zero regulators, so the search must climb through
every partially-regulated mechanism before it reaches one that binds them all. When
the effectors are already known — phosphofructokinase and its five regulators, say —
that lower shelf is pure waste.

`EnzymeRates.seed_mechanisms(reaction, required_allosteric, required_competitive)`
starts the search where it belongs. It grows every mechanism of
`init_mechanisms`, variants included, through the regulator-binding moves alone —
go allosteric, add an allosteric regulator, add a competitive inhibitor — under
three constraints that keep the set small: each
required regulator binds at its own single-ligand site, every allosteric state stays
cheap (`:OnlyA`/`:OnlyI`, never `:NonequalAI`), and no partially-regulated mechanism
survives. The result is every fully-regulated mechanism those moves build, and
nothing that binds fewer regulators; like the starting set, it holds mixed
parameter counts. The beam then refines these seeds with the detail moves —
steady-state flips, splits, and `:NonequalAI` relaxations — as usual.

By default every declared regulator is required. `identify_rate_equation`'s
`optional_allosteric_regulators` and `optional_competitive_inhibitors` keywords move
named regulators back to optional, so the beam adds them as refinements rather than
forcing them into every seed; listing every regulator as optional recovers the
`init_mechanisms` starting set. When no mechanism binds every required regulator —
a uni-uni reaction whose substrate is also declared a competitive inhibitor, say —
`seed_mechanisms` stops with an error that names the regulators and these two
keywords. Declaring a regulator's type in the reaction —
`::Activator` or `::Inhibitor` — pins it to one allosteric state and shrinks the set
further. The [Identify tutorial](@ref) works a concrete example.

## The seven expansion moves

`EnzymeRates.expand_mechanisms(mechs, rxn)` applies all seven moves to every input
mechanism and returns the resulting child mechanisms pooled into one list. Each
move is applied in every applicable way, so one input mechanism yields many
children — every set of rapid-equilibrium groups that can flip to steady state,
every way a group can be split by binding context, every way an inhibitor can
bind, and so on. Every child is checked to conserve atoms before it is returned.

The seven moves:

### 1. Flip rapid-equilibrium groups to steady state

Turns whole kinetic groups from rapid equilibrium (RE) to steady state (SS): every
step in a group converts together, keeping one pair of rate constants per group.
The move emits one child per smallest set of groups whose joint flip divides a
rapid-equilibrium segment in two. On a seed every group cuts on its own, so each
RE group the move offers gives one child. Once splits have separated a metabolite's
binding steps, a single group may be bridged by an RE route through the others;
the move then flips the bridging groups together, because a steady-state step
whose two ends stay in one equilibrated segment never reaches the rate equation.

Competitive-inhibitor binding never flips: it stays at rapid equilibrium by modeling
choice. A steady-state group must carry net flux. A group none of whose steps would
carry flux with every step at steady state is never offered for flipping, and a set of
flips that leaves one of its groups without flux in the child, because an equilibrated
route around it carries the turnover, counts as failed and is extended like a set that
divides no segment. Such a group's two constants enter the equation only as their ratio,
which its equilibrium form already has.

A `Mechanism` never flips a flank of a steady-state chain. A *qualifying chain*
is an isomerization alone in its group between two forms that each have exactly
one other step, a binding into that form alone in its group; those two bindings
are the chain's *flanks*. When the isomerization is steady state, the rest of the
mechanism sees the chain only through its flux and the enzyme it holds, which the
chain with both flanks at rapid equilibrium already covers, so flipping a flank
adds one or two parameters and no rate law (`_chain_flank_groups`; see
[Apparent constants on a steady-state chain](@ref)). The uni-uni seed is such a
chain and has no flip child. An `AllostericMechanism` keeps flipping flanks: over
more than one catalytic subunit a flank's flip can show through the two
conformations.

An allosteric parent with more than one catalytic subunit also drops every child
whose catalytic scheme would carry a concentration to a power of its own — the
random-order steady-state pattern in which a metabolite binds on two steps that
one King–Altman tree can hold, or an equilibrated segment already carries the
metabolite once. See
**Conformational mechanisms carry hyperbolic catalytic schemes** under
[Modeling choices](@ref).

**Parameter delta:** +1 per flipped group in most cases. A flip can give an
equation that is the parent's up to renaming the constants even though it divides
a segment, as a flank of a steady-state chain would. The move skips those flanks;
any other such child is fit once and loses to its parent on parsimony.

### 2. Split a kinetic group by binding context

Divides one kinetic group into two, so the two parts have separate rate
constants. The division is always by context: the steps whose enzyme form
already carries some other ligand Y, sits in a given conformation, or carries
a given covalent residual, against the rest. It encodes the hypothesis that
the affinity for a metabolite depends on what is bound next to it, on the
enzyme's conformation, or on its covalent state; with Y a competitive
inhibitor it separates a catalytic step from its inhibitor-bound mirror.

A single split often frees no parameter. In random-order binding, thermodynamics
ties the new constant straight back to the old one, and the equation is
unchanged. The move therefore emits the smallest sets of simultaneous splits, at
most one per group, whose combined effect raises the number of independent
parameters; the thermodynamic constraint solver decides. Random-order bi-bi needs
the A group and the B group split together before either constant is free.

A part of a steady-state group that carries no flux, such as an abortive binding
separated from the catalytic binding whose constants it shared, is emitted at
rapid equilibrium: its two rates would enter the equation only as their ratio, and
the equilibrium form is the same family with one constant fewer. The count test
runs on that form, so a constant the thermodynamic ties pull back is still
rejected. No split child holds a redundant copy (see move 3). The move judges
each child whole: splitting one group can complete a copy's dwell gauge, and
splitting a group that forms the productive complexes the copy duplicates can
break it, so whether a part of a copy group is redundant depends on every group
the child splits.

**Parameter delta:** at least +1 by construction. A child that would add nothing
is never emitted.

### 3. Add a competitive inhibitor binding site

Adds binding steps for a `CompetitiveInhibitor` declared in the reaction, as a
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
  binding's source form, and for a Theorell–Chance step `E(A) + B ⇌ E(Q) + P`
  the form `E(A)` for `B` and `E(Q)` for `P` — and never at a form that holds it
  productively, while a form that carries only its copy is a site. A form that
  already carries a competing ligand never receives the inhibitor,
  so the binding step of a competing substrate is never mirrored and the
  inhibitor-bound branch can never complete the net reaction. In a ping-pong
  mechanism it can carry out the one half-reaction whose ligands the inhibitor
  does not compete with.
- **No redundant copy.** A dead-end complex duplicates a productive complex, one
  that carries no competitive inhibitor, when a rapid-equilibrium route joins the
  two and binds exactly the copied ligand along it, so their weights are
  proportional; a complex of the same composition formed by a steady-state binding
  is not a duplicate. A substrate or product declared as its own competitive
  inhibitor then binds in a second orientation of a complex the mechanism already
  has. A copy is
  not emitted when every complex it forms duplicates a productive complex and a
  dwell gauge absorbs its constant into the existing binding's, so that the data
  cannot separate the two; a copy whose complexes all
  duplicate productive complexes but whose gauge fails, as when a shared kinetic
  group pins the existing binding, is emitted, since its constant may be
  identifiable. A complex that matches only an inhibitor-bound complex counts as
  new: the sites that pin each inhibitor may keep both constants apart. In an
  allosteric mechanism the test is made in every conformation the copy binds, with
  a kinetic group shared by both conformations rescaled alike in each, and every
  conformation holds its free enzyme: a copy that duplicates a complex in the
  active conformation but is the only such complex in the inactive one keeps a
  visible constant and is emitted. On a uni-uni seed with one conformation every
  placement is redundant, and `seed_mechanisms` reports the unsatisfiable
  requirement.

The inhibitor's own binding steps form one fresh kinetic group (one new
dissociation constant `K_R`).

**Parameter delta:** +1 per competitive-inhibitor group added.

### 4. Promote a non-allosteric to allosteric mechanism

Converts a `Mechanism` to an `AllostericMechanism` variant set. An `:OnlyA`
catalytic binding asserts `K_I → ∞`: the inactive conformation cannot bind
that metabolite, so it cannot complete the catalytic cycle. The MWC reading is a
**catalytically-dead** inactive conformation — every chemistry group `:OnlyA` —
that binds ligands but runs no chemistry. A chemistry group holds an
isomerization, a fused binding or a Theorell–Chance step (`_is_chemistry`), so a
merged mechanism's variants tag its fused binding `:OnlyA` with the rest of its
chemistry and never tag that binding on its own. The engine emits, per
multiplicity: every non-empty subset of binding groups `:OnlyA`, each with every
chemistry group `:OnlyA` (a K-type mechanism, emitted bare); plus the empty subset
(every chemistry group `:OnlyA`, every binding group `:EqualAI`) paired with a
declared allosteric regulator, one variant per `(regulator, tag)` with
`tag ∈ {:OnlyA, :OnlyI}` (a V-type mechanism, where the regulator makes `L`
identifiable). The all-`:EqualAI` baseline is never emitted (the conformations
would be identical and `L` would cancel). A binding subset that leaves a
binding-only Wegscheider cycle unsatisfiable is dropped (see
[Thermodynamic constraints of MWC equations](@ref)). Enumeration runs over
`allowed_catalytic_multiplicities`. No-op on an already allosteric input.

A K-type variant does not always reveal `L`. Over more than one catalytic
subunit the conformational equilibrium enters the rate with the multiplicity's
power. Over one subunit the inactive conformation is a dead-end branch of the
free enzyme, and `L` is a phantom in some K-type variants and not in others. In
every K-type variant of the rapid-equilibrium ordered bi-bi seed a rescaling of
the active state's constants absorbs it, whatever that conformation binds. In
the merged ordered bi-bi whose B group mixes a fused and a plain binding it is a
phantom only in the variants whose inactive conformation binds nothing: no
rescaling absorbs an A or Q binding kept there. Binding nothing is not enough:
in a random-order merged bi-bi, whose A group takes A up at E and at E(B) under
one rate constant, `L` shows in every K-type variant. The engine emits the
variants where `L` is a phantom all the same, so a fit can carry an `L` the data
cannot determine.

A parent whose catalytic scheme already carries concentration powers, such as a
random-order scheme with steady-state binding or a substrate that traps a
steady-state intermediate in an abortive complex, is promoted at catalytic
multiplicity 1 only; see
**Conformational mechanisms carry hyperbolic catalytic schemes** under
[Modeling choices](@ref).

In MWC terminology the A-state corresponds to the R-state and the I-state
to the T-state of the original Monod–Wyman–Changeux nomenclature; this
package uses A/I throughout.

**Parameter delta:** +1 for a K-type variant — the conformational equilibrium
constant `L` is the sole new parameter, a phantom in some variants over one
subunit. +2 for a V-type variant — `L` plus the paired regulator's binding
constant. `:OnlyA` zeroes out the I-state and adds no parameter of its own.

### 5. Add an allosteric ligand

Adds one `AllostericRegulator` at a new or an existing regulatory site, with
an allosteric-state tag drawn from `{:OnlyA, :OnlyI, :NonequalAI}` (plus
`:EqualAI` only at a mixed existing site). No-op on a non-allosteric input.

**Parameter delta:**
- `:OnlyA` or `:OnlyI` tag: **+1** (one binding constant `K_A` or `K_I`).
- `:NonequalAI` tag: **+2** (independent `K_A` and `K_I`).

### 6. Relax an `:EqualAI` group to independent A/I states

Changes one `:EqualAI` or `:OnlyA`/`:OnlyI` allosteric-state tag to
`:NonequalAI`, giving the active-state and inactive-state versions of the
group independent parameters. No-op on a non-allosteric input.

**Parameter delta:**
- RE binding group: **+1** (one shared K splits into `K_A` and `K_I`).
- SS rate group: **+2** (each SS rate `kf` and `kr` splits into A/I pairs,
  adding two new independent constants).

### 7. Merge two regulatory sites

Combines two regulatory sites into one shared site, so their ligands compete for it
rather than binding independently. For each merged pair the move also enumerates the
antagonist forms: one ligand may switch to `:EqualAI` — binding both conformations
equally, with no allosteric effect of its own — so that it acts purely by displacing
the other ligand from the shared site. An activator that displaces an inhibitor still
reads as activation, and an inhibitor that displaces an activator still reads as
inhibition, so the observable effect is preserved. No-op on a non-allosteric input or a
single-site mechanism.

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
Theorell–Chance steps as well.

**Conformational mechanisms carry hyperbolic catalytic schemes.** An MWC
conformational equilibrium and a random-order steady-state catalytic scheme
each put a metabolite concentration to a power in the rate equation, and
sigmoidal data cannot tell the two sources apart. The moves therefore never
combine them: promotion to allosteric gives a parent whose catalytic scheme
already carries a power no variant above one catalytic subunit, and the flip
move on an allosteric parent with more than one subunit drops a child that
would introduce one. With one catalytic subunit the conformational equilibrium
only reweights each enzyme form and adds no power of its own, so such a
mechanism may carry any catalytic scheme. The test is structural
(`_hyperbolic_catalysis`): the equation's degree in a metabolite is read from
the rapid-equilibrium segment graph without deriving it. Every binding of a
substrate or product at its catalytic site counts, abortive complexes included,
so an ordered scheme whose substrate traps a steady-state intermediate, as
pyruvate does with E·NAD⁺ in lactate dehydrogenase, stays out of allosteric
mechanisms. Binding of a declared competitive inhibitor does not count, because
an inhibitor binds a site of its own; declaring a substrate as a dead-end
inhibitor is how substrate inhibition enters an allosteric mechanism, provided the
copy is not redundant (move 3).
Hand-written mechanisms are not subject to the rule: an `@allosteric_mechanism`
with random-order steady-state binding still derives and fits.

**Every step but a plain binding is chemistry.** A binding takes up exactly one
metabolite and gives off none. It is plain when its product form is its source
form with that metabolite added, the conformation free to change, and fused
otherwise. The moves tell catalysis from binding by `_is_chemistry`, true for
every step but a plain binding: isomerizations, fused bindings and
Theorell–Chance steps. Kinetic groups go by stoichiometry instead: a group's
steps take up and give off the same metabolites, so a fused and a plain binding
of one metabolite may share its constant, as the B group of a merged seed with a
dead-end complex does (`E(A) + B ⇌ E(P, Q)` and `E(Q) + B ⇌ E(B, Q)`). A step
that only gives off one metabolite is stored as the binding it reverses, so a
hand-written release that carries chemistry (`E(A, B) <--> E(Q) + P`) is a fused
binding of `P` into `E(A, B)`: a merged complex on the substrate side, which the
enumerator never builds. `expand_mechanisms` accepts such a parent, but its
children are not the ones a merged seed would give.

**Merging happens in the starting set only.** No move merges a complex or
eliminates a form, and no move rewrites a child into another mechanism of its
family. Merged and Theorell–Chance mechanisms descend from the starting
variants. The flip rule judges each flip when it is made, so a chain that a later
split makes qualifying keeps its steady-state flank; such a mechanism stays as a
tolerated phantom.

**Competitive-inhibitor binding stays at equilibrium.** An inhibitor bound to
two forms that a catalytic step connects carries flux through its mirror step,
so keeping its binding at rapid equilibrium is a modeling choice ("inhibitor
binding is fast"). What competition does decide is turnover: the inhibitor
competes with at least one substrate and one product, and the binding step of
a competing ligand is never mirrored onto the inhibitor-bound branch, so that
branch never completes the net reaction. Ping-pong is the one enumerated family
whose catalytic cycle has more than one chemistry step; there an inhibitor that
competes with the first
half-reaction's ligands can still let the modified enzyme run the second half
with the inhibitor bound, which is the two-site picture, and the inhibitor
must leave before the next cycle.

**Splits follow binding context.** A group is divided only by whether another
ligand is already bound, by conformation, or by covalent residual. Arbitrary
partitions are not tried, a group is never split three ways in one move, and
groups never merge. A partition that no sequence of context splits produces is
unreachable.

**A child is never a provable copy of its parent.** The moves reject a child
whose equation can be shown to equal the parent's: a split the constraint solver
ties back, a flip that leaves the segment count unchanged, a flip that leaves a
steady-state group without flux, a flip of a steady-state chain's flank, a
redundant dead-end copy. A few children whose equation is the parent's up to
renaming the constants survive; they cost one fit each and never win selection.
