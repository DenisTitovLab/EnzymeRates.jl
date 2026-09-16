# The enumeration engine

The enumeration engine is the core of EnzymeRates.jl: it defines the *universe*
of enzyme mechanisms that [`identify_rate_equation`](@ref) can sample. The
engine starts from the simplest mechanisms and applies a fixed set of moves that
each turn one mechanism into a slightly more complex one, so the search visits
candidates in order of increasing complexity. A rate equation the engine cannot
construct is never fit, and so can never be selected.

At a high level the engine has two halves. `EnzymeRates.init_mechanisms` builds
the starting set: the minimum-parameter catalytic mechanisms for a reaction —
each binding order, with optional dead-end substrate and product inhibition — at
their lowest parameter count. `EnzymeRates.expand_mechanisms` then grows the set
through a fixed set of single moves, each taking a mechanism and returning
slightly more complex children: splitting a shared rate constant, flipping a
rapid-equilibrium group to steady state, adding a regulator, making the enzyme
allosteric, and so on. After every step the candidates are deduplicated, since
different move sequences can sometimes reach the same mechanism. Both functions
are internal, but understanding them explains exactly which equations the search
can reach.

## `EnzymeRates.init_mechanisms`

`EnzymeRates.init_mechanisms(reaction)` returns the minimum-parameter mechanisms
for a reaction. For each catalytic topology it enumerates all substrate/product
dead-end subsets, assigns one steady-state catalytic step, and collapses binding
steps that share the same `(metabolite, RE/SS)` class into a single kinetic
group — the lowest-parameter starting point for the rate-equation search. When a
reaction's atom inventory allows a covalent fragment to persist between
half-reactions, it also builds ping-pong mechanisms: the modified enzyme stays
on conformation `:E` carrying a residual rather than a separate conformation
label, and a step that would return the enzyme to free `:E` with an empty
residual mid-cycle is rejected, since it would split the reaction into two
disconnected half-cycles. See [Ping-pong mechanisms](@ref) for detail on these
mechanisms.

## `EnzymeRates.seed_mechanisms`

For a reaction that declares regulators, `init_mechanisms` starts far below the
useful region: every seed binds zero regulators, so the search must climb through
every partially-regulated mechanism before it reaches one that binds them all. When
the effectors are already known — phosphofructokinase and its five regulators, say —
that lower shelf is pure waste.

`EnzymeRates.seed_mechanisms(reaction, required_allosteric, required_competitive)`
starts the search where it belongs. It grows `init_mechanisms` through the
regulator-binding moves alone — go allosteric, add an allosteric regulator, add a
competitive inhibitor — under three constraints that keep the set small: each
required regulator binds at its own single-ligand site, every allosteric state stays
cheap (`:OnlyA`/`:OnlyI`, never `:NonequalAI`), and no partially-regulated mechanism
survives. The result is every fully-regulated mechanism at its minimum parameter
count, and nothing beneath it. The beam then refines these seeds with the detail
moves — steady-state flips, splits, and `:NonequalAI` relaxations — as usual.

By default every declared regulator is required. `identify_rate_equation`'s
`optional_allosteric_regulators` and `optional_competitive_inhibitors` keywords move
named regulators back to optional, so the beam adds them as refinements rather than
forcing them into every seed; listing every regulator as optional recovers the
`init_mechanisms` starting set. Declaring a regulator's type in the reaction —
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
rapid-equilibrium segment in two. On a starting mechanism every group cuts on its
own, so each RE group gives one child. Once splits have separated a metabolite's
binding steps, a single group may be bridged by an RE route through the others;
the move then flips the bridging groups together, because a steady-state step
whose two ends stay in one equilibrated segment never reaches the rate equation.

Two kinds of group never flip. Competitive-inhibitor binding stays at rapid
equilibrium by modeling choice. A group whose steps all lie on dead-end branches
carries no net flux at steady state, so the equation could only ever see its
equilibrium ratio; flipping it would add a parameter the data cannot determine.

**Parameter delta:** +1 per flipped group in most cases. A few flips give an
equation that is the parent's up to renaming the constants, even though they
divide a segment; uni-uni is the classic case, where the steady-state and
rapid-equilibrium laws have the same form. Those children are fit once and lose
to their parent on parsimony.

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

**Parameter delta:** at least +1 by construction. A child that would add nothing
is never emitted.

### 3. Add a competitive inhibitor binding site

Adds binding steps for a `CompetitiveInhibitor` declared in the reaction, as a
dead-end complex with the enzyme. The move enumerates the combinations of enzyme
species the inhibitor can bind, emitting one child mechanism per combination,
subject to two rules:

- **Binding capacity.** An enzyme form holds at most as many metabolites as the
  larger of the substrate and product counts, `max(#substrates, #products)`, so
  the inhibitor is not added to a form already at capacity.
- **Mirror steps.** If the inhibitor binds two enzyme forms that a catalytic
  step already connects, a mirror step is added between the two inhibitor-bound
  forms, so the inhibitor-bound branch stays connected to the cycle. Each mirror
  inherits its counterpart's kinetic group and adds no parameter. Because the
  inhibitor competes with at least one substrate, the form that binds that
  substrate never carries the inhibitor, so the inhibitor-bound branch can
  never complete the net reaction. In a ping-pong mechanism it can carry out
  the one half-reaction whose ligands the inhibitor does not compete with.

The inhibitor's own binding steps form one fresh kinetic group (one new
dissociation constant `K_R`).

**Parameter delta:** +1 per competitive-inhibitor group added.

### 4. Promote a non-allosteric to allosteric mechanism

Converts a `Mechanism` to an `AllostericMechanism` variant set. An `:OnlyA`
catalytic binding asserts `K_I → ∞`: the inactive conformation cannot bind
that metabolite, so it cannot complete the catalytic cycle. The MWC reading is a
**catalytically-dead** inactive conformation — every chemistry step `:OnlyA`
— that binds ligands but runs no chemistry. The engine emits, per
multiplicity: every non-empty subset of binding groups `:OnlyA`, each with all
chemical steps `:OnlyA` (a K-type mechanism, emitted bare — the bound
metabolite's concentration reveals `L`); plus the empty subset (all chemical
steps `:OnlyA`, all `:EqualAI` binding) paired with a declared allosteric regulator,
one variant per `(regulator, tag)` with `tag ∈ {:OnlyA, :OnlyI}` (a V-type
mechanism, where the regulator makes `L` identifiable). The all-`:EqualAI`
baseline is never emitted (the conformations would be identical and `L` would
cancel). A binding subset that leaves a binding-only Wegscheider cycle
unsatisfiable is dropped (see [Thermodynamic constraints of MWC equations](@ref)).
Enumeration runs over `allowed_catalytic_multiplicities`. No-op on an already
allosteric input.

In MWC terminology the A-state corresponds to the R-state and the I-state
to the T-state of the original Monod–Wyman–Changeux nomenclature; this
package uses A/I throughout.

**Parameter delta:** +1 for a K-type variant — the conformational equilibrium
constant `L` is the sole new parameter. +2 for a V-type variant — `L` plus the
paired regulator's binding constant. `:OnlyA` zeroes out the I-state and adds no
parameter of its own.

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

**Dead-end branches stay at equilibrium.** A group whose every step lies off the
catalytic cycle carries no net flux at steady state. Its forward and reverse rates
enter the equation only as their ratio, which the equilibrium form already has,
so the group never flips.

**Competitive-inhibitor binding stays at equilibrium.** An inhibitor bound to
two forms that a catalytic step connects carries flux through its mirror step,
so keeping its binding at rapid equilibrium is a modeling choice ("inhibitor
binding is fast"). What competition does decide is turnover: the inhibitor
competes with at least one substrate and one product, so the form that binds
the competing substrate never carries it, and no inhibitor-bound branch
completes the net reaction. A ping-pong mechanism is the one case with more
than one chemistry step; there an inhibitor that competes with the first
half-reaction's ligands can still let the modified enzyme run the second half
with the inhibitor bound, which is the two-site picture, and the inhibitor
must leave before the next cycle.

**Splits follow binding context.** A group is divided only by whether another
ligand is already bound, by conformation, or by covalent residual. Arbitrary
partitions are not tried, a group is never split three ways in one move, and
groups never merge. A partition that no sequence of context splits produces is
unreachable.

**A child is never a provable copy of its parent.** Both moves reject a child
whose equation can be shown to equal the parent's: a split the constraint solver
ties back, a flip that leaves the segment count unchanged, a flip of a dead-end
group. A few children whose equation is the parent's up to renaming the
constants survive, uni-uni flips among them; they cost one fit each and never
win selection.
