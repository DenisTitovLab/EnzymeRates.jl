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
rapid-equilibrium step to steady state, adding a regulator, making the enzyme
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
children — every rapid-equilibrium group that can flip to steady state, every
step that can be split into its own group, every way an inhibitor can bind, and
so on. Every child is checked to conserve atoms before it is returned.

The seven moves:

### 1. Flip a rapid-equilibrium group to steady state

Flips one entire kinetic group from rapid equilibrium (RE) to steady state
(SS), atomically — all steps in the group convert together. The move is a
no-op for a group already in SS.

**Parameter delta:** +1 for most groups (the SS reverse rate is a new
independent parameter). A Haldane/Wegscheider constraint can make the reverse
rate dependent on existing parameters, giving a net **+0** — which is why the
search counts each mechanism's actual fitted parameters rather than assuming a
fixed +1 per move.

### 2. Give one step its own kinetic group

For each kinetic group with two or more steps, carves one step into a fresh
singleton group. The split step then has an independent rate constant rather
than sharing one with its former group members.

**Parameter delta:** +1 in the simplest case — the split-off constant becomes
independent. As with the RE→SS move, a Haldane/Wegscheider constraint can make
that constant dependent instead, giving a net **+0**; the search counts each
mechanism's actual fitted parameters rather than assuming +1.

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
  inherits its catalytic counterpart's kinetic group and adds no parameter.

The inhibitor's own binding steps form one fresh kinetic group (one new
dissociation constant `K_R`).

**Parameter delta:** +1 per competitive-inhibitor group added.

### 4. Promote a non-allosteric to allosteric mechanism

Converts a `Mechanism` to an `AllostericMechanism` variant set. An `:OnlyA`
catalytic binding asserts `K_I → ∞`: the inactive conformation cannot bind
that metabolite, so it cannot complete the catalytic cycle. The MWC reading is a
**catalytically-dead** inactive conformation — every isomerization (chemical)
step `:OnlyA` — that binds ligands but runs no chemistry. The engine emits, per
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
