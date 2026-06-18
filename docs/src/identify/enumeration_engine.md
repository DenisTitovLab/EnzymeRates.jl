# The enumeration engine

The enumeration engine is the core of EnzymeRates.jl: it defines the *universe*
of enzyme mechanisms that [`identify_rate_equation`](@ref) can sample. The
engine starts from the simplest mechanisms and applies a fixed set of moves that
each turn one mechanism into a slightly more complex one, so the search visits
candidates in order of increasing complexity. This reachability is decisive — a
rate equation the engine cannot construct is never fit, and so can never be
selected, however well it would describe the data. Knowing what the engine can
and cannot build is knowing what the package can find.

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
group — the lowest-parameter starting point for the beam search. When a
reaction's atom inventory allows a covalent fragment to persist between
half-reactions, it also builds ping-pong mechanisms: the modified enzyme stays
on conformation `:E` carrying a residual rather than a separate conformation
label, and a step that would return the enzyme to free `:E` with an empty
residual mid-cycle is rejected, since it would split the reaction into two
disconnected half-cycles. See [Ping-pong mechanisms](@ref) for detail on these
mechanisms.

## The six expansion moves

`EnzymeRates.expand_mechanisms(mechs, rxn)` applies all six moves to every input
mechanism and returns the resulting child mechanisms pooled into one list. Each
move is applied in every applicable way, so one input mechanism yields many
children — every rapid-equilibrium group that can flip to steady state, every
step that can be split into its own group, every way an inhibitor can bind, and
so on. Every child is checked to conserve atoms before it is returned.

The six moves:

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

**Parameter delta:** +1 — a single previously-shared constant becomes two
independent ones.

### 3. Add a competitive inhibitor binding site

Adds binding steps for a `CompetitiveInhibitor` declared in the reaction, as a
dead-end complex with the enzyme. The move enumerates the combinations of enzyme
species the inhibitor can bind, emitting one child mechanism per combination.
The new steps form one fresh kinetic group (one new dissociation constant
`K_R`); mirror steps share their catalytic counterpart's kinetic group.

**Parameter delta:** +1 per competitive-inhibitor group added.

### 4. Promote a non-allosteric to allosteric mechanism

Converts a `Mechanism` to an `AllostericMechanism` variant set. The baseline
variant uses the all-`:EqualAI` state (all groups have the same binding
constants in the active A-state and inactive I-state); one additional variant
per group uses `:OnlyA` for that group (the group is zeroed in the I-state).
Enumeration runs over `allowed_catalytic_multiplicities`. No-op on an already
allosteric input.

In MWC terminology the A-state corresponds to the R-state and the I-state
to the T-state of the original Monod–Wyman–Changeux nomenclature; this
package uses A/I throughout.

**Parameter delta:** +1 — the conformational equilibrium constant `L` is the
sole new parameter. `:OnlyA` variants zero out the I-state, not adding
parameters.

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
