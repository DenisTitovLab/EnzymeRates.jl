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
different move sequences often reach the same mechanism. Both functions are
internal — not exported, so you never call them directly
([`identify_rate_equation`](@ref) drives the loop) — but understanding them
explains exactly which equations the search can reach.

Deduplication is `unique!` on the raw `Vector` — no separate canonicalization
pass. The `Mechanism` and `AllostericMechanism` constructors canonicalize step
order, group order, and regulatory-site order at construction time, so two
mechanistically identical structs built in any order compare equal and hash
identically. Structural duplicates are removed by exact equality.

## `EnzymeRates.init_mechanisms`

`EnzymeRates.init_mechanisms(reaction)` returns a `Vector{Mechanism}` at minimum
parameter count. For each catalytic topology, it enumerates all
substrate/product dead-end subsets, assigns one steady-state catalytic step,
and collapses binding steps that share the same `(metabolite, RE/SS)` class
into a single kinetic group. The result is the lowest-parameter starting point
for the beam search.

When a reaction's atom inventory allows a covalent fragment to persist between
half-reactions, `EnzymeRates.init_mechanisms` also builds ping-pong mechanisms. The modified
enzyme stays on conformation `:E` carrying a residual rather than a separate
conformation label, and a step that would return the enzyme to free `:E` with an
empty residual mid-cycle is rejected — it would split the reaction into two
disconnected half-cycles. See [Ping-pong mechanisms](@ref) for the resulting
rate-equation form.

## The six expansion moves

`EnzymeRates.expand_mechanisms(mechs, rxn)` applies all six moves to every input
mechanism and returns the children as a flat `Vector{Union{Mechanism,
AllostericMechanism}}`. Bucketing by parameter count is the caller's job.

Every child is asserted atom-conserving before being returned.

The six moves:

### 1. Flip a rapid-equilibrium group to steady state

Flips one entire kinetic group from rapid equilibrium (RE) to steady state
(SS), atomically — all steps in the group convert together. The move is a
no-op for a group already in SS.

**Parameter delta:** +1 for most groups (the SS reverse rate is a new
independent parameter). A Haldane/Wegscheider constraint can make the reverse
rate dependent on existing parameters, giving a net **+0**. This is precisely
why the search buckets by actual fitted-param count rather than assuming +1.

### 2. Give one step its own kinetic group

For each kinetic group with two or more steps, carves one step into a fresh
singleton group. The split step then has an independent rate constant rather
than sharing one with its former group members.

**Parameter delta:** +1 — a single previously-shared constant becomes two
independent ones.

### 3. Add a competitive inhibitor binding site

Adds binding steps for a `CompetitiveInhibitor` declared in the reaction.
The new steps form one fresh kinetic group (one new dissociation constant
`K_R`); mirror steps share their catalytic counterpart's kinetic group.

**Parameter delta:** +1 per competitive-inhibitor group added.

### 4. Promote a non-allosteric mechanism to MWC

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

## Actual-count bucketing and `max_param_count`

Exact fitted-parameter counts come from
`length(fitted_params(compile_mechanism(m)))`, evaluated per child during
fitting. The beam search buckets by this actual count, not by a structural
estimate. A child whose actual count exceeds `max_param_count` is dropped
before fitting — this caps search depth without bounding per-mechanism compile
cost.

## Loud failures and artifacts

A mechanism that throws during compilation or fitting becomes a `FitFailure`
carrying the exception text. Failures are never silently discarded: they appear
in the CSV rows with `retcode` and `error` columns populated. If every
mechanism in the base tier fails, the search re-raises the first exception —
an unsupported optimizer kwarg or a memory overflow surfaces immediately.

`save_dir` defaults to a dated results directory; the search writes:
- `initial_mechanisms.csv` — all base-tier fits plus any failures.
- `equation_search_iteration_N.csv` — one file per expansion iteration that
  produced at least one row.
- `progress.log` — one line per master-level stage, appended and
  flushed after each write so cluster job logs stay current.
