# The enumeration engine

The mechanism search is built from three composable functions, not a
monolithic pipeline. The caller owns the expansion loop and the beam; the
enumeration functions produce the candidates.

```
init_mechanisms(rxn)           → Vector{Mechanism}
expand_mechanisms(mechs, rxn)  → Vector{Union{Mechanism, AllostericMechanism}}
unique!(children)              → same vector, structurally deduped
```

Deduplication is `unique!` on the raw `Vector` — no separate canonicalization
pass. The `Mechanism` and `AllostericMechanism` constructors canonicalize step
order, group order, and regulatory-site order at construction time, so two
mechanistically identical structs built in any order compare equal and hash
identically. `_dedup_flat!` mentioned in some internal specs is `unique!`:
the name was proposed but never added to the source; the live code calls
`unique!` directly.

## `init_mechanisms`

`init_mechanisms(reaction)` returns a `Vector{Mechanism}` at minimum
parameter count. For each catalytic topology from `_catalytic_topologies`, it
enumerates all substrate/product dead-end subsets, assigns one steady-state
catalytic step, and collapses binding steps that share the same `(metabolite,
RE/SS)` class into a single kinetic group. The result is the lowest-parameter
starting point for the beam search.

## The six expansion moves

`expand_mechanisms(mechs, rxn)` applies all six moves to every input
mechanism and returns the children as a flat `Vector{Union{Mechanism,
AllostericMechanism}}`. Bucketing by parameter count is the caller's job.

Every child is asserted atom-conserving before being returned.

The six moves, dispatched in `_add_expansions_mech!` order:

### 1. `_expand_re_to_ss` — flip a rapid-equilibrium group to steady state

Flips one entire kinetic group from rapid equilibrium (RE) to steady state
(SS), atomically — all steps in the group convert together. The move is a
no-op for a group already in SS.

**Parameter delta:** +1 for most groups (the SS reverse rate is a new
independent parameter). A Haldane/Wegscheider constraint can make the reverse
rate dependent on existing parameters, giving a net **+0**. This is precisely
why the search buckets by actual fitted-param count rather than assuming +1.

### 2. `_expand_split_kinetic_group` — give one step its own kinetic group

For each kinetic group with two or more steps, carves one step into a fresh
singleton group. The split step then has an independent rate constant rather
than sharing one with its former group members.

**Parameter delta:** +1 — a single previously-shared constant becomes two
independent ones.

### 3. `_expand_add_dead_end_regulator` — add a competitive inhibitor binding site

Adds binding steps for a `CompetitiveInhibitor` declared in the reaction.
The new steps form one fresh kinetic group (one new dissociation constant
`K_R`); mirror steps share their catalytic counterpart's kinetic group.

**Parameter delta:** +1 per competitive-inhibitor group added.

### 4. `_expand_to_allosteric` — promote a non-allosteric mechanism to MWC

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

### 5. `_expand_add_allosteric_regulator` — add an allosteric ligand

Adds one `AllostericRegulator` at a new or an existing regulatory site, with
an allosteric-state tag drawn from `{:OnlyA, :OnlyI, :NonequalAI}` (plus
`:EqualAI` only at a mixed existing site). No-op on a non-allosteric input.

**Parameter delta:**
- `:OnlyA` or `:OnlyI` tag: **+1** (one binding constant `K_A` or `K_I`).
- `:NonequalAI` tag: **+2** (independent `K_A` and `K_I`).

### 6. `_expand_change_allo_state` — relax an `:EqualAI` group to independent A/I states

Changes one `:EqualAI` or `:OnlyA`/`:OnlyI` allosteric-state tag to
`:NonequalAI`, giving the active-state and inactive-state versions of the
group independent parameters. No-op on a non-allosteric input.

**Parameter delta:**
- RE binding group: **+1** (one shared K splits into `K_A` and `K_I`).
- SS rate group: **+2** (each SS rate `kf` and `kr` splits into A/I pairs,
  adding two new independent constants).

## Actual-count bucketing and `max_param_count`

Exact fitted-parameter counts come from
`length(fitted_params(compile_mechanism(m)))`, called per child inside
`_process_batch`. The beam search buckets by this actual count, not by a
structural estimate. A child whose actual count exceeds `max_param_count` is
dropped before fitting — this caps search depth without bounding per-mechanism
compile cost.

## Loud failures and artifacts

A mechanism that throws during compilation or fitting becomes a `FitFailure`
carrying the exception text. Failures are never silently discarded: they appear
in the CSV rows with `retcode` and `error` columns populated. If every
mechanism in the base tier fails, the search re-raises the first exception —
an unsupported optimizer kwarg or a memory overflow surfaces immediately.

`save_dir` is mandatory. The search writes:
- `initial_mechanisms.csv` — all base-tier fits plus any failures.
- `equation_search_iteration_N.csv` — one file per expansion iteration that
  produced at least one row.
- `progress.log` — one line per master-level stage, appended and
  flushed after each write so cluster job logs stay current.

## Enumeration in practice

`init_mechanisms` produces the minimum-parameter mechanisms for a reaction:

```@example enum
using EnzymeRates

rxn = @enzyme_reaction begin
    substrates: S[C]
    products:   P[C]
end

mechs = EnzymeRates.init_mechanisms(rxn)
(count = length(mechs), eltype = eltype(mechs), nonempty = !isempty(mechs))
```

Applying the expansion moves grows the candidate set; `expand_mechanisms`
returns a flat vector of `Mechanism` / `AllostericMechanism`, which `unique!`
collapses to the structurally distinct ones:

```@example enum
children = EnzymeRates.expand_mechanisms(mechs, rxn)
unique!(children)
(eltype = eltype(children), nonempty = !isempty(children))
```

The result is a non-empty `Union{Mechanism, AllostericMechanism}` vector. The
exact counts vary with the enumeration rules; the invariants — non-empty,
correct eltype, no structural duplicates — do not.
