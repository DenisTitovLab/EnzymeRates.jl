# Refactor EnzymeRates to concrete types instead of Symbols — design

**Branch:** `refactor-to-concrete-types-instead-of-symbols`
**Status:** spec; pending implementation plan
**Date:** 2026-05-20

## 1. Motivation

The current EnzymeRates codebase has two structural problems that compound
each other:

1. **Two parallel data worlds.** Enumeration uses regular structs
   (`MechanismSpec`, `AllostericMechanismSpec`, `StepSpec`) with
   `Vector`/`Dict` fields. Derivation uses singleton-typed
   `EnzymeMechanism{Metabolites, Reactions}` whose data lives entirely as
   nested tuples of `Symbol` and `Int` in type parameters. A converter
   bridges them. The same domain concepts (a step, a form, a regulator,
   a parameter) are expressed twice in code.

2. **Pervasive Symbol-string dispatch.** Enzyme form names are built by
   string concatenation (`_form_name`) and parsed back by string
   splitting (`_parse_bound`). Parameter names are built positionally
   (`Symbol("K$idx")`, `Symbol("k$(idx)f")`). Parameters are classified
   by string prefix (`is_k_parameter` = `startswith(string(sym), "k") &&
   isdigit(s[2])`). T-state mirrors are built by string append
   (`* "_T"`). Regulator K names by concatenation
   (`K_<lig>_reg<i>`). The canonical-equation hash in
   `identify_rate_equation.jl` regex-walks the printed rate-equation
   string with multiple passes and an ad-hoc multiplicative-run sorter.
   Every one of these is a fragile workaround for not having structured
   data.

This refactor replaces both with **one concrete struct family** shared
between enumeration and derivation. Every Symbol that appears in the
public API is rendered deterministically from struct fields via a single
chokepoint accessor; nothing parses Symbols back.

A previous attempt (`refactor-to-use-structs-throughout`) tried the same
goal but went off the rails: 24k LOC added, parallel NonSingleton +
Singleton representations everywhere, never merged. The lessons that
shaped this spec:

- **No parallel representations** — replace, don't augment. Each commit
  switches one chunk; old types do not coexist with new types for the
  same purpose.
- **Tests adapt synchronously** — same commit. Never "tests broken on
  this branch, will fix next stage".
- **Code size is the success metric** — see §3.

## 2. Test integrity (NON-NEGOTIABLE)

Under NO circumstances may an existing test be deleted, commented out,
or have its assertion logic changed. This is how the previous attempt
failed — tests got dropped piecemeal and the package lost its safety
net. This rule applies to every commit on every stage.

The ONLY change permitted to a test is **mechanical syntax adaptation**
to the new struct surface:

| Old | New (mechanical adaptation) |
|---|---|
| `(:E, (:S,), :ES, false, 1)` tuple literal | `Step(Species([], :E), Species([Substrate(:S)], :E), Substrate(:S), false)` constructor call |
| `m.reactions[1][1]` (tuple indexing) | `m.steps[1][1].from_species` (accessor) |
| `:K1` Symbol literal compared to `parameters(m)` output | unchanged — `parameters(m)` still returns the same Symbol |
| `:K_ATP_reg1` Symbol built by hand in a test fixture | `name(Kreg(site, ligand, :None), m)` — only IF the test was asserting the rendered name |

**Forbidden:**

- Removing a `@testset` block
- Removing a `@test` line
- Changing a hardcoded numerical value (e.g., a kcat or rate in a fixture)
- Changing a hardcoded Symbol value the test compares against (e.g.,
  `:K1`, `:E_total`, the strings inside `MECHANISM_TEST_SPECS`)
- "Skipping" a test with `@test_skip` or by commenting it out
- Adding `@test_broken`

If a test fails after a mechanical adaptation, the test is telling you
the refactor regressed. STOP and fix the src. Do not change the test.

Stage CI passes ⇔ all tests green, no `@test_broken`, no `@test_skip`,
no commented-out tests, full assertion strength preserved.

### § 2.1 NARROW EXCEPTION — tests of deleted helpers

The ONE permitted form of test deletion is when a `src/` helper function
is itself deleted by the refactor AND that helper had its own dedicated
unit tests that no longer make sense without the helper. Example:
`_is_ss_rate_constant(sym::Symbol)` is a Symbol-string classifier;
when Parameter struct dispatch replaces it (Stage 3.4), the function
is deleted and the `@testset "_is_ss_rate_constant"` block in
`test/test_rate_eq_derivation.jl` (~lines 799–805) has no surviving
referent.

**Conditions for invoking this exception:**

1. The helper being deleted is private (`_` prefix) AND only consumed
   by the test in question + the src code being rewritten in the same
   stage.
2. The behavior being tested is preserved by the replacement — i.e.,
   if `_is_ss_rate_constant` answered "is this Symbol a steady-state
   rate constant?", the new code answers the equivalent structural
   question via Parameter dispatch (`p isa Kon || p isa Koff || p isa
   Kfor || p isa Krev`). The downstream behavior the helper enabled
   is still tested via the integration tests that depend on it
   (rate_equation correctness, parameters output, etc.).
3. **Every deletion is recorded in the per-PR deleted-tests log
   (below) with reason and replacement evidence.**

**Deleted-tests log — REQUIRED.**

The refactor PR must include a top-level file
`docs/superpowers/refactor-deleted-tests.md` listing every test
deleted under this exception, with for each:

- Test file + testset name + commit SHA that deleted it
- Deleted helper(s) the testset covered
- The replacement code path that subsumes the same behavior
- Integration-test names that exercise the replacement

The log is reviewable evidence that no test deletion was a silent
weakening. Reviewers check the log at PR review time alongside the
diff. Format:

```markdown
## Stage 3.4 (commit abc1234)

### test_rate_eq_derivation.jl `@testset "_is_ss_rate_constant"`
- Lines: 799–805 (also referenced 664, 677)
- Helper deleted: `_is_ss_rate_constant(sym::Symbol)`
- Replacement: `p isa Union{Kon,Koff,Kfor,Krev}` struct dispatch
  in `_kcat_groups_from_polys` and downstream.
- Integration coverage:
  - `test_rate_eq_derivation.jl` `@testset "_kcat_forward correctness"`
  - `test_rate_eq_derivation.jl` `test_rate_equation_performance`
  (the kcat computation that depended on the classifier still works)
```

This exception is NARROW. If a deletion doesn't fit the three conditions
above, it falls back under the unconditional no-deletion rule and is
forbidden.

## 3. Code-size goal (LOAD-BEARING — the refactor's primary objective)

Reducing src LOC is the **measure of success** for this refactor,
alongside test correctness and perf gates. Per
[[feedback-simplification-means-less-code]]: in this codebase, less code
IS simpler code.

**Baseline (main, `wc -l src/*.jl`):**

```
types.jl                                       744
dsl.jl                                         679
rate_eq_derivation.jl                        1,460
mechanism_enumeration.jl                     2,434
thermodynamic_constr_for_rate_eq_derivation    365
sym_poly_for_rate_eq_derivation                322
fitting.jl                                     211
identify_rate_equation.jl                      884
EnzymeRates.jl                                  37
                                            ──────
                                             7,136
```

**Target: ≤ 3,600 src LOC (≥ 50% reduction).** This is the headline
goal. If we land with <30% reduction we have failed the simplification
mission regardless of how clean the type hierarchy looks.

**Per-stage LOC tracking — mandatory:**

- Every stage commit message ends with one line:
  `src delta: -X / +Y net Z, cumulative: -W`.
- Computed as `wc -l src/*.jl` before and after the commit.
- **Negative Z is the default expectation.** Stages 1–2 (foundations +
  DSL) may carry small positive Z while the new type defs land; from
  Stage 3 onward Z must be negative.
- Positive Z requires explicit justification in the commit body
  (`"+150 because new struct definitions; offset by deletions in stage 4"`).
- **Mid-refactor checkpoint:** after Stage 4 the cumulative delta must
  be at least −500. If not net-negative by half-way, STOP and redesign.
- **Final gate:** Stage 7 cleanup must drive cumulative to ≤ −3,500.

**Implementor discipline:**

After every commit, re-read each changed `src/*.jl` file end-to-end.
Hunt for:

- Dead helpers (functions never called after this stage)
- Redundant accessors (two functions computing the same thing)
- Symbol-juggling that snuck back in (any `Symbol("X$idx")`, any
  `startswith(string(sym), ...)`, any string concat building a parameter
  name — should appear ONLY inside `name(::Parameter, m)` body)
- Single-use helpers (inline them — except the chokepoint accessors
  protected by [[feedback-chokepoint-accessors-for-future-migrations]])
- Comments explaining "what" rather than "why" (delete; well-named
  identifiers self-document)

**Where the LOC budget comes from (the kill list):**

| Region | Now | After | Δ |
|---|---|---|---|
| `mechanism_enumeration.jl` Symbol/form-name helpers + canonicalizer regex | ~1,100 | ~150 | −950 |
| `rate_eq_derivation.jl` Symbol classifiers + rename helpers + T-state renamers + reg-name builders | ~600 | ~100 | −500 |
| `thermodynamic_constr_for_rate_eq_derivation.jl` `Symbol("K$idx")` builders | ~150 | ~50 | −100 |
| Bridging code (`MechanismSpec` ↔ `EnzymeMechanism` converters) | ~200 | ~80 | −120 |
| `identify_rate_equation.jl` regex canonicalizer + factor sort | ~200 | ~50 | −150 |
| New struct definitions in `types.jl` | n/a | +400 | +400 |
| **Net** | **~7,136** | **~3,600** | **−3,500** |

Rough but honest sizing.

## 4. Goals and non-goals

**Goals:**

- One concrete struct family shared between enumeration and derivation.
- All operations go through accessors; no Symbol-string parsing in any
  src file.
- Code-size reduction per §3.
- `rate_equation` perf invariant preserved (`allocs == 0`, `t < 100ns`
  per call for every mechanism in `MECHANISM_TEST_SPECS`).
- `loss!` perf preserved (no regression vs main baseline).
- Compile-time regression gate added so future changes can't silently
  blow up compilation (trace-compile budget + macro wall-clock).
- DSL macros adopt the prior branch's cleaner grammar while emitting
  the new concrete structs directly.

**Non-goals:**

- New fitting algorithms or enumeration strategies.
- Performance changes other than preserving the existing invariants.
- Path-based parameter names (`:K_ATP_from_E`). Today's positional names
  (`:K1`, `:k1f`, `:K1_T`, `:K_<lig>_reg<i>`) are preserved, BUT all
  Symbol production routes through the `name(p::Parameter, m)`
  chokepoint accessor. A future refactor can switch to path-based names
  by changing that one function body.

## 5. Type hierarchy

All types are concrete structs (immutable). Two parametric types
(`EnzymeMechanism{Sig}`, `AllostericEnzymeMechanism{Sig}`) are the
**only** parametric structs in the design — everything else is a plain
concrete struct used both in enumeration's `Vector{T}` storage and at
@generated-body-build time when derivation lifts the Sig.

### 5.1 Metabolite hierarchy

```julia
abstract type Metabolite end

# Reactant = things that participate in the net reaction.
abstract type Reactant <: Metabolite end
struct Substrate <: Reactant; name::Symbol end
struct Product   <: Reactant; name::Symbol end

# Regulator = things that bind but don't appear in net stoichiometry.
abstract type Regulator <: Metabolite end
struct AllostericRegulator  <: Regulator; name::Symbol end
struct CompetitiveInhibitor <: Regulator; name::Symbol end
```

### 5.2 Residual

```julia
struct Residual
    added::Vector{Substrate}      # molecules added to residual
    subtracted::Vector{Product}   # molecules released from residual
end
Residual() = Residual(Substrate[], Product[])
```

Empty `Residual()` means no covalent adduct on the enzyme. Residual is
orthogonal to conformation — the `:Estar` vs `:E` label lives on
`Species.conformation`.

### 5.3 Species

```julia
struct Species
    bound::Vector{Metabolite}     # canonically sorted by name
    conformation::Symbol          # :E or :Estar (extensible)
    residual::Residual            # empty == no covalent adduct
end
Species(bound, conformation) = Species(bound, conformation, Residual())
```

### 5.4 Step

```julia
struct Step
    from_species::Species
    to_species::Species
    bound_metabolite::Union{Metabolite,Nothing}   # nothing for iso steps
    is_equilibrium::Bool
end
```

No `kinetic_group::Int` field — kinetic-group sharing is **structural**,
expressed by grouping `Step`s into inner `Vector{Step}` at the mechanism
level (§5.7).

**Rep-idx convention (matches today's behavior exactly):** today's
`_canonicalize!` (`src/mechanism_enumeration.jl:2199-2208`) sorts steps
by structural key and re-numbers `kinetic_group` so groups are
contiguous 1..N in sorted order. Then today's
`rep = first(steps_in_group(cm, g))` (`src/rate_eq_derivation.jl:125`)
uses the post-canonicalized step POSITION within the group as `rep_idx`.

The new `name(p, m)` computes rep_idx the same way structurally:

```julia
function _rep_idx_for_step(step::Step, m::Mechanism)
    pos = 0
    for group in steps(m)
        if step in group
            return pos + findfirst(==(step), group)  # or just pos + 1 since
                                                     # rep = first step
        end
        pos += length(group)
    end
    error("Step not found in mechanism")
end
```

This matches today's `Symbol("K$rep")` output for every mechanism — no
behavior change, no `source_idx` field, no propagation burden through
every `_expand_*` move. The Sig encoding is also smaller because the
step encoding doesn't carry source_idx (one fewer leaf per step).

Canonical direction (enforced by the constructor):

- Binding/release steps (`bound_metabolite !== nothing`): metabolite is
  on the `from_species` side (i.e., the metabolite is bound during the
  step direction LHS→RHS). The constructor swaps if the user authored
  the other direction.
- Iso steps (`bound_metabolite === nothing`): direction picked by
  substrate-content → product-content → lex on `name(from_species)`.

### 5.5 RegulatorySite

```julia
struct RegulatorySite
    ligands::Vector{AllostericRegulator}
    multiplicity::Int
    allo_states::Vector{Symbol}    # one per ligand, parallel to ligands
end
```

`allo_states[i]` is the allosteric state (`:OnlyR`, `:OnlyT`,
`:EqualRT`, `:NonequalRT`) of `ligands[i]` at this site.

### 5.6 Parameter

```julia
abstract type Parameter end

# Step-bound RE parameters.
struct Kd   <: Parameter; step::Step; state::Symbol end   # RE binding
struct Kiso <: Parameter; step::Step; state::Symbol end   # RE iso

# Step-bound SS parameters.
struct Kon  <: Parameter; step::Step; state::Symbol end   # SS binding forward
struct Koff <: Parameter; step::Step; state::Symbol end   # SS binding reverse
struct Kfor <: Parameter; step::Step; state::Symbol end   # SS iso forward
struct Krev <: Parameter; step::Step; state::Symbol end   # SS iso reverse

# Regulator-site parameter (needs both site and ligand because the
# same ligand can appear at multiple sites and the same site can
# have multiple ligands).
struct Kreg <: Parameter
    site::RegulatorySite
    ligand::AllostericRegulator
    state::Symbol
end

# Mechanism-level scalars.
struct Keq   <: Parameter end
struct Etot  <: Parameter end
struct Lallo <: Parameter end
```

`state::Symbol` is one of `:R`, `:T`, `:None`:

- `:None` — non-allosteric step (`:NoneState`) or `:EqualRT` (shared
  parameter)
- `:R` — R-state branch of an `:OnlyR` or `:NonequalRT` step
- `:T` — T-state branch of an `:OnlyT` or `:NonequalRT` step

### 5.7 EnzymeReaction (bundled per-reactant and per-regulator data)

```julia
struct ReactantAtoms
    metabolite::Reactant                            # Substrate or Product
    atoms::Vector{Pair{Symbol,Int}}
end

struct RegulatorMults
    regulator::Regulator
    allowed_multiplicities::Vector{Int}
end

struct EnzymeReaction
    reactants::Vector{ReactantAtoms}                # substrates + products unified
    regulators::Vector{RegulatorMults}
    allowed_catalytic_multiplicities::Vector{Int}
end
```

Substrates and products live in the same `reactants` field;
`substrates(r)` / `products(r)` are derived filters by metabolite
subtype.

### 5.8 Mechanism (non-parametric) and EnzymeMechanism (parametric)

**Non-parametric form — what enumeration produces.** Single Julia type;
instantiating millions of these costs nothing in compile time.

```julia
struct Mechanism
    reaction::EnzymeReaction
    steps::Vector{Vector{Step}}             # structural kinetic groups
end

struct AllostericMechanism
    reaction::EnzymeReaction
    cat_steps::Vector{Vector{Step}}         # catalytic-cycle steps
    cat_allo_states::Vector{Symbol}         # one per kinetic group; parallel to cat_steps;
                                            # values in (:OnlyR, :EqualRT, :NonequalRT) —
                                            # :OnlyT rejected by constructor (R-state-active)
    catalytic_multiplicity::Int
    regulatory_sites::Vector{RegulatorySite}
end
```

**Parametric form — what derivation dispatches on.** Each unique `Sig`
triggers one `@generated rate_equation` body-build. Only created for
mechanisms we actually need fast `rate_equation` evaluation for.

```julia
abstract type AbstractEnzymeMechanism end
struct EnzymeMechanism{Sig}            <: AbstractEnzymeMechanism end
struct AllostericEnzymeMechanism{Sig}  <: AbstractEnzymeMechanism end
```

`Sig` is a minimal hashable tuple-of-tuples-of-`Symbol`s and `Int`s
encoding the `Mechanism`'s data. The shape is internal — users never
touch it. `_sig_of(m::Mechanism)` and `_mechanism_from_sig(Sig::Tuple)`
are the conversion functions.

**Sig shape (concrete form):** initially `Sig === (reaction_sig::Tuple,
steps_sig::Tuple)` where `reaction_sig` and `steps_sig` are
tuple-of-tuples-of-Symbols/Ints encodings of `EnzymeReaction` and
`Vector{Vector{Step}}` respectively. The exact internal layout may be
tightened in Stage 7 cleanup if profiling reveals waste; the
`_mechanism_from_sig` / `_sig_of` pair is the boundary of this
internal-only contract.

**Repacking the type-parameter signature.** Today's `EnzymeMechanism`
takes two type parameters (`{Metabolites, Reactions}`). Stage 1
collapses that to one (`{Sig}`) where `Sig === (Metabolites, Reactions)`
content-equivalent — this is a mechanical change across every src
where-clause (`where {M, R}` → `where {Sig}` plus
`(reaction_sig, steps_sig) = Sig` destructuring inside body). Zero
behavior change; pure spelling refactor that gives later stages a
single-Sig API to evolve.

**Bidirectional conversion at the boundary:**

```julia
EnzymeMechanism(m::Mechanism)              = EnzymeMechanism{_sig_of(m)}()
Mechanism(::EnzymeMechanism{Sig}) where Sig = _mechanism_from_sig(Sig)
```

`_mechanism_from_sig` is a **regular function** (not `@generated`).
It's called from inside the `@generated rate_equation` body at
body-build time. Nested `@generated` would risk world-age issues; a
plain function called from inside a `@generated` body is fine because
body-build runs in normal eval scope.

```julia
@generated function rate_equation(
    ::EnzymeMechanism{Sig}, concs::NamedTuple, params::NamedTuple,
) where {Sig}
    mech = _mechanism_from_sig(Sig)        # plain call at body-build
    _build_rate_body(mech)                 # returns arithmetic Expr
end
```

## 6. Accessor surface

The "no Symbol-string parsing anywhere" rule means every operation goes
through these accessors. Per
[[feedback-chokepoint-accessors-for-future-migrations]], the
`name(p::Parameter, m)` function is the single chokepoint for parameter
Symbol production — a future refactor to path-based names is a
single-function edit.

```julia
# Metabolite / Reactant / Regulator
name(::Substrate)             :: Symbol
name(::Product)               :: Symbol
name(::AllostericRegulator)   :: Symbol
name(::CompetitiveInhibitor)  :: Symbol

# Residual
added(::Residual)             :: Vector{Substrate}
subtracted(::Residual)        :: Vector{Product}
isempty(::Residual)           :: Bool

# Species
bound(s::Species)             :: Vector{Metabolite}
conformation(s::Species)      :: Symbol
residual(s::Species)          :: Residual
has_residual(s::Species)      :: Bool             # !isempty(residual(s))
name(s::Species)              :: Symbol           # rendered from fields

# Step
from_species(::Step)          :: Species
to_species(::Step)            :: Species
bound_metabolite(::Step)      :: Union{Metabolite,Nothing}
is_equilibrium(::Step)        :: Bool
is_binding(::Step)            :: Bool             # bound_metabolite !== nothing
is_iso(::Step)                :: Bool             # bound_metabolite === nothing
direction(::Step)             :: Symbol           # :binding or :iso

# Parameter — THE CHOKEPOINT
name(p::Parameter, m::AbstractEnzymeMechanism) :: Symbol
governing_step(p::Parameter)  :: Step             # for step-bound params
is_t_state(p::Parameter)      :: Bool             # state === :T

# RegulatorySite
ligands(::RegulatorySite)             :: Vector{AllostericRegulator}
multiplicity(::RegulatorySite)        :: Int
allo_states(::RegulatorySite)         :: Vector{Symbol}

# EnzymeReaction
reactants(::EnzymeReaction)                       :: Vector{ReactantAtoms}
substrates(::EnzymeReaction)                      :: Vector{Substrate}
products(::EnzymeReaction)                        :: Vector{Product}
regulators(::EnzymeReaction)                      :: Vector{RegulatorMults}
allowed_catalytic_multiplicities(::EnzymeReaction) :: Vector{Int}

# Bundling-struct accessors
metabolite(::ReactantAtoms)               :: Reactant
atoms(::ReactantAtoms)                    :: Vector{Pair{Symbol,Int}}
regulator(::RegulatorMults)               :: Regulator
allowed_multiplicities(::RegulatorMults)  :: Vector{Int}

# Mechanism (both parametric and non-parametric forms)
reaction(::Mechanism)            :: EnzymeReaction
reaction(::EnzymeMechanism)      :: EnzymeReaction        # via Mechanism(m)
steps(::Mechanism)               :: Vector{Vector{Step}}
steps(::EnzymeMechanism)         :: Vector{Vector{Step}}  # via Mechanism(m)
kinetic_groups(m)                :: UnitRange{Int}        # 1:length(steps(m))
rep_step(m, g::Int)              :: Step                  # first(steps(m)[g])
n_steps(m)                       :: Int                   # sum(length, steps(m))

# AllostericMechanism / AllostericEnzymeMechanism — similar
catalytic_mechanism(m::AllostericMechanism)   :: Mechanism
catalytic_multiplicity(m)                     :: Int
regulatory_sites(m)                           :: Vector{RegulatorySite}
cat_allo_state(m::AllostericMechanism, g::Int) :: Symbol   # per kinetic group
allosteric_regulators(m)                      :: Vector{AllostericRegulator}
competitive_inhibitors(m)                     :: Vector{CompetitiveInhibitor}
```

**Chokepoint scope (explicit).** `name(p::Parameter, m)` is THE Symbol-rendering site for
parameters in the rate equation Expr leaves, `rate_equation_string` display, and the
`parameters(m)`/`fitted_params(m)` user-facing outputs. The following sites legitimately
bypass it:

- **DSL parse-time** — `@enzyme_mechanism` etc. emit Symbol literals for declared
  metabolite names (`:S`, `:P`) directly from the user's source. These are concrete-struct
  field values, not Parameter renderings.
- **Canonical-hash internals** — `_canonical_rate_eq_hash` (spec §8.4) keys on the
  Parameter STRUCT (via `_parameter_canonical_key`) and NEVER calls `name(p, m)`. This is
  intentional: the canonical hash must be invariant to positional rep-idx, so it cannot
  use the positional-renamer output.
- **Metabolite name(::Substrate)** etc. — these are the metabolite-side chokepoint, not
  the parameter side. They share the function name `name` via dispatch but answer a
  different question.

Stage 7 cleanup audit must confirm NO `Symbol("K…")` / `Symbol("k…")` / `Symbol("k_…")`
string-concat building of parameter Symbols outside `name(p::Parameter, m)`.

`name(p::Parameter, m)` initial implementation (Stage 1):

```julia
# Today's positional naming routed through one function. Future path-based
# refactor: change these bodies only.
function name(p::Kd, m::AbstractEnzymeMechanism)
    rep = _rep_idx_for_step(p.step, m)
    p.state === :T ? Symbol("K$(rep)_T") : Symbol("K$rep")
end
function name(p::Kon, m::AbstractEnzymeMechanism)
    rep = _rep_idx_for_step(p.step, m)
    p.state === :T ? Symbol("k$(rep)f_T") : Symbol("k$(rep)f")
end
# ... Koff, Kfor, Krev, Kiso similarly
function name(p::Kreg, m::AllostericEnzymeMechanism)
    site_idx = _site_idx_of(p.site, m)
    lig_name = name(p.ligand)
    p.state === :T ? Symbol("K_$(lig_name)_T_reg$site_idx") :
                     Symbol("K_$(lig_name)_reg$site_idx")
end
name(::Keq, _)   = :Keq
name(::Etot, _)  = :E_total
name(::Lallo, _) = :L
```

`_rep_idx_for_step` and `_site_idx_of` look up position in
`steps(m)`/`regulatory_sites(m)`. Linear scan is fine — N is small;
called only at @generated body-build time, not per `rate_equation` call.

## 7. DSL grammar (adopted from prior branch)

The `@enzyme_reaction`, `@enzyme_mechanism`, and `@allosteric_mechanism`
macros adopt the prior branch's cleaner grammar but emit constructors
of the new concrete struct family directly (no NonSingleton↔Singleton
indirection).

### 7.1 `@enzyme_reaction`

```julia
@enzyme_reaction begin
    substrates: S[C6H12O6], ATP[C10H16N5O13P3]
    products:   P[C6H13O9P]
    competitive_inhibitors: I            # bare OK; mults match catalytic
    allosteric_regulators: A(1, 2, 4)    # per-reg mults required
    allowed_catalytic_multiplicities: (1, 2, 4)
    # or shorthand:
    # oligomeric_state: 2
end
# Emits: EnzymeReaction(
#   [ReactantAtoms(Substrate(:S), [:C=>6, :H=>12, :O=>6]), ...],
#   [RegulatorMults(CompetitiveInhibitor(:I), [1,2,4]),
#    RegulatorMults(AllostericRegulator(:A), [1,2,4])],
#   [1, 2, 4],
# )
```

### 7.2 `@enzyme_mechanism`

```julia
@enzyme_mechanism begin
    substrates: S
    products:   P
    regulators: I                        # all treated as CompetitiveInhibitor
    steps: begin
        E + S ⇌ E(S)                     # RE binding
        (E(S) ⇌ E(S), E(S) <--> E(S))    # parenthesized = shared kinetic group
        E(P) ⇌ E + P                     # RE release
    end
end
# Emits: Mechanism(reaction, [
#   [Step(E, ES, S, true)],
#   [Step(ES, ES, nothing, true), Step(ES, ES, nothing, false)],   # shared K
#   [Step(EP, E, P, true)],
# ])
```

Species notation:

- Bare `E` is `E()` (no bound, no residual)
- `E(S)`: E with S bound
- `Estar`: empty Estar
- `Estar(B)`: Estar with B bound
- `Estar(; residual = A - P)`: Estar with empty bound, residual = +A, −P
- `Estar(B; residual = A - P)`: Estar with B bound, residual = +A, −P

Conformation labels (`:E`, `:Estar`) come from the species call's
function-name part; cannot shadow declared metabolite names.

### 7.3 `@allosteric_mechanism`

```julia
@allosteric_mechanism begin
    substrates: F6P
    products:   F16BP
    catalytic_multiplicity: 2
    allosteric_regulators: A::OnlyR, I::OnlyT

    catalytic_steps: begin
        E + F6P ⇌ E(F6P)        :: EqualRT
        E(F6P) <--> E(F16BP)    :: EqualRT
        (E(F16BP) ⇌ E + F16BP)  :: EqualRT
    end

    regulatory_site(multiplicity = 4): begin
        ligands: A
    end
    regulatory_site(multiplicity = 4): begin
        ligands: I
    end
end
# Emits: AllostericMechanism(reaction, cat_steps, 2,
#   [RegulatorySite([AllostericRegulator(:A)], 4, [:OnlyR]),
#    RegulatorySite([AllostericRegulator(:I)], 4, [:OnlyT])])
```

## 8. Pipeline architecture

### 8.1 Enumeration

Today's `init_mechanisms` / `expand_mechanisms` / `dedup!` operate on
`MechanismSpec` / `StepSpec` (Vector-of-Symbol structs). After the
refactor they operate on `Mechanism` and `Vector{Vector{Step}}` using
the same Step/Species/Residual structs that derivation uses.

```julia
init_mechanisms(reaction::EnzymeReaction) :: Vector{Mechanism}
expand_mechanisms(specs::Vector{Mechanism}, reaction::EnzymeReaction) ::
    Dict{Int, Vector{Mechanism}}    # keyed by n_fit_params_estimate
dedup!(cache::Dict{Int, Vector{Mechanism}})
```

**Dedup is structural.** `Mechanism` gets `==`/`hash` via field equality
on `steps::Vector{Vector{Step}}` after canonical sort within and across
groups. No string canonicalization, no regex.

**Compile-time cost is bounded.** `Vector{Mechanism}` is one Julia type;
enumeration of 10⁴ mechanisms creates 10⁴ instances of that one type at
runtime — zero new compilations.

### 8.2 Conversion to derivation form

After enumeration + dedup + complexity filtering, each surviving
`Mechanism` is converted to `EnzymeMechanism{Sig}` ONLY when it needs
fast `rate_equation` evaluation (i.e., during fitting). This is where
the `@generated rate_equation` body builds the per-mechanism arithmetic
Expr.

```julia
# Called from identify_rate_equation per candidate:
em = EnzymeMechanism(m)                   # one @generated build
fp = FittingProblem(em, data; Keq=...)    # stores parametric form
loss!(x, fp)                              # hot loop; uses @generated rate_equation
```

Compile-time cost = number of mechanisms actually fitted (post-dedup
post-filter). Matches today.

### 8.3 Derivation

`@generated rate_equation`, `parameters`, `_dependent_param_exprs`,
`rate_equation_string`, `_kcat_forward`, and the AllostericEnzymeMechanism
counterparts all lift `Sig` once via `_mechanism_from_sig` at body-build
time, then walk the resulting `Mechanism`'s `steps::Vector{Vector{Step}}`
and access `Step.from_species`, `Step.to_species`, `Step.bound_metabolite`
directly. The body Expr they emit is unchanged in shape (still pure
arithmetic, still 0-alloc / <100ns per call); only the body-build code
changes from "walk Symbol tuples" to "walk Step/Species structs".

**Deleted Symbol-juggling helpers** (representative — exhaustive list
emerges during implementation):

- `_form_name`, `_parse_bound`, `_dead_end_form_name`, `_atoms_dict`,
  `_is_estar_form`, `_can_pingpong`, `_subtract_atoms` (form-name
  string building / parsing)
- `_split_reaction_side` (Symbol-based enzyme vs metabolite split)
- `is_k_parameter`, `_is_ss_rate_constant` (Symbol classifiers)
- `_rename_params_T`, `_reg_param_name`, `_T_rename`, `_onlyR_syms`,
  `_all_t_state_names`, `_group_param_symbols` (Symbol building /
  renaming for allosteric)
- `_build_kinetic_rename_map` (or vastly simplified)
- `_canonicalize_rate_eq_with_map`, `_sort_run_factors`,
  `_factor_sort_key`, the regex pipeline in
  `identify_rate_equation.jl`

### 8.4 Identify_rate_equation canonicalizer

Replace the printed-string regex pipeline with a structural hash. **The
new hasher MUST preserve the existing `(UInt64, name_map)` return shape**
because `src/identify_rate_equation.jl:117,424-425,465-476,492-494`
threads `name_map::Dict{String,String}` through `_CachedFitResult` and
uses it in `_project_cached_params(cached.params, canon_to_rep,
c.name_map, c.fitted_keys)` to relabel cached fit values across
equivalent mechanisms. Returning only `UInt64` would silently break
that param-projection pathway — the dedup cache would still hit but
return wrong-slot values (correctness regression) or force re-fitting
every equivalent mechanism (perf regression).

```julia
# Old: build rate_equation_string → regex over the string → (hash, hex, name_map).
# New: walk Mechanism's steps + parameters → normalize relabeling of
#      kinetic-group rep indices structurally → emit (canon, name_map);
#      derive the 3-tuple contract.
function _canonical_rate_eq_hash_data(m::AbstractEnzymeMechanism)
    mech = _to_mechanism(m)
    canon, name_map = _canonicalize_for_hash(mech)
    h = hash(canon)
    (h, string(h, base=16, pad=16), name_map)
end

function _canonical_rate_eq_hash(m::AbstractEnzymeMechanism)
    first(_canonical_rate_eq_hash_data(m))
end
```

**Return-arity contract: `(UInt64, String, Dict{String,String})`.** The
existing `_canonical_rate_eq_hash_data` (`src/mechanism_enumeration.jl:2422`)
returns 3-tuple and `src/identify_rate_equation.jl:424` destructures all
three (`h_full, h_short, name_map = ...`) into `_CachedFitResult` fields
where `h_short::String` (the 16-char hex display) surfaces as the
user-facing `eq_hash` column in the results DataFrame (line 503). The
new struct hash MUST preserve all three return slots — a 2-tuple
rewrite would either `BoundsError` or silently shred values at Stage 6
integration.

`name_map::Dict{String,String}` maps per-mechanism parameter Symbol names
(as Strings — e.g., `"K1"`) to canonical tokens (e.g., `"p_3"`) that are
stable across equivalent mechanisms. The Plan Task 6.1 implementation
proof obligation: two structurally-equivalent mechanisms must produce
`name_map`s assigning the same canonical token to bijection-matched
parameters, so `_project_cached_params` continues to project cached fit
values correctly across the equivalence class.

No regex. No string canonicalization. Dedup-equivalent mechanisms hash
to the same UInt64.

## 9. Perf gates (CI-enforced)

Three hard gates run in CI for every commit on the branch:

1. **`rate_equation` per-call**: `allocs == 0` AND `t < 100e-9` per
   call, for every mechanism in `MECHANISM_TEST_SPECS`. (Today's test;
   must stay green.)
2. **`loss!` per-call**: no regression vs main baseline. (New test;
   establishes baseline in the pre-PR.)
3. **Compile-time budget**: trace-compile primary + macro wall-clock
   secondary. Sub-budgets for `init_mechanisms`, `expand_mechanisms`,
   and `using EnzymeRates`. (Ships in the standalone pre-PR.)

Trace-compile test shape (cribbed from prior branch's working version):

```julia
@testset "trace-compile budget for init_mechanisms" begin
    cmd = Cmd([julia_exe, "--trace-compile=$trace_file",
               "--project=.", "-e", runner_script])
    run(cmd; wait=true)
    relevant = filter(occursin(EnzymeRates_pattern),
                      collect(eachline(trace_file)))
    @test length(relevant) <= INIT_BUDGET
end

@testset "wall-clock for init_mechanisms" begin
    init_mechanisms(uni_uni)   # JIT warmup
    t = @elapsed init_mechanisms(uni_uni)
    @test t < 100e-3   # 100ms (calibrated against main + headroom)
end
```

Budgets are calibrated against main during the pre-PR (Stage 0).

## 10. Behavior changes (user-visible)

These are intentional, but flagged so the PR description can call them out:

1. **DSL grammar change**: `@enzyme_reaction` regulators now declare
   per-regulator multiplicities (`R(1, 2, 4)`). The
   `competitive_inhibitors:` / `allosteric_regulators:` labels replace
   the old `dead_end_inhibitors:` / `allosteric_regulators:` split
   (already partly done on main). `@enzyme_mechanism` uses bare
   metabolite names (no atom brackets — atoms belong in
   `@enzyme_reaction`). `@allosteric_mechanism` uses
   `regulatory_site(multiplicity = N): begin ligands: ... end` blocks.
   - Migration path: existing user code using the old DSL needs
     rewriting. README + docstrings updated in Stage 7.
2. **Structural kinetic groups**: today's `kinetic_group::Int` field
   on `StepSpec` is gone; grouping is via `Vector{Vector{Step}}`. The
   user-facing `rate_equation` / `parameters` / `fitted_params`
   outputs are unchanged.
3. **Parameter naming**: positional (`:K1`, `:k1f`, etc.) preserved
   — same Symbols emitted. Chokepoint accessor enables a future
   refactor to path-based names without re-touching the rest of the
   codebase.

No other user-visible changes. `rate_equation(m, concs, params)`
returns the same numerical values; `fitted_params(m)` returns the same
tuple of Symbols; `parameters(m, Full)` / `parameters(m, Reduced)`
return the same tuples.

## 11. Staging plan

One large refactor PR on `refactor-to-concrete-types-instead-of-symbols`
with atomic stage commits. The compile-gate test infrastructure ships
**separately to main** as a small pre-PR.

### Stage 0 (separate PR to main, ~50 LOC)

- Add `test/test_compile_budget.jl` with trace-compile subprocess +
  wall-clock tests for `init_mechanisms`, `expand_mechanisms`, and
  `using EnzymeRates`.
- Calibrate budgets against current main; set thresholds with ~20%
  headroom over current values.
- Commit message: `Add compile-time regression test infrastructure`.
- Merge to main before opening the refactor PR.

### Stage 1 (branch, TDD-driven — see §12)

**Goal:** add the new concrete struct family + lift function + `name(p,
m)` chokepoint, with full TDD coverage. Existing src consumers unchanged
(they still use the old types — old types coexist with new types
**only** for this one transitional stage; replaced in Stage 3+).

**Files touched:**

- `src/types.jl` — add all structs from §5 (Substrate, Product,
  AllostericRegulator, CompetitiveInhibitor, Residual, Species, Step,
  RegulatorySite, Kd/Kiso/Kon/Koff/Kfor/Krev/Kreg/Keq/Etot/Lallo,
  ReactantAtoms, RegulatorMults, EnzymeReaction, Mechanism,
  AllostericMechanism) and the accessors from §6. The existing
  `EnzymeReaction` parametric singleton gets REPLACED by the new
  concrete struct (Stage 1 deletes the old parametric form, since
  every consumer that reads it can be mechanically adapted to the new
  field accessors in this same commit).
- `src/types.jl` — **repack `EnzymeMechanism{Metabolites, Reactions}`
  → `EnzymeMechanism{Sig}`** with `Sig === (Metabolites, Reactions)`
  content-equivalent. Every `where {M, R}` in src becomes `where {Sig}`
  plus `(metabolites_sig, reactions_sig) = Sig` destructuring inside
  body. Zero behavior change; pure spelling refactor. Same for
  `AllostericEnzymeMechanism`.
- `src/types.jl` — add `name(p::Parameter, m)` implementations
  computing today's positional Symbols (the chokepoint accessor).
- `src/types.jl` — add `_sig_of(::Mechanism)`,
  `_mechanism_from_sig(::Tuple)`, `EnzymeMechanism(::Mechanism)`,
  `Mechanism(::EnzymeMechanism)`. These replace the existing
  `EnzymeMechanism(spec::MechanismSpec)` converter from
  `mechanism_enumeration.jl` (which depends on the new Mechanism
  struct — the old MechanismSpec stays alive until Stage 5 enumeration
  rewrite, but its converter to EnzymeMechanism now goes through
  `Mechanism` as intermediate: `MechanismSpec → Mechanism →
  EnzymeMechanism{Sig}`).
- `test/test_types.jl` — add tests for every new struct + accessor,
  written FIRST per TDD.
- Existing tests untouched at the assertion level. Test files that
  consumed the old `EnzymeReaction{S,P,R,N}` parametric form get
  mechanically adapted to the new concrete form (e.g.,
  `substrates(r)` returns `Vector{Substrate}` instead of `Tuple{...}`)
  in this same commit.

**Exit criteria:**

- All tests green (including the new tests for new types).
- Compile-budget tests stay within main's calibrated budget +20%.
- src delta tracked. Expect +500 to +700 (foundation; recouped in
  Stage 3+).

### Stage 2 (branch — DSL)

**Goal:** rewrite the three DSL macros to emit constructors of the new
struct family using the prior branch's grammar. Old macros deleted.

**Files touched:**

- `src/dsl.jl` — full rewrite per §7.
- `test/test_dsl.jl` — mechanical syntax adaptation (DSL output shape
  changes from old tuple/struct mix to new concrete structs). No test
  deleted.
- Existing tests for DSL macros that fixture-use the old grammar:
  these get the grammar adapted in their input expression (but the
  output assertions stay strong; if the new macro produces a different
  Mechanism, the test fails and src is fixed).

**Exit criteria:**

- All DSL tests green.
- Compile-budget tests pass.
- src delta: expect roughly flat or slight negative as
  prior-branch-style DSL is more compact than current grammar.

### Stage 3 (branch — `EnzymeMechanism` derivation switch-over)

**Goal:** rewrite the `EnzymeMechanism` (non-allosteric) derivation
pipeline to consume `_mechanism_from_sig` output. Delete the Symbol
classifiers.

**Files touched:**

- `src/rate_eq_derivation.jl` — rewrite `rate_equation`,
  `parameters`, `rate_equation_string`, `_dependent_param_exprs`,
  `_kcat_forward` for `EnzymeMechanism` and friends. All consume
  `Mechanism(em)` output and walk `Step`/`Species` structs. Delete
  `_build_kinetic_rename_map` Symbol-based code, `_split_reaction_side`,
  `_compute_re_groups` Symbol-form versions.
- `src/thermodynamic_constr_for_rate_eq_derivation.jl` — rewrite to
  use Parameter struct family instead of `Symbol("K$idx")` builders.
- `src/sym_poly_for_rate_eq_derivation.jl` — replace `is_k_parameter`
  Symbol-regex logic with struct-typed dispatch.
- `test/test_rate_eq_derivation.jl` — mechanical syntax adaptation
  for assertions that build Mechanism / Parameter objects in tests.
  Numerical assertions + perf gates unchanged.

**Exit criteria:**

- All rate-eq tests green; `rate_equation` 0-alloc/<100ns gate green.
- Compile-budget tests pass (this is the highest-risk stage; if
  body-build cost balloons, fix here).
- src delta: target −800 to −1200.

### Stage 4 (branch — AllostericEnzymeMechanism derivation switch-over)

**Goal:** rewrite the allosteric derivation pipeline.

**Files touched:**

- `src/rate_eq_derivation.jl` (allosteric section) — rewrite
  AllostericEnzymeMechanism `rate_equation`, `parameters`, T-state
  building. Delete `_T_rename`, `_onlyR_syms`, `_all_t_state_names`,
  `_reg_param_name`, `_group_param_symbols` Symbol-based code.
- `test/test_rate_eq_derivation.jl` (allosteric tests) — mechanical
  adaptation.

**Exit criteria:**

- All allosteric tests green; perf gates green.
- **Mid-refactor checkpoint:** review cumulative src delta. Must be at
  least −500. Re-read spec; revise if anything's wrong.
- src delta: target −300 to −500.

### Stage 5 (branch — enumeration)

**Goal:** rewrite the enumeration pipeline to produce
`Vector{Mechanism}` directly using new structs. Delete
`MechanismSpec`/`StepSpec`/`AllostericMechanismSpec` and all the
Symbol-form-name helpers.

**Files touched:**

- `src/mechanism_enumeration.jl` — full rewrite. `init_mechanisms`,
  `expand_mechanisms`, `dedup!` all operate on `Mechanism`. Each
  `_expand_*` move produces new `Mechanism` instances. Delete
  `_form_name`, `_parse_bound`, `_dead_end_form_name`, `_atoms_dict`,
  `_is_estar_form`, all the per-form Symbol building.
- Delete old `EnzymeMechanism(spec::MechanismSpec)` converter (replaced
  by Stage 1's `EnzymeMechanism(::Mechanism)`).
- `test/test_mechanism_enumeration.jl` — mechanical adaptation.

**Exit criteria:**

- All enumeration tests green.
- Compile-budget tests pass (low risk — Mechanism is one type, not
  parametric, so no body-build explosion).
- src delta: target −1200 to −1500.

### Stage 6 (branch — identify_rate_equation canonicalizer)

**Goal:** replace the regex canonicalizer with a struct-based hash.

**Files touched:**

- `src/identify_rate_equation.jl` — replace
  `_canonicalize_rate_eq_with_map`, delete `_sort_run_factors`,
  `_factor_sort_key`, the regex pattern construction. New
  `_canonical_rate_eq_hash` walks `Mechanism(em)`'s steps + Parameter
  structs.
- `src/mechanism_enumeration.jl` — `_canonicalize!`/`_dedup_key` may
  also simplify; check for further dead code.
- `test/test_identify_rate_equation.jl` — mechanical adaptation; cache
  hit/miss tests adjust to the new hash. Numerical results unchanged.

**Exit criteria:**

- All identify tests green.
- src delta: target −300 to −500.

### Stage 7 (branch — fitting touch-up + cleanup + docs)

**Goal:** final dead-code sweep + documentation update.

**Files touched:**

- `src/fitting.jl` — minor adaptation if any (FittingProblem stores
  AbstractEnzymeMechanism; loss! uses fitted_params, metabolites,
  rate_equation — all unchanged interfaces).
- All src files — re-read end-to-end, delete dead code, inline
  single-use helpers, prune comments per §3 implementor discipline.
- `README.md` + docstrings — update DSL examples and any
  out-of-date pointers.
- `CLAUDE.md` — update architecture section to reflect the new
  hierarchy.
- Final `wc -l src/*.jl` check: cumulative delta must be ≤ −3,500.

**Exit criteria:**

- All tests green; all perf gates green.
- Code-size goal hit: ≤ 3,600 src LOC.
- README + CLAUDE.md current.
- PR ready for review.

## 12. TDD discipline (Stage 1 in particular)

Per CLAUDE.md Rule §TDD: every new struct and accessor in Stage 1
follows strict test-first development.

For each new type, the loop is:

1. **Write a failing test** in `test/test_types.jl` (or appropriate
   file) that exercises the constructor/accessor with concrete
   inputs and asserts the expected return.
2. **Run it**; confirm the test fails because the symbol doesn't exist
   yet (not because of an unrelated error).
3. **Write the minimum code** in `src/types.jl` to make the test pass.
4. **Run again**; confirm green.
5. **Refactor** if obvious simplification applies; tests stay green.

Per struct: ~5–10 tests covering:

- Constructor validation (rejection of invalid inputs)
- Accessor return types and values
- Equality / hashing semantics
- Pretty printing (where applicable)
- Canonicalization (Step's iso/binding direction, Species's `bound`
  sorting)

The DSL-driven round-trip tests come in Stage 2; Stage 1's tests
exercise structs via direct constructor calls.

In subsequent stages (3–6), TDD applies whenever new behavior is
introduced. Pure "switch from Symbol path A to struct path B" changes
follow the existing tests as the regression suite — write new tests
when the new code path adds behavior or changes invariants.

## 13. Risks and mitigations

| Risk | Mitigation |
|---|---|
| `@generated rate_equation` body-build cost explodes (prior branch's failure mode) | Compile-budget gate from Stage 0; trace-compile counts checked per commit; mid-refactor checkpoint reviews cumulative trend. |
| `rate_equation` 0-alloc/<100ns invariant breaks | Existing perf test fails the commit; fix before proceeding. Body-build code (not per-call code) is where struct lifts happen — per-call body is unchanged shape. |
| `loss!` regression | New baseline-comparison test added in Stage 0; flagged per commit. |
| Code-size goal not hit | Per-stage cumulative tracking; mid-refactor STOP-AND-REDESIGN gate after Stage 4. |
| Stages get out of sequence; tests fail mid-commit | "Tests green at every commit" CI; commits that fail can't land. Atomic commits per stage. |
| Long-lived branch accumulates merge conflicts with main | Weekly rebase onto main; small main during this window (no parallel large refactors). |
| Design defect discovered mid-refactor | Mid-refactor checkpoint after Stage 4 is a planned revision opportunity; the spec can amend. |
| Test deletions creep in under stress | NON-NEGOTIABLE banner (§2); per-stage review confirms no `@test_skip`/`@test_broken`/commented tests. |
| DSL grammar change breaks user code | Behavior change called out in PR description; README has migration examples; Stage 7 updates docstrings. |

## 14. Open questions

None at spec-time. Any open implementation questions surface during
execution and are resolved either inline (if narrow) or via spec
amendment + commit before continuing (if structural).

## 15. References

- Prior attempt branch: `refactor-to-use-structs-throughout` (1383-line
  spec, 3557-line plan, 201 commits, ~24k LOC delta — never merged).
- Memory note [[project-structs-throughout-refactor]] — historical
  context.
- Memory note [[feedback-simplification-means-less-code]] — refactor
  success metric.
- Memory note [[feedback-chokepoint-accessors-for-future-migrations]] —
  rationale for the `name(p::Parameter, m)` chokepoint design.
- Memory note [[feedback-subagent-driven-workflow]] — preferred
  execution pattern for substantial refactors.
- CLAUDE.md `rate_equation` perf invariant section — non-negotiable
  constraint.
