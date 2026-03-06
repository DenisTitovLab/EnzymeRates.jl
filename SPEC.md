# EnzymeRates.jl API Specification

## Purpose

EnzymeRates.jl identifies the best enzyme rate equation from kinetic data.
Given a reaction definition and experimental measurements at varying
substrate/product/regulator concentrations, the package enumerates all
biochemically valid mechanisms, fits each to the data, and selects the
mechanism with the fewest parameters that adequately describes the data
based on cross-validation.

Users can also use the package to derive rate equations from manually
defined mechanisms and to fit a known mechanism to data, but these are
secondary use cases.

---

## Exported API

### Types

| Type | Description |
|------|-------------|
| `EnzymeReaction` | Overall reaction specification (substrates, products, regulators). Created via `@enzyme_reaction`. |
| `EnzymeMechanism` | Full mechanism with species, steps, RE/SS flags, and constraints. Created via `@enzyme_mechanism` or from selection results. |
| `OligomericEnzymeMechanism` | Multi-subunit MWC allosteric enzyme. Created via `@enzyme_mechanism` with `site()`/`conformations:` DSL. |
| `FittingProblem` | Single-mechanism fitting problem (mechanism + data + Keq). |
| `IdentifyRateEquationProblem` | Multi-mechanism selection problem (reaction + data + Keq + search config). |
| `IdentifyRateEquationResults` | Results from `identify_rate_equation`: all fitted candidates with CV scores. |

### Macros

| Macro | Description |
|-------|-------------|
| `@enzyme_reaction` | Create an `EnzymeReaction` from a DSL block (substrates, products, regulators). |
| `@enzyme_mechanism` | Create an `EnzymeMechanism` or `OligomericEnzymeMechanism` from DSL blocks. |

### Constants

| Constant | Description |
|----------|-------------|
| `Full` | Rate equation mode: all 2N raw rate constants + E_total. |
| `Reduced` | Rate equation mode (default): thermodynamically independent k's + Keq + E_total. |

### Functions

| Function | Description |
|----------|-------------|
| `identify_rate_equation` | Select the best mechanism from all candidates for a reaction. Primary entry point. |
| `fit_rate_equation` | Fit a single mechanism to data. |
| `rate_equation` | Compute the steady-state rate: `rate_equation(m, concs, params, [mode])`. Default mode is `Reduced`. |
| `rate_equation_string` | Human-readable string representation of a mechanism's rate equation. |
| `parameters` | Parameter names required by a mechanism: `parameters(m, [mode])`. Default mode is `Reduced`. |
| `metabolites` | Distinct metabolite names as a tuple of Symbols: `metabolites(m) → (:S, :P)`. |
| `structural_identifiability_deficit` | Structural identifiability deficit (non-positive = identifiable). |
| `rescale_parameter_values` | Rescale SS rate constants so kcat equals target (default 1.0). K's, Keq, E_total unchanged. |
| `compile_mechanism` | Convert a `MechanismSpec` to `EnzymeMechanism` or `OligomericEnzymeMechanism`. |

---

## Core Workflow

```julia
using EnzymeRates

# 1. Define the reaction
rxn = @enzyme_reaction begin
    substrates: A[C], B[N]
    products:   P[C], Q[N]
end

# 2. Construct the selection problem
prob = IdentifyRateEquationProblem(rxn, data; Keq=5.0)

# 3. Identify the best rate equation
results = identify_rate_equation(prob)

# 4. Inspect results
best = results.best
best.mechanism        # EnzymeMechanism
best.params           # NamedTuple of fitted parameter values
best.cv_score         # cross-validation score
best.loss             # training loss

rate_equation_string(best.mechanism)

# 5. Use the identified equation
v = rate_equation(best.mechanism, (A=1.0, B=0.5, P=0.1, Q=0.05), best.params)
```

---

## Data Format

The data table must be any Tables.jl-compatible object with the following
columns:

| Column | Type | Description |
|--------|------|-------------|
| `group` | Any | Identifies a group of measurements sharing the same E_total (e.g., same experiment, same enzyme preparation). Each unique `group` value defines one CV fold for leave-one-group-out cross-validation. Must have >= 2 unique values. |
| `Rate` | Real | Measured reaction rate. Must be nonzero (zero rates produce -Inf in log-space). |
| One column per metabolite | Real | Concentration of each metabolite. Column names must match `metabolites(mechanism)`. |

Example:

```julia
data = (
    group = ["exp1", "exp1", "exp1", "exp2", "exp2", "exp2"],
    Rate  = [0.5, 1.2, 2.1, 0.4, 1.1, 1.9],
    A     = [0.1, 0.5, 1.0, 0.1, 0.5, 1.0],
    B     = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
    P     = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    Q     = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
)
```

### Migration from Current Format (NOT YET DONE)

The current `FittingProblem` implementation still requires `Article` and
`Fig` columns. This spec calls for migrating to a single `group` column.
The `group` column will replace the combination of `Article` and `Fig` —
users who had both should combine them (e.g., `group = "ArticleName_Fig1"`).

---

## `IdentifyRateEquationProblem`

### Constructor

```julia
IdentifyRateEquationProblem(
    reaction::EnzymeReaction,
    data;
    Keq::Real,
    max_re_groups::Int = 7,
    catalytic_n::Int = 0,
)
```

**Construction is eager**: the `MechanismIterator` is created at
construction time (enumeration stages 1-4 are materialized, stage 5 is
lazy). This lets users inspect `prob.mechanism_iterator` and
`length(prob.mechanism_iterator)` before calling `identify_rate_equation`.

### Fields (User-Accessible)

| Field | Type | Description |
|-------|------|-------------|
| `reaction` | `EnzymeReaction` | The reaction being analyzed. |
| `data` | NamedTuple | The data table (columnar). |
| `Keq` | `Float64` | Thermodynamic equilibrium constant. |
| `mechanism_iterator` | `MechanismIterator` | Lazy iterator over all candidate mechanisms. |

---

## `identify_rate_equation`

### Signature

```julia
identify_rate_equation(
    prob::IdentifyRateEquationProblem;
    # Search strategy
    strategy::Symbol = :exhaustive,     # :exhaustive or :beam (future)
    # Fitting parameters
    optimizer = LBFGSB(),               # Optimization.jl optimizer
    n_restarts::Int = 10,               # multi-start restarts per mechanism
    maxtime::Real = 60.0,               # max time per mechanism fit (seconds)
    # Cross-validation
    cv_fraction::Float64 = 0.2,         # fraction of groups held out per fold
    # Output
    save_path::Union{Nothing,String} = nothing,  # CSV save path (incremental)
    # Progress
    show_progress::Bool = true,
) → IdentifyRateEquationResults
```

### Algorithm

1. **Group mechanisms by parameter count** using the lazy iterator.
2. **For each parameter-count group** (ascending order):
   a. Compile each mechanism (`MechanismSpec` → `EnzymeMechanism`).
   b. Fit to full training data via `fit_rate_equation`.
   c. Compute leave-one-group-out CV score.
   d. Record results.
   e. If `save_path` is set, append results for this group to CSV.
3. **Show progress** within each parameter-count group.
4. **Select best**: mechanism with fewest parameters whose CV score is
   adequate (details of the selection criterion are a future design
   decision — initially just return the full ranking and let the user
   decide).
5. Return `IdentifyRateEquationResults`.

### Future: Beam Search Strategy

When `strategy = :beam` (not yet implemented):
1. Fit all mechanisms with the minimum parameter count.
2. Take the top-N by CV score.
3. For each, expand to mechanisms with one additional parameter.
4. Repeat until no improvement or max parameters reached.

---

## `IdentifyRateEquationResults`

### Fields

| Field | Type | Description |
|-------|------|-------------|
| `best` | `CandidateResult` | Best mechanism by selection criterion. |
| `candidates` | `Vector{CandidateResult}` | All fitted candidates, sorted by (n_params, cv_score). |
| `problem` | `IdentifyRateEquationProblem` | Reference to the problem that produced these results. |

### `CandidateResult` (Internal Struct, Accessed via Results)

| Field | Type | Description |
|-------|------|-------------|
| `mechanism` | `EnzymeMechanism` | The compiled mechanism. |
| `params` | `NamedTuple` | Fitted parameter values (independent k's + Keq + E_total). |
| `n_params` | `Int` | Number of fitted parameters. |
| `cv_score` | `Float64` | Leave-one-group-out cross-validation score. |
| `loss` | `Float64` | Training loss on full data. |

### CSV Output Format

When `save_path` is provided, results are written incrementally as each
parameter-count group completes. One CSV file at the specified path.

Columns:

| Column | Description |
|--------|-------------|
| `n_params` | Number of fitted parameters |
| `cv_score` | Cross-validation score |
| `loss` | Training loss |
| `mechanism_type` | Eval-able string that reconstructs the `EnzymeMechanism` type (e.g., the full type signature) |
| `rate_equation` | Human-readable rate equation string |
| `param_1`, `param_2`, ... | Fitted parameter values |

The `mechanism_type` column contains a string that can be `eval`-ed to
produce the `EnzymeMechanism` instance:

```julia
row = CSV.Row(...)
m = eval(Meta.parse(row.mechanism_type))
```

---

## `FittingProblem` (Updated)

### Constructor

```julia
FittingProblem(mechanism::EnzymeMechanism, data; Keq::Real)
```

The data table uses the same format as `IdentifyRateEquationProblem`:
a `group` column identifies measurement groups sharing the same E_total.
The `Rate` column and one column per metabolite are required.

### Usage

```julia
m = @enzyme_mechanism begin
    species: begin
        substrates: S[C]
        products:   P[C]
        enzymes:    E, ES[C]
    end
    steps: begin
        [E, S] <--> [ES]
        [ES] <--> [E, P]
    end
end

fp = FittingProblem(m, data; Keq=5.0)
result = fit_rate_equation(fp, LBFGSB(); n_restarts=10, maxtime=60.0)
```

---

## `fit_rate_equation`

### Signature

```julia
fit_rate_equation(
    fp::FittingProblem,
    optimizer;
    n_restarts::Int = 10,
    maxtime::Real = 60.0,
    lb = fill(-15.0, length(fitted_params(fp.mechanism))),
    ub = fill(15.0, length(fitted_params(fp.mechanism))),
    kwargs...,
) → NamedTuple
```

Returns a NamedTuple of fitted parameter values (log-space internally,
converted to natural scale in the result).

---

## Cross-Validation

### Strategy: Leave-One-Group-Out

Each unique value of the `group` column defines one fold. In each fold:

1. Hold out all data points from one group.
2. Fit the mechanism to the remaining groups.
3. Predict the held-out group's rates using fitted parameters.
4. Compute the held-out loss.

The CV score is the mean held-out loss across all folds.

### Rationale

Data points within a group share systematic errors (e.g., same enzyme
preparation, same E_total scaling). Leave-one-group-out respects this
structure and gives an honest estimate of generalization to new
experimental conditions.

### Requirements

- At least 2 unique `group` values are required.
- The `cv_fraction` parameter in `identify_rate_equation` controls what
  fraction of groups are held out per fold (default 0.2). If there are
  5 groups, each fold holds out 1 group. If there are 20 groups, each
  fold holds out 4 groups.

---

## Secondary Use Cases

### Derive Rate Equations from Manual Mechanisms

```julia
m = @enzyme_mechanism begin
    species: begin
        substrates: A[CX], B[N]
        products:   P[C], Q[NX]
        enzymes:    E, EA[CX], FP[CX], F[X], FB[NX], EQ[NX]
    end
    steps: begin
        [E, A] <--> [EA]
        [EA] <--> [FP]
        [FP] <--> [F, P]
        [F, B] <--> [FB]
        [FB] <--> [EQ]
        [EQ] <--> [E, Q]
    end
end

rate_equation_string(m)
parameters(m)
structural_identifiability_deficit(m)
```

### Fit a Known Mechanism to Data

```julia
fp = FittingProblem(m, data; Keq=5.0)
result = fit_rate_equation(fp, LBFGSB())
```

---

## Complete Exported Symbol List

```julia
# Types
export EnzymeReaction, EnzymeMechanism, OligomericEnzymeMechanism
export FittingProblem
export IdentifyRateEquationProblem, IdentifyRateEquationResults

# Macros
export @enzyme_reaction, @enzyme_mechanism

# Rate equation modes
export Full, Reduced

# Core: model selection
export identify_rate_equation

# Core: single-mechanism fitting
export fit_rate_equation

# Core: rate equation evaluation
export rate_equation, rate_equation_string

# Parameters & metabolites
export parameters, metabolites

# Identifiability
export structural_identifiability_deficit

# Parameter rescaling
export rescale_parameter_values

# Mechanism compilation
export compile_mechanism
```

Total: 6 types + 2 macros + 2 constants + 9 functions = **19 exported symbols**.

---

## Known Limitations

### `rate_equation` compilation for large mechanisms

`rate_equation(m, concs, params)` uses `@generated` functions that perform
symbolic King-Altman/Cha rate equation derivation at compile time. Because
each unique `EnzymeMechanism` is a distinct type (species, reactions, RE/SS
flags, and constraints are all encoded in type parameters), every new
mechanism triggers a full symbolic derivation during compilation.

For mechanisms with many enzyme forms and elementary steps (e.g., Bi-Bi
reactions with regulators and dead-end complexes), this compilation can:

- Take a very long time (minutes)
- Exhaust available memory
- Cause a StackOverflow in the compiler

This is an inherent consequence of the type-parameter architecture that
enables zero-allocation rate evaluation at runtime. The
`identify_rate_equation` pipeline should mitigate this by processing
candidates in ascending order of `param_count_estimate` and enforcing
time/memory budgets per mechanism.
