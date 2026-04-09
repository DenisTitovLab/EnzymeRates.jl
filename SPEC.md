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
| `AllostericEnzymeMechanism` | Multi-subunit MWC allosteric enzyme. Created via `@enzyme_mechanism` with `site()` DSL. |
| `FittingProblem` | Single-mechanism fitting problem (mechanism + data + Keq). |
| `IdentifyRateEquationProblem` | Multi-mechanism selection problem (reaction + data + Keq + search config). |
| `IdentifyRateEquationResults` | Results from `identify_rate_equation`: all fitted candidates with CV scores. |

### Macros

| Macro | Description |
|-------|-------------|
| `@enzyme_reaction` | Create an `EnzymeReaction` from a DSL block (substrates, products, regulators). |
| `@enzyme_mechanism` | Create an `EnzymeMechanism` or `AllostericEnzymeMechanism` from DSL blocks. |

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
| `compile_mechanism` | Convert a `MechanismSpec` to `EnzymeMechanism` or `AllostericEnzymeMechanism`. |

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
results.best              # AbstractEnzymeMechanism — best mechanism
results.cv_results        # DataFrame with LOOCV results for top candidates

rate_equation_string(results.best)

# 5. Use the identified equation
params = results.cv_results[1, :]  # get params from CV results DataFrame
v = rate_equation(results.best, (A=1.0, B=0.5, P=0.1, Q=0.05), params)
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
)
```

Holds the reaction, data, and equilibrium constant. The beam search
pipeline is driven by `identify_rate_equation`, not by a pre-built iterator.

### Fields (User-Accessible)

| Field | Type | Description |
|-------|------|-------------|
| `reaction` | `EnzymeReaction` | The reaction being analyzed. |
| `data` | NamedTuple | The data table (columnar, with `:group` column). |
| `Keq` | `Float64` | Thermodynamic equilibrium constant. |

---

## `identify_rate_equation`

### Signature

```julia
identify_rate_equation(
    prob::IdentifyRateEquationProblem;
    # Beam search
    min_beam_width::Int = 200,          # minimum mechanisms to keep per level
    beam_fraction::Float64 = 0.1,       # fraction of mechanisms to keep
    max_param_count::Int = 20,          # stop expanding beyond this
    # Fitting parameters
    optimizer = PyCMAOpt(),             # Optimization.jl optimizer
    n_restarts::Int = 10,               # multi-start restarts per mechanism
    maxtime::Real = 60.0,               # max time per mechanism fit (seconds)
    maxiters::Int = 10_000_000,         # max iterations per optimizer run
    popsize::Int = 200,                 # population size for optimizer
    # Model selection
    n_cv_candidates::Int = 5,           # LOOCV top N per param count
    # Output
    save_dir::Union{Nothing,String} = nothing,  # directory for per-level CSVs
    # Parallelism
    pmap_function::Function = map,      # e.g. Distributed.pmap for HPC
) → IdentifyRateEquationResults
```

### Algorithm

Two-phase beam search pipeline:

**Phase 1 — Beam search** (finds candidate mechanisms):
1. `init_mechanisms(reaction)` → fit all on full data → rank by loss.
2. Keep top N by loss (beam width = `max(beam_fraction * n, min_beam_width)`).
3. `expand_mechanisms` on beam → `dedup!` → fit → rank by loss.
4. Repeat until no new mechanisms or `max_param_count` reached.
5. Save all fitted mechanisms to per-param-count CSV files if `save_dir` set.

**Phase 2 — Model selection** (finds optimal complexity):
1. Take top `n_cv_candidates` per param count (by loss from Phase 1).
2. Leave-one-group-out CV each candidate.
3. Best mechanism = lowest training loss at the param count with best CV score.

### Note on `optimizer`

`optimizer` is passed as a keyword argument to `identify_rate_equation` but
forwarded as a **positional** argument to `fit_rate_equation` internally.

---

## `IdentifyRateEquationResults`

### Fields

| Field | Type | Description |
|-------|------|-------------|
| `best` | `AbstractEnzymeMechanism` | Best mechanism (lowest loss at optimal param count). |
| `cv_results` | `DataFrame` | LOOCV results for top candidates per param count. |

### CSV Output Format

When `save_dir` is provided, one CSV file per parameter count level:
`save_dir/params_7.csv`, `save_dir/params_8.csv`, etc.

Columns:

| Column | Description |
|--------|-------------|
| `n_params` | Number of fitted parameters |
| `loss` | Training loss |
| `mechanism_type` | Eval-able string that reconstructs the mechanism type |
| `rate_equation` | Human-readable rate equation string |
| `K_S`, `k1f`, ... | Fitted parameter values (vary per mechanism) |

The `mechanism_type` column contains a string that can be `eval`-ed to
produce the mechanism instance:

```julia
row = CSV.Row(...)
m = eval(Meta.parse(row.mechanism_type))
```

---

## `FittingProblem` (Updated)

### Constructor

```julia
FittingProblem(mechanism::AbstractEnzymeMechanism, data; Keq::Real)
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
- LOOCV is used: each fold holds out exactly one group.

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
export EnzymeReaction, EnzymeMechanism, AllostericEnzymeMechanism
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
