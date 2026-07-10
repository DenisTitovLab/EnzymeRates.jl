# Getting Started

This page walks through the three pillars of EnzymeRates.jl in one arc:
**define** a mechanism, **derive** its rate equation, **fit** it to data, and
**identify** the best mechanism from data automatically. Each section links
forward to the full pillar page for details.

---

## Define a reaction and a mechanism

A reaction names the substrates, products, and any regulators. A mechanism
adds the elementary steps.

```@example getting-started
using EnzymeRates

rxn = @enzyme_reaction begin
    substrates: S[C]
    products:   P[C]
end

m = @enzyme_mechanism begin
    substrates: S
    products:   P
    steps: begin
        E + S ⇌ E(S)
        E(S) <--> E(P)
        E(P) ⇌ E + P
    end
end

m isa EnzymeMechanism
```

`⇌` marks a rapid-equilibrium binding step, described by one equilibrium
constant K. `<-->` marks a steady-state interconversion, described by
independent forward and reverse rate constants. The result is an
`EnzymeMechanism` value. The expression `m isa EnzymeMechanism` (rather than
`typeof(m)`) is shown deliberately: the `Sig` type parameter that encodes the
mechanism structure is an unreadable string, so `isa` is the useful check.

---

## Derive the rate equation

[Derivation tutorial](@ref) and the other deriving pages cover all
mechanism forms in depth.

```@example getting-started
parameters(m)
```

```@example getting-started
metabolites(m)
```

```@example getting-started
print(rate_equation_string(m))
```

`parameters(m)` returns the independent parameter symbols. `metabolites(m)`
returns the concentration symbols the rate equation expects. `Keq` is the
**equilibrium constant** of the overall reaction — the ratio of product to
substrate concentrations once the reaction has reached equilibrium, fixed by
the reaction's thermodynamics. It is always user-supplied; the package never
estimates it from data.

---

## Evaluate the rate equation numerically

Supply a `NamedTuple` of parameters and a `NamedTuple` of concentrations:

```@example getting-started
params = (
    K_S_E = 1.0,
    K_P_E = 1.0,
    k_ES_to_EP = 10.0,
    Keq = 2.0,
    E_total = 1.0,
)
concs = (S = 2.0, P = 0.5)

rate_equation(m, concs, params)
```

---

## Estimate kinetic parameters from data

[Fitting tutorial & data format](@ref) covers the data format, loss function,
and optimizer options in depth.

Generate noiseless synthetic data from the mechanism above:

```@example getting-started
using OptimizationCMAEvolutionStrategy, Random
Random.seed!(1)

groups = String[]; Rate = Float64[]; Svals = Float64[]; Pvals = Float64[]
for g in 1:3, _ in 1:8
    s = 0.1 + 9.9 * rand()
    p = 0.1 + 9.9 * rand()
    push!(groups, "G$g")
    push!(Rate, rate_equation(m, (S = s, P = p), params))
    push!(Svals, s); push!(Pvals, p)
end
data = (group = groups, Rate = Rate, S = Svals, P = Pvals)

fp = FittingProblem(m, data; Keq = 2.0)
fitted_params, loss, retcode = fit_rate_equation(fp, CMAEvolutionStrategyOpt(); n_restarts = 1, maxtime = 2.0)
retcode
```

The data table has a `group` column identifying measurement batches that share
one `E_total`, a `Rate` column with measured rates, and one column per
metabolite or regulator name (here `S` and `P`). `Keq` is a required keyword,
always user-supplied.

---

## Identify the enzyme mechanism from data

[Identify tutorial](@ref) covers the full rate-equation search, cross-validation, and
production settings.

Build an `IdentifyRateEquationProblem` from the reaction and data, then run
`identify_rate_equation` with the width-1 beam settings so the search finishes
in seconds:

```@example getting-started
prob = IdentifyRateEquationProblem(rxn, data; Keq = 2.0)

results = identify_rate_equation(prob;
    min_beam_width = 1,
    loss_rel_threshold = 1.0,
    loss_abs_threshold = 0.0,
    max_param_count = 6,
    n_cv_candidates = 1,
    save_dir = mktempdir(),
    optimizer = CMAEvolutionStrategyOpt(),
    n_restarts = 1, maxtime = 1.0,
    show_progress = false)

print(rate_equation_string(results.best))
```

`min_beam_width=1` and `loss_rel_threshold=1.0` collapse the beam to exactly
one survivor per parameter-count level, making the search fast and
deterministic. The full production search uses the wider defaults
(`min_beam_width=50`, `loss_rel_threshold=2.0`, `loss_abs_threshold=0.01`,
`loss_parsimony_threshold=1.01`, `max_param_count=20`, `eq_complexity_filter=337`)
and would often run for many hours and require a High
Performance Compute cluster (see [Running in parallel](identify/parallel.md)).
`save_dir` is mandatory; the search writes its progress and results there:
`progress.log`, `initial_mechanisms.csv`, one `equation_search_iteration_N.csv`
per beam iteration, and — once cross-validation finishes —
`loocv_results.csv` (the full leave-one-group-out table for every candidate that
entered CV) and `best_equation.csv` (the single selected equation with its
fitted parameters). The last two persist the model-selection outcome, so a
cluster run's result is available without re-deriving it from the iteration
files.
