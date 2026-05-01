# EnzymeRates.jl

Identify the best enzyme rate equation from kinetic data.

Given a reaction definition and experimental rate measurements at varying
substrate, product, and regulator concentrations, EnzymeRates enumerates
all biochemically valid mechanisms, fits each to the data, and selects
the simplest mechanism that adequately describes the data based on
leave-one-group-out cross-validation.

The package has first-class support for MWC allostery and for mechanisms
that mix steady-state and rapid-equilibrium elementary steps.
Thermodynamic constraints (Haldane, Wegscheider) are derived
automatically from the cycle structure of the mechanism, so users supply
only the independent rate constants plus a measured equilibrium
constant.

## Installation

```julia
# README-SKIP-IN-TEST
using Pkg
Pkg.add(url="https://github.com/DenisTitovLab/EnzymeRates.jl")
```

## Define a mechanism, derive its rate equation, fit data

The running example is a uni-uni reaction `S ⇌ P` catalyzed by an MWC
homodimer in which substrate, product, the catalytic interconversion,
and an allosteric activator `A` all operate exclusively in the R
conformation. The T conformation is catalytically silent — a textbook
K-type allosteric activator.

```julia
using EnzymeRates

m = @allosteric_mechanism begin
    substrates: S
    products:   P
    allosteric_regulators: A::OnlyR

    site(:catalytic, 2): begin
        steps: begin
            [E, S] ⇌ [ES]      :: OnlyR
            [ES] <--> [EP]     :: OnlyR
            [EP] ⇌ [E, P]      :: OnlyR
        end
    end
end
```

The two `⇌` steps mark binding events at rapid equilibrium (one
binding constant `K` per step); the `<-->` step marks a steady-state
catalytic interconversion (independent forward and reverse rate
constants `kf`, `kr`). The `::OnlyR` annotation tells the framework
that those steps fire only in the R conformation, so the T-state
contribution to the rate equation is identically zero.

`parameters(m)` lists the names the framework needs at evaluation
time, `rate_equation_string(m)` prints the symbolic rate equation, and
`rate_equation(m, concs, params)` evaluates the rate numerically:

```julia
parameters(m)
rate_equation_string(m)
```

We generate synthetic data by evaluating `rate_equation` on a grid of
concentrations and adding multiplicative log-normal noise. Multiple
`group` values represent independent experimental batches that share
the same `E_total`; the framework's loss function is invariant to a
per-group `E_total` rescaling.

```julia
using OptimizationPyCMA, Random
Random.seed!(42)

true_params = (
    K1 = 1.0,         # S binding K (R-state)
    k2f = 5.0,        # catalytic SS forward rate (R-state)
    K3 = 0.5,         # P binding K (R-state)
    K_A_reg1 = 2.0,   # activator binding K (R-state)
    L = 0.1,          # conformational [T]/[R] for free enzyme
    Keq = 2.0,
    E_total = 1.0,
)

data_rows = NamedTuple[]
for grp in 1:5
    for _ in 1:10
        S = exp(randn() * 0.8)         # ~lognormal around 1
        A = 0.05 + 5.0 * rand()        # uniform 0.05..5.05
        P = rand() < 0.5 ? 0.05 : 0.5  # two product levels
        v_true = rate_equation(m, (S=S, P=P, A=A), true_params)
        v_obs = v_true * exp(0.05 * randn())   # 5% log-normal noise
        push!(data_rows, (group="G$grp", Rate=v_obs, S=S, P=P, A=A))
    end
end
data = (
    group = [r.group for r in data_rows],
    Rate  = [r.Rate  for r in data_rows],
    S     = [r.S     for r in data_rows],
    P     = [r.P     for r in data_rows],
    A     = [r.A     for r in data_rows],
)
```

The fit runs `fit_rate_equation` on a `FittingProblem`, using the PyCMA
optimizer (multi-start CMA-ES) recommended for rate-equation fitting.
Fitted rate constants are returned with kcat normalized to 1.0 by
default — the absolute scale is recovered by multiplying with a
separately measured kcat.

```julia
fp = FittingProblem(m, data; Keq=2.0)
result = fit_rate_equation(fp, PyCMAOpt();
    n_restarts=3, maxtime=5.0, popsize=50)
result.params       # K1, k2f, K3 recover near true; K_A_reg1 and L are
                    # not jointly identifiable from rates alone — see
                    # structural_identifiability_deficit(m)
result.loss         # final loss value (~5% noise floor)
```

## Recover the mechanism with `identify_rate_equation`

If the mechanism is unknown — only the overall reaction and its
regulators are — `identify_rate_equation` enumerates biochemically
valid mechanisms, fits each to the data, and returns the simplest that
generalizes (judged by leave-one-group-out cross-validation). The same
chemistry from Section 2, declared as a *reaction*:

```julia
rxn = @enzyme_reaction begin
    substrates: S
    products:   P
    regulators: A
    oligomeric_state: 2
end
```

`regulators: A` declares `A` with an unspecified role; the search
enumerates dead-end-inhibitor and allosteric variants and selects
between them on cross-validation score. (If you already know `A` is
allosteric, declare it with `allosteric_regulators: A` instead and the
search skips dead-end variants.)

```julia
prob = IdentifyRateEquationProblem(rxn, data; Keq=2.0)
```

The actual search runs `fit_rate_equation` on each candidate and is
slow on a laptop (minutes); skip the next block if you just want to
read along, or run it when you have time.

```julia
# README-SKIP-IN-TEST
results = identify_rate_equation(prob;
    optimizer=PyCMAOpt(),
    max_param_count=10,
    pmap_function=map,            # serial; pass `pmap` for distributed
)
results.best                       # the recovered mechanism
rate_equation_string(results.best) # printed rate equation
first(results.cv_results, 5)       # top rows of the CV-score DataFrame
```

`results.best` is the mechanism with minimum training loss at the
parameter-count level whose CV score is lowest — i.e., the simplest
mechanism that generalizes. For the synthetic data we generated, the
recovered mechanism agrees with the one we used to generate it.
