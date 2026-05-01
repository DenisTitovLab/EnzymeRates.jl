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

`results.best` is the simplest mechanism that generalizes:
`identify_rate_equation` first picks the parameter-count with the
lowest CV score, then within that level picks the candidate with the
lowest training loss. For the synthetic data we generated, the
recovered mechanism reproduces the rate equation up to the
unidentifiable parameters flagged earlier (`K_A_reg1`, `L`) — typically
the simplest variant of the generating mechanism with non-identifiable
parameters dropped.

## How rate-equation derivation works

### Steady-state vs rapid equilibrium

`<-->` denotes a steady-state (QSSA) elementary step — both the forward
and reverse rate constants enter the rate equation as independent
parameters, and the King-Altman determinant assembles them into the
denominator polynomial. `⇌` denotes a rapid-equilibrium step — only the
binding constant `K` matters, because the framework collapses the
forward and reverse rates into a single equilibrium relation. A typical
mechanism mixes both, and the framework handles the mixed Cha-style
derivation automatically. `parameters(m)` reflects this: each RE step
contributes one `K`; each SS step contributes a forward `kf` and a
reverse `kr`.

### Haldane and Wegscheider relationships

When the mechanism contains thermodynamic cycles — any closed loop of
binding and catalytic steps — the rate constants around the cycle are
constrained by the equilibrium constant of the overall reaction. The
framework detects these cycles automatically (via the null space of the
mechanism's enzyme-form incidence matrix), declares one rate constant per
cycle as *dependent*, and computes it from the rest plus a user-supplied
`Keq`. You fit the *independent* rate constants; dependent constants
are derived. `structural_identifiability_deficit(m)` reports the deficit
of the mechanism's parameter map: zero means every independent
parameter can in principle be identified from the rate equation.

### Allostery: the MWC R/T model

For multi-subunit enzymes, the framework uses the Monod-Wyman-Changeux
two-state model: the enzyme exists in an active R conformation and an
inactive T conformation, with `L = [T]/[R]` the conformational
equilibrium for the bare enzyme and the same `L` propagating to all
ligand-bound species. Each kinetic group (binding step, catalytic
interconversion) can be `:OnlyR`, `:EqualRT`, or `:NonequalRT`; each
regulatory ligand can additionally be `:OnlyT`. The four tags:

- `:OnlyR` — the symbol exists in R only; T-state contributions are
  zero. A `:OnlyR` activator binds R preferentially and shifts the
  population toward R.
- `:OnlyT` — symbol exists in T only. A `:OnlyT` regulator binds T
  preferentially and shifts the population toward T (a typical
  allosteric inhibitor).
- `:EqualRT` — same `K` (or `kf`, `kr`) in both conformations. Useful
  for ligands that bind without conformational preference.
- `:NonequalRT` — independent R and T parameters (`K_R`, `K_T`).

The full rate equation is then the sum of R-state and T-state numerator
terms, weighted by the partition function `(R-state polynomial)^n +
L*(T-state polynomial)^n`, where `n` is the oligomeric state. The
example mechanism above uses `:OnlyR` everywhere — the T-state
contributions vanish and the printed rate equation simplifies
accordingly.
