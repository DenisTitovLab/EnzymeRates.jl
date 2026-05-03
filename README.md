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
homodimer with V-type allosteric activation. Substrate and product
bind both conformations with the same affinity (`:EqualRT`), but only
the R conformation catalyzes (`:OnlyR`); an allosteric activator `A`
binds R preferentially and shifts the population from the catalytically
silent T toward R.

```julia
using EnzymeRates

m = @allosteric_mechanism begin
    substrates: S
    products:   P
    allosteric_regulators: A::OnlyR

    site(:catalytic, 2): begin
        steps: begin
            E + S ⇌ ES      :: EqualRT
            ES <--> EP     :: OnlyR
            EP ⇌ E + P      :: EqualRT
        end
    end
end
```

The two `⇌` steps mark binding events at rapid equilibrium (one
binding constant `K` per step); the `<-->` step marks a steady-state
catalytic interconversion (independent forward and reverse rate
constants `kf`, `kr`). The `::EqualRT` annotation says the
corresponding K is shared between R and T conformations; `::OnlyR`
says the catalytic step fires only in R; and the `A::OnlyR` regulator
means the activator binds R only.

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
    K1 = 1e-4,           # S binding K (shared R/T)
    k2f = 100.0,         # catalytic SS forward rate (R only)
    K3 = 1e-3,           # P binding K (shared R/T)
    K_A_reg1 = 1e-5,     # activator binding K (R only)
    L = 10000.0,         # conformational [T]/[R] for free enzyme
    Keq = 2.0,
    E_total = 1.0,
)

# Sample concentrations log-uniformly across 100x below to 100x above
# each metabolite's K — the regime where the rate equation is informative.
logu(K) = K * 10.0 ^ (rand() * 4 - 2)

data_rows = NamedTuple[]
for grp in 1:5
    for _ in 1:10
        S = logu(true_params.K1)
        P = logu(true_params.K3)
        A = logu(true_params.K_A_reg1)
        v_true = rate_equation(m, (S=S, P=P, A=A), true_params)
        v_obs = v_true * exp(0.05 * randn())   # 5% multiplicative log-normal noise
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
result.params       # K1, K3, K_A_reg1, L recover near true.
                    # k2f is normalized so kcat = 1.0 (its true
                    # value 100.0 is the kcat scale).
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
    substrates: S[C]
    products:   P[C]
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
slow on a laptop (~1 hour on a single core with default settings;
faster with `pmap` distributed across workers). Skip the next block
if you just want to read along, or run it when you have time.

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
recovered mechanism agrees with the one we used to generate it.

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
are derived.

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

## How `identify_rate_equation` works

### Enumeration as composable building blocks

Mechanism enumeration is built from three small functions, not a
monolithic pipeline:

- `init_mechanisms(reaction)` produces the biochemically minimal
  mechanisms for a reaction by combining catalytic topologies (orderings
  of substrate binding, catalytic interconversion, and product release —
  random-order, ordered, ping-pong) with subsets of dead-end inhibition
  steps. Steps that bind the same metabolite share a kinetic group, so
  the parameter count starts at the smallest physically meaningful
  value.
- `expand_mechanisms(specs, reaction)` applies a fixed set of
  single-move expansions to each spec — converting an RE step to SS,
  splitting a kinetic group, adding a dead-end regulator, becoming
  allosteric, changing an allosteric state — and returns the expanded
  candidates keyed by their estimated parameter count.
- `dedup!(cache)` canonicalizes specs (sorted steps; renumbered kinetic
  groups by first occurrence) and removes structural duplicates.

The enumeration is grounded in chemical reasoning rather than blind
combinatorics: a step is "elementary" only if it changes one binding
site by one event with atom balance preserved, and only catalytic
topologies that satisfy bounds on bound-metabolite count, isomerization
size, and substrate participation are emitted.

### Beam search across parameter counts

`identify_rate_equation` walks parameter counts in ascending order:

1. Fit all candidates at the smallest parameter count on the full
   data; record training loss.
2. Keep the top fraction by training loss (at least
   `min_beam_width` candidates, or all of them if there are
   fewer than that).
3. Apply `expand_mechanisms` to surviving specs to produce candidates
   at the next parameter-count level.
4. `dedup!` and fit; rank by training loss.
5. Repeat until no new candidates appear or `max_param_count` is
   reached.

The beam width balances coverage (more candidates explored) against
runtime (every kept candidate gets a multi-restart fit).

### Model selection by leave-one-group-out cross-validation

After beam search, the top `n_cv_candidates` mechanisms per parameter
count enter LOOCV. Each unique value of the `group` column defines one
fold: the mechanism is fit on every group except one, then evaluated on
the held-out group. The CV score is the mean held-out loss across
folds. The "best" mechanism is the one with minimum training loss at
the parameter count whose CV score is lowest — *the simplest mechanism
that generalizes*. The `group` column reflects experimental batches
that share an `E_total`; LOOCV respects this structure and gives an
honest estimate of generalization to new conditions.
