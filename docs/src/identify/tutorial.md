# Identify tutorial

`identify_rate_equation` takes an `EnzymeReaction` and rate data, enumerates
biochemically valid mechanisms, and fits the most promising candidates at each
parameter count — fitting every mechanism is far too expensive — returning the
simplest one that generalizes by leave-one-group-out
cross-validation.

This page walks through a fully runnable example: an MWC allosteric enzyme,
recovered from noiseless data in about a minute. The full production search
widens the beam to the defaults (`min_beam_width=50`, `loss_rel_threshold=2.0`,
`loss_abs_threshold=0.01`, `loss_parsimony_threshold=1.01`, `max_param_count=20`,
`eq_complexity_filter=337`) and would often run for many
hours and require a High Performance Compute cluster (see
[Running in parallel](@ref)).

```@setup identify_fast
using EnzymeRates
```

## A reaction and a generating mechanism

The example is a reversible uni-uni reaction `S ⇌ P` run by a dimeric enzyme
with an allosteric activator `A`. The reaction declares the regulator and the
oligomeric state; from these [The enumeration engine](@ref) builds the MWC
variants the search fits.

```@example identify_fast
using EnzymeRates

rxn = @enzyme_reaction begin
    substrates: S[C]
    products:   P[C]
    allosteric_regulators: A
    oligomeric_state: 2
end

# A concrete MWC mechanism to generate the data. A binds only the active
# conformation (:OnlyA), which makes it an allosteric activator.
generator = @allosteric_mechanism begin
    substrates: S
    products:   P
    catalytic_multiplicity: 2
    allosteric_regulators: A::OnlyA
    catalytic_steps: begin
        E + S ⇌ E(S)         :: EqualAI
        E(S) <--> E(P)       :: OnlyA
        E(P) ⇌ E + P         :: EqualAI
    end
end

println("fitted params: ", EnzymeRates.fitted_params(generator))
println("metabolites:   ", metabolites(generator))
```

The mechanism has five independent parameters: the binding constants `K_S_E`
and `K_P_E` (shared by both conformations, `:EqualAI`), the active-state
catalytic constant `k_A_ES_to_EP` (`:OnlyA`), the activator binding constant
`K_A_Areg`, and the conformational equilibrium `L = [T]/[R]`. `Keq` is
user-supplied and `E_total` is absorbed into the rate scale.
[Mechanisms with allosteric regulators](@ref) covers the allosteric-state tags
and the partition-function structure.

Each substrate and product declares its **atom inventory** in the bracket:
`S[C]` is one carbon, and `A[C2N1]` would be two carbons and one nitrogen.
These counts let the enumerator enforce atom conservation across steps and
recognise covalent (ping-pong) intermediates; this transfer needs only a single
`[C]` placeholder.

## Simulate noiseless data

We evaluate the generator on a concentration grid. The activator `A` spans well
below to well above its binding constant, so the grid fully samples the
allosteric response, and `Keq = 10` keeps every net rate strictly positive — a
rate of exactly zero has no logarithm, and the loss works in log space.

```@example identify_fast
Keq = 10.0
true_params = (K_S_E = 1.0, K_P_E = 1.0, k_A_ES_to_EP = 5.0,
               K_A_Areg = 1.0, L = 100.0, Keq = Keq, E_total = 1.0)

concs = [(S = s, P = p, A = a)
         for s in (0.3, 1.0, 3.0, 10.0) for p in (0.1, 0.5)
         for a in (0.03, 0.3, 3.0, 30.0)]
groups = [i ≤ length(concs) ÷ 2 ? "G1" : "G2" for i in eachindex(concs)]

data = (group = groups,
        Rate  = [rate_equation(generator, c, true_params) for c in concs],
        S = [c.S for c in concs],
        P = [c.P for c in concs],
        A = [c.A for c in concs])
nothing # hide
```

The `data` table has a `:group` column, a `:Rate` column, and one column per
metabolite — substrate, product, and the regulator `A`. Each unique `group`
value becomes one cross-validation fold, so at least two groups are required.

## The data contract

`IdentifyRateEquationProblem` validates the table at construction:

- `:group` and `:Rate` columns must be present.
- One column per substrate, product, and regulator (names match
  `metabolites(mechanism)` exactly).
- Every `Rate` must be nonzero — the loss function works in log space.
- At least two distinct `group` values are required for cross-validation.

`Keq` is a required keyword argument, always user-supplied; the package never
estimates it from data. Most enzyme reactions have a known `Keq` — measure it
directly, or compute it from a resource such as
[eQuilibrator](https://equilibrator.weizmann.ac.il).

## Run the search

```@example identify_fast
using OptimizationCMAEvolutionStrategy, Random
Random.seed!(123)   # reproducible multi-start fits

prob = IdentifyRateEquationProblem(rxn, data; Keq = Keq)

results = identify_rate_equation(prob;
    optimizer       = CMAEvolutionStrategyOpt(),
    min_beam_width  = 10,
    max_param_count = 5,
    n_restarts      = 5,
    maxtime         = 4.0,
)
nothing # hide
```

The progress lines above trace the search, which is a *beam search*: it walks
parameter counts in ascending order, and at each count it fits the candidates,
keeps the most promising (the beam — here `min_beam_width = 10`), and expands
the survivors into the next count. Which mechanisms stay in the beam is set by
`min_beam_width` and the `loss_rel_threshold` / `loss_abs_threshold` /
`loss_parsimony_threshold` cutoff — see [Best mechanism selection](@ref).
Two filters bound the search: `max_param_count` (here `5`, the generating
mechanism's size, to keep the example quick) caps it by fitted-parameter count,
and `eq_complexity_filter` caps it by rate-equation complexity — roughly the
number of terms in the equation's denominator, its default passing a fully
steady-state random-order bi-bi. Each drops a candidate before fitting it. [Best mechanism selection](@ref) also covers the cross-validation
rule that picks the winner. The
production search widens the beam to 50 and the cap to 20, and distributes the
fits across workers (see
[Running in parallel](@ref)).

### Seeding from the regulator

Because the reaction declares `A`, the search seeds from mechanisms that already bind
it — every fully-regulated mechanism at its minimum parameter count — and never fits
the non-allosteric mechanisms beneath. This is the default: every declared regulator
is required, which is what lets the search reach the generating MWC mechanism so
quickly here.

Two controls tune it. To let the search decide whether `A` matters, and fit the
non-regulated mechanisms too, mark it optional:

```julia
identify_rate_equation(prob; optimizer = CMAEvolutionStrategyOpt(),
                       optional_allosteric_regulators = [:A])
```

To go the other way and shrink the seed set, declare `A`'s type. An activator binds
the active conformation, so `A::Activator` pins it to `:OnlyA` and halves the seeds an
undesignated regulator would otherwise produce:

```julia
allosteric_regulators: A::Activator
```

[The enumeration engine](@ref) explains how the seed set is built, and the
[Roadmap](@ref) tracks the moves that refine it.

## Read the result

`IdentifyRateEquationResults` has exactly two fields: `best` and `cv_results`.

```@example identify_fast
results.best
```

`results.best` is an `AbstractEnzymeMechanism`. Pass it to
[`rate_equation_string`](@ref) to see the symbolic rate equation:

```@example identify_fast
print(rate_equation_string(results.best))
```

On noiseless data the search recovers the generating MWC mechanism. The active
and inactive conformations enter as a partition function,
`(active polynomial)² + L·(inactive polynomial)²` — squared because the enzyme
is a dimer — the activator contributes `(1 + A/K_A_Areg)²` to the active state,
and the Haldane constraint fixes the dependent reverse rate `k_A_EP_to_ES` from
`Keq` and the independent constants.

`results.cv_results` is a `DataFrame` with one row per candidate equation that
entered cross-validation, scored as detailed on the [Best mechanism selection](@ref)
page:

```@example identify_fast
first(results.cv_results, 5)
```

`results.cv_results` and the selected `results.best` are also written to
`save_dir`, as `loocv_results.csv` (this whole table) and `best_equation.csv`
(the winning row: equation string plus fitted parameters). They sit alongside
the per-iteration `equation_search_iteration_N.csv` files, so a cluster run's
model-selection outcome is saved without re-running cross-validation.

## Loud failures

A mechanism that throws during compilation or fitting becomes a `FitFailure`
carrying the exception text. Failures are never silently discarded — they appear
in `cv_results` (and the saved CSVs) with the `retcode` and `error` columns
populated. If every mechanism in the base tier fails, the search re-raises the
first exception, so an unsupported optimizer keyword or a memory overflow
surfaces immediately rather than being swallowed.
