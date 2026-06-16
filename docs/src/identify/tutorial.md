# Identify tutorial

`identify_rate_equation` takes an `EnzymeReaction` and rate data, enumerates
biochemically valid mechanisms, and fits a beam of the most promising
candidates at each parameter count — fitting every mechanism is far too
expensive — returning the simplest one that generalizes by leave-one-group-out
cross-validation.

This page walks through a fully runnable example: an MWC allosteric enzyme,
recovered from noiseless data in about a minute. The full production search
widens the beam to the defaults (`min_beam_width=50`, `loss_rel_threshold=2.0`,
`loss_abs_threshold=0.01`, `max_param_count=20`) and would often run for many
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

The progress lines above trace the search. It walks parameter counts in
ascending order: at each count it fits the candidates, keeps a beam of the
most promising candidates (here `min_beam_width = 10`), and expands the
survivors into the next count. `max_param_count = 5` stops at the generating mechanism's size to
keep the example quick. [Model selection](@ref) details the beam cutoff and the
cross-validation rule that picks the winner. The production search widens the
beam to 50 and the cap to 20, and distributes the fits across workers (see
[Running in parallel](@ref)).

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
entered cross-validation, scored as detailed on the [Model selection](@ref)
page:

```@example identify_fast
first(results.cv_results, 5)
```
