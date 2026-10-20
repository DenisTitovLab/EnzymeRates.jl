# Mechanisms with allosteric regulators

Allosteric regulation is one of the main ways metabolic enzymes are controlled.
An **allosteric regulator** binds a regulatory site that is distinct from the
catalytic site and, from there, either inhibits or activates the enzyme. A
classical description of how this works is the Monod–Wyman–Changeux (MWC) model
[Monod1965](@cite), and EnzymeRates can derive rate equations under that model.

In the MWC picture the enzyme exists in two conformations — an **active** state
(A) and an **inactive** state (I) — in equilibrium with one another. (In the
classic MWC notation A is the R state and I is the T state; the package uses A/I
throughout, so `K_A_…` names correspond to the literature's R-state constants.) A
coupling constant `L` sets the inactive-to-active population ratio of the free
enzyme. What distinguishes an activator from an inhibitor is simply which
conformation it prefers: an **activator binds the active state more tightly** and
shifts the population toward A, while an **inhibitor binds the inactive state**
and shifts it toward I. That preference is expressed entirely through the
regulator's binding constants for the two states — its `K_A` versus its `K_I`.

## Setting up an MWC mechanism

Use [`@allosteric_mechanism`](@ref), tagging each catalytic step and each
regulator with its allosteric-state behavior. The printed rate equation shows the
result:

```@example mwc
using EnzymeRates
allo = @allosteric_mechanism begin
    substrates: S
    products:   P
    catalytic_multiplicity: 2
    allosteric_regulators: A, I
    steps: begin
        E + S ⇌ E(S)        :: EqualAI
        E(S) <--> E(P)       :: NonequalAI
        E(P) ⇌ E + P        :: EqualAI
    end
    regulatory_site(2): A::OnlyA
    regulatory_site(2): I::OnlyI
end
print(rate_equation_string(allo))
```

The independent parameters the mechanism fits:

```@example mwc
parameters(allo)
```

Each catalytic kinetic group and each regulator carries one of four tags that
controls how its A- and I-state symbols are rendered:

| Tag | Meaning | Example symbol |
|-----|---------|----------------|
| `:EqualAI` | One shared symbol for both states | `K_ES_to_E_S` |
| `:NonequalAI` | Independent A and I symbols | `k_A_ES_to_EP`, `k_I_ES_to_EP` |
| `:OnlyA` | Present in the active state only; zeroed in the inactive polynomial | `K_A_Areg` |
| `:OnlyI` | Present in the inactive state only; zeroed in the active polynomial | `K_I_Ireg` |

In the example the two binding steps are `:EqualAI`, so they share `K_ES_to_E_S`
and `K_EP_to_E_P` with no A/I token; the isomerization step is `:NonequalAI`, giving
independent `k_A_ES_to_EP` and `k_I_ES_to_EP`; the activator `A` is `:OnlyA`
(`K_A_Areg`, active state only) and the inhibitor `I` is `:OnlyI` (`K_I_Ireg`,
inactive state only). An `:EqualAI` regulator `R` would take one bare constant, `K_Rreg`,
in both states. `L` is the free-enzyme coupling constant.
A catalytic group takes `:OnlyI` only when its steps bind modifiers and nothing else (a
name under `modifiers:` or a reactant's `::Mod` copy), such as an inhibitor that binds only
the inactive conformation. The active conformation must still catalyze, so the
constructor rejects `:OnlyI` on a step that takes up or gives off a substrate or product,
or changes what else its form holds.

The mechanism also declares **multiplicities**, which capture the enzyme's
oligomeric state. `catalytic_multiplicity` is the number of identical catalytic
subunits; it becomes the exponent on each conformation's catalytic partition
function. Each `regulatory_site(m)` says how many copies of that
site the oligomer carries, raising the site's binding factor to the power `m`. In
the example both are `2` — a homodimer with two copies of each regulatory site.

!!! note
    The derivation reports the rate **per catalytic site**: `catalytic_multiplicity`
    enters only through the partition-function exponents, never as an overall
    multiplier on the numerator.

## Derivation of the MWC rate equation

The MWC rate equation is assembled from two ordinary single-conformation rate
equations — one for the active state, one for the inactive state — combined by a
partition function. The derivation rests on two rapid-equilibrium assumptions.
Allosteric regulators bind fast relative to catalysis. And **the enzyme changes
conformation only while free**: it settles into the active or inactive state
before it binds anything — the free-enzyme inactive-to-active ratio is `L` — and
keeps that conformation for a full catalytic cycle. The overall rate is then each
conformation's catalytic cycle weighted by its free-enzyme population.
The conformations interconvert only at free E, the form `E` with nothing bound, and `L` is
its inactive-to-active ratio. A mechanism may hold other free forms, such as a
conformation `F` with `F <--> E`; they change conformation only by way of E. When another
free form sits at rapid equilibrium with E, its own inactive-to-active ratio is `L·r_I/r_A`,
where `r_c` is the equilibrium ratio `[F]/[E]` in conformation c across the step that joins
F to E, so a fit reports `L` at E whatever the free forms are named.

An equivalent statement lets every form interconvert between the A and I states,
but with a thermodynamically fixed ratio: `L` for the free enzyme, and `L` scaled
by the form's active-to-inactive binding-constant ratio for a bound form (so `L`
itself for an `:EqualAI` ligand). The two statements agree except when catalysis
differs between the conformations (`:NonequalAI`); there, committing to one
conformation per cycle — the model here — routes each molecule's turnover through
a single conformation rather than blending the two, a difference of a few percent.

Treating regulator binding as steady state instead would not help — a regulator
is not consumed, so a non-RE binding step contributes only the ratio of its on-
and off-rate constants, which is exactly the equilibrium constant the RE
treatment already uses, and the two rates would not be separately identifiable
from rate data.

For a mechanism with catalytic multiplicity `cat_n`, the package normalizes each
conformation to its own free-enzyme population and combines the two:

```
v = E_total * num / den
num = (N_A/D_A) * (Q_A/D_A)^(cat_n-1) * W_A + L * (N_I/D_I) * (Q_I/D_I)^(cat_n-1) * W_I
den = (Q_A/D_A)^cat_n * W_A                 + L * (Q_I/D_I)^cat_n                 * W_I
```

`N_A`, `Q_A` are the active state's catalytic numerator and partition function
(the King–Altman/Cha polynomials for the catalytic cycle, on the active-state
graph with its `:OnlyI` groups pruned), and `D_A` is free E's
weight in `Q_A`: the spanning-tree weight of its free-enzyme segment, times the
concentration monomial with which the derivation clears reciprocal concentrations.
`N_I`, `Q_I`, `D_I` are the
inactive-state counterparts, derived the same way from the inactive-state graph —
`:OnlyA` groups pruned, `:EqualAI` groups keeping the shared active-state
constants, and `:NonequalAI` groups carrying their own I-state names. `W_A`,
`W_I` are products of the regulatory-site factors `(1 + lig/K)^m`, and `L`
weights the inactive branch throughout.

The free-enzyme weight `D` is `1` whenever the catalytic graph is a single
rapid-equilibrium segment without a ping-pong covalent form — the common case,
where `Q/D = Q` and the equation is the textbook MWC form. `D` differs from `1`
when an `:OnlyA`/`:OnlyI` binding fragments the graph, a steady-state binding
splits it into segments, or a covalent form shares free E's segment: such a form
sits at a reciprocal concentration relative to free E (`Q/B`, say), and clearing
it gives free E the weight `B`. The normalization then keeps the two
conformations on a common free-enzyme basis,
which is what stops a bare rate constant from leaking into the `L`-term. The
package divides `Q/D` directly when `D` is a single rate constant — giving the
readable `1 + [S]/K + …` form — and clears the fraction to a polynomial when `D`
carries a metabolite.
When the `:OnlyA` groups stop every catalytic cycle of the inactive conformation, as the
tagging move's bases do, its numerator `N_I` is 0: the inactive branch still contributes
to the denominator as enzyme mass but carries no forward flux.

In the printed equation above `cat_n = 2`, so each partition function is squared:
the active branch carries `(1 + A / K_A_Areg) ^ 2` and a trivial `1 ^ 2` from the
inhibitor site (absent in the active state), and the inactive branch mirrors it
with `(1 + I / K_I_Ireg) ^ 2`.

## Tags as limits, and admission

An `:OnlyA` tag says what the inactive conformation lacks, as the limit of ordinary
constants; `:OnlyI` says the same of the active conformation:

- a plain binding loses its affinity: its dissociation constant grows without bound, or, for
  a steady-state binding, its binding rate constant falls to zero with its release rate
  fixed;
- any other step (an isomerization, a fused binding, a Theorell–Chance step) is blocked:
  both rate constants fall to zero with their ratio held by the thermodynamic cycles, and a
  rapid-equilibrium step keeps the constant the cycles set while its speed vanishes;
- each step of a tagged group vanishes at its own rate, and a group may mix a plain and a
  fused binding of one metabolite.

The derivation reads each conformation's law by deleting the steps it lacks (the `:OnlyA`
groups in the inactive conformation, the `:OnlyI` groups in the active one) and the forms
that the deletion cuts off from free E. **Admission** decides when that is the limit,
whatever the rates at which the tagged constants vanish:

1. **Existence.** Join the forms by the steps the conformation keeps (outside the groups it
   lacks), and draw each tagged plain binding as an arrow from its unbound side's component
   to its bound side's. A path to the limit exists exactly when no arrow stays inside one
   component and the arrows close no cycle.
2. **Uniqueness.** Every form that free E reaches only through tagged steps lies behind a
   tagged plain binding, on its bound side, and the limit is the same on every face where
   data sit (below).
3. **Turnover.** With `:OnlyI` groups, the active conformation still turns the reaction
   over.

A tag set that fails is rejected with an error that names the step. Tagging `E + S ⇌ E(S)`
`:OnlyA` while the catalytic step and the product binding stay `:EqualAI` fails existence:
in the inactive conformation those steps still join E to E(S), so no path takes that
affinity to zero. Admission is sound but not exact. It also rejects two rare shapes whose
limit exists, equilibrium bridges and stranded basins fed only through vanishing steps,
which the repository's `ROADMAP.md` lists.

### What the package rejects

A mechanism is well defined when its rate law is the unique fast-equilibrium limit of
mass-action kinetics wherever data sit: at interior concentrations and on every **face**,
where every product is zero, every substrate, each reactant alone on a side of two or more,
or any of these with every modifier-only name zero. The macros throw an error that names
the step and the face; `init_mechanisms` and `expand_mechanisms` drop the mechanism. The
package rejects:

1. no catalytic cycle: no cycle of steps turns substrates into products, so the rate is
   zero;
2. a cycle of rapid-equilibrium steps that turns the reaction over, so the rate is
   infinite;
3. a face where the law is not the unique limit: a rapid-equilibrium segment with no form
   whose weight keeps the segment's lowest power of each of the face's metabolites, a
   segment whose lowest forms are joined only through other forms, a trap that the face
   feeds and its fast exits cannot drain, or a second closed set of forms;
4. tags with no limit (admission, above);
5. a law with no value (0/0) on a face, which the macros reject; the derivation checks the
   derived denominator once more.

The faces leave out the corners that zero one substrate with every product. A law the
package accepts can be 0/0 there, where a well-defined law is 0, so `rate_equation`
returns `NaN` at such a point; a fit scores the NaN as it scores the 0 a well-defined law
predicts there.

## Thermodynamic constraints of MWC equations

A single conformation already carries thermodynamic constraints: every cycle of
steps must agree with the overall equilibrium constant `Keq` (a **Haldane**
relation), and every cycle of binding steps must be path-independent (a
**Wegscheider** relation). The derivation reports both and expresses the
dependent rate constants through the independent ones. An MWC equation solves the
relations of both conformations together, adding constraints neither has alone. Here
we review several examples of such additional constraints.

Tag a lone substrate binding `:NonequalAI` while catalysis and product release
stay `:EqualAI`, and the split collapses:

```@example mwc
collapse = @allosteric_mechanism begin
    substrates: S
    products:   P
    catalytic_multiplicity: 2
    allosteric_regulators: R
    steps: begin
        E + S ⇌ E(S)      :: NonequalAI
        E(S) <--> E(P)    :: EqualAI
        E(P) ⇌ E + P      :: EqualAI
    end
    regulatory_site(2): R::OnlyI
end
print(rate_equation_string(collapse))
```

The equation reports `K_I_ES_to_E_S = K_A_ES_to_E_S`. The Haldane relation fixes `Keq` from
the catalytic rate constants and the two affinities in each conformation; with
catalysis, product release, and `Keq` all shared, the substrate affinity is
pinned as well. The tag asked for a difference thermodynamics forbids, so the
substrate stops distinguishing the states.

Tag **both** the substrate and product bindings `:NonequalAI`, and one degree of
freedom survives:

```@example mwc
coupled = @allosteric_mechanism begin
    substrates: S
    products:   P
    catalytic_multiplicity: 2
    allosteric_regulators: R
    steps: begin
        E + S ⇌ E(S)      :: NonequalAI
        E(S) <--> E(P)    :: EqualAI
        E(P) ⇌ E + P      :: NonequalAI
    end
    regulatory_site(2): R::OnlyI
end
print(rate_equation_string(coupled))
```

Now `K_I_ES_to_E_S = K_A_ES_to_E_S · K_I_EP_to_E_P / K_A_EP_to_E_P`, i.e.
`K_S^I / K_S^A = K_P^I / K_P^A`:
the affinities differ freely as long as they differ *together*, holding their
ratio fixed. This is the thermodynamically consistent form of a **K-system**, and
it needs two coupled `:NonequalAI` bindings, not one.

Often a K-system uses exclusive `:OnlyA` binding: the substrate binds the active
state alone, and the inactive cycle is pruned (`N_I = 0`). Its catalytic step
**must** be `:OnlyA` too. This is a thermodynamic requirement, not a convention: with
`:EqualAI` or `:NonequalAI` catalysis and an `:EqualAI` product binding, the inactive
conformation still joins E to E(S) through catalysis and product release, so no path takes
the substrate affinity to zero while every cycle keeps `Keq`. Admission's existence test
rejects the combination.

The exception is a balanced K-system. When `:OnlyA` binds both a substrate and a
product, the two `K_I → ∞` limits cancel in the Haldane relation: the affinities
diverge together and their ratio stays free, so `:EqualAI` catalysis is legal.

A steady-state binding splits the constraint further, carrying an affinity (a
ratio of its two rate constants) the cycles constrain and a speed (their
product) they do not: a forbidden affinity collapses while the speed stays
free, so the two conformations bind with the same `Kd` but different kinetics.
Tagging the substrate binding steady-state and `:NonequalAI` derives its
reverse rate (`k_I_ES_to_E_S = k_A_ES_to_E_S · k_I_E_S_to_ES / k_A_E_S_to_ES`)
and leaves the forward rate a free parameter:

```@example mwc
ss = @allosteric_mechanism begin
    substrates: S
    products:   P
    catalytic_multiplicity: 2
    allosteric_regulators: R
    steps: begin
        E + S <--> E(S)    :: NonequalAI
        E(S) <--> E(P)     :: EqualAI
        E(P) ⇌ E + P       :: EqualAI
    end
    regulatory_site(2): R::OnlyI
end
parameters(ss)
```

## See also

- [`@allosteric_mechanism`](@ref) — DSL reference
- [`parameters`](@ref) — list independent parameter names
- [`rate_equation_string`](@ref) — inspect the derived equation
