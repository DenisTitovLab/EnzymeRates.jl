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
    allosteric_regulators: A::OnlyA, I::OnlyI
    catalytic_steps: begin
        E + S ⇌ E(S)        :: EqualAI
        E(S) <--> E(P)       :: NonequalAI
        E(P) ⇌ E + P        :: EqualAI
    end
    regulatory_site(multiplicity = 2): begin
        ligands: A
    end
    regulatory_site(multiplicity = 2): begin
        ligands: I
    end
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
| `:EqualAI` | One shared symbol for both states | `K_S_E` |
| `:NonequalAI` | Independent A and I symbols | `k_A_ES_to_EP`, `k_I_ES_to_EP` |
| `:OnlyA` | Present in the active state only; zeroed in the inactive polynomial | `K_A_Areg` |
| `:OnlyI` | Present in the inactive state only; zeroed in the active polynomial | `K_I_Ireg` |

In the example the two binding steps are `:EqualAI`, so they share `K_S_E` and
`K_P_E` with no A/I token; the isomerization step is `:NonequalAI`, giving
independent `k_A_ES_to_EP` and `k_I_ES_to_EP`; the activator `A` is `:OnlyA`
(`K_A_Areg`, active state only) and the inhibitor `I` is `:OnlyI` (`K_I_Ireg`,
inactive state only). `L` is the free-enzyme coupling constant. A catalytic step
can never be `:OnlyI` — the active state must always be able to catalyze — so
tagging one that way is a hard constructor error.

The mechanism also declares **multiplicities**, which capture the enzyme's
oligomeric state. `catalytic_multiplicity` is the number of identical catalytic
subunits; it becomes the exponent on each conformation's catalytic partition
function. Each `regulatory_site(multiplicity = m)` says how many copies of that
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
(the King–Altman/Cha polynomials for the catalytic cycle), and `D_A` is the
spanning-tree weight of its free-enzyme segment. `N_I`, `Q_I`, `D_I` are the
inactive-state counterparts, derived the same way from the inactive-state graph —
`:OnlyA` groups pruned, `:EqualAI` groups keeping the shared active-state
constants, and `:NonequalAI` groups carrying their own I-state names. `W_A`,
`W_I` are products of the regulatory-site factors `(1 + lig/K)^m`, and `L`
weights the inactive branch throughout.

The free-enzyme weight `D` is `1` whenever the catalytic graph is a single
rapid-equilibrium segment — the common case, where `Q/D = Q` and the equation is
the textbook MWC form. `D` differs from `1` only when an `:OnlyA`/`:OnlyI`
binding fragments the graph or a steady-state binding splits it into segments;
the normalization then keeps the two conformations on a common free-enzyme basis,
which is what stops a bare rate constant from leaking into the `L`-term. The
package divides `Q/D` directly when `D` is a single rate constant — giving the
readable `1 + [S]/K + …` form — and clears the fraction to a polynomial when `D`
carries a metabolite. When a catalytic group is `:OnlyA`, the inactive cycle
cannot close, so the package sets `N_I = 0`: the inactive branch still
contributes to the denominator as enzyme mass but carries no forward flux.

In the printed equation above `cat_n = 2`, so each partition function is squared:
the active branch carries `(1 + A / K_A_Areg) ^ 2` and a trivial `1 ^ 2` from the
inhibitor site (absent in the active state), and the inactive branch mirrors it
with `(1 + I / K_I_Ireg) ^ 2`.

## Thermodynamic constraints of MWC equations

A single conformation already carries thermodynamic constraints: every cycle of
steps must agree with the overall equilibrium constant `Keq` (a **Haldane**
relation), and every cycle of binding steps must be path-independent (a
**Wegscheider** relation). The derivation reports both and expresses the
dependent rate constants through the independent ones. An MWC equation stacks two
such equations, one per conformation, adding constraints neither has alone. Here
we review several examples of such additional constraints.

Tag a lone substrate binding `:NonequalAI` while catalysis and product release
stay `:EqualAI`, and the split collapses:

```@example mwc
collapse = @allosteric_mechanism begin
    substrates: S
    products:   P
    catalytic_multiplicity: 2
    allosteric_regulators: R::OnlyI
    catalytic_steps: begin
        E + S ⇌ E(S)      :: NonequalAI
        E(S) <--> E(P)    :: EqualAI
        E(P) ⇌ E + P      :: EqualAI
    end
    regulatory_site(multiplicity = 2): begin
        ligands: R
    end
end
print(rate_equation_string(collapse))
```

The equation reports `K_I_S_E = K_A_S_E`. The Haldane relation fixes `Keq` from
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
    allosteric_regulators: R::OnlyI
    catalytic_steps: begin
        E + S ⇌ E(S)      :: NonequalAI
        E(S) <--> E(P)    :: EqualAI
        E(P) ⇌ E + P      :: NonequalAI
    end
    regulatory_site(multiplicity = 2): begin
        ligands: R
    end
end
print(rate_equation_string(coupled))
```

Now `K_I_P_E = K_A_P_E · K_I_S_E / K_A_S_E`, i.e. `K_P^I / K_S^I = K_P^A / K_S^A`:
the affinities differ freely as long as they differ *together*, holding their
ratio fixed. This is the thermodynamically consistent form of a **K-system**, and
it needs two coupled `:NonequalAI` bindings, not one.

Often a K-system uses exclusive `:OnlyA` binding: the substrate binds the active
state alone, and the inactive cycle is pruned (`N_I = 0`). Its catalytic step
**must** be `:OnlyA` too. This is a thermodynamic requirement, not a convention:
the two conformations share `Keq`, so a ligand that binds only the active state
(`K_I → ∞`) forces the inactive catalytic rate to zero, and an `:EqualAI` or
`:NonequalAI` catalytic tag would assert a rate the Haldane relation forbids. The
constructor rejects the inconsistent combination.

The exception is a balanced K-system. When `:OnlyA` binds both a substrate and a
product, the two `K_I → ∞` limits cancel in the Haldane relation: the affinities
diverge together and their ratio stays free, so `:EqualAI` catalysis is legal.

A steady-state binding splits the constraint further, carrying an affinity
(`kon/koff`) the cycles constrain and a speed (`kon·koff`) they do not: a
forbidden affinity collapses while the speed stays free, so the two conformations
bind with the same `Kd` but different kinetics. Tagging the substrate binding
steady-state and `:NonequalAI` derives its reverse rate
(`koff_I_S_E = koff_A_S_E · kon_I_S_E / kon_A_S_E`) and leaves the forward rate a
free parameter:

```@example mwc
ss = @allosteric_mechanism begin
    substrates: S
    products:   P
    catalytic_multiplicity: 2
    allosteric_regulators: R::OnlyI
    catalytic_steps: begin
        E + S <--> E(S)    :: NonequalAI
        E(S) <--> E(P)     :: EqualAI
        E(P) ⇌ E + P       :: EqualAI
    end
    regulatory_site(multiplicity = 2): begin
        ligands: R
    end
end
parameters(ss)
```

## See also

- [`@allosteric_mechanism`](@ref) — DSL reference
- [`parameters`](@ref) — list independent parameter names
- [`rate_equation_string`](@ref) — inspect the derived equation
