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
coupling constant `L` sets the inactive-to-active population ratio of the bare
enzyme. What distinguishes an activator from an inhibitor is simply which
conformation it prefers: an **activator binds the active state more tightly** and
shifts the population toward A, while an **inhibitor binds the inactive state**
and shifts it toward I. That preference is expressed entirely through the
regulator's binding constants for the two states — its `K_A` versus its `K_I`.

## Setting up an MWC mechanism

Use [`@allosteric_mechanism`](@ref), tagging each catalytic step and each
regulator with its allosteric-state behavior:

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
parameters(allo)
```

## Allosteric-state tags and multiplicities

Each catalytic kinetic group and each regulator carries one of four tags that
controls how its A- and I-state symbols are rendered:

| Tag | Meaning | Example symbol |
|-----|---------|----------------|
| `:EqualAI` | One shared symbol for both states | `K_S_E` |
| `:NonequalAI` | Independent A and I symbols | `k_A_ES_to_EP`, `k_I_ES_to_EP` |
| `:OnlyA` | Present in the active state only; zeroed in the inactive polynomial | `K_A_Areg` |
| `:OnlyI` | Present in the inactive state only; zeroed in the active polynomial | `K_I_Ireg` |

In the example above the two binding steps are `:EqualAI`, so they share `K_S_E`
and `K_P_E` with no A/I token; the isomerization step is `:NonequalAI`, giving
independent `k_A_ES_to_EP` and `k_I_ES_to_EP`; the activator `A` is `:OnlyA`
(`K_A_Areg`, active state only) and the inhibitor `I` is `:OnlyI` (`K_I_Ireg`,
inactive state only). `L` is the bare-enzyme coupling constant.

A catalytic step can never be `:OnlyI` — the active state must always be able to
catalyze — so tagging one that way is a hard constructor error:

```@example mwc
try
    @allosteric_mechanism begin
        substrates: S
        products:   P
        catalytic_steps: begin
            E + S ⇌ E(S)   :: OnlyI
            E(S) <--> E(P)  :: EqualAI
            E(P) ⇌ E + P   :: EqualAI
        end
    end
catch err
    showerror(stdout, err)
end
```

The mechanism also declares **multiplicities**, which capture the enzyme's
oligomeric state. `catalytic_multiplicity` is the number of identical catalytic
subunits; it becomes the exponent on each conformation's catalytic partition
function. Each `regulatory_site(multiplicity = m)` says how many copies of that
site the oligomer carries, raising the site's binding factor to the power `m`. In
the example both are `2` — a homodimer with two copies of each regulatory site.

## Derivation of the MWC rate equation

The MWC rate equation is assembled from two ordinary single-conformation rate
equations — one for the active state, one for the inactive state — combined by a
partition function. The derivation rests on the MWC assumption that **both the
binding of allosteric regulators and the transitions between the A and I states
are rapid-equilibrium**: each conformation equilibrates quickly relative to
catalysis, so the overall rate is a population-weighted average of the two
conformations' catalytic rates. Treating regulator binding as steady state
instead would not help — a regulator is not consumed, so a non-RE binding step
contributes only the ratio of its on- and off-rate constants, which is exactly
the equilibrium constant the RE treatment already uses, and the two rates would
not be separately identifiable from rate data.

For a mechanism with catalytic multiplicity `cat_n`, the package builds

```
v = E_total * num / den
num = cat_n * (N_A * Q_A^(cat_n-1) * W_A + L * N_I * Q_I^(cat_n-1) * W_I)
den =          Q_A^cat_n * W_A           + L * Q_I^cat_n           * W_I
```

where `N_A`, `Q_A` are the active-state catalytic numerator and partition
function (the King–Altman/Cha polynomials for the catalytic cycle); `N_I`, `Q_I`
are the inactive-state counterparts, obtained by zeroing `:OnlyA` symbols and
renaming `:NonequalAI` symbols to their I-state names; `W_A`, `W_I` are products
of the regulatory-site factors `(1 + lig/K)^m`; and `L` weights the inactive
branch throughout. When a catalytic group is `:OnlyA`, the inactive cycle cannot
close, so the package sets `N_I = 0`: the inactive branch still contributes to
the denominator as enzyme mass but carries no forward flux.

```@example mwc
print(rate_equation_string(allo))
```

Here `cat_n = 2`, so each partition function is squared. The active branch
carries `(1 + A / K_A_Areg) ^ 2` and a trivial `1 ^ 2` from the inhibitor site
(absent in the active state); the inactive branch mirrors it with
`(1 + I / K_I_Ireg) ^ 2`.

!!! warning "Known rendering bug"
    For an allosteric mechanism, `rate_equation_string` currently renders the
    **active-state** Haldane constraint lines with un-prefixed names — for
    example `k_EP_to_ES = (1 / Keq) * K_P_E * (1 / K_S_E) * k_ES_to_EP` — instead
    of the active-state names `k_A_EP_to_ES` and `k_A_ES_to_EP`. Those names
    appear in neither the `params` list nor the `v` expression, so the printed
    allosteric string is not directly runnable. The compiled
    [`rate_equation`](@ref) is unaffected — it computes the dependent constants
    correctly — so this is purely a display bug in `rate_equation_string`, and it
    will be fixed.

## See also

- [`@allosteric_mechanism`](@ref) — DSL reference
- [`parameters`](@ref) — list independent parameter names
- [`rate_equation_string`](@ref) — inspect the derived equation
