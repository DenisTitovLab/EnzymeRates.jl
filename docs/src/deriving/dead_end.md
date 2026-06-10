# Dead-end inhibitor binding

A dead-end (competitive) inhibitor binds an enzyme form and prevents catalysis
without contributing to the reaction. Because the inhibitor-bound form has no
catalytic exit, it appears only in the rate-equation denominator — it consumes
enzyme but produces no flux.

## How dead-end binding enters the rate equation

A dead-end inhibitor binding step is always rapid equilibrium. The enumerator
creates it as `Step(free_form, inhibited_form, CompetitiveInhibitor(name), true)`
— the `true` marks rapid equilibrium. Because the inhibitor-bound form is a
dead end, the King–Altman/Cha derivation adds a term to the denominator for each
such form, reducing the fraction of enzyme available for the catalytic cycle.

For a single competitive inhibitor `I` on free enzyme `E`, the effect is one
additional term in the denominator:

```
+ I / K_Iinh_E
```

leaving the numerator unchanged.

## Declaring inhibitors

In [`@enzyme_reaction`](@ref), list inhibitors under `dead_end_inhibitors:` or
`competitive_inhibitors:` — both map to role `:competitive`. In a hand-written
[`@enzyme_mechanism`](@ref), declare the inhibitor under `regulators:` and tag
the inhibitor-bound form with `::Inh` in the step:

```@example deadend
using EnzymeRates
de = @enzyme_mechanism begin
    substrates: S
    products:   P
    regulators: I
    steps: begin
        E + S ⇌ E(S)
        E(S) <--> E(P)
        E(P) ⇌ E + P
        E + I ⇌ E(I::Inh)
    end
end
print(rate_equation_string(de))
```

The `::Inh` tag binds `I` in its `CompetitiveInhibitor` role while keeping the
real metabolite name, so `concs.I` drives it at evaluation time. The rendered
parameter name carries an `inh` marker — `K_Iinh_E` — keeping the
inhibitor-bound form `:EIinh` distinct from a hypothetical product-bound form
`:EI`.

## The parameter list

```@example deadend
parameters(de)
```

`parameters(de)` returns `(:K_Iinh_E, :K_P_E, :K_S_E, :k_ES_to_EP, :Keq,
:E_total)`. The only new parameter relative to the inhibitor-free textbook
mechanism is `K_Iinh_E` — the dissociation constant for `I` on free enzyme.
The numerator and all catalytic parameters are identical to the uninhibited case.

## Multiple inhibitors on the same site

When multiple regulators bind the same site, the denominator factor takes the
form `(1 + R1/K_R1 + R2/K_R2)^m`, where `m` is the site multiplicity. For
competitive inhibitors the multiplicity is 1, so each inhibitor contributes an
independent additive term in the denominator.

## Enumeration

When searching for the best rate equation, the enumeration engine can add a
dead-end regulator to any enzyme form as one of its expansion moves. See the
mechanism enumeration page for details on how `expand_mechanisms` generates
inhibitor variants.

See also: [`@enzyme_mechanism`](@ref), [`rate_equation_string`](@ref).
