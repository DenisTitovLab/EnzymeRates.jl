# Mechanisms with regulators

A **regulator** is a species that binds the enzyme and changes its rate without
being consumed by the reaction. The derivation treats regulators no differently
from substrates and products: it walks whatever enzyme-form graph the steps
define and reads the rate law off it. A regulator's *kind* — competitive
inhibitor, essential activator, general modifier — is therefore not a setting you
choose but an emergent consequence of **where** it binds.

This makes the derivation more general than the automatic search. When
[`identify_rate_equation`](@ref) enumerates mechanisms it adds only competitive
(dead-end) inhibitors and allosteric regulators. A hand-written
[`@enzyme_mechanism`](@ref), though, can express the full range — competitive,
non-competitive, and uncompetitive inhibition, essential and non-essential
activation — simply by choosing which enzyme forms the regulator binds. Declare
the regulator with `regulators: <name>` and write its binding steps.

## A competitive inhibitor

A competitive inhibitor binds the **free enzyme** only, competing with substrate.
The inhibitor-bound form has no catalytic exit, so it appears only in the
denominator — it ties up enzyme but produces no flux:

```@example reg
using EnzymeRates
competitive = @enzyme_mechanism begin
    substrates: S
    products:   P
    regulators: I
    steps: begin
        E + S ⇌ E(S)
        E(S) <--> E(P)
        E(P) ⇌ E + P
        E + I ⇌ E(I)
    end
end
print(rate_equation_string(competitive))
```

The inhibitor adds a single `I / K_Iinh_E` term to the denominator and leaves the
numerator untouched — the textbook competitive form. (`I` renders with an `inh`
marker, `K_Iinh_E`, because it occupies a regulator binding site.)

## An essential activator

An essential activator must be bound for catalysis to occur at all: there is no
catalytic cycle on the activator-free enzyme. The activator `R` binds the free
enzyme, and every catalytic form carries it:

```@example reg
essential = @enzyme_mechanism begin
    substrates: S
    products:   P
    regulators: R
    steps: begin
        E + R ⇌ E(R)
        E(R) + S ⇌ E(S, R)
        E(S, R) <--> E(P, R)
        E(R) + P ⇌ E(P, R)
    end
end
print(rate_equation_string(essential))
```

Now the activator concentration `R` appears in **every numerator term**, so the
rate vanishes as `R → 0` — the defining behaviour of an essential activator,
opposite to the competitive inhibitor whose regulator stayed in the denominator.
Nothing in the DSL marks `R` as an activator; its role follows entirely from
binding before the substrate and riding along through catalysis.

See also: [`@enzyme_mechanism`](@ref), [`rate_equation_string`](@ref).
