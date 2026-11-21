# Kinetic groups

A **kinetic group** is a set of steps that share one kinetic parameter. The
point is parsimony: when several steps are chemically the same — the same
metabolite binding the enzyme, or the same chemical conversion happening in more
than one place — there is no reason for each to carry its own constant. Putting
them in one group makes the model say they behave identically, and the rate
equation then carries a single parameter for the whole set instead of one per
step.

By default `@enzyme_mechanism` keeps each step in its own group. To declare that
several steps share a parameter, wrap them in **parentheses** — a parenthesized
step-group binds the steps into one kinetic group:

```julia
steps: begin
    (E + Q ⇌ EQ, EA + Q ⇌ EAQ)   # both Q-binding steps share one constant
    ...
end
```

The steps of a group must take up and give off the same metabolites and carry
one flag, all rapid-equilibrium or all steady-state. Every binding of one
metabolite is the same kind of step, whether or not it runs chemistry: a fused
binding of `B` (`EA + B ⇌ EPQ`) may share a group with a plain one
(`EQ + B ⇌ EBQ`). Every isomerization is one kind, and a Theorell–Chance
step cannot share a group with a binding. Each reaction appears once, in one group.
Building a mechanism that breaks these rules raises an error.

The isomerizations of one group share their constants in the direction each is written.
Write a group of conformational changes, such as `F <--> E` and its mirror `FI <--> EI` at
the inhibitor-bound forms, with every step the same way round, from the forms of one
conformation to those of the other. Turning one of them round writes another mechanism,
in which the forward constant of one step is the reverse constant of the other; the
package cannot tell which you meant, since a form's name means nothing to it.

This is a modeling choice, distinct from the thermodynamic reduction
([Thermodynamic constraints](@ref)). The reduction also removes parameters, but
only when a Haldane or Wegscheider relation *forces* two constants to be equal.
Kinetic grouping is a decision you make: it ties together steps the
thermodynamics leaves independent.

## An example: a two-substrate, two-product mechanism

Consider a random-order reaction `A + B ⇌ P + Q` with two abortive complexes:
`EAQ`, where the enzyme binds substrate `A` together with the product `Q` of
the other half-reaction, and `EBP`, its mirror. Like the catalytic
complexes, each abortive complex forms either way — `EAQ` by `Q` binding
`EA` or `A` binding `EQ`. Every metabolite therefore binds the enzyme in
several places: `A` binds free `E`, the `EB` complex, and `EQ`, and the same
holds for `B`, `P`, and `Q`.

With every step in its own group, each of those bindings gets its own constant:

```@example kingroups
using EnzymeRates
ungrouped = @enzyme_mechanism begin
    substrates: A, B
    products:   P, Q
    steps: begin
        E + A ⇌ EA
        E + B ⇌ EB
        EA + B ⇌ EAB
        EB + A ⇌ EAB
        EAB <--> EPQ
        EPQ ⇌ EP + Q
        EPQ ⇌ EQ + P
        EP ⇌ E + P
        EQ ⇌ E + Q
        EA + Q ⇌ EAQ
        EQ + A ⇌ EAQ
        EB + P ⇌ EBP
        EP + B ⇌ EBP
    end
end
print(rate_equation_string(ungrouped))
```

There are nine independent constants, with form-specific names: `K_EA_to_E_A` for `A`
on free E, `K_EAB_to_EA_B` and `K_EAQ_to_EA_Q` for `B` and `Q` on `EA`, and similar
families for the others. Four more bindings are not fit at all but fixed by Wegscheider
relations (`K_EAB_to_EB_A`, `K_EAQ_to_EQ_A`, `K_EBP_to_EP_B`, `K_EPQ_to_EQ_P`), since
each catalytic and abortive loop closes a thermodynamic cycle.

Now group every binding of a given metabolite together — all `A`-binding steps
in one group, all `B`-binding in another, and likewise for `P` and `Q`. Each
metabolite then has a single binding constant:

```@example kingroups
grouped = @enzyme_mechanism begin
    substrates: A, B
    products:   P, Q
    steps: begin
        (E + A ⇌ EA, EB + A ⇌ EAB, EQ + A ⇌ EAQ)
        (E + B ⇌ EB, EA + B ⇌ EAB, EP + B ⇌ EBP)
        EAB <--> EPQ
        (EPQ ⇌ EQ + P, EP ⇌ E + P, EB + P ⇌ EBP)
        (EPQ ⇌ EP + Q, EQ ⇌ E + Q, EA + Q ⇌ EAQ)
    end
end
print(rate_equation_string(grouped))
```

The parameter list collapses from nine constants to five — one binding constant
per metabolite (`K_EA_to_E_A`, `K_EB_to_E_B`, `K_EP_to_E_P`, `K_EQ_to_E_Q`) plus
the catalytic `k_EAB_to_EPQ`:

```@example kingroups
(ungrouped = parameters(ungrouped), grouped = parameters(grouped))
```

The denominator becomes symmetric: the two abortive complexes read
`A * Q / (K_EA_to_E_A * K_EQ_to_E_Q)` and `B * P / (K_EB_to_E_B * K_EP_to_E_P)`,
with no form-specific suffixes. The Wegscheider section is gone, too — once
each metabolite has a single binding constant, the loop-closing relations
become identities and drop out. Grouping did the collapsing the thermodynamics
could not: it is what keeps the parameter count at the lowest physically
meaningful value, and the mechanism enumeration starts there, splitting groups
back apart only as the data warrant (see [The enumeration engine](@ref)).
