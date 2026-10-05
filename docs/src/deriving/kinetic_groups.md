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
    (E + Q ⇌ E(Q), E(A) + Q ⇌ E(A, Q))   # both Q-binding steps share one constant
    ...
end
```

The steps of a group must take up and give off the same metabolites and carry
one flag, all rapid-equilibrium or all steady-state. Every binding of one
metabolite is the same kind of step, whether or not it runs chemistry: a fused
binding of `B` (`E(A) + B ⇌ E(P, Q)`) may share a group with a plain one
(`E(Q) + B ⇌ E(B, Q)`). Every isomerization is one kind, and a Theorell–Chance
step cannot share a group with a binding. Each reaction appears once, in one
group. Building a mechanism that breaks these rules raises an error.

This is a modeling choice, distinct from the thermodynamic reduction
([Thermodynamic constraints](@ref)). The reduction also removes parameters, but
only when a Haldane or Wegscheider relation *forces* two constants to be equal.
Kinetic grouping is a decision you make: it ties together steps the
thermodynamics leaves independent.

## An example: a two-substrate, two-product mechanism

Consider a random-order reaction `A + B ⇌ P + Q` with two abortive complexes:
`E(A, Q)`, where the enzyme binds substrate `A` together with the product `Q` of
the other half-reaction, and `E(B, P)`, its mirror. Like the catalytic
complexes, each abortive complex forms either way — `E(A, Q)` by `Q` binding
`E(A)` or `A` binding `E(Q)`. Every metabolite therefore binds the enzyme in
several places: `A` binds free `E`, the `E(B)` complex, and `E(Q)`, and the same
holds for `B`, `P`, and `Q`.

With every step in its own group, each of those bindings gets its own constant:

```@example kingroups
using EnzymeRates
ungrouped = @enzyme_mechanism begin
    substrates: A, B
    products:   P, Q
    steps: begin
        E + A ⇌ E(A)
        E + B ⇌ E(B)
        E(A) + B ⇌ E(A, B)
        E(B) + A ⇌ E(A, B)
        E(A, B) <--> E(P, Q)
        E(P, Q) ⇌ E(P) + Q
        E(P, Q) ⇌ E(Q) + P
        E(P) ⇌ E + P
        E(Q) ⇌ E + Q
        E(A) + Q ⇌ E(A, Q)
        E(Q) + A ⇌ E(A, Q)
        E(B) + P ⇌ E(B, P)
        E(P) + B ⇌ E(B, P)
    end
end
print(rate_equation_string(ungrouped))
```

There are nine independent constants, with form-specific names: `K_EA_to_E_A`,
`K_EAB_to_EB_A`, and `K_EAQ_to_EQ_A` for `A` on three different forms, and
similar families for the others. Four more bindings are not fit at all but
fixed by Wegscheider relations (`K_EAB_to_EA_B`, `K_EBP_to_EB_P`,
`K_EAQ_to_EA_Q`, `K_EPQ_to_EP_Q`), since each catalytic and abortive loop
closes a thermodynamic cycle.

Now group every binding of a given metabolite together — all `A`-binding steps
in one group, all `B`-binding in another, and likewise for `P` and `Q`. Each
metabolite then has a single binding constant:

```@example kingroups
grouped = @enzyme_mechanism begin
    substrates: A, B
    products:   P, Q
    steps: begin
        (E + A ⇌ E(A), E(B) + A ⇌ E(A, B), E(Q) + A ⇌ E(A, Q))
        (E + B ⇌ E(B), E(A) + B ⇌ E(A, B), E(P) + B ⇌ E(B, P))
        E(A, B) <--> E(P, Q)
        (E(P, Q) ⇌ E(Q) + P, E(P) ⇌ E + P, E(B) + P ⇌ E(B, P))
        (E(P, Q) ⇌ E(P) + Q, E(Q) ⇌ E + Q, E(A) + Q ⇌ E(A, Q))
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
