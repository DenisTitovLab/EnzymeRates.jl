# Kinetic groups

A **kinetic group** is a set of steps that share one kinetic parameter. The
point is parsimony: when several steps are chemically the same — the same
metabolite binding the enzyme, or the same chemical conversion happening in
more than one place in the mechanism — there is no reason for each to carry its
own constant. Putting them in one group makes the model say they behave
identically, and the rate equation then carries a single parameter for the
whole set instead of one per step.

By default `@enzyme_mechanism` keeps each step in its own group. To declare that
several steps share a parameter, wrap them in **parentheses** — a parenthesized
step-group binds the steps into one kinetic group:

```julia
steps: begin
    E + S ⇌ E(S)
    (E + I ⇌ E(I), E(P) + I ⇌ E(P, I))   # these two share one constant
    ...
end
```

The thermodynamic reduction ([Thermodynamic constraints](@ref)) also removes
parameters, but it does so only when a Haldane or Wegscheider relation *forces*
two constants to be equal. Kinetic grouping is a separate, modeling choice: it
ties steps together that the thermodynamics leaves independent.

## An example: a dead-end inhibitor on two enzyme forms

Take a reversible uni-uni mechanism with a dead-end inhibitor `I` that can bind
both the free enzyme `E` and the `E(P)` complex. Each `I`-binding step is in its
own group, so each gets its own dissociation constant:

```@example kingroups
using EnzymeRates
ungrouped = @enzyme_mechanism begin
    substrates: S
    products:   P
    regulators: I
    steps: begin
        E + S ⇌ E(S)
        E(S) <--> E(P)
        E(P) ⇌ E + P
        E + I ⇌ E(I)
        E(P) + I ⇌ E(P, I)
    end
end
print(rate_equation_string(ungrouped))
```

The denominator carries two inhibitor constants: `K_Iinh_E` in the `I / K_Iinh_E`
term (`I` on free `E`) and a distinct `K_Iinh_EP` in the
`I * P / (K_Iinh_EP * K_P_E)` term (`I` on `E(P)`). Nothing forces these equal —
the two are independent dead-end branches, and the thermodynamic reduction
leaves both in the parameter list.

Now group the two `I`-binding steps by wrapping them in parentheses. This
asserts that `I` binds `E` and `E(P)` with the *same* affinity:

```@example kingroups
grouped = @enzyme_mechanism begin
    substrates: S
    products:   P
    regulators: I
    steps: begin
        E + S ⇌ E(S)
        E(S) <--> E(P)
        E(P) ⇌ E + P
        (E + I ⇌ E(I), E(P) + I ⇌ E(P, I))
    end
end
print(rate_equation_string(grouped))
```

`K_Iinh_EP` is gone. Both inhibitor terms now use the single shared constant
`K_Iinh_E`, so the `E(P)` term reads `I * P / (K_Iinh_E * K_P_E)`. The grouped
mechanism has one fewer parameter:

```@example kingroups
(ungrouped = parameters(ungrouped), grouped = parameters(grouped))
```

Grouping is how the mechanism enumeration keeps the parameter count at its
lowest physically meaningful value and adds parameters only as the data warrant;
see [The enumeration engine](@ref).
