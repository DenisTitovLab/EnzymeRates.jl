# Thermodynamic constraints

Every reversible mechanism contains cycles: sequences of steps that return the
enzyme to its starting form. Thermodynamics constrains the rate constants around
each cycle — their product must equal a fixed function of the overall equilibrium
constant. The package identifies these constraints automatically and removes one
rate constant per independent cycle, reducing the number of fitted parameters and
guaranteeing that every fitted equation is thermodynamically consistent.

## Haldane and Wegscheider relations

Two kinds of constraint arise [Haldane1930](@cite) [Wegscheider1901](@cite):

- A **Haldane constraint** appears when a cycle's net metabolite change is
  proportional to the overall reaction. The dependent rate constant is expressed
  as a power product involving `Keq`. For example, in the textbook reversible
  Michaelis–Menten mechanism the reverse SS rate is fully determined by `Keq`,
  the two binding constants, and the forward SS rate:

  ```
  k_EP_to_ES = (1 / Keq) * K_P_E * (1 / K_S_E) * k_ES_to_EP
  ```

- A **Wegscheider constraint** appears when a cycle is closed internally (zero
  net metabolite change). The dependent rate constant is expressed as a pure
  power product of other rate constants, with no `Keq` factor. Wegscheider
  constraints arise in mechanisms where two thermodynamically equivalent binding
  paths exist — for example, when the same enzyme form is reachable by binding
  either substrate first. When only two binding constants are linked, the
  constraint collapses to a direct equality (`K_A_EQ = K_A_E`), and the
  absorbed symbol is substituted silently into the rate equation.

## How the constraints are found

The package represents each step as a column in the **enzyme incidence matrix**,
where rows are enzyme forms. Each column has `+1` at the to-form and `-1` at the
from-form. The **exact-integer null space** of this matrix (`_integer_nullspace`,
`src/thermodynamic_constr_for_rate_eq_derivation.jl`) gives a basis of
independent cycles. The computation uses `Rational{BigInt}` throughout — no
floating point — so each null-space vector is reduced to a primitive integer
vector with an exact sign convention.

Each cycle is then classified by its stoichiometry: if its stoichiometry vector is
proportional to the net reaction's stoichiometry, it is a Haldane cycle; if its
stoichiometry is zero, it is a Wegscheider cycle. Any other result is a hard
error — the mechanism is thermodynamically contradictory and the package raises
rather than producing a silent or wrong equation.

## Which constant is made dependent

Gaussian elimination selects one rate constant per cycle to eliminate. The pivot
priority (`_step_priority`) places internal isomerizations first (most
eliminable), then metabolite steps on non-free enzyme forms, then free-enzyme SS
binding steps, and finally free-enzyme RE binding steps (least eliminable, since
they are the structurally primary parameters). Because step order is
canonicalized in the `Mechanism` constructor, the dependent-parameter choice is
deterministic: two mechanisms with the same structure but written in different
step order produce the identical reduced equation.

A thermodynamically contradictory mechanism — where elimination reduces a
constraint row to `0 = c * log(Keq)` with `c ≠ 0` — raises an error
immediately.

## A concrete example

```@example thermo
using EnzymeRates
m = @enzyme_mechanism begin
    substrates: S
    products:   P
    steps: begin
        E + S ⇌ E(S)
        E(S) <--> E(P)
        E(P) ⇌ E + P
    end
end
print(rate_equation_string(m))
```

The `# Haldane constraints:` section shows that `k_EP_to_ES` is eliminated: the
reverse isomerization rate is determined by `Keq`, the two binding constants,
and the forward rate `k_ES_to_EP`. The parameter list therefore omits
`k_EP_to_ES` — it is not a fitted parameter. `parameters(m)` confirms:

```@example thermo
parameters(m)
```

This is the Haldane relation [Haldane1930](@cite): for a reversible enzyme, the
ratio of forward to reverse catalytic rates must equal the overall equilibrium
constant scaled by the binding affinity ratio. Any mechanism where this relation
is violated cannot be in thermodynamic equilibrium.

## Wegscheider constraints in expanded mechanisms

Mechanisms with multiple substrate- or product-binding orders can contain
Wegscheider cycles in addition to Haldane cycles. The `# Wegscheider
constraints:` section of [`rate_equation_string`](@ref) lists each such tie.
Entries marked `(substituted into v)` have been absorbed directly into the rate
equation; only entries with a multi-parameter right-hand side appear as explicit
assignments.

See also: [`parameters`](@ref), [`rate_equation_string`](@ref).
