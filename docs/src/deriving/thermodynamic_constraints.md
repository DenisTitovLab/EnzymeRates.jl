# Thermodynamic constraints

Every reversible mechanism contains cycles: sequences of steps that return the
enzyme to its starting form. Thermodynamics constrains the rate constants around
each cycle — their product must equal a fixed function of the overall equilibrium
constant. The package identifies these constraints automatically and removes one
rate constant per independent cycle, reducing the number of fitted parameters and
guaranteeing that every fitted equation is thermodynamically consistent.

## Haldane and Wegscheider relations

Two kinds of constraint arise [Haldane1930](@cite) [Wegscheider1901](@cite), one
for each kind of cycle.

A **Haldane constraint** comes from a cycle whose net change matches the overall
reaction — there is one per independent cycle that turns substrate into product.
Its dependent rate constant is a power product that **carries `Keq`**, so the
forward and reverse rates around the cycle cannot be chosen independently of the
equilibrium constant. For the reversible Michaelis–Menten mechanism the reverse
catalytic rate is fixed by `Keq`, the two binding constants, and the forward
rate:

```
k_EP_to_ES = (1 / Keq) * K_P_E * (1 / K_S_E) * k_ES_to_EP
```

A **Wegscheider constraint** comes from a cycle that closes on itself with zero
net metabolite change — one per independent internal loop. Its dependent rate
constant is a power product of the others **with no `Keq` factor**. These appear
when the same enzyme form can be reached by more than one binding order:
detailed balance forces the binding constants around the loop to multiply
consistently. For a random-order mechanism in which substrates `A` and `B` can
bind in either order, the four binding constants are tied:

```
K_A_E * K_B_EA = K_B_E * K_A_EB
```

## How the constraints are found

The package builds the enzyme-form **incidence matrix** — one column per step,
with `+1` at the step's to-form and `-1` at its from-form — and computes its
exact-integer null space. Each null-space basis vector is one independent cycle,
so the constraints come out **linearly independent by construction**; there is
no redundant set to thin down. The whole computation uses exact rational
arithmetic, never floating point, so the same cycles are reproduced exactly every
time.

Each cycle is then classified by its net metabolite change: proportional to the
overall reaction makes it a Haldane cycle, zero makes it a Wegscheider cycle. If
a mechanism is thermodynamically infeasible — a cycle whose constraint reduces to
`0 = c · log(Keq)` with `c ≠ 0` — the package raises an error rather than
emitting a silent or wrong equation.

## The choice of independent vs dependent parameters

Each cycle costs one rate constant, and the package chooses which constant to
keep and which to express in terms of the others by a fixed priority, designed to
keep the **biochemically meaningful** parameters independent:

1. **Free-enzyme binding constants** — the affinity of a substrate or product for
   the free enzyme, such as `K_S_E` — are kept independent whenever possible.
   These are the quantities an experimentalist measures and reports.
2. Binding steps on already-occupied enzyme forms are eliminated next, when a
   cycle needs them.
3. **Internal isomerization rates** are made dependent first: they are the most
   eliminable and the least directly measurable.

Within a steady-state step, the reverse rate is preferred dependent over the
forward rate. Because step and group order are canonicalized when a mechanism is
constructed, this choice is deterministic: two mechanisms with the same structure
written in a different step order produce the identical reduced equation.

## A concrete example

A random-order mechanism — two substrates that can bind in either order — carries
**both** kinds of constraint at once: a Wegscheider tie among its binding
constants and a Haldane tie on its catalytic step.

```@example thermo
using EnzymeRates
m = @enzyme_mechanism begin
    substrates: A, B
    products:   P
    steps: begin
        E + A ⇌ E(A)
        E + B ⇌ E(B)
        E(A) + B ⇌ E(A, B)
        E(B) + A ⇌ E(A, B)
        E(A, B) <--> E(P)
        E(P) ⇌ E + P
    end
end
print(rate_equation_string(m))
```

The `# Wegscheider constraints:` section ties the binding constants of the two
orders together (no `Keq`), and the `# Haldane constraints:` section fixes the
reverse catalytic rate from `Keq` and the forward constants. Two rate constants
are therefore dependent and do not appear in the fitted-parameter list:

```@example thermo
parameters(m)
```

Both eliminated constants are absent from `parameters(m)`, so the fit explores
only the thermodynamically independent directions — every candidate equation the
fitter evaluates already satisfies detailed balance.

See also: [`parameters`](@ref), [`rate_equation_string`](@ref).
