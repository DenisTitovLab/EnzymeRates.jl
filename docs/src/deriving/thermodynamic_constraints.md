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
k_EP_to_ES = K_EP_to_E_P * k_ES_to_EP / (K_ES_to_E_S * Keq)
```

A **Wegscheider constraint** comes from a cycle that closes on itself with zero
net metabolite change — one per independent internal loop. Its dependent rate
constant is a power product of the others **with no `Keq` factor**. These appear
when the same enzyme form can be reached by more than one binding order:
detailed balance forces the binding constants around the loop to multiply
consistently. For a random-order mechanism in which substrates `A` and `B` can
bind in either order, the four binding constants are tied:

```
K_EA_to_E_A * K_EAB_to_EA_B = K_EB_to_E_B * K_EAB_to_EB_A
```

When the other constants of a loop cancel, the loop fixes the remaining one at 1
and its line reads `K_E_to_F = 1`. An inhibitor that binds conformations `E` and
`F` into one complex `EI` (`E + I ⇌ EI` and `F + I ⇌ EI`) under one shared
dissociation constant does this: the loop `E → EI → F → E` forces `E` and `F`
to equal weight.

## How the constraints are found

The package takes the fundamental cycles of a spanning forest of the enzyme-form graph,
one per step outside the forest, so the cycle rows come out **linearly independent by
construction**. Each cycle is classified by its net metabolite change: proportional to
the overall reaction makes it a Haldane cycle, zero makes it a Wegscheider cycle. An
allosteric mechanism writes every cycle once per conformation, and the steps of a kinetic
group add rows that tie their constants together. One exact rational elimination over
all of these rows chooses the dependent constants (integer arithmetic that switches to
big integers if a value grows too large), so the same constraints come out every time.
If a mechanism is thermodynamically infeasible — its rows combine to
`0 = c · log(Keq)` with `c ≠ 0` — the package raises an error rather than emitting a
silent or wrong equation.

## The choice of independent vs dependent parameters

Each independent cycle costs one constant. The package decides which constant becomes
dependent by a fixed priority, read in turn; the earlier a constant ranks, the sooner it
becomes dependent:

1. **Conformation.** In an allosteric mechanism, the inactive conformation's own
   constants go before the active ones.
2. **Kind of step.** Isomerizations (steps that take up and give off nothing, on-enzyme
   chemistry included) go first, then steps that move more than one metabolite (a
   Theorell–Chance step), then bindings, which take up or give off one metabolite,
   whether or not they also convert.
3. **Occupancy.** Steps neither of whose enzyme forms is free E go before steps that
   touch free E.
4. **Speed.** Steady-state constants go before rapid-equilibrium ones.
5. **Direction.** Within a steady-state step, the rate constant that runs against the
   reaction goes first: the reverse rate, except for a step stored against the reaction,
   such as the release of a product, which is stored as its uptake. So in
   `ES <--> E + P`, stored as the binding `E + P → ES`, the product uptake is dependent
   and the catalytic rate stays fitted.

The free-enzyme binding constants a biochemist measures therefore stay independent
whenever a cycle leaves them free. Ties go to the names of each step's forms and
metabolites, never to the order in which steps were written, so two mechanisms that
differ only in step order, step direction or how the kinetic groups were labelled have
the same dependents. A kinetic group's constant takes the name of the step the solve
keeps, so splitting a group whose parts the cycles already tie changes no fitted name.
The order reads names, though: renaming a metabolite or an enzyme form can change which
of two tied constants is dependent.

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
        E + A ⇌ EA
        E + B ⇌ EB
        EA + B ⇌ EAB
        EB + A ⇌ EAB
        EAB <--> EP
        EP ⇌ E + P
    end
end
print(rate_equation_string(m))
```

The `# Wegscheider constraints:` section ties the binding constants of the two
orders together (no `Keq`), and the `# Haldane constraints:` section fixes the
reverse catalytic rate from `Keq` and the forward constants. Two constants are
therefore dependent and do not appear in the fitted-parameter list:

```@example thermo
parameters(m)
```

Both eliminated constants are absent from `parameters(m)`, so the fit explores
only the thermodynamically independent directions — every candidate equation the
fitter evaluates already satisfies detailed balance.

See also: [`parameters`](@ref), [`rate_equation_string`](@ref).
