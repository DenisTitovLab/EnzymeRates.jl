# Rapid equilibrium vs steady state

Every step in a mechanism is either **rapid equilibrium** (RE) or
**steady state** (SS).
That single flag drives parameter count: RE steps cost one parameter each,
while SS steps cost two.

## The flag in the DSL


In [`@enzyme_mechanism`](@ref), `⇌` marks a step as rapid equilibrium and
`<-->` marks it as steady state.
The distinction is stored on each `Step` as the `is_equilibrium` field.

## Parameter count per step


- An **RE binding step** contributes one parameter: a dissociation constant
  `Kd`, rendered as `K_<metabolite>_<form>` (for example, `K_S_E`).
- An **RE isomerization step** contributes one parameter: an isomerization
  constant `Kiso`, rendered as `Kiso_<from>_to_<to>` (for example,
  `Kiso_ES_to_EP`).
- An **SS binding step** contributes two rate constants: `kon_<met>_<form>`
  and `koff_<met>_<form>`.
- An **SS isomerization step** contributes two directed rate constants:
  `k_<from>_to_<to>` and `k_<to>_to_<from>`.

All parameter names flow through the `name(p, m)` chokepoint, which uses
metabolite and form names rather than arbitrary indices.

## A concrete comparison

The textbook reversible Michaelis–Menten mechanism has RE binding steps and
one SS catalytic step:

```@example revss
using EnzymeRates
re = @enzyme_mechanism begin
    substrates: S
    products:   P
    steps: begin
        E + S ⇌ E(S)
        E(S) <--> E(P)
        E(P) ⇌ E + P
    end
end
parameters(re)
```

The same skeleton with all steps made steady state:

```@example revss
ss = @enzyme_mechanism begin
    substrates: S
    products:   P
    steps: begin
        E + S <--> E(S)
        E(S) <--> E(P)
        E(P) <--> E + P
    end
end
parameters(ss)
```

The RE form fits five parameters: `(:K_P_E, :K_S_E, :k_ES_to_EP, :Keq,
:E_total)`.
The SS form fits seven: `(:k_ES_to_EP, :koff_P_E, :koff_S_E, :kon_P_E,
:kon_S_E, :Keq, :E_total)`.
Each binding step that switches from RE to SS gains a separate on-rate and
off-rate in place of the single equilibrium constant.

After thermodynamic reduction, one rate constant per independent cycle is
eliminated by the Haldane/Wegscheider constraint — but SS steps still add
more *independent* parameters before that reduction, because they start with
two rate constants each.

## The RE assumption

Rapid equilibrium assumes the binding step relaxes to equilibrium on a time
scale much faster than the catalytic step.
When that assumption holds, the single `Kd` is sufficient.
When it does not hold, the full on/off pair is needed.

The [The Cha / King–Altman algorithm](@ref) solves the full Cha steady state
regardless of whether individual steps are RE or SS; RE steps simply factor
out of the rate matrix as pre-equilibrium segments, giving the familiar
`K_met_form` notation in the denominator.

## The RE→SS expansion move

When the package searches for the best rate equation, it can promote a whole
kinetic group from RE to SS atomically — every step sharing that group
converts together.
This expansion is one of the moves in the enumeration engine.
See the enumeration engine page for details.
