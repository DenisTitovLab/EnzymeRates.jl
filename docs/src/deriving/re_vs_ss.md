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

## A concrete comparison

A one-substrate reaction barely distinguishes the two assumptions — the
rapid-equilibrium and steady-state rate laws differ only in how three or four
constants are named. The difference becomes obvious for even the simplest
mechanisms with more than one substrate or product.
Take an ordered bi-uni reaction — `A` binds, then `B`, the complex isomerizes,
and `P` leaves — first as a rapid-equilibrium mechanism with steady-state
catalysis:

```@example revss
using EnzymeRates
re = @enzyme_mechanism begin
    substrates: A, B
    products:   P
    steps: begin
        E + A ⇌ E(A)
        E(A) + B ⇌ E(A, B)
        E(A, B) <--> E(P)
        E(P) ⇌ E + P
    end
end
parameters(re)
```

Six parameters, and a compact rate law:

```@example revss
print(rate_equation_string(re))
```

Now the same skeleton with every step made steady state:

```@example revss
ss = @enzyme_mechanism begin
    substrates: A, B
    products:   P
    steps: begin
        E + A <--> E(A)
        E(A) + B <--> E(A, B)
        E(A, B) <--> E(P)
        E(P) <--> E + P
    end
end
parameters(ss)
```

Nine parameters — each binding step trades its single `K` for an independent
`kon`/`koff` pair — and the rate law is far larger:

```@example revss
print(rate_equation_string(ss))
```

The contrast is structural, not cosmetic. The rapid-equilibrium denominator has
one term per reachable enzyme form — `1`, `A`, `A·B`, and `P`, four terms —
because every binding step factors out as a pre-equilibrium segment. The
steady-state denominator keeps those same four terms but adds new ones in `B`
and `B·P`: nothing factors out, so the King–Altman treatment carries a term for
every enzyme-form pattern, with products of on/off rates throughout. Those extra
terms give the steady-state equation qualitatively different behaviour — not the
same form with renamed constants — and after thermodynamic reduction it still
keeps more independent parameters, because each binding step starts with two
rate constants instead of one. Mechanisms with **random-order** binding diverge
even more drastically: their steady-state rate equations can carry *squared*
concentration terms, which users can confirm for themselves by building a
random-order mechanism and printing its rate equation.

## The RE assumption

Rapid equilibrium assumes the binding step relaxes to equilibrium on a time
scale much faster than the catalytic step.
When that assumption holds, the single `Kd` is sufficient.
When it does not hold, the full on/off pair is needed.

The [The Cha / King–Altman algorithm](@ref) solves the full Cha steady state
regardless of whether individual steps are RE or SS; RE steps simply factor
out of the rate matrix as pre-equilibrium segments, giving the familiar
`K_met_form` notation in the denominator.
