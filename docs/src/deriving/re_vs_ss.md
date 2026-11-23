# Rapid equilibrium vs steady state

Every step in a mechanism is either **rapid equilibrium** (RE) or
**steady state** (SS).
That single flag drives parameter count: RE steps cost one parameter each,
while SS steps cost two.

## The flag in the DSL


In [`@enzyme_mechanism`](@ref), `⇌` marks a step as rapid equilibrium and
`<-->` marks it as steady state.
The mechanism stores it as each step's `rapid_equilibrium` flag.

## Parameter count per step


- An **RE binding step** contributes one parameter: a dissociation constant
  `Kd`, named in the release direction — the bound side to the free side —
  as `K_<bound>_to_<free>` (for example, `K_ES_to_E_S`). A binding is any step
  that takes up exactly one metabolite and gives off none, so a fused binding,
  whose bound form is not its free form with the metabolite added, is named the
  same way: `EA + B ⇌ EPQ` gives `K_EPQ_to_EA_B`. A step that gives off
  one metabolite and takes up none is the binding it reverses: the fused
  release `EAB ⇌ EQ + P` binds `P` into `EAB`, with dissociation
  constant `K_EAB_to_EQ_P` = [EQ]·[P]/[EAB].
- An **RE isomerization step** contributes one parameter: an equilibrium
  constant, named in its stored direction as `K_<from>_to_<to>` (for
  example, `K_ES_to_EP`).
- An **SS binding step** contributes two rate constants, one per direction:
  `k_<free>_to_<bound>` and `k_<bound>_to_<free>` (for example,
  `k_E_S_to_ES` and `k_ES_to_E_S`).
- An **SS isomerization step** contributes two directed rate constants:
  `k_<from>_to_<to>` and `k_<to>_to_<from>`.
- Every other step — a Theorell–Chance step (`EA + B <--> EQ + P`), or one
  with several metabolites on a side — is named by its two sides, as an
  isomerization is. A side is a step's enzyme form followed by its free
  metabolites. The mechanism, not the direction the step is written in, decides
  which side is `<from>`. A step that only takes up or only gives off metabolites is
  stored as the uptake, whatever their number: the release `EPQ ⇌ E + P + Q` is stored
  as `E + P + Q ⇌ EPQ` and gives `K_E_P_Q_to_EPQ`, named in the uptake direction, unlike a
  binding's dissociation constant. A step with metabolites on both sides is stored the
  way the reaction runs: the way it takes up more substrates and gives off more products,
  net. When that count does not decide, as for an exchange of one modifier for another,
  and for an isomerization, the step runs the way the shortest cycle through it that
  turns the reaction over runs; with no such cycle it runs away from free `E`, and with
  both ends equally far, from the form whose name comes first. So writing the step
  backwards gives the same names: at rapid equilibrium the Theorell–Chance step above
  gives `K_EA_B_to_EQ_P`. Such a constant can carry concentration units, as
  `K_E_P_Q_to_EPQ` does.

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
        E + A ⇌ EA
        EA + B ⇌ EAB
        EAB <--> EP
        EP ⇌ E + P
    end
end
parameters(re)
```

Four parameters, and a compact rate law:

```@example revss
print(rate_equation_string(re))
```

Now the same skeleton with every step made steady state:

```@example revss
ss = @enzyme_mechanism begin
    substrates: A, B
    products:   P
    steps: begin
        E + A <--> EA
        EA + B <--> EAB
        EAB <--> EP
        EP <--> E + P
    end
end
parameters(ss)
```

Seven parameters — each binding step trades its single `K` for an independent
forward/reverse rate-constant pair — and the rate law is far larger:

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
`K_ES_to_E_S`-style notation in the denominator.

## Apparent constants on a steady-state chain

An RE label does not always claim that a binding is fast. Take a *chain*: an
isomerization alone in its kinetic group, between two forms that each have
exactly one other step, a binding into that form alone in its group. Those two
bindings are the chain's *flanks*. The uni-uni mechanism `E + S ⇌ ES`, `ES <--> EP`,
`EP ⇌ E + P` is one chain. When the isomerization is steady state, the rest of
the mechanism sees the chain only through its flux and the enzyme it holds, and
rapid-equilibrium flanks already give every rate law that steady-state flanks
give. An RE flank's constant is then an apparent constant, a Michaelis-type
combination of the chain's rate constants, not the dissociation constant of its
complex: here `K_ES_to_E_S` plays the part of the Michaelis constant for `S`. The
mechanism gives the same rate laws as the chain with all three steps steady
state, with three parameters instead of five. The enumeration therefore never
flips such a flank to steady state (see [The enumeration engine](@ref)).
