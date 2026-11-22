# Mechanisms with allosteric regulators

Allosteric regulation is one of the main ways metabolic enzymes are controlled.
An **allosteric regulator** binds a regulatory site that is distinct from the
catalytic site and, from there, either inhibits or activates the enzyme. A
classical description of how this works is the Monod–Wyman–Changeux (MWC) model
[Monod1965](@cite), and EnzymeRates can derive rate equations under that model.

In the MWC picture the enzyme exists in two conformations — an **active** state
(A) and an **inactive** state (I) — in equilibrium with one another. (In the
classic MWC notation A is the R state and I is the T state; the package uses A/I
throughout, so `K_A_…` names correspond to the literature's R-state constants.) A
coupling constant `L` sets the inactive-to-active population ratio of the free
enzyme. What distinguishes an activator from an inhibitor is simply which
conformation it prefers: an **activator binds the active state more tightly** and
shifts the population toward A, while an **inhibitor binds the inactive state**
and shifts it toward I. That preference is expressed entirely through the
regulator's binding constants for the two states — its `K_A` versus its `K_I`.

## Setting up an MWC mechanism

Write it with [`@enzyme_mechanism`](@ref), ending each catalytic statement and each
regulator with its allosteric state. The printed rate equation shows the result:

```@example mwc
using EnzymeRates
allo = @enzyme_mechanism begin
    substrates: S
    products:   P
    oligomeric_state: 2
    steps: begin
        E + S ⇌ ES     :: EqualAI
        ES <--> EP     :: NonequalAI
        EP ⇌ E + P     :: EqualAI
    end
    regulatory_site(2): A::OnlyA
    regulatory_site(2): I::OnlyI
end
print(rate_equation_string(allo))
```

The independent parameters the mechanism fits:

```@example mwc
parameters(allo)
```

Each catalytic kinetic group and each regulator carries one of four allosteric states,
which controls how its A- and I-state symbols are rendered:

| State | Meaning | Example symbol |
|-----|---------|----------------|
| `:EqualAI` | One shared symbol for both states | `K_ES_to_E_S` |
| `:NonequalAI` | Independent A and I symbols | `k_A_ES_to_EP`, `k_I_ES_to_EP` |
| `:OnlyA` | Present in the active state only; zeroed in the inactive polynomial | `K_A_Areg` |
| `:OnlyI` | Present in the inactive state only; zeroed in the active polynomial | `K_I_Ireg` |

In the example the two binding steps are `:EqualAI`, so they share `K_ES_to_E_S`
and `K_EP_to_E_P` with no A/I token; the isomerization step is `:NonequalAI`, giving
independent `k_A_ES_to_EP` and `k_I_ES_to_EP`; the activator `A` is `:OnlyA`
(`K_A_Areg`, active state only) and the inhibitor `I` is `:OnlyI` (`K_I_Ireg`,
inactive state only). An `:EqualAI` regulator `R` would take one bare constant, `K_Rreg`,
in both states. `L` is the free-enzyme coupling constant.
Any catalytic group may be `:OnlyI`, such as an inhibitor's binding that only the inactive
conformation has, or one of two reaction paths, as long as the active conformation still
turns substrate into product.

The mechanism also declares its **oligomeric state** and the **multiplicity** of each
regulatory site. `oligomeric_state` is the number of identical catalytic subunits, which
change conformation together; it becomes the exponent on each conformation's catalytic
partition function. Each `regulatory_site(m)` says how many copies of that
site the oligomer carries, raising the site's binding factor to the power `m`. In
the example both are `2` — a homodimer with two copies of each regulatory site.
A regulator is named only at its site: it binds the catalytic site too only where a step
takes it up, and is then declared under `modifiers:` as well.

!!! note
    The rate equation gives the rate **per catalytic site**: `oligomeric_state`
    enters only through the partition-function exponents, never as an overall
    multiplier on the numerator.

## Derivation of the MWC rate equation

The MWC rate equation is assembled from two ordinary single-conformation rate
equations — one for the active state, one for the inactive state — combined by a
partition function. The derivation rests on two rapid-equilibrium assumptions.
Allosteric regulators bind fast relative to catalysis. And **the enzyme changes
conformation only while free**: it settles into the active or inactive state
before it binds anything — the free-enzyme inactive-to-active ratio is `L` — and
keeps that conformation for a full catalytic cycle. The overall rate is then each
conformation's catalytic cycle weighted by its free-enzyme population.
The conformations interconvert only at free E, the form `E` with nothing bound, and `L` is
its inactive-to-active ratio. A mechanism may hold other free forms, such as a
conformation `F` with `F <--> E`; they change conformation only by way of E. When another
free form sits at rapid equilibrium with E, its own inactive-to-active ratio is `L·r_I/r_A`,
where `r_c` is the equilibrium ratio `[F]/[E]` in conformation c across the step that joins
F to E, so a fit reports `L` at E whatever the free forms are named.

An equivalent statement lets every form interconvert between the A and I states,
but with a thermodynamically fixed ratio: `L` for the free enzyme, and `L` scaled
by the form's active-to-inactive binding-constant ratio for a bound form (so `L`
itself for an `:EqualAI` ligand). The two statements agree except when catalysis
differs between the conformations (`:NonequalAI`); there, committing to one
conformation per cycle — the model here — routes each molecule's turnover through
a single conformation rather than blending the two, a difference of a few percent.

Treating regulator binding as steady state instead would not help — a regulator
is not consumed, so a non-RE binding step contributes only the ratio of its on-
and off-rate constants, which is exactly the equilibrium constant the RE
treatment already uses, and the two rates would not be separately identifiable
from rate data.

For a mechanism with oligomeric state `n`, the package derives each conformation's
law on the graph it keeps (the active conformation without its `:OnlyI` groups, the
inactive one without its `:OnlyA` groups, each without the forms that dropping those groups
cuts off from free E) and combines the two into the rate per catalytic site:

```
v = (H) * (F_I^n * Q_A * D_A^(n-1) * W_A + L * F_A^n * Q_I * D_I^(n-1) * W_I) /
          (F_I^n * D_A^n * W_A + L * F_A^n * D_I^n * W_I)
```

`H = Π[S] − Π[P]/Keq` is the Haldane factor both conformations share. `Q_A` is the active
conformation's numerator over `H` and `D_A` its denominator, the King–Altman/Cha
polynomials of its catalytic graph; `F_A` is free E's weight in `D_A`, the spanning-tree
sum of free E's segment times the concentration monomial that brings the law to lowest
terms. `Q_I`, `D_I` and `F_I` are the inactive conformation's, with `:EqualAI` groups
keeping the shared constants and `:NonequalAI` groups their own `I` names. `W_A` and `W_I`
multiply the regulatory-site factors `(1 + X/K)^m`, and `L` weighs the inactive branch.
Each branch carries the other conformation's free-E weight, so both stand on free E, where
the conformations interconvert; that keeps a bare rate constant out of the `L` term.

The printed law lists these polynomials by name. Each conformation's `Q`, `D` and `F` are
divided by that conformation's own normalizer, `N_A` or `N_I` ([Normalization](@ref)),
which the combination does not see. Every term of `D_A` and `D_I` is then dimensionless,
except in a conformation with no normalizer, which prints its polynomials as derived.
`F_A` and `F_I` print only when they differ. They are 1 whenever each conformation's
denominator has a constant term and the graph it keeps is one rapid-equilibrium segment,
with E its only free form and no covalent form. That is the common case, and the law is
then the textbook MWC form. When the `:OnlyA` groups stop every
catalytic cycle of the inactive conformation, as the search's allosteric bases do, the law
has no `Q_I`: the inactive branch adds enzyme to the denominator and no flux. When the
inactive conformation is free E alone, the
law has no `D_I` either, and its branch of the denominator is `L` times `F_A^n` and `W_I`,
each left out when it is 1. With `coefficients = true` each polynomial reads in named
coefficients (`dA_S` for the coefficient of `S` in `D_A`, `nI_0` for the constant term of
`Q_I`).

In the printed equation above `n = 2`, so each `D` is squared in the denominator, the
activator's site gives `W_A = (1 + A / K_A_Areg)^2`, and the inhibitor's
`W_I = (1 + I / K_I_Ireg)^2`.

The term cap counts both conformations: a law whose two conformations hold more than 5,000
spanning-tree products together is refused, and identify's `eq_complexity_filter` compares
the same sum.

## Allosteric states as limits, and admission

An `:OnlyA` state excludes a group's steps from the inactive conformation, and `:OnlyI`
excludes them from the active one. The derivation reads each conformation's law by deleting
the steps it excludes and the forms that the deletion cuts off from free E. **Admission**
decides whether that law is the limit of ordinary mass-action kinetics as the excluded steps
vanish, each at its own rate. When it cannot show that it is, it rejects the states with an
error that says why and names the steps at fault where it has them.

### How an excluded step is switched off

The package reads an excluded step as nearly absent, not absent: in the conformation that
excludes it, the step has a very high dissociation constant or very low rate constants, and
every thermodynamic constraint, each Haldane and Wegscheider relation, still holds in each
conformation.

A step keeps its kind in both conformations: a rapid-equilibrium step stays at rapid
equilibrium where it is excluded, and a steady-state step stays at steady state. Rapid
equilibrium and steady state differ more than 100-fold in speed, and turning a fast step
into the slowest step of the other conformation is not a reasonable approximation. A
mechanism that needs it is written with that step steady-state.

The picture is one of energies. In each conformation the steps it keeps join the forms into
**islands**: sets of forms joined to each other by those steps. An island is neither a
kinetic group nor a rapid-equilibrium segment. Free E's island holds the forms of that
conformation's law. Every other island is cut off: it lies higher in energy than free E's
island and holds no enzyme. An excluded step is switched off according to where it lies:

| Excluded step | How it is switched off |
|:--|:--|
| Between two islands | It loses its rate uphill, into the higher island, and keeps its rate downhill. A binding whose bound side lies higher loses its affinity: its dissociation constant grows without bound. A rapid-equilibrium step stays fast. |
| Between two islands, steady state | It may stop instead: both rate constants fall to zero. |
| Inside one island, steady state | It stops: both rate constants fall to zero, their ratio fixed by thermodynamics. |
| Inside one island, rapid equilibrium | Error. Its two sides lie at one energy, so thermodynamics fixes its equilibrium constant, and it could be switched off only by turning slow. |

The mechanism fixes the islands but not which lies above which. The package looks for a
choice of energies, and of the steady-state steps that stop, under which it can show that
the law is the limit, and accepts the mechanism when it finds one. When the excluded steps
allow too many choices to search, it tries only three, and rejects the mechanism when it
can show the limit under none of them. The choice does not change the law, which leaves the
excluded steps out. The package takes the law to be the limit under a choice when:

1. **Existence.** The steps that lose a rate, each drawn uphill, close no cycle and never
   point into free E's island, so the energies exist and no island lies below free E's.
2. **Uniqueness.** They lead from free E's island to every cut-off island, so each cut-off
   island empties. On every face where data sit (below), the steady state does not depend
   on where the enzyme starts, no cut-off form fills faster than it empties, and no
   excluded step opens a faster route between two forms of the law, unless the cut-off
   forms on that route cannot change the limit.
3. **Turnover.** With `:OnlyI` groups, the active conformation still turns the reaction
   over.

Four examples follow, each excluding steps from the inactive conformation. In the first,
the substrate binding and the catalytic step are excluded:

```@example mwc
both_excluded = @enzyme_mechanism begin
    substrates: S
    products:   P
    steps: begin
        E + S ⇌ ES     :: OnlyA
        ES <--> EP     :: OnlyA
        EP ⇌ E + P     :: EqualAI
    end
end
print(rate_equation_string(both_excluded))
```

In the inactive conformation only the product binding runs, so E and EP form free E's
island and ES is cut off. The substrate binding loses its affinity, and the catalytic step
stops or loses its rate into ES. The inactive conformation binds P and turns nothing
over: the law has `D_I` and no `Q_I`.

With the substrate binding excluded alone, the package rejects the mechanism:

```jldoctest; filter = r"\s+" => " "
julia> using EnzymeRates

julia> @enzyme_mechanism begin
           substrates: S
           products:   P
           steps: begin
               E + S ⇌ ES     :: OnlyA
               ES <--> EP     :: EqualAI
               EP ⇌ E + P     :: EqualAI
           end
       end
ERROR: Mechanism: in the inactive conformation the rapid-equilibrium step E + S ⇌ ES is
tagged :OnlyA, but the steps that conformation keeps still join its two sides.
Thermodynamics then fixes its equilibrium constant, and a step keeps its kind in both
conformations, so it cannot be switched off. Make its kinetic group steady-state, or tag
:OnlyA the steps that join its sides: E + P ⇌ EP; ES <--> EP.
```

The catalytic step and the product binding still join ES to free E, so ES stays in free E's
island, made from P, and the substrate binding lies inside one island. `Keq` and the other
steps' constants fix its dissociation constant, which cannot grow.

The same binding at steady state is accepted:

```@example mwc
ss_excluded = @enzyme_mechanism begin
    substrates: S
    products:   P
    steps: begin
        E + S <--> ES  :: OnlyA
        ES <--> EP     :: EqualAI
        EP ⇌ E + P     :: EqualAI
    end
end;
```

The substrate binding again lies inside free E's island, and at steady state it stops. The
inactive conformation keeps ES, made from P, and turns nothing over.
`identify_rate_equation` never proposes such a mechanism, in which a conformation keeps a
form that only excluded bindings enter, but you may write one.

The last example is a ping-pong whose inactive conformation runs neither half-reaction.
Write F for the modified enzyme `E_res_plusA_minusP`, and FP and FB for its forms with P
and B bound, `EP_res_plusA_minusP` and `EB_res_plusA_minusP`:

```@example mwc
pingpong = @enzyme_mechanism begin
    substrates: A, B
    products:   P, Q
    steps: begin
        E + A ⇌ EA                                     :: EqualAI
        EA <--> EP_res_plusA_minusP                    :: OnlyA
        EP_res_plusA_minusP ⇌ E_res_plusA_minusP + P   :: EqualAI
        E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP   :: EqualAI
        EB_res_plusA_minusP ⇌ EQ                       :: OnlyA
        EQ ⇌ E + Q                                     :: EqualAI
    end
end;
```

In the inactive conformation E, EA and EQ form free E's island, and F, FP and FB form
one cut-off island. Each half-reaction runs between the two islands: the first loses its
rate from EA into FP, the second its rate from EQ into FB, and each keeps its rate back.
The modified enzyme then empties through both halves, whatever the ratio of the two
vanishing rates, and the inactive conformation binds A and Q and turns nothing over.

Two of the face checks, that no cut-off form fills faster than it empties and that the
cut-off forms on a faster route cannot change the limit, rest on numerical tests, not on a
proof. The package was run on 78,554 mechanisms. Among those it accepts only through its
search for energies or through cut-off forms that cannot change the limit, two numerical
oracles agree with it on every mechanism they checked; [Well-posedness guard](@ref), on the
Developer page, gives the evidence. Admission is not exact: it rejects some mechanisms whose
law is the limit, which the repository's `ROADMAP.md` lists.

When the steps a conformation keeps join the two sides of an excluded rapid-equilibrium
step, as in the second example, the error says so and offers each fix that the package
accepts: make the step's kinetic group steady-state, so that it stops (the third example),
or exclude the steps that keep its two sides joined, here the product binding and the
catalytic step, so that one side is cut off. In the active conformation that fix makes
them `:OnlyI`.

### What the package rejects

A mechanism is well defined when its rate law is the unique fast-equilibrium limit of
mass-action kinetics wherever data sit: at interior concentrations and on every **face**,
where every product is zero, every substrate, each reactant alone on a side of two or more,
or any of these with every modifier-only name zero. `@enzyme_mechanism` throws an error
that states the fault and names the steps at fault and the face where it has them, and the
mechanism search drops the mechanism. The package rejects:

1. an unbalanced cycle: a cycle of steps that takes up or gives off a metabolite in another
   proportion than the reaction, as when a step leaves out a metabolite;
2. no catalytic cycle: no cycle of steps turns substrates into products, so the rate is
   zero;
3. a cycle of rapid-equilibrium steps that turns the reaction over, so the rate is
   infinite;
4. a face where the law is not the unique limit: a rapid-equilibrium segment with no form
   whose weight keeps the segment's lowest power of each of the face's metabolites, a
   segment whose lowest forms are joined only through other forms, a trap that the face
   feeds and its fast exits cannot drain, a second closed set of forms, or a law with no
   value (0/0) there;
5. allosteric states for which the package finds no choice of energies under which it can
   show that the law is the limit. When the states allow too many choices to search, it
   tries only three (admission, above).

The faces leave out the corners that zero one substrate with every product. A law the
package accepts can be 0/0 there, where a well-defined law is 0, so the compiled rate
equation returns `NaN` at such a point; a fit scores the NaN as it scores the 0 a
well-defined law predicts there.

## Thermodynamic constraints of MWC equations

A single conformation already carries thermodynamic constraints: every cycle of
steps must agree with the overall equilibrium constant `Keq` (a **Haldane**
relation), and every cycle of binding steps must be path-independent (a
**Wegscheider** relation). The derivation reports both and expresses the
dependent rate constants through the independent ones. An MWC equation solves the
relations of both conformations together, adding constraints neither has alone. Here
we review several examples of such additional constraints.

Make a lone substrate binding `:NonequalAI` while catalysis and product release
stay `:EqualAI`, and the split collapses:

```@example mwc
collapse = @enzyme_mechanism begin
    substrates: S
    products:   P
    oligomeric_state: 2
    steps: begin
        E + S ⇌ ES     :: NonequalAI
        ES <--> EP     :: EqualAI
        EP ⇌ E + P     :: EqualAI
    end
    regulatory_site(2): R::OnlyI
end
print(rate_equation_string(collapse))
```

The equation reports `K_I_ES_to_E_S = K_A_ES_to_E_S`. The Haldane relation fixes `Keq` from
the catalytic rate constants and the two affinities in each conformation; with
catalysis, product release, and `Keq` all shared, the substrate affinity is
pinned as well. The state asked for a difference thermodynamics forbids, so the
substrate stops distinguishing the conformations.

Make **both** the substrate and product bindings `:NonequalAI`, and one degree of
freedom survives:

```@example mwc
coupled = @enzyme_mechanism begin
    substrates: S
    products:   P
    oligomeric_state: 2
    steps: begin
        E + S ⇌ ES     :: NonequalAI
        ES <--> EP     :: EqualAI
        EP ⇌ E + P     :: NonequalAI
    end
    regulatory_site(2): R::OnlyI
end
print(rate_equation_string(coupled))
```

Now `K_I_ES_to_E_S = K_A_ES_to_E_S · K_I_EP_to_E_P / K_A_EP_to_E_P`, i.e.
`K_S^I / K_S^A = K_P^I / K_P^A`:
the affinities differ freely as long as they differ *together*, holding their
ratio fixed. This is the thermodynamically consistent form of a **K-system**, and
it needs two coupled `:NonequalAI` bindings, not one.

Often a K-system uses exclusive `:OnlyA` binding: the substrate binds the active
state alone, and the inactive cycle is pruned (no `Q_I`). When that binding is at rapid
equilibrium, its catalytic step **must** be `:OnlyA` too. This is a thermodynamic
requirement, not a convention: with `:EqualAI` or `:NonequalAI` catalysis and an `:EqualAI`
product binding, the inactive conformation still joins E to ES through catalysis and
product release, so thermodynamics fixes the substrate affinity there, and the package
rejects the combination ([How an excluded step is switched off](@ref)). A steady-state
substrate binding is stopped instead: the inactive conformation keeps ES, made from the
product, and turns nothing over.

The exception is a balanced K-system. When `:OnlyA` binds both a substrate and a
product, ES and EP form one cut-off island in the inactive conformation, which rises as
a whole: the two `K_I → ∞` limits cancel in the Haldane relation, the affinities diverge
together and their ratio stays free, so `:EqualAI` catalysis is legal.

A steady-state binding splits the constraint further, carrying an affinity (a
ratio of its two rate constants) the cycles constrain and a speed (their
product) they do not: a forbidden affinity collapses while the speed stays
free, so the two conformations bind with the same `Kd` but different kinetics.
Making the substrate binding steady-state and `:NonequalAI` derives its
reverse rate (`k_I_ES_to_E_S = k_A_ES_to_E_S · k_I_E_S_to_ES / k_A_E_S_to_ES`)
and leaves the forward rate a free parameter:

```@example mwc
ss = @enzyme_mechanism begin
    substrates: S
    products:   P
    oligomeric_state: 2
    steps: begin
        E + S <--> ES  :: NonequalAI
        ES <--> EP     :: EqualAI
        EP ⇌ E + P     :: EqualAI
    end
    regulatory_site(2): R::OnlyI
end
parameters(ss)
```

## See also

- [`@enzyme_mechanism`](@ref) — DSL reference
- [`parameters`](@ref) — list independent parameter names
- [`rate_equation_string`](@ref) — inspect the derived equation
