# ABOUTME: The DSL macros @enzyme_reaction and @enzyme_mechanism, which build a reaction
# ABOUTME: or a mechanism, plain or allosteric, and check a mechanism with the guard.

"""
    @enzyme_reaction begin
        substrates: S[C6H12O6], ...
        products: P[C6H13O9P], ...
        competitive_inhibitors: I, ...
        allosteric_regulators: A, ...
        oligomeric_state: 2
        shared_catalytic_site: (S, P), ...
    end

Build the [`EnzymeReaction`](@ref EnzymeRates.Mechanisms.EnzymeReaction) that the
mechanism search starts from: its reactants with their atoms, and the regulators and
oligomeric states the search may use.

- `substrates:` and `products:` (required) list each reactant with its atoms, one
  chemical formula in brackets: `S[C6H12O6]`. The atoms of the two sides balance, and a
  name is a substrate or a product, not both.
- `competitive_inhibitors:` lists the competitive inhibitors, names that bind the
  catalytic site without reacting. A substrate or product listed here may also bind
  where it does not react, as an inhibitor.
- `allosteric_regulators:` lists the allosteric regulators, names that bind regulatory
  sites. A type tag, `A::Activator` or `A::Inhibitor`, pins the regulator to the
  allosteric state that activates or inhibits; an untagged one (`:unspecified`) may take
  either. Only these entries take a type tag. A name may be both a competitive
  inhibitor and an allosteric regulator.
- `oligomeric_state:` gives the oligomeric states the search may assign: one
  positive Int, or a tuple of distinct ones such as `(1, 2, 4)`. It defaults to 1 and is
  given at most once.
- `shared_catalytic_site:` lists `(substrate, product)` pairs, in either order, that bind
  one catalytic site, so never both at once.
- The list labels may repeat; their entries accumulate.

No metabolite is named `Keq`, `L`, `N`, `params`, `concs` or `v`, or holds `_`: a
rate equation keeps those names, and its constant names join names with `_`. Nor
is one named `E`, the name of the free enzyme form; a name is an identifier that is no
keyword. The mechanism search names an enzyme form by joining its bound metabolites'
names without a separator, so the names may not spell one string two ways: with `Ac`,
`CoA` and `AcCoA`, E with Ac and CoA bound and E with AcCoA bound would both be
`EAcCoA`. The search rejects such a reaction when it starts and names both spellings
(`Ac + CoA and AcCoA both spell AcCoA`); rename one (`Acetyl`). A reactant bound where it
does not react is written `Xmod` in a form's name, which counts as a name too.

A reaction prints as the block that builds it, its names sorted:

```jldoctest
julia> using EnzymeRates

julia> @enzyme_reaction begin
           substrates: S[C6H12O6], ATP[C10H16N5O13P3]
           products: P[C6H13O9P], ADP[C10H15N5O10P2]
           competitive_inhibitors: P
           allosteric_regulators: A::Activator, I::Inhibitor
           oligomeric_state: (1, 2)
       end
@enzyme_reaction begin
    substrates: ATP[C10H16N5O13P3], S[C6H12O6]
    products: ADP[C10H15N5O10P2], P[C6H13O9P]
    competitive_inhibitors: P
    allosteric_regulators: A::Activator, I::Inhibitor
    oligomeric_state: (1, 2)
end
```
"""
macro enzyme_reaction(block)
    return reaction_expr(block)
end

"""
    @enzyme_mechanism begin
        substrates: S, ...
        products: P, ...
        modifiers: I, ...
        oligomeric_state: 2              # an allosteric mechanism only
        steps: begin
            E + S ⇌ ES                   # :: EqualAI on every statement, if allosteric
            ...
        end
        regulatory_site(2): X::OnlyA     # an allosteric mechanism only
    end

Build a mechanism from its elementary steps. The value is a [`Mechanism`](@ref
EnzymeRates.Mechanisms.Mechanism), which
[`compile_rate_equation`](@ref EnzymeRates.Derivation.compile_rate_equation) turns into a
rate law that evaluates without allocating. Two writings of one mechanism give equal
values. It prints as the block that builds it, each list sorted and one line per kinetic
group in canonical order.

Labels:
- `substrates:` and `products:` (required) name the reactants, without atoms; atoms
  belong in [`@enzyme_reaction`](@ref).
- `modifiers:` names the metabolites that bind the catalytic site without reacting, such
  as inhibitors and activators. A name is declared once: a substrate or product that
  also binds where it does not react is not listed here, since a step that binds it
  simply takes it up (`ES + S ⇌ ESS`).
- `steps:` (required) holds the steps in one `begin ... end` block.
- `substrates:`, `products:` and `modifiers:` may repeat; their names accumulate.
- `oligomeric_state:` and `regulatory_site(n):` belong to an allosteric mechanism (below).

Steps, one per line:
- `lhs ⇌ rhs` is a rapid-equilibrium step and `lhs <--> rhs` a steady-state one; either
  direction may be written, but the isomerizations of one kinetic group share their
  constants in the direction each is written, so write them the same way round.
- Each side holds one enzyme form and the metabolites the step takes up (left) or gives
  off (right): `E + S ⇌ ES`, `EA + B <--> EQ + P`. A name on a side that no label
  declares is an enzyme form, and each side holds exactly one.
- An enzyme form's name is any identifier that names no metabolite, such as `ES`, `EAB`
  or `Fstar`; `E` is the free enzyme. The name says nothing to the package, but it is
  part of its constants' names, so a name that tells what the form holds keeps the rate
  law readable: the mechanism search writes `EAB` for E with A and B bound, and
  `E_res_plusA_minusP` for a ping-pong's covalent intermediate after A bound and P left.
- A name is one form wherever it is written, so a misspelt name is another form:
  `EI + S ⇌ EIS` beside `ES + I ⇌ ESI` makes `EIS` and `ESI` two forms. The constructor
  rejects a form that no step joins to `E`, and the guard some other misspellings, such
  as one that leaves no cycle turning substrates into products, but others pass and
  write another mechanism.
- A parenthesized tuple of steps is one kinetic group, whose steps share their
  constants: `(E + I ⇌ EI, ES + I ⇌ EIS)`.
- A substrate or product binds wherever a step takes it up. Substrate inhibition is
  `ES + S ⇌ ESS`; the form's name tells the two bindings apart.

Every mechanism holds the free enzyme `E`, and its steps join every form to it. Each
declared name takes part in some step: a step takes up or gives off each substrate and
product, and binds each modifier.

Constants are named after the forms and metabolites a step joins: `K_ES_to_E_S` for
`E + S ⇌ ES`, `k_ES_to_EP` and `k_EP_to_ES` for `ES <--> EP`, `K_ESSmod_to_ES_S` for
`ES + S ⇌ ESSmod`.

An allosteric (MWC) mechanism: the enzyme is active (A) or inactive (I), the free
enzyme's inactive-to-active ratio is `L`, and allosteric regulators bind regulatory sites
apart from the catalytic one. A mechanism is allosteric when its statements carry
allosteric states: every statement of `steps:` ends with one, a step after its right-hand
side, a parenthesized group after its closing parenthesis, and `(step,) :: State` for a
group of one. Some statements with a state and others without is an error. Only an
allosteric mechanism takes:
- `oligomeric_state:` (default 1), the number of catalytic subunits, which change
  conformation together.
- `regulatory_site(n): X::State, Y::State`, one regulatory site, present `n` times, and
  its ligands in any order, each with its state. A ligand sits on one site line. It binds
  only its regulatory site, unless `modifiers:` lists it too, which lets it bind the
  catalytic site as well.

A state says how a statement's constants, or a ligand's binding, differ between the
conformations: `EqualAI` shares one constant, `NonequalAI` takes one for each (`K_A_…`
and `K_I_…`), `OnlyA` acts in the active conformation only and `OnlyI` in the inactive
one only. A catalytic statement or a ligand takes any of the four, as long as the active
conformation still turns substrate into product. The conformations interconvert only at
free E.

When the macro's code runs, it rejects a mechanism whose rate law is not defined where data
sit, with an error that names the steps at fault: no cycle of steps turns substrates into
products; a cycle of rapid-equilibrium steps does, so the rate has no finite value; or, with
every product, every substrate, one reactant or every modifier absent, the steady state is
not unique or is not the limit the rate law gives as the rapid-equilibrium steps speed up.
It also rejects a law with no value there (0/0).

A step with a state is switched off in the conformation that lacks it, and it keeps its
kind, rapid equilibrium or steady state, in both. The steps that conformation keeps join
the forms into islands, and every island but free E's lies higher in energy and holds no
enzyme: a switched-off step between two islands loses its rate uphill and keeps its rate
downhill, or, at steady state, may stop; one inside an island stops at steady state, and
at rapid equilibrium cannot be switched off, so the macro rejects it. The macro looks for
energies of the islands under which the rate law is the limit, and rejects the states when
it finds none (under each, a face's steady state depends on where the enzyme starts, or
the macro cannot show that the limit does not depend on how fast the switched-off
constants vanish), or when the states allow more choices than it searches and none of the
few it tries first gives the limit (see [How an excluded step is switched off](@ref)).

Below, the inhibitor I binds E and ES with one constant, and a second S binds ES, which
inhibits at high substrate:

```jldoctest
julia> using EnzymeRates

julia> m = @enzyme_mechanism begin
           substrates: S
           products: P
           modifiers: I
           steps: begin
               E + S ⇌ ES
               ES <--> EP
               EP ⇌ E + P
               (E + I ⇌ EI, ES + I ⇌ EIS)
               ES + S ⇌ ESSmod
           end
       end;

julia> foreach(println, parameters(m))
K_EI_to_E_I
K_EP_to_E_P
K_ESSmod_to_ES_S
K_ES_to_E_S
k_ES_to_EP
```

An allosteric mechanism with two subunits and a regulatory site for an activator X and an
inhibitor Y:

```jldoctest
julia> using EnzymeRates

julia> am = @enzyme_mechanism begin
           substrates: A, B
           products: P, Q
           modifiers: I
           oligomeric_state: 2
           steps: begin
               E + A ⇌ EA                   :: EqualAI
               EA + B ⇌ EAB                 :: NonequalAI
               EAB <--> EPQ                 :: OnlyA
               EPQ ⇌ EQ + P                 :: EqualAI
               EQ ⇌ E + Q                   :: EqualAI
               (E + I ⇌ EI, EA + I ⇌ EAI)   :: EqualAI
               EQ + P ⇌ EPmodQ              :: EqualAI
           end
           regulatory_site(2): X::OnlyA, Y::OnlyI
       end;

julia> foreach(println, parameters(am))
K_A_EAB_to_EA_B
K_EA_to_E_A
K_EI_to_E_I
K_EPQ_to_EQ_P
K_EPmodQ_to_EQ_P
K_EQ_to_E_Q
K_I_EAB_to_EA_B
k_A_EAB_to_EPQ
K_A_Xreg
K_I_Yreg
L
```
"""
macro enzyme_mechanism(block)
    return :($_checked_mechanism($(raw_mechanism(block))))
end
