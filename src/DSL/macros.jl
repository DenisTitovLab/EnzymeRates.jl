# ABOUTME: The DSL macros @enzyme_reaction, @enzyme_mechanism and @allosteric_mechanism,
# ABOUTME: which check a mechanism with the guard and compile it.

"""
    @enzyme_reaction begin
        substrates: S[C6H12O6], ...
        products: P[C6H13O9P], ...
        competitive_inhibitors: I, ...
        allosteric_regulators: A, ...
        oligomeric_state: 2
        shared_catalytic_site: (S, P), ...
    end

Build the [`EnzymeReaction`](@ref) that the mechanism search starts from: its reactants
with their atoms, and the regulators and oligomeric states the search may use.

- `substrates:` and `products:` (required) list each reactant with its atoms, one
  chemical formula in brackets: `S[C6H12O6]`. The atoms of the two sides balance, and a
  name is a substrate or a product, not both.
- `competitive_inhibitors:` lists the competitive inhibitors, names that bind the
  catalytic site without reacting. A substrate or product listed here may also bind as
  its modifier copy, which a mechanism writes `X::Mod`.
- `allosteric_regulators:` lists the allosteric regulators, names that bind regulatory
  sites. A type tag, `A::Activator` or `A::Inhibitor`, pins the regulator to the
  allosteric state that activates or inhibits; an untagged one (`:unspecified`) may take
  either. Only these entries take a type tag. A name may be both a competitive
  inhibitor and an allosteric regulator.
- `oligomeric_state:` gives the catalytic multiplicities the search may assign: one
  positive Int, or a tuple of distinct ones such as `(1, 2, 4)`. It defaults to 1 and is
  given at most once.
- `shared_catalytic_site:` lists `(substrate, product)` pairs, in either order, that bind
  one catalytic site, so never both at once.
- The list labels may repeat; their entries accumulate.

No metabolite is named `E_total`, `Keq`, `L`, `N`, `params`, `concs` or `v`, or holds
`_`: a rate equation keeps those names, and its constant names join names with `_`. Nor
is one named `E`, the label of the free enzyme form; a name is an identifier that is no
keyword. An enzyme form's name joins its bound metabolites' names without a separator,
so the names may not spell one string two ways: with `Ac`, `CoA` and `AcCoA`, E with Ac
and CoA bound and E with AcCoA bound would both be `EAcCoA`. `@enzyme_reaction` rejects
such a reaction and names both spellings (`Ac + CoA and AcCoA both spell AcCoA`); rename
one (`Acetyl`). A reactant's modifier copy is named `Xmod`, which counts as a name too.

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
        steps: begin
            E + S ⇌ E(S)
            ...
        end
    end

Build a mechanism from its elementary steps. The value is an `EnzymeMechanism`, the
compiled form that [`rate_equation`](@ref) evaluates without allocating. It prints as the
block that builds it, each list sorted and one line per kinetic group in canonical order.

Labels:
- `substrates:` and `products:` (required) name the reactants, without atoms; atoms
  belong in [`@enzyme_reaction`](@ref).
- `modifiers:` names the metabolites that bind the enzyme without reacting, such as
  inhibitors and activators. A substrate or product listed here may bind as a modifier
  as well.
- `steps:` (required) holds the steps in one `begin ... end` block.
- `substrates:`, `products:` and `modifiers:` may repeat; their names accumulate. The
  allosteric labels and step tags belong in [`@allosteric_mechanism`](@ref).

Steps, one per line:
- `lhs ⇌ rhs` is a rapid-equilibrium step and `lhs <--> rhs` a steady-state one; either
  direction may be written.
- Each side holds one enzyme form and the metabolites the step takes up (left) or gives
  off (right): `E + S ⇌ E(S)`, `E(A) + B <--> E(Q) + P`.
- An enzyme form is a conformation label with its bound metabolites in parentheses:
  `E`, `E(S)`, `E(A, B)`. `E` is the free enzyme; another conformation takes another
  label, a capital letter followed by lowercase letters and digits such as `Estar` or
  `F2`, never `A` or `I` (constant names tag the active and inactive conformations with
  `A_` and `I_`) and never a metabolite's name. A label such as `ES` reads as an opaque
  name and is rejected; write `E(S)`.
- `residual =` records the atoms a covalent intermediate keeps, the substrates it took
  up minus the products it gave off (ping-pong): `E(; residual = A - P)` is the enzyme
  without bound metabolites after A bound and P left.
- A parenthesized tuple of steps is one kinetic group, whose steps share their
  constants: `(E + I ⇌ E(I), E(S) + I ⇌ E(I, S))`.
- A name binds in its own role; a name that is only a modifier is written bare. A
  substrate or product that `modifiers:` lists binds as a modifier where it is written
  `X::Mod`: `E(S) + S::Mod ⇌ E(S, S::Mod)`.

Every mechanism holds the free enzyme `E`, with nothing bound and no residual, and its
steps join every form to it. Each declared name takes part in some step: a step takes up
or gives off each substrate and product, and binds each modifier.

Constants are named after the forms a step joins, each form written as its label and
its bound names without a separator: `K_ES_to_E_S` for `E + S ⇌ E(S)`, `k_ES_to_EP` and
`k_EP_to_ES` for `E(S) <--> E(P)`. A modifier copy writes `Smod`, and a residual
`_res_plusA_minusP`.

When the macro's code runs, it rejects a mechanism whose rate law is not defined where data
sit, with an error that names the steps at fault: no cycle of steps turns substrates into
products; a cycle of rapid-equilibrium steps does, so the rate has no finite value; or, with
every product, every substrate, one reactant or every modifier absent, the steady state is
not unique or is not the limit the rate law gives as the rapid-equilibrium steps speed up.
It also rejects a law with no value there (0/0).

Below, the inhibitor I binds E and E(S) with one constant, and a second S binds E(S) as
a modifier, which inhibits at high substrate:

```jldoctest
julia> using EnzymeRates

julia> m = @enzyme_mechanism begin
           substrates: S
           products: P
           modifiers: I, S
           steps: begin
               E + S ⇌ E(S)
               E(S) <--> E(P)
               E(P) ⇌ E + P
               (E + I ⇌ E(I), E(S) + I ⇌ E(I, S))
               E(S) + S::Mod ⇌ E(S, S::Mod)
           end
       end;

julia> foreach(println, parameters(m))
K_EI_to_E_I
K_EP_to_E_P
K_ESSmod_to_ES_Smod
K_ES_to_E_S
k_ES_to_EP
Keq
E_total
```
"""
macro enzyme_mechanism(block)
    return :($_compiled_mechanism($(raw_mechanism(block, false))))
end

"""
    @allosteric_mechanism begin
        substrates: S, ...
        products: P, ...
        modifiers: I, ...
        allosteric_regulators: X, ...
        catalytic_multiplicity: 2
        steps: begin
            E + S ⇌ E(S) :: EqualAI
            ...
        end
        regulatory_site(2): X::OnlyA, ...
    end

Build a two-conformation (MWC) mechanism: the enzyme is active (A) or inactive (I), the
free enzyme's inactive-to-active ratio is `L`, and allosteric regulators bind regulatory
sites apart from the catalytic one. The value is an `AllostericEnzymeMechanism`, the
compiled form that [`rate_equation`](@ref) evaluates without allocating. It prints as the
block that builds it, each list sorted, `catalytic_multiplicity:` written out and one line
per kinetic group in canonical order.

It takes the labels and step syntax of [`@enzyme_mechanism`](@ref), and:
- `allosteric_regulators:` names the metabolites that bind regulatory sites. An
  allosteric regulator binds only its regulatory site, unless `modifiers:` lists it too,
  which lets it bind the catalytic site as well.
- `catalytic_multiplicity:` (default 1) is the number of catalytic subunits.
- `regulatory_site(n): X::Tag, Y::Tag` is one regulatory site, present `n` times, and
  its ligands in any order, each with its tag. Each ligand is an allosteric regulator,
  and each allosteric regulator sits on exactly one site line.
- Every statement of `steps:` ends with its tag: a step after its right-hand side, a
  parenthesized group after its closing parenthesis, and `(step,) :: Tag` for a group of
  one.

A tag says how a statement's constants, or a ligand's binding, differ between the
conformations: `EqualAI` shares one constant, `NonequalAI` takes one for each (`K_A_…`
and `K_I_…`), `OnlyA` acts in the active conformation only and `OnlyI` in the inactive
one only. A catalytic statement takes `EqualAI`, `NonequalAI` or `OnlyA`, or `OnlyI` when
its steps bind only modifiers (`modifiers:` names or `::Mod` copies); a ligand takes any of
the four. The conformations interconvert only at free E.

A tagged step's constants vanish in the conformation that lacks it, and the macro rejects a
tag set whose limit it cannot show to exist and be unique: an `:OnlyA` or `:OnlyI` binding
whose two sides the steps that conformation keeps still join; tagged bindings that close a
cycle; forms that free E reaches only through tagged steps other than a binding toward
them; `:OnlyI` groups without which the active conformation turns nothing over; a face
where the steady state depends on where the enzyme starts; or a face where it cannot show
that the limit does not depend on how fast the tagged constants vanish. The checks of
[`@enzyme_mechanism`](@ref) apply as well.

```jldoctest
julia> using EnzymeRates

julia> am = @allosteric_mechanism begin
           substrates: A, B
           products: P, Q
           modifiers: I, P
           allosteric_regulators: X, Y
           catalytic_multiplicity: 2
           steps: begin
               E + A ⇌ E(A)                          :: EqualAI
               E(A) + B ⇌ E(A, B)                    :: NonequalAI
               E(A, B) <--> E(P, Q)                  :: OnlyA
               E(P, Q) ⇌ E(Q) + P                    :: EqualAI
               E(Q) ⇌ E + Q                          :: EqualAI
               (E + I ⇌ E(I), E(A) + I ⇌ E(A, I))    :: EqualAI
               E(Q) + P::Mod ⇌ E(P::Mod, Q)          :: EqualAI
           end
           regulatory_site(2): X::OnlyA, Y::OnlyI
       end;

julia> foreach(println, parameters(am))
K_A_EAB_to_EA_B
K_EA_to_E_A
K_EI_to_E_I
K_EPQ_to_EQ_P
K_EPmodQ_to_EQ_Pmod
K_EQ_to_E_Q
K_I_EAB_to_EA_B
k_A_EAB_to_EPQ
K_A_Xreg
K_I_Yreg
L
Keq
E_total
```
"""
macro allosteric_mechanism(block)
    return :($_compiled_mechanism($(raw_mechanism(block, true))))
end

"""
    _compiled_mechanism(raw) -> EnzymeMechanism or AllostericEnzymeMechanism

The compiled mechanism of `raw` (`_checked_mechanism`): the macros' expansion calls it,
so its checks run when the expanded code runs.
"""
_compiled_mechanism(raw::Union{Mechanisms.RawMechanism,
                               Mechanisms.RawAllostericMechanism}) =
    Mechanisms.compile_mechanism(_checked_mechanism(raw))
