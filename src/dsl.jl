# ABOUTME: The DSL macros @enzyme_reaction, @enzyme_mechanism and @allosteric_mechanism,
# ABOUTME: which read a block into plain data, and a mechanism printed as its macro text.

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
- `competitive_inhibitors:` lists the `CompetitiveInhibitor`s, names that bind the
  catalytic site without reacting. A substrate or product listed here may also bind as
  its modifier copy, which a mechanism writes `X::Mod`.
- `allosteric_regulators:` lists the `AllostericRegulator`s, names that bind regulatory
  sites. A type tag, `A::Activator` or `A::Inhibitor`, pins the regulator to the
  allosteric state that activates or inhibits; an untagged one (`:unspecified`) may take
  either. Only these entries take a type tag. A name may be both a competitive
  inhibitor and an allosteric regulator.
- A regulator may list its multiplicities, each once, `I(1, 2)` or `A(4)`; a bare name
  takes the `oligomeric_state:` ones. The search reads only the mechanism's catalytic
  multiplicity, which it also gives every regulatory site; the per-regulator lists are
  stored but not read.
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
           allosteric_regulators: A::Activator, I(2)::Inhibitor
           oligomeric_state: (1, 2)
       end
@enzyme_reaction begin
    substrates: ATP[C10H16N5O13P3], S[C6H12O6]
    products: ADP[C10H15N5O10P2], P[C6H13O9P]
    competitive_inhibitors: P
    allosteric_regulators: A::Activator, I(2)::Inhibitor
    oligomeric_state: (1, 2)
end
```
"""
macro enzyme_reaction(block)
    (; reactants, regs, mults, shared) = _parse_reaction_block(block)
    return _reaction_expr(reactants, regs, mults, shared)
end

"""
Parse the `@enzyme_reaction` body. Returns a `NamedTuple`:
- `reactants ::Vector{Tuple{Symbol, Symbol, Vector{Pair{Symbol,Int}}}}` —
  `(type, name, atoms)` with `type ∈ (:Substrate, :Product)` and `atoms` the
  `element => count` pairs.
- `regs  ::Vector{Tuple{Symbol, Symbol, Vector{Int}, Symbol}}` —
  `(type, name, mults, reg_type)` where `type ∈ (:CompetitiveInhibitor,
  :AllostericRegulator)`, a bare entry takes the catalytic multiplicities, and
  `reg_type ∈ (:activator, :inhibitor, :unspecified)`.
- `mults ::Vector{Int}` — allowed catalytic multiplicities (default `[1]`).
- `shared ::Vector{Tuple{Symbol,Symbol}}` — `shared_catalytic_site:` pairs.
"""
function _parse_reaction_block(block)
    block isa Expr && block.head === :block ||
        error("@enzyme_reaction: expected a `begin ... end` block, got $block")
    reactants = Tuple{Symbol, Symbol, Vector{Pair{Symbol,Int}}}[]
    regs  = Tuple{Symbol, Symbol, Union{Nothing, Vector{Int}}, Symbol}[]
    mults::Union{Nothing, Vector{Int}} = nothing
    shared = Tuple{Symbol,Symbol}[]

    for arg in block.args
        arg isa LineNumberNode && continue
        label, values = _parse_labeled_line(arg)
        if label === :substrates || label === :products
            type = label === :substrates ? :Substrate : :Product
            for (name, atoms) in _parse_atom_bracket_entries(values, label)
                push!(reactants, (type, name, atoms))
            end
        elseif label === :competitive_inhibitors
            append!(regs, _parse_regulator_entries(values, :CompetitiveInhibitor))
        elseif label === :allosteric_regulators
            append!(regs, _parse_regulator_entries(values, :AllostericRegulator))
        elseif label === :oligomeric_state
            mults === nothing ||
                error("@enzyme_reaction: `oligomeric_state:` given more than once.")
            mults = _parse_multiplicity_tuple(values, label)
        elseif label === :shared_catalytic_site
            append!(shared, _parse_shared_site_pairs(values))
        else
            error("@enzyme_reaction: unknown label `$label:`. Valid labels: " *
                  "substrates, products, competitive_inhibitors, " *
                  "allosteric_regulators, oligomeric_state, shared_catalytic_site.")
        end
    end

    any(r -> r[1] === :Substrate, reactants) ||
        error("@enzyme_reaction: `substrates:` not specified.")
    any(r -> r[1] === :Product, reactants) ||
        error("@enzyme_reaction: `products:` not specified.")
    catalytic_mults = something(mults, [1])
    regs = [(t, n, something(ms, catalytic_mults), rt) for (t, n, ms, rt) in regs]
    (; reactants, regs, mults = catalytic_mults, shared)
end

"""
Parse `S[C6H12O6]` entries, each a name and one chemical formula in brackets, into
`(name, atoms)` pairs, `atoms` a `Vector{Pair{Symbol,Int}}` of `element => count`.
"""
function _parse_atom_bracket_entries(values, label)
    map(values) do v
        v isa Expr && v.head === :ref && length(v.args) == 2 && v.args[1] isa Symbol ||
            error("@enzyme_reaction: `$label:` lists each reactant with one chemical " *
                  "formula in brackets, such as `S[C6H12O6]`; got `$v`.")
        (v.args[1]::Symbol, _parse_chemical_formula(string(v.args[2])))
    end
end

"""
Parse `shared_catalytic_site:` entries. Each value is a two-name tuple
`(sub, prod)`; roles are resolved by the `EnzymeReaction` constructor.
Returns `Vector{Tuple{Symbol,Symbol}}`.
"""
function _parse_shared_site_pairs(values)
    map(values) do v
        v isa Expr && v.head === :tuple && length(v.args) == 2 &&
            all(a -> a isa Symbol, v.args) || error(
                "@enzyme_reaction: `shared_catalytic_site:` each entry must be " *
                "a `(substrate, product)` pair of two names; got $v.")
        (v.args[1]::Symbol, v.args[2]::Symbol)
    end
end

"""
    _parse_chemical_formula(s::String)

Parse a chemical formula string like `"C6H12O6"` into `element => count` pairs.
Elements are identified by an uppercase letter optionally followed by lowercase
letters, then an optional integer count (defaults to 1).
"""
function _parse_chemical_formula(s::String)
    occursin(r"^([A-Z][a-z]*\d*)*$", s) || error("Invalid chemical formula: \"$s\"")
    [Symbol(m[1]) => (isempty(m[2]) ? 1 : parse(Int, m[2]::SubString))
     for m in eachmatch(r"([A-Z][a-z]*)(\d*)", s)]
end

"""`v` as an `Int`; errors unless it is a positive integer literal, which `true` is not.
`what` names it."""
_positive_int(v, what) = v isa Integer && !(v isa Bool) && v >= 1 ? Int(v) :
    error("$what must be a positive Int, got $v.")

"""
Parse `R`, `R(1, 2)`, or `R(4)` entries, optionally tagged with a regulator
type: `R::Activator`, `R::Inhibitor`, or `R(1, 2)::Inhibitor`, into
`(type, name, mults, reg_type)`. Bare `R` produces `nothing` mults (filled from the
`oligomeric_state:` multiplicities, default `[1]`) and reg_type `:unspecified`.
`RegulatorMults` rejects a type tag on a competitive inhibitor.
"""
function _parse_regulator_entries(values, type::Symbol)
    map(values) do v
        reg_type, inner = :unspecified, v
        if v isa Expr && v.head === :(::)
            inner, tag = v.args
            reg_type = tag === :Activator ? :activator :
                   tag === :Inhibitor ? :inhibitor :
                   error("@enzyme_reaction: unknown regulator type tag ::$tag " *
                         "on $v; use ::Activator or ::Inhibitor.")
        end
        inner isa Symbol && return (type, inner, nothing, reg_type)
        inner isa Expr && inner.head === :call && length(inner.args) >= 2 &&
            inner.args[1] isa Symbol ||
            error("@enzyme_reaction: cannot parse regulator entry $v. " *
                  "Use `R` or `R(1, 2)`.")
        name = inner.args[1]::Symbol
        what = "@enzyme_reaction: regulator $name multiplicity"
        (type, name, Int[_positive_int(a, what) for a in inner.args[2:end]], reg_type)
    end
end

"""Parse `(1, 2, 4)` or `4` into `Vector{Int}`, each multiplicity listed once: a bare
regulator takes the list as its own."""
function _parse_multiplicity_tuple(values, label)
    v = values[1]
    length(values) == 1 && v != Expr(:tuple) ||
        error("@enzyme_reaction: `$label:` takes one Int or a tuple of one or more " *
              "Ints, such as `(1, 2, 4)`; got $(join(values, ", ")).")
    ms = Int[_positive_int(a, "@enzyme_reaction: `$label:` entry")
             for a in (v isa Expr && v.head === :tuple ? v.args : (v,))]
    allunique(ms) ||
        error("@enzyme_reaction: `$label:` lists a multiplicity more than once; got $v.")
    ms
end

"""
Build the `EnzymeReaction` `Expr` from parsed declarations: `reactants` as
`(type, name, atoms)` and `regulators` as `(type, name, mults, reg_type)`, each
`type` naming the `Metabolite` subtype to construct. The macro returns it unescaped, so
hygiene resolves its names in EnzymeRates itself and the macro works in a module that
imports the package under another name or imports only the macro.
"""
function _reaction_expr(reactants, regulators, mults, shared)
    reactant_exprs = [
        :(EnzymeRates.ReactantAtoms(EnzymeRates.$type($(QuoteNode(name))),
            Pair{Symbol,Int}[$((:($(QuoteNode(e)) => $c) for (e, c) in atoms)...)]))
        for (type, name, atoms) in reactants]
    regulator_exprs = [
        :(EnzymeRates.RegulatorMults(EnzymeRates.$type($(QuoteNode(name))),
                                     Int[$(ms...)], $(QuoteNode(reg_type))))
        for (type, name, ms, reg_type) in regulators]
    shared_exprs = [:(($(QuoteNode(a)), $(QuoteNode(b)))) for (a, b) in shared]
    :(EnzymeRates.EnzymeReaction(
        EnzymeRates.ReactantAtoms[$(reactant_exprs...)],
        EnzymeRates.RegulatorMults[$(regulator_exprs...)], Int[$(mults...)];
        shared_catalytic_site = Tuple{Symbol,Symbol}[$(shared_exprs...)]))
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
compiled form that [`rate_equation`](@ref) evaluates without allocating.

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
    return :($_compiled_mechanism($(_raw_mechanism(block, false))))
end

"""
Parse a labeled-line, returning `(label, values_vector)`. `values_vector` is the
list of args after the label, in source order. Values are returned as parsed
(bare Symbols, `name::Tag`, `Sym[atoms]`, `begin ... end` blocks, literals, ...);
each caller checks the shape it expects. Handles both Julia parse shapes:
  - `Expr(:call, :(:), label, value)` — single labeled value.
  - `Expr(:tuple, Expr(:call, :(:), label, first), rest...)` — multi-element labeled.
"""
function _parse_labeled_line(arg)
    head, rest = arg isa Expr && arg.head == :tuple && !isempty(arg.args) ?
                 (arg.args[1], arg.args[2:end]) : (arg, Any[])
    head isa Expr && head.head == :call && head.args[1] == :(:) ||
        error("Expected `label: value` or `label: v1, v2, ...`; got $arg")
    head.args[2], Any[head.args[3], rest...]
end

"""
Coerce labeled-line values to bare Symbols. Reject any other value, with a hint where an
atom bracket (`S[C]`) or a role tag (`P::Mod`) goes.
"""
function _bare_symbols_from_values(values, label)
    for v in values
        v isa Symbol && continue
        name, untagged = _bare_name(v), _untagged(v)
        hint = !(v isa Expr) ? "" :
               untagged isa Expr && untagged.head === :ref ?
                   " Atom brackets belong in @enzyme_reaction." :
               v.head === :(::) && name !== nothing ?
                   " List `$name` bare; a role tag goes on the metabolite in a step." : ""
        error("`$label:` expects bare names; got $(join(values, ", ")).$hint")
    end
    Symbol[values...]
end

"""`x` without the `::` tags written after it, however many: `X` for `X::A::B`."""
_untagged(x) = x isa Expr && x.head === :(::) ? _untagged(x.args[1]) : x

"""
The name `x` writes under its `::` tags (`X` for `X::A::B`), or `nothing` when no name
stands there. It unwraps only a multiplicity call: a Symbol head that `Base.isoperator`
does not flag, applied to one or more Integers (`X` for `X(2)::A` and `X(1, 2)::A`). That
test passes the word operators, so `in(2)::A` gives `in`. Every other call gives `nothing`
(`(A + B)::A`, `(A in B)::A`, `(-(2))::A`), as does a number (`1::A`).
"""
function _bare_name(x)
    x = _untagged(x)
    x isa Expr && x.head === :call && length(x.args) >= 2 && x.args[1] isa Symbol &&
        !Base.isoperator(x.args[1]) && all(a -> a isa Integer, x.args[2:end]) &&
        return x.args[1]
    x isa Symbol ? x : nothing
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
compiled form that [`rate_equation`](@ref) evaluates without allocating.

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
one only. A catalytic statement takes `EqualAI`, `NonequalAI` or `OnlyA`; a ligand takes
any of the four.

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
    return :($_compiled_mechanism($(_raw_mechanism(block, true))))
end

# ─── Mechanism macro blocks as written ─────────────────────────────────

"""
    _raw_mechanism(block, allosteric::Bool) -> RawMechanism or RawAllostericMechanism

The mechanism a mechanism macro's `begin ... end` block writes, as raw input:
`@allosteric_mechanism`'s grammar when `allosteric`, else `@enzyme_mechanism`'s. It reads
the block's expression without evaluating it or building code, and throws the DSL's
syntax errors; the constructors (`Mechanism(raw)`, `AllostericMechanism(raw)`) and the
boundary check (`_check_well_posed`) check the mechanism.

`substrates:`, `products:` and `modifiers:` list bare names and may repeat; `steps:`
holds one `begin ... end` block (`_raw_steps`). `@allosteric_mechanism` adds
`allosteric_regulators:` (bare names), `catalytic_multiplicity:` (one positive Int, 1
when omitted) and `regulatory_site(n):` lines, each listing its ligands as `X::Tag`
(`_raw_site`); every allosteric regulator sits on exactly one of them.
"""
function _raw_mechanism(block, allosteric::Bool)
    macro_name = allosteric ? "@allosteric_mechanism" : "@enzyme_mechanism"
    block isa Expr && block.head === :block ||
        error("$macro_name: expected a `begin ... end` block, got $block")
    substrates, products, modifiers, regulators = Symbol[], Symbol[], Symbol[], Symbol[]
    steps_block, multiplicity, sites = nothing, nothing, RawSite[]
    for arg in block.args
        arg isa LineNumberNode && continue
        label, values = _parse_labeled_line(arg)
        if label === :substrates
            append!(substrates, _bare_symbols_from_values(values, label))
        elseif label === :products
            append!(products, _bare_symbols_from_values(values, label))
        elseif label === :modifiers
            append!(modifiers, _bare_symbols_from_values(values, label))
        elseif label === :steps
            steps_block === nothing ||
                error("$macro_name: `steps:` given more than once.")
            length(values) == 1 || error("$macro_name: `steps:` takes one " *
                                         "`begin ... end` block; got $(length(values)) " *
                                         "values.")
            steps_block = only(values)
            steps_block isa Expr && steps_block.head === :block ||
                error("$macro_name: `steps:` takes a `begin ... end` block of steps; " *
                      "got $steps_block.")
        elseif allosteric && label === :allosteric_regulators
            tagged = findfirst(v -> v isa Expr && v.head === :(::), values)
            if tagged !== nothing
                name = _bare_name(values[tagged])
                error("$macro_name: `allosteric_regulators:` lists names; got " *
                      "`$(values[tagged])`. A ligand's tag goes on its " *
                      "`regulatory_site(n):` line" * (name === nothing ? "." :
                      ", as in `regulatory_site(2): $name::OnlyA`."))
            end
            append!(regulators, _bare_symbols_from_values(values, label))
        elseif allosteric && label === :catalytic_multiplicity
            multiplicity === nothing ||
                error("$macro_name: `catalytic_multiplicity:` given more than once.")
            length(values) == 1 ||
                error("$macro_name: `catalytic_multiplicity:` takes one positive Int; " *
                      "got $(join(values, ", ")).")
            multiplicity = _positive_int(only(values),
                                         "$macro_name: `catalytic_multiplicity:`")
        elseif allosteric && label isa Expr && label.head === :call &&
               label.args[1] === :regulatory_site
            push!(sites, _raw_site(label, values))
        else
            error("$macro_name: unknown label `$label:`. Valid labels: " * (allosteric ?
                  "substrates, products, modifiers, allosteric_regulators, " *
                  "catalytic_multiplicity, steps, regulatory_site(n)." :
                  "substrates, products, modifiers, steps. Allosteric declarations " *
                  "(`allosteric_regulators:`, `catalytic_multiplicity:`, " *
                  "`regulatory_site(n):`) belong in @allosteric_mechanism."))
        end
    end
    isempty(substrates) && error("$macro_name: `substrates:` not specified.")
    isempty(products) && error("$macro_name: `products:` not specified.")
    steps_block === nothing && error("$macro_name: `steps:` not specified.")
    names = (; substrates, products, modifiers, regulators)
    steps, group, tags = _raw_steps(steps_block, names, macro_name, allosteric)
    raw = RawMechanism(substrates, products, modifiers, steps, group)
    allosteric || return raw
    _check_site_ligands(regulators, sites)
    RawAllostericMechanism(raw, tags, something(multiplicity, 1), sites)
end

"""
A `regulatory_site(n): X::Tag, ...` line as written: `label` is its `regulatory_site(n)`
call, whose one argument is the site's multiplicity, and `values` its ligand entries,
each a name with its tag.
"""
function _raw_site(label, values)
    length(label.args) == 2 && !(label.args[2] isa Expr) ||
        error("@allosteric_mechanism: `regulatory_site(n):` takes one positive Int " *
              "multiplicity, as in `regulatory_site(2): A::OnlyA, I::OnlyI`; got " *
              "`$label`.")
    multiplicity = _positive_int(label.args[2],
                                 "@allosteric_mechanism: `regulatory_site` multiplicity")
    ligands, tags = Symbol[], Symbol[]
    for v in values
        v isa Expr && v.head === :(::) && length(v.args) == 2 &&
            v.args[1] isa Symbol && v.args[2] isa Symbol ||
            error("@allosteric_mechanism: `regulatory_site(n):` lists each ligand with " *
                  "its tag, as in `regulatory_site(2): A::OnlyA, I::OnlyI`; got " *
                  (v isa Expr && v.head === :block ? "a `begin ... end` block." : "`$v`."))
        push!(ligands, v.args[1])
        push!(tags, v.args[2])
    end
    RawSite(ligands, multiplicity, tags)
end

"""
Error unless the `regulatory_site(n):` lines `sites` place each allosteric regulator of
`regulators`, the names `allosteric_regulators:` lists once each, exactly once, and list
no other ligand.
"""
function _check_site_ligands(regulators::Vector{Symbol}, sites::Vector{RawSite})
    repeated = unique(n for (i, n) in enumerate(regulators)
                      if n in view(regulators, 1:i-1))
    isempty(repeated) ||
        error("@allosteric_mechanism: `allosteric_regulators:` lists " *
              "$(join(("`$n`" for n in repeated), ", ")) more than once.")
    for site in sites, l in site.ligands
        l in regulators ||
            error("@allosteric_mechanism: ligand `$l` on a `regulatory_site(n):` line is " *
                  "not declared in `allosteric_regulators:`.")
    end
    for r in regulators
        n = sum(count(==(r), site.ligands) for site in sites; init = 0)
        n == 0 && error("@allosteric_mechanism: allosteric regulator `$r` sits at no " *
                        "regulatory site; list it on one `regulatory_site(n):` line, as " *
                        "in `regulatory_site(2): $r::OnlyA`.")
        n == 1 || error("@allosteric_mechanism: allosteric regulator `$r` is listed $n " *
                        "times on `regulatory_site(n):` lines; an allosteric regulator " *
                        "sits at exactly one regulatory site.")
    end
end

"""
The steps of a `steps:` block as written, each statement's group label per step and,
for an allosteric block, each statement's tag: `(steps, group, tags)`. A statement is
one kinetic group: a step (`lhs ⇌ rhs` rapid equilibrium, `lhs <--> rhs` steady state)
or a parenthesized tuple of steps (`(step,)` for one). An allosteric statement ends in
its tag (`step :: Tag`, `(step, step) :: Tag`); a plain one carries none. `names` holds
the block's declared names.
"""
function _raw_steps(steps_block, names, macro_name, allosteric::Bool)
    steps, group, tags = RawStep[], Int[], Symbol[]
    missing_tag = "its tag, `:: EqualAI`, `:: NonequalAI` or `:: OnlyA`; add it after"
    n = 0
    for arg in steps_block.args
        arg isa LineNumberNode && continue
        n += 1
        if arg isa Expr && arg.head === :(::) && length(arg.args) == 2 &&
           arg.args[1] isa Expr && arg.args[1].head in (:tuple, :call)
            inner, tag = arg.args
            allosteric || error("$macro_name: tag annotation `$arg` is not allowed")
            inner isa Expr && inner.head === :tuple ||
                error("$macro_name: a single step carries its tag after its right-hand " *
                      "side, as in `E + S ⇌ E(S) :: EqualAI`. Got `$arg`.")
            tag isa Symbol ||
                error("$macro_name: step-group tag must be a Symbol; got $tag")
            push!(tags, tag)
            members = inner.args
        elseif arg isa Expr && arg.head === :tuple
            allosteric && error("$macro_name: parenthesized step group `$arg` is missing " *
                                "$missing_tag the closing paren.")
            members = arg.args
        elseif arg isa Expr && arg.head === :call && length(arg.args) == 3
            op, lhs, rhs = _step_parts(arg, macro_name)
            tag = _peel_tag!(rhs)
            if tag !== nothing
                tag isa Symbol || error("$macro_name: step tag must be a Symbol; got $tag")
                allosteric || error("$macro_name: tag annotation on `$arg` is not allowed")
                push!(tags, tag)
            elseif allosteric
                error("$macro_name: step `$arg` is missing $missing_tag the step " *
                      "expression.")
            end
            push!(steps, _raw_step(op, lhs, rhs, names, macro_name))
            push!(group, n)
            continue
        else
            error("$macro_name: expected step or step-group; got $arg")
        end
        isempty(members) && error("$macro_name: a parenthesized group holds no step; a " *
                                  "kinetic group holds at least one step.")
        for s in members
            push!(steps, _raw_step(_step_parts(s, macro_name)..., names, macro_name))
            push!(group, n)
        end
    end
    steps, group, tags
end

"""A step expression's operator, `⇌` or `<-->`, and the terms of its two sides; a side
is a `+` chain of terms or a single term."""
function _step_parts(expr, macro_name)
    expr isa Expr && expr.head === :call ||
        error("$macro_name: expected lhs ⇌ rhs or lhs <--> rhs; got $expr")
    op = expr.args[1]
    op === :⇌ || op === :(<-->) ||
        error("$macro_name: expected ⇌ or <--> step operator; got $op")
    length(expr.args) == 3 ||
        error("$macro_name: expected lhs ⇌ rhs or lhs <--> rhs; got $expr")
    op, _side_terms(expr.args[2]), _side_terms(expr.args[3])
end

"""The terms of a step side: the operands of a `+` chain, or the side itself."""
_side_terms(side) = side isa Expr && side.head === :call && side.args[1] === :+ ?
                    side.args[2:end] : Any[side]

"""
The tag a step's right-hand terms `rhs` end in, taken off the last term, or `nothing`.
`::` binds tighter than `+`, so `E(P) ⇌ E + P :: EqualAI` reads `P :: EqualAI` as its
last term, while `E(P) ⇌ (E + P) :: EqualAI` tags the whole side, whose terms then
replace it. The role tag `::Mod` marks a metabolite, never a step. `::Inh` is read the
same way, so `_raw_met` reports it as an unknown role tag rather than as a step's state.
"""
function _peel_tag!(rhs::Vector{Any})
    t = rhs[end]
    t isa Expr && t.head === :(::) && length(t.args) == 2 &&
        t.args[2] ∉ (:Mod, :Inh) || return nothing
    rhs[end] = t.args[1]
    length(rhs) == 1 && splice!(rhs, 1, _side_terms(t.args[1]))
    t.args[2]
end

"""
One step as written: `op` is `⇌` (rapid equilibrium) or `<-->` (steady state), `lhs`
and `rhs` the terms of its two sides. Each side holds exactly one enzyme form and any
number of metabolites (`_raw_met`): the left-hand metabolites are consumed, the
right-hand ones released. Every term of both sides is read before either side's enzyme
forms are counted, so a misspelled bound metabolite is reported as undeclared rather
than as a second enzyme form. Each form's conformation label is checked
(`_check_conformation_label`).
"""
function _raw_step(op, lhs, rhs, names, macro_name)
    is_met(t) = t isa Symbol ? _is_declared(t, names) : t isa Expr && t.head === :(::)
    read_term(t) = is_met(t) ? _raw_met(t, names, macro_name) :
                   _raw_form(t, names, macro_name)
    term_name(t) = t isa Expr && !isempty(t.args) ? t.args[1] : t
    parsed = (Any[read_term(t) for t in lhs], Any[read_term(t) for t in rhs])
    (from, consumed), (to, released) = map((lhs, rhs), parsed) do terms, side
        i = findall(x -> x isa RawForm, side)
        isempty(i) && error("$macro_name: step side has no enzyme-form term (terms: " *
                            "$(join(term_name.(terms), ", "))).")
        length(i) == 1 ||
            error("$macro_name: step side has more than one enzyme-form term " *
                  "($(term_name(terms[i[1]])), $(term_name(terms[i[2]]))); each " *
                  "elementary step has exactly one enzyme form per side.")
        _check_conformation_label(side[only(i)].conformation, macro_name)
        side[only(i)], Met[x for x in side if x isa Met]
    end
    RawStep(from, to, consumed, released, op === :⇌)
end

"""Whether `x` is a name the block declares (`names`: substrates, products, modifiers
and allosteric regulators)."""
_is_declared(x, names) = x in names.substrates || x in names.products ||
                         x in names.modifiers || x in names.regulators

"""The block's declared names (`names`), each once, in name order, for error messages."""
_declared_names(names) =
    sort!(unique!([names.substrates; names.products; names.modifiers; names.regulators]))

"""
A metabolite term as written: a declared name `X`, or `X::Mod`. A substrate or product
written bare binds as itself, a name under `modifiers:` that is no reactant as a modifier,
and `X::Mod` as the modifier copy of a reactant `X` that `modifiers:` also lists. An
allosteric regulator binds only at its regulatory site. `species` is the enzyme-form term
that binds the metabolite, or `nothing` for a free metabolite of a step.
"""
function _raw_met(t, names, macro_name, species = nothing)
    name, tag = t isa Expr && length(t.args) == 2 ? (t.args[1], t.args[2]) : (t, nothing)
    in_species = species === nothing ? "" : " in species `$species`"
    name isa Symbol || error("$macro_name: expected `name::Mod`$in_species; got $t")
    tag === nothing || tag === :Mod ||
        error("$macro_name: unknown role tag ::$tag on $name; the one role tag is " *
              "::Mod, which marks a substrate or product bound as a modifier.")
    _is_declared(name, names) ||
        error("$macro_name: " * (species === nothing ? "tagged metabolite `$name` in `$t`" :
                                 "bound metabolite `$name` in species `$species`") *
              " is not declared. Declared: $(_declared_names(names)).")
    reactant = name in names.substrates || name in names.products
    if tag === :Mod
        reactant ||
            error("$macro_name: `$name::Mod` marks a substrate or product bound as a " *
                  "modifier, and `$name` is neither; write a name under `modifiers:` " *
                  "bare.")
        name in names.modifiers ||
            error("$macro_name: `$name::Mod` names `$name`, which `modifiers:` does not " *
                  "declare; add `$name` to `modifiers:` to let it bind the catalytic " *
                  "site as a modifier.")
        return Met(name, true)
    end
    reactant && return Met(name, false)
    name in names.modifiers && return Met(name, true)
    error("$macro_name: `$name` is an allosteric regulator; it binds only at its " *
          "regulatory site.")
end

"""
An enzyme-form term as written: a bare conformation label `E`, or a call
`E(S, X::Mod; residual = A - P)` whose head is the conformation label, whose positional
arguments are the bound metabolites (`_raw_met`) and whose `residual` lists the
substrates added to and the products subtracted from the enzyme (`_read_residual!`). A
conformation label may not be a declared name.
"""
function _raw_form(t, names, macro_name)
    t isa Symbol || t isa Expr && t.head === :call && t.args[1] isa Symbol ||
        error("$macro_name: expected metabolite Symbol or species expression on step " *
              "side; got $t")
    conformation, args = t isa Symbol ? (t, Any[]) : (t.args[1], t.args[2:end])
    _is_declared(conformation, names) &&
        _label_collision_error(conformation, macro_name, "declared metabolite")
    bound, added, subtracted = Met[], Symbol[], Symbol[]
    for a in args
        if a isa Expr && a.head === :parameters
            for kw in a.args
                kw isa Expr && kw.head === :kw && kw.args[1] === :residual ||
                    error("$macro_name: unknown keyword in species `$t`: $kw. Only " *
                          "`residual = ...` is allowed.")
                _read_residual!(added, subtracted, kw.args[2], true, names, macro_name)
            end
        else
            a isa Symbol || a isa Expr && a.head === :(::) ||
                error("$macro_name: invalid entry in species `$t`: $a")
            push!(bound, _raw_met(a, names, macro_name, t))
        end
    end
    RawForm(conformation, bound, added, subtracted)
end

"""
Read a residual expression (`A`, `A - P`, `S1 + S2 - P1 - P3`, ...) into `added` (its
positive names) and `subtracted` (its negative ones); `positive` is the sign of `e`. A
residual adds substrates and subtracts products.
"""
function _read_residual!(added, subtracted, e, positive::Bool, names, macro_name)
    if e isa Symbol
        _is_declared(e, names) ||
            error("$macro_name: residual entry `$e` is not a declared metabolite. " *
                  "Declared: $(_declared_names(names)).")
        e in (positive ? names.substrates : names.products) ||
            error("$macro_name: a residual adds substrates and subtracts products; got " *
                  "`$e`.")
        push!(positive ? added : subtracted, e)
    elseif e isa Expr && e.head === :call && e.args[1] === :+ && length(e.args) >= 3
        for a in e.args[2:end]
            _read_residual!(added, subtracted, a, positive, names, macro_name)
        end
    elseif e isa Expr && e.head === :call && e.args[1] === :- && length(e.args) in (2, 3)
        length(e.args) == 3 &&
            _read_residual!(added, subtracted, e.args[2], positive, names, macro_name)
        _read_residual!(added, subtracted, e.args[end], !positive, names, macro_name)
    else
        error("$macro_name: invalid residual expression: $e")
    end
end

"""
    _compiled_mechanism(raw) -> EnzymeMechanism or AllostericEnzymeMechanism

The compiled mechanism of `raw`, the walker's value for a mechanism macro's block
(`_raw_mechanism`): the macros' expansion calls it, so its checks run when the expanded
code runs. It checks that `raw`'s names spell every string one way only
(`_check_uniquely_decodable`, as a reaction's are checked), then the constructors
(`Mechanism(raw)`, `AllostericMechanism(raw)`) and the boundary check
(`_check_well_posed`) check the mechanism.
"""
function _compiled_mechanism(raw::RawMechanism)
    _check_uniquely_decodable([raw.substrates; raw.products], raw.modifiers, "Mechanism")
    EnzymeMechanism(_well_posed(Mechanism(raw)))
end
function _compiled_mechanism(raw::RawAllostericMechanism)
    c = raw.catalytic
    _check_uniquely_decodable([c.substrates; c.products], c.modifiers, "Mechanism")
    AllostericEnzymeMechanism(_well_posed(AllostericMechanism(raw)))
end

# ─── Mechanisms as text ─────────────────────────────────────────────────

"""
Print `m` as the `@enzyme_mechanism` text that builds it: its substrates, products and
modifiers (the names its steps bind as modifiers), each list sorted and an empty one left
out, then its steps (`_print_steps`). The text is a function of `m` alone, so equal
mechanisms print one text.
"""
function Base.show(io::IO, m::Mechanism)
    print(io, "@enzyme_mechanism begin")
    _print_names(io, m)
    _print_steps(io, m, Symbol[])
    print(io, "\nend")
    nothing
end

"""
Print `am` as the `@allosteric_mechanism` text that builds it: its catalytic graph's
names, its allosteric regulators (left out when it has no site), its catalytic
multiplicity, its steps with their tags (`_print_steps`), then one `regulatory_site(n):`
line per site.
"""
function Base.show(io::IO, am::AllostericMechanism)
    print(io, "@allosteric_mechanism begin")
    _print_names(io, am.catalytic)
    _print_line(io, "allosteric_regulators", sort(allosteric_regulators(am)))
    _print_line(io, "catalytic_multiplicity", [am.catalytic_multiplicity])
    _print_steps(io, am.catalytic, am.cat_allo_states)
    for s in am.regulatory_sites
        _print_line(io, "regulatory_site($(s.multiplicity))",
                    ["$l::$t" for (l, t) in zip(s.ligands, s.allo_states)])
    end
    print(io, "\nend")
    nothing
end

# A compiled mechanism prints its mechanism's text. Each method compiles once, whatever
# the Sig: it reads the Sig at run time (`Mechanism(em)`), and `@nospecializeinfer` keeps
# a caller compiled for one Sig, such as `sprint`, from inferring a show for that Sig.
Base.@nospecializeinfer Base.show(io::IO, @nospecialize(em::EnzymeMechanism)) =
    show(io, Mechanism(em))
Base.@nospecializeinfer Base.show(io::IO, @nospecialize(aem::AllostericEnzymeMechanism)) =
    show(io, AllostericMechanism(aem))

"""Print the line `label: item, item, ...` of a mechanism's block, or nothing when
`items` is empty."""
function _print_line(io::IO, label::String, items::Vector)
    isempty(items) || print(io, "\n    ", label, ": ", join(items, ", "))
    nothing
end

"""Print the substrates:, products: and modifiers: lines of `m`."""
function _print_names(io::IO, m::Mechanism)
    _print_line(io, "substrates", m.substrates)
    _print_line(io, "products", m.products)
    _print_line(io, "modifiers", modifiers(m))
end

"""
Print the `steps:` block of `m`: one line per kinetic group in group order, a group of
one step as the step (`_step_text`), a larger group as the tuple of its steps in edge
order, and, when `tags` holds one tag per group, each line ending in its group's tag.
"""
function _print_steps(io::IO, m::Mechanism, tags::Vector{Symbol})
    print(io, "\n    steps: begin")
    for (g, members) in enumerate(_group_members(m))
        text = join((_step_text(m, k) for k in members), ", ")
        print(io, "\n        ", length(members) == 1 ? text : "($text)")
        isempty(tags) || print(io, " :: ", tags[g])
    end
    print(io, "\n    end")
    nothing
end
