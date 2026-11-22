# ABOUTME: @enzyme_mechanism and parse_mechanism, which build a mechanism from its text
# ABOUTME: and check it with the guard, and the reader of a mechanism block as written.

# ─── @enzyme_mechanism ──────────────────────────────────────────────────

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
  rejects a form that no step joins to `E`, and the guard some other slips, such as a
  name that leaves no cycle turning substrates into products, or a step that forgets a
  metabolite (`ESI <--> ES` for `ESI <--> ES + I`). Others pass and write another
  mechanism.
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
sit, with an error that names the steps at fault, when:
- a cycle of steps takes up or gives off a metabolite in another proportion than the
  reaction, as when a step forgets one;
- no cycle of steps turns substrates into products;
- a cycle of rapid-equilibrium steps does, so the rate has no finite value;
- with every product, every substrate, one reactant or every modifier absent, the steady
  state is not unique, or is not the limit the rate law gives as the rapid-equilibrium
  steps speed up, or the law has no value (0/0).

A step with a state is switched off in the conformation that lacks it, and it keeps its
kind, rapid equilibrium or steady state, in both. The steps that conformation keeps join
the forms into islands; every island but free E's lies higher in energy and holds no
enzyme. A switched-off step between two islands loses its rate uphill and keeps its rate
downhill, or, at steady state, may stop. One inside an island stops at steady state; at
rapid equilibrium it cannot be switched off, and the macro rejects it. The macro also
rejects the states when no choice of the islands' energies makes the rate law the limit,
or when the states allow more choices than it searches and none of the few it tries first
gives the limit. [How an excluded step is switched off](@ref) gives the rules.

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

# ─── Bare names of a block ──────────────────────────────────────────────

"""
Coerce labeled-line values to bare Symbols. Reject any other value, with a hint where an
atom bracket (`S[C]`) goes, or that a name tagged (`P::Mod`) is listed bare.
"""
function _bare_symbols_from_values(values, label)
    for v in values
        v isa Symbol && continue
        name, untagged = _bare_name(v), _untagged(v)
        hint = !(v isa Expr) ? "" :
               untagged isa Expr && untagged.head === :ref ?
                   " Atom brackets belong in @enzyme_reaction." :
               v.head === :(::) && name !== nothing ?
                   " List `$name` bare." : ""
        error("`$label:` expects bare names; got $(join(values, ", ")).$hint")
    end
    Symbol[values...]
end

"""`x` without the `::` tags written after it, however many: `X` for `X::A::B`."""
_untagged(x) =
    x isa Expr && x.head === :(::) && length(x.args) == 2 ? _untagged(x.args[1]) : x

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

# ─── Mechanism macro blocks as written ─────────────────────────────────

"""
    raw_mechanism(block) -> NamedTuple

The mechanism an `@enzyme_mechanism` block writes, as written. It reads the block's
expression without evaluating it or building code, and throws the DSL's syntax errors;
the constructors (`_construct`) and the boundary check (`check_valid_rate_equation`) check
the mechanism. The block's grammar is the one `@enzyme_mechanism` documents. The value
holds:
- the declared `substrates`, `products` and `modifiers`;
- the steps as `forms` and `steps` (`_raw_steps`): each step joins two forms of `forms` in
  its written direction, takes up and gives off metabolites by their place among
  `[substrates; products; modifiers]`, and its `kinetic_group` is the label of its
  statement;
- `allosteric_states`, the state of each group label, empty when no statement carries one;
- the `oligomeric_state`, 1 when not given;
- the `regulatory_sites`, each as written (`_raw_site`).
An `ArgumentError` says when some statements carry a state and others none, or when a
mechanism without states gives `oligomeric_state:` or a `regulatory_site(n):` line.
"""
function raw_mechanism(block)
    block isa Expr && block.head === :block ||
        error("@enzyme_mechanism: expected a `begin ... end` block, got $block")
    substrates, products, modifiers = Symbol[], Symbol[], Symbol[]
    steps_block, multiplicity, sites = nothing, nothing, _WrittenSite[]
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
                error("@enzyme_mechanism: `steps:` given more than once.")
            length(values) == 1 || error("@enzyme_mechanism: `steps:` takes one " *
                                         "`begin ... end` block; got $(length(values)) " *
                                         "values.")
            steps_block = only(values)
            steps_block isa Expr && steps_block.head === :block ||
                error("@enzyme_mechanism: `steps:` takes a `begin ... end` block of " *
                      "steps; got $steps_block.")
        elseif label === :oligomeric_state
            multiplicity === nothing ||
                error("@enzyme_mechanism: `oligomeric_state:` given more than once.")
            length(values) == 1 ||
                error("@enzyme_mechanism: `oligomeric_state:` takes one positive Int; " *
                      "got $(join(values, ", ")).")
            multiplicity = _positive_int(only(values),
                                         "@enzyme_mechanism: `oligomeric_state:`")
        elseif label isa Expr && label.head === :call &&
               label.args[1] === :regulatory_site
            push!(sites, _raw_site(label, values))
        else
            error("@enzyme_mechanism: unknown label `$label:`. Valid labels: substrates, " *
                  "products, modifiers, oligomeric_state, steps, regulatory_site(n).")
        end
    end
    isempty(substrates) && error("@enzyme_mechanism: `substrates:` not specified.")
    isempty(products) && error("@enzyme_mechanism: `products:` not specified.")
    steps_block === nothing && error("@enzyme_mechanism: `steps:` not specified.")
    regulators = unique!([l for site in sites for l in site.ligands])
    names = (; substrates, products, modifiers, regulators)
    forms, steps, states = _raw_steps(steps_block, names)
    if isempty(states)
        allosteric_only = "is read only for an allosteric mechanism, whose statements " *
                          "each carry an allosteric state (`:: EqualAI`, " *
                          "`:: NonequalAI`, `:: OnlyA` or `:: OnlyI`)"
        multiplicity === nothing ||
            throw(ArgumentError("@enzyme_mechanism: `oligomeric_state:` " *
                                "$allosteric_only; no statement carries one."))
        isempty(sites) ||
            throw(ArgumentError("@enzyme_mechanism: `regulatory_site(n):` " *
                                "$allosteric_only; no statement carries one."))
    end
    (; substrates, products, modifiers, forms, steps, allosteric_states = states,
     oligomeric_state = something(multiplicity, 1), regulatory_sites = sites)
end

"""A regulatory site as written: its ligands in written order, its multiplicity and each
ligand's tag."""
const _WrittenSite = NamedTuple{(:ligands, :multiplicity, :allosteric_states),
                                Tuple{Vector{Symbol}, Int, Vector{Symbol}}}

"""
A `regulatory_site(n): X::Tag, ...` line as written: `label` is its `regulatory_site(n)`
call, whose one argument is the site's multiplicity, and `values` its ligand entries,
each a name with its tag.
"""
function _raw_site(label, values)
    length(label.args) == 2 && !(label.args[2] isa Expr) ||
        error("@enzyme_mechanism: `regulatory_site(n):` takes one positive Int " *
              "multiplicity, as in `regulatory_site(2): A::OnlyA, I::OnlyI`; got " *
              "`$label`.")
    multiplicity = _positive_int(label.args[2],
                                 "@enzyme_mechanism: `regulatory_site` multiplicity")
    ligands, tags = Symbol[], Symbol[]
    for v in values
        v isa Expr && v.head === :(::) && length(v.args) == 2 &&
            v.args[1] isa Symbol && v.args[2] isa Symbol ||
            error("@enzyme_mechanism: `regulatory_site(n):` lists each ligand with " *
                  "its tag, as in `regulatory_site(2): A::OnlyA, I::OnlyI`; got " *
                  (v isa Expr && v.head === :block ? "a `begin ... end` block." : "`$v`."))
        push!(ligands, v.args[1])
        push!(tags, v.args[2])
    end
    _WrittenSite((ligands, multiplicity, tags))
end

"""
The steps of a `steps:` block as written and each statement's allosteric state, as
`(forms, steps, states)`; `names` holds the block's declared names. `forms` names each
enzyme form a step side writes, once, in the order the sides first write them. Each step
joins two of them and takes up and gives off metabolites by their place among
`[names.substrates; names.products; names.modifiers]` (`_raw_step!`). Each statement, a
step or a parenthesized group, is one kinetic group; its steps' `kinetic_group` is its
label, the statements numbered from 1. Either every statement carries a state, and
`states` holds them in statement order, or none does and `states` is empty; an
`ArgumentError` names a statement without a state when others carry one.
"""
function _raw_steps(steps_block, names)
    forms, steps = Symbol[], Mechanisms.Step[]
    # Each statement as written, and its state or `nothing`.
    statements, states = Any[], Union{Symbol, Nothing}[]
    n = 0
    for arg in steps_block.args
        arg isa LineNumberNode && continue
        n += 1
        push!(statements, arg)
        if arg isa Expr && arg.head === :(::) && length(arg.args) == 2 &&
           arg.args[1] isa Expr && arg.args[1].head in (:tuple, :call)
            inner, state = arg.args
            inner isa Expr && inner.head === :tuple ||
                error("@enzyme_mechanism: a single step carries its tag after its " *
                      "right-hand side, as in `E + S ⇌ ES :: EqualAI`. Got `$arg`.")
            state isa Symbol ||
                error("@enzyme_mechanism: step-group tag must be a Symbol; got $state")
            push!(states, state)
            members = inner.args
        elseif arg isa Expr && arg.head === :tuple
            push!(states, nothing)
            members = arg.args
        elseif arg isa Expr && arg.head === :call && length(arg.args) == 3
            op, lhs, rhs = _step_parts(arg)
            state = _peel_tag!(rhs)
            state === nothing || state isa Symbol ||
                error("@enzyme_mechanism: step tag must be a Symbol; got $state")
            push!(states, state)
            push!(steps, _raw_step!(forms, op, lhs, rhs, n, names))
            continue
        else
            error("@enzyme_mechanism: expected step or step-group; got $arg")
        end
        isempty(members) && error("@enzyme_mechanism: a parenthesized group holds no " *
                                  "step; a kinetic group holds at least one step.")
        for s in members
            op, lhs, rhs = _step_parts(s)
            _peel_tag!(rhs) === nothing ||
                error("@enzyme_mechanism: the step `$s` carries a state inside a " *
                      "parenthesized group; a group's state is written after its closing " *
                      "parenthesis, as in `(E + I ⇌ EI, ES + I ⇌ EIS) :: EqualAI`.")
            push!(steps, _raw_step!(forms, op, lhs, rhs, n, names))
        end
    end
    with = findfirst(!isnothing, states)
    without = findfirst(isnothing, states)
    with === nothing || without === nothing ||
        throw(ArgumentError(
            "@enzyme_mechanism: the statement `$(statements[without])` carries no " *
            "allosteric state, while `$(statements[with])` carries one; a mechanism is " *
            "allosteric when every statement ends in its state (`:: EqualAI`, " *
            "`:: NonequalAI`, `:: OnlyA` or `:: OnlyI`), and plain when none does."))
    forms, steps, Symbol[s for s in states if s !== nothing]
end

"""A step expression's operator, `⇌` or `<-->`, and the terms of its two sides; a side
is a `+` chain of terms or a single term."""
function _step_parts(expr)
    expr isa Expr && expr.head === :call ||
        error("@enzyme_mechanism: expected lhs ⇌ rhs or lhs <--> rhs; got $expr")
    op = expr.args[1]
    op === :⇌ || op === :(<-->) ||
        error("@enzyme_mechanism: expected ⇌ or <--> step operator; got $op")
    length(expr.args) == 3 ||
        error("@enzyme_mechanism: expected lhs ⇌ rhs or lhs <--> rhs; got $expr")
    op, _side_terms(expr.args[2]), _side_terms(expr.args[3])
end

_side_terms(side) = side isa Expr && side.head === :call && side.args[1] === :+ ?
                    side.args[2:end] : Any[side]

"""
The tag a step's right-hand terms `rhs` end in, taken off the last term, or `nothing`.
`::` binds tighter than `+`, so `EP ⇌ E + P :: EqualAI` reads `P :: EqualAI` as its last
term, while `EP ⇌ (E + P) :: EqualAI` tags the whole side, whose terms then replace it.
`::Mod` is left on its term, which `_raw_term` rejects with the text to write instead.
"""
function _peel_tag!(rhs::Vector{Any})
    t = rhs[end]
    t isa Expr && t.head === :(::) && length(t.args) == 2 && t.args[2] !== :Mod ||
        return nothing
    rhs[end] = t.args[1]
    length(rhs) == 1 && splice!(rhs, 1, _side_terms(t.args[1]))
    t.args[2]
end

"""
One step as written, in kinetic group `group`: `op` is `⇌` (rapid equilibrium) or
`<-->` (steady state), `lhs` and `rhs` the terms of its two sides (`_raw_term`). Each side
holds exactly one enzyme form and any number of declared metabolites: the left-hand
metabolites are taken up, the right-hand ones given off, each by its place among
`[names.substrates; names.products; names.modifiers]`. A form not yet in `forms` is added
to it, and the step joins the two forms by their places there.
"""
function _raw_step!(forms::Vector{Symbol}, op, lhs, rhs, group::Int, names)
    declared = [names.substrates; names.products; names.modifiers]
    (from, taken_up), (to, given_off) = map((lhs, rhs)) do terms
        side = [_raw_term(t, names) for t in terms]
        found = findall(first, side)
        isempty(found) &&
            error("@enzyme_mechanism: step side has no enzyme-form term (terms: " *
                  "$(join(last.(side), ", "))); a name that is no declared metabolite " *
                  "names an enzyme form.")
        length(found) == 1 ||
            error("@enzyme_mechanism: step side has more than one enzyme-form term " *
                  "($(join(last.(side[found]), ", "))); each elementary step has exactly " *
                  "one enzyme form per side, and a name that is no declared metabolite " *
                  "names an enzyme form.")
        form = last(side[only(found)])
        f = findfirst(==(form), forms)
        f === nothing && (push!(forms, form); f = length(forms))
        f, sort!(Int[findfirst(==(x), declared) for (is_form, x) in side if !is_form])
    end
    Mechanisms.Step(from, to, taken_up, given_off, op === :⇌, group)
end

_is_declared(x, names) = x in names.substrates || x in names.products ||
                         x in names.modifiers

"""
A term of a step side as written, `(is_form, name)`: a declared metabolite (`_is_declared`)
is `(false, name)`, and any other name `(true, name)`, the name of an enzyme form, used as
written. An allosteric regulator binds only at its regulatory site, unless `modifiers:`
lists it too. A reactant's modifier copy `X::Mod`, a form in call notation (`E(A, B)`) and
a count (`2S`) are errors that say what to write instead.
"""
function _raw_term(t, names)
    if t isa Symbol
        !_is_declared(t, names) && t in names.regulators &&
            error("@enzyme_mechanism: `$t` is an allosteric regulator; it binds only at " *
                  "its regulatory site, unless `modifiers:` lists it too.")
        return (!_is_declared(t, names), t)
    end
    if t isa Expr && t.head === :(::) && length(t.args) == 2 && t.args[1] isa Symbol
        x, tag = t.args
        tag === :Mod &&
            error("@enzyme_mechanism: `$t` writes a modifier copy, which a mechanism " *
                  "does not hold: a step that binds $x as a modifier takes $x up; write " *
                  "`$x`.")
        error("@enzyme_mechanism: unknown tag ::$tag on $x; a step ends in its " *
              "allosteric state, written after its right-hand side.")
    end
    t isa Expr && t.head === :call && t.args[1] === :* &&
        error("@enzyme_mechanism: `$t` writes a count; a step names a metabolite once " *
              "for each copy it takes up or gives off, as in `E + S + S ⇌ ESS`.")
    if t isa Expr && t.head === :call && t.args[1] isa Symbol
        name = _called_form_name(t)
        error("@enzyme_mechanism: `$t` is call notation; an enzyme form is a name, so " *
              (name === nothing ? "write one, such as `EAB` for `E(A, B)`." :
               "write `$name`."))
    end
    error("@enzyme_mechanism: expected a metabolite or an enzyme form's name on a step " *
          "side; got `$t`.")
end

"""The name of the form that the call `t` writes in call notation, for the error that asks
for one: its head, then its arguments sorted by name (`E(B, A)` gives `EAB`), when every
argument is a bare name; else `nothing`."""
function _called_form_name(t::Expr)
    args = t.args[2:end]
    all(a -> a isa Symbol, args) || return nothing
    Symbol(t.args[1], sort!(Symbol[args...])...)
end

"""
    _checked_mechanism(raw) -> Mechanism

The mechanism that `raw`, the walker's value for a mechanism block (`raw_mechanism`),
writes, checked: the constructors (`_construct`) and the boundary check
(`check_valid_rate_equation`) check the mechanism.
The macro's expansion and the parser (`parse_mechanism`) build every mechanism through
it, so they accept the same mechanisms, and its checks run when the expanded code runs.
"""
function _checked_mechanism(raw::NamedTuple)
    m = _construct(raw)
    Derivation.check_valid_rate_equation(m)
    m
end

"""
    _construct(raw) -> Mechanism

The mechanism that `raw` (`raw_mechanism`) writes, built and checked by the constructors
from its declared names, forms and steps, the allosteric state of each group label, its
oligomeric state and its regulatory sites.
"""
function _construct(raw::NamedTuple)
    ns, np = length(raw.substrates), length(raw.products)
    Mechanisms.Mechanism(
        Mechanisms.CatalyticScheme([raw.substrates; raw.products; raw.modifiers], raw.forms,
                                   collect(1:ns), collect((ns + 1):(ns + np)), raw.steps),
        Mechanisms.Allostery(raw.allosteric_states, raw.oligomeric_state,
                             [Mechanisms.RegulatorySite(s.ligands, s.multiplicity,
                                                        s.allosteric_states)
                              for s in raw.regulatory_sites]))
end

"""
    parse_mechanism(text) -> Mechanism

The mechanism `text` writes: one unqualified `@enzyme_mechanism` block, such as `show`
prints for a mechanism, with any whitespace around it. It reads the parsed expression with
the macro's walker (`raw_mechanism`), never evaluating it, and checks the mechanism as the
macro does (`_checked_mechanism`). A block the macro rejects fails with the macro's
message and text with a syntax error with `Meta.parse`'s `ParseError`; any other text,
incomplete text included, fails with this function's own. The text `show` prints for a
mechanism the macro accepts reads back as that mechanism.
"""
function parse_mechanism(text::AbstractString)
    # The parser runs in a fixed world age, where a SubString of a string type from a
    # package (CSV.jl's DataString) cannot be iterated.
    ex = Meta.parse(strip(String(text)))
    ex isa Expr && ex.head === :macrocall && length(ex.args) == 3 &&
        ex.args[1] === Symbol("@enzyme_mechanism") ||
        error("parse_mechanism: expected the text of one @enzyme_mechanism block; got $ex")
    _checked_mechanism(raw_mechanism(ex.args[3]))
end
