# ABOUTME: Reads @enzyme_reaction and mechanism macro blocks into plain data without
# ABOUTME: checking them, and reads the text of a mechanism back as that mechanism.

"""
Parse the `@enzyme_reaction` body. Returns a `NamedTuple`:
- `reactants ::Vector{Tuple{Symbol, Symbol, Vector{Pair{Symbol,Int}}}}` —
  `(role, name, atoms)` with `role ∈ (:substrate, :product)` and `atoms` the
  `element => count` pairs.
- `regs  ::Vector{Tuple{Symbol, Symbol, Symbol}}` —
  `(kind, name, effect)` where `kind ∈ (:competitive, :allosteric)` and
  `effect ∈ (:activator, :inhibitor, :unspecified)`.
- `mults ::Vector{Int}` — allowed oligomeric states (default `[1]`).
- `shared ::Vector{Tuple{Symbol,Symbol}}` — `shared_catalytic_site:` pairs.
"""
function _parse_reaction_block(block)
    block isa Expr && block.head === :block ||
        error("@enzyme_reaction: expected a `begin ... end` block, got $block")
    reactants = Tuple{Symbol, Symbol, Vector{Pair{Symbol,Int}}}[]
    regs  = Tuple{Symbol, Symbol, Symbol}[]
    mults::Union{Nothing, Vector{Int}} = nothing
    shared = Tuple{Symbol,Symbol}[]

    for arg in block.args
        arg isa LineNumberNode && continue
        label, values = _parse_labeled_line(arg)
        if label === :substrates || label === :products
            role = label === :substrates ? :substrate : :product
            for (name, atoms) in _parse_atom_bracket_entries(values, label)
                push!(reactants, (role, name, atoms))
            end
        elseif label === :competitive_inhibitors
            append!(regs, _parse_regulator_entries(values, :competitive))
        elseif label === :allosteric_regulators
            append!(regs, _parse_regulator_entries(values, :allosteric))
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
    (; reactants, regs, mults = something(mults, [1]), shared)
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
Parse `R` entries, optionally tagged with a regulator type, `R::Activator` or
`R::Inhibitor`, into `(kind, name, effect)`; an untagged `R` has effect
`:unspecified`. A multiplicity list such as `R(1, 2)` is an error: a regulator takes the
oligomeric states that `oligomeric_state:` lists. `Regulator` rejects a type tag on a
competitive inhibitor.
"""
function _parse_regulator_entries(values, kind::Symbol)
    map(values) do v
        effect, inner = :unspecified, v
        if v isa Expr && v.head === :(::)
            inner, tag = v.args
            effect = tag === :Activator ? :activator :
                   tag === :Inhibitor ? :inhibitor :
                   error("@enzyme_reaction: unknown regulator type tag ::$tag " *
                         "on $v; use ::Activator or ::Inhibitor.")
        end
        inner isa Symbol && return (kind, inner, effect)
        inner isa Expr && inner.head === :call && length(inner.args) >= 2 &&
            inner.args[1] isa Symbol && all(a -> a isa Integer, inner.args[2:end]) &&
            error("@enzyme_reaction: regulator $(inner.args[1]) lists multiplicities, " *
                  "got `$v`; a regulator takes the oligomeric states that " *
                  "`oligomeric_state:` lists, so list $(inner.args[1]) bare.")
        error("@enzyme_reaction: cannot parse regulator entry $v. Use `R`.")
    end
end

"""Parse `(1, 2, 4)` or `4` into `Vector{Int}`."""
function _parse_multiplicity_tuple(values, label)
    v = values[1]
    length(values) == 1 && v != Expr(:tuple) ||
        error("@enzyme_reaction: `$label:` takes one Int or a tuple of one or more " *
              "Ints, such as `(1, 2, 4)`; got $(join(values, ", ")).")
    Int[_positive_int(a, "@enzyme_reaction: `$label:` entry")
        for a in (v isa Expr && v.head === :tuple ? v.args : (v,))]
end

"""
    reaction_expr(block) -> Expr

The `EnzymeReaction` `Expr` that `@enzyme_reaction`'s `block` writes, from its parsed
declarations (`_parse_reaction_block`): `reactants` as `(role, name, atoms)` and
`regs` as `(kind, name, effect)`. The macro returns it unescaped, so hygiene resolves its
names in `DSL` itself and the macro works in a module that imports the package under
another name or imports only the macro.
"""
function reaction_expr(block)
    (; reactants, regs, mults, shared) = _parse_reaction_block(block)
    reactant_exprs(role) = [
        :(Mechanisms.Reactant($(QuoteNode(name)),
            Pair{Symbol,Int}[$((:($(QuoteNode(e)) => $c) for (e, c) in atoms)...)]))
        for (r, name, atoms) in reactants if r === role]
    regulator_exprs = [
        :(Mechanisms.Regulator(
            $(QuoteNode(name)), $(QuoteNode(kind)), $(QuoteNode(effect))))
        for (kind, name, effect) in regs]
    shared_exprs = [:(($(QuoteNode(a)), $(QuoteNode(b)))) for (a, b) in shared]
    :(Mechanisms.EnzymeReaction(
        Mechanisms.Reactant[$(reactant_exprs(:substrate)...)],
        Mechanisms.Reactant[$(reactant_exprs(:product)...)],
        Mechanisms.Regulator[$(regulator_exprs...)], Int[$(mults...)];
        shared_catalytic_site = Tuple{Symbol,Symbol}[$(shared_exprs...)]))
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
the mechanism. The value holds the declared `substrates`, `products` and `modifiers`; the
steps as `forms` and `steps` (`_raw_steps`): each step joins two forms of `forms` in its
written direction, takes up and gives off metabolites by their place among
`[substrates; products; modifiers]`, and its `kinetic_group` is the label of its
statement; the allosteric state of each group label (`allosteric_states`, empty when no
statement carries one); the `oligomeric_state`, 1 when not given; and the
`regulatory_sites`, each as written (`_raw_site`).

`substrates:`, `products:` and `modifiers:` list bare names and may repeat; `steps:`
holds one `begin ... end` block (`_raw_steps`). A mechanism is allosteric when its
statements carry allosteric states, every statement one; only then does it take
`oligomeric_state:` (one positive Int) and `regulatory_site(n):` lines, each listing its
ligands as `X::State`, each ligand on one line. An `ArgumentError` says when some
statements carry a state and others none, or when a mechanism without states gives
`oligomeric_state:` or a `regulatory_site(n):` line.
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
The steps of a `steps:` block as written and each statement's allosteric state: `(forms,
steps, states)`. `forms` names each enzyme form a step side writes, once, in the order the
sides first write them; each step joins two of them and takes up and gives off metabolites
by their place among the declared names, `[names.substrates; names.products;
names.modifiers]` (`_raw_step!`), and its `kinetic_group` is its statement's label, the
statements numbered from 1. A statement is one kinetic group: a step (`lhs ⇌ rhs` rapid
equilibrium, `lhs <--> rhs` steady state) or a parenthesized tuple of steps (`(step,)` for
one), and it may end in its state (`step :: State`, `(step, step) :: State`), a step inside
a group carrying none of its own. Either every statement carries a state, and `states`
holds them in statement order, or none does and `states` is empty; an `ArgumentError`
names a statement without a state when others carry one. `names` holds the block's
declared names.
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

"""The terms of a step side: the operands of a `+` chain, or the side itself."""
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

"""Whether `x` is a metabolite the block declares (`names`: substrates, products and
modifiers)."""
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
