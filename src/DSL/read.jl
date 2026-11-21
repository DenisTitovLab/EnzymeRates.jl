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

    any(r -> r[1] === :substrate, reactants) ||
        error("@enzyme_reaction: `substrates:` not specified.")
    any(r -> r[1] === :product, reactants) ||
        error("@enzyme_reaction: `products:` not specified.")
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
`oligomeric_state:` multiplicities. `Regulator` rejects a type tag on a
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
                  "got `$v`; a regulator takes the `oligomeric_state:` multiplicities, " *
                  "so list $(inner.args[1]) bare.")
        error("@enzyme_reaction: cannot parse regulator entry $v. Use `R`.")
    end
end

"""Parse `(1, 2, 4)` or `4` into `Vector{Int}`, each multiplicity listed once."""
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
steps as `nodes` and `steps` (`_raw_steps`): step `j` joins nodes `2j - 1` and `2j`, its
written ends, in its written direction, and its `group` is the label of its statement;
the allosteric state of each group label (`allosteric_states`, empty when no statement
carries one); the `oligomeric_state`, 1 when not given; and the `regulatory_sites`, each
as written (`_raw_site`).

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
    nodes, steps, states = _raw_steps(steps_block, names)
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
    _check_site_ligands(sites)
    (; substrates, products, modifiers, nodes, steps, allosteric_states = states,
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

"""Error unless the `regulatory_site(n):` lines `sites` list each ligand once: a ligand,
an allosteric regulator, sits at exactly one regulatory site."""
function _check_site_ligands(sites::Vector{_WrittenSite})
    for r in unique(l for site in sites for l in site.ligands)
        n = sum(count(==(r), site.ligands) for site in sites)
        n == 1 || error("@enzyme_mechanism: allosteric regulator `$r` is listed $n " *
                        "times on `regulatory_site(n):` lines; an allosteric regulator " *
                        "sits at exactly one regulatory site.")
    end
end

"""
The steps of a `steps:` block as written and each statement's allosteric state: `(nodes,
steps, states)`. Each step joins its two written ends, which it adds to `nodes`
(`_raw_step!`), and its `group` is its statement's label, the statements numbered from 1.
A statement is one kinetic group: a step (`lhs ⇌ rhs` rapid equilibrium, `lhs <--> rhs`
steady state) or a parenthesized tuple of steps (`(step,)` for one), and it may end in its
state (`step :: State`, `(step, step) :: State`). Either every statement carries a state,
and `states` holds them in statement order, or none does and `states` is empty; an
`ArgumentError` names a statement without a state when others carry one. `names` holds
the block's declared names.
"""
function _raw_steps(steps_block, names)
    nodes, steps = Mechanisms.EnzymeForm[], Mechanisms.Step[]
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
                      "right-hand side, as in `E + S ⇌ E(S) :: EqualAI`. Got `$arg`.")
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
            push!(steps, _raw_step!(nodes, op, lhs, rhs, n, names))
            continue
        else
            error("@enzyme_mechanism: expected step or step-group; got $arg")
        end
        isempty(members) && error("@enzyme_mechanism: a parenthesized group holds no " *
                                  "step; a kinetic group holds at least one step.")
        for s in members
            push!(steps, _raw_step!(nodes, _step_parts(s)..., n, names))
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
    nodes, steps, Symbol[s for s in states if s !== nothing]
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
One step as written, in kinetic group `group`: `op` is `⇌` (rapid equilibrium) or
`<-->` (steady state), `lhs` and `rhs` the terms of its two sides. Each side holds
exactly one enzyme form and any number of metabolites (`_raw_met`): the left-hand
metabolites are taken up, the right-hand ones given off. The left-hand form and then the
right-hand one are added to `nodes`, and the step joins them. Every term of both sides is
read before either side's enzyme forms are counted, so a misspelled bound metabolite is
reported as undeclared rather than as a second enzyme form. Each form's conformation
label is checked (`check_conformation_label`).
"""
function _raw_step!(nodes::Vector{Mechanisms.EnzymeForm}, op, lhs, rhs, group::Int, names)
    is_met(t) = t isa Symbol ? _is_declared(t, names) : t isa Expr && t.head === :(::)
    read_term(t) = is_met(t) ? _raw_met(t, names) : _raw_form(t, names)
    term_name(t) = t isa Expr && !isempty(t.args) ? t.args[1] : t
    parsed = (Any[read_term(t) for t in lhs], Any[read_term(t) for t in rhs])
    (from, taken_up), (to, given_off) = map((lhs, rhs), parsed) do terms, side
        i = findall(x -> x isa Mechanisms.EnzymeForm, side)
        isempty(i) && error("@enzyme_mechanism: step side has no enzyme-form term " *
                            "(terms: $(join(term_name.(terms), ", "))).")
        length(i) == 1 ||
            error("@enzyme_mechanism: step side has more than one enzyme-form term " *
                  "($(term_name(terms[i[1]])), $(term_name(terms[i[2]]))); each " *
                  "elementary step has exactly one enzyme form per side.")
        Mechanisms.check_conformation_label(side[only(i)].conformation,
                                            "@enzyme_mechanism")
        side[only(i)],
        sort!(Mechanisms.Metabolite[x for x in side if x isa Mechanisms.Metabolite])
    end
    push!(nodes, from, to)
    Mechanisms.Step(length(nodes) - 1, length(nodes), taken_up, given_off, op === :⇌,
                    group)
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
function _raw_met(t, names, species = nothing)
    name, tag = t isa Expr && length(t.args) == 2 ? (t.args[1], t.args[2]) : (t, nothing)
    in_species = species === nothing ? "" : " in species `$species`"
    name isa Symbol || error("@enzyme_mechanism: expected `name::Mod`$in_species; got $t")
    tag === nothing || tag === :Mod ||
        error("@enzyme_mechanism: unknown role tag ::$tag on $name; the one role tag is " *
              "::Mod, which marks a substrate or product bound as a modifier.")
    _is_declared(name, names) ||
        error("@enzyme_mechanism: " *
              (species === nothing ? "tagged metabolite `$name` in `$t`" :
               "bound metabolite `$name` in species `$species`") *
              " is not declared. Declared: $(_declared_names(names)).")
    reactant = name in names.substrates || name in names.products
    if tag === :Mod
        reactant ||
            error("@enzyme_mechanism: `$name::Mod` marks a substrate or product bound " *
                  "as a modifier, and `$name` is neither; write a name under " *
                  "`modifiers:` bare.")
        name in names.modifiers ||
            error("@enzyme_mechanism: `$name::Mod` names `$name`, which `modifiers:` " *
                  "does not declare; add `$name` to `modifiers:` to let it bind the " *
                  "catalytic site as a modifier.")
        return Mechanisms.Metabolite(name, true)
    end
    reactant && return Mechanisms.Metabolite(name, false)
    name in names.modifiers && return Mechanisms.Metabolite(name, true)
    error("@enzyme_mechanism: `$name` is an allosteric regulator; it binds only at its " *
          "regulatory site.")
end

"""
An enzyme-form term as written, a node with each list sorted: a bare conformation label
`E`, or a call `E(S, X::Mod; residual = A - P)` whose head is the conformation label,
whose positional arguments are the bound metabolites (`_raw_met`) and whose `residual`
lists the substrates added to and the products subtracted from the enzyme
(`_read_residual!`). The constructors reject a conformation label that is a declared name.
"""
function _raw_form(t, names)
    t isa Symbol || t isa Expr && t.head === :call && t.args[1] isa Symbol ||
        error("@enzyme_mechanism: expected metabolite Symbol or species expression on " *
              "step side; got $t")
    conformation, args = t isa Symbol ? (t, Any[]) : (t.args[1], t.args[2:end])
    bound, added, subtracted = Mechanisms.Metabolite[], Symbol[], Symbol[]
    for a in args
        if a isa Expr && a.head === :parameters
            for kw in a.args
                kw isa Expr && kw.head === :kw && kw.args[1] === :residual ||
                    error("@enzyme_mechanism: unknown keyword in species `$t`: " *
                          "$(Meta.isexpr(kw, :kw) ? Expr(:(=), kw.args...) : kw). Only " *
                          "`residual = ...` is allowed.")
                _read_residual!(added, subtracted, kw.args[2], true, names)
            end
        else
            a isa Symbol || a isa Expr && a.head === :(::) ||
                error("@enzyme_mechanism: invalid entry in species `$t`: $a")
            push!(bound, _raw_met(a, names, t))
        end
    end
    Mechanisms.EnzymeForm(conformation, sort!(bound), sort!(added), sort!(subtracted))
end

"""
Read a residual expression (`A`, `A - P`, `S1 + S2 - P1 - P3`, ...) into `added` (its
positive names) and `subtracted` (its negative ones); `positive` is the sign of `e`. A
residual adds substrates and subtracts products.
"""
function _read_residual!(added, subtracted, e, positive::Bool, names)
    if e isa Symbol
        _is_declared(e, names) ||
            error("@enzyme_mechanism: residual entry `$e` is not a declared metabolite. " *
                  "Declared: $(_declared_names(names)).")
        e in (positive ? names.substrates : names.products) ||
            error("@enzyme_mechanism: a residual adds substrates and subtracts products; " *
                  "got `$e`.")
        push!(positive ? added : subtracted, e)
    elseif e isa Expr && e.head === :call && e.args[1] === :+ && length(e.args) >= 3
        for a in e.args[2:end]
            _read_residual!(added, subtracted, a, positive, names)
        end
    elseif e isa Expr && e.head === :call && e.args[1] === :- && length(e.args) in (2, 3)
        length(e.args) == 3 &&
            _read_residual!(added, subtracted, e.args[2], positive, names)
        _read_residual!(added, subtracted, e.args[end], !positive, names)
    else
        error("@enzyme_mechanism: invalid residual expression: $e")
    end
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
from its nodes and steps, the allosteric state of each group label, its oligomeric state
and its regulatory sites. The names under `modifiers:` are checked here, since the
constructors take every name bound as a modifier as one: each is listed once, and some
step binds each as a modifier.
"""
function _construct(raw::NamedTuple)
    allunique(raw.modifiers) || error("Mechanism: duplicate modifier names")
    m = Mechanisms.Mechanism(
        Mechanisms.CatalyticScheme(raw.substrates, raw.products, raw.nodes, raw.steps),
        Mechanisms.Allostery(raw.allosteric_states, raw.oligomeric_state,
                             [Mechanisms.RegulatorySite(s.ligands, s.multiplicity,
                                                        s.allosteric_states)
                              for s in raw.regulatory_sites]))
    unused = findfirst(∉(m.catalytic.modifiers), raw.modifiers)
    unused === nothing ||
        error("Mechanism: modifier $(raw.modifiers[unused]) is declared, but no step " *
              "binds it")
    m
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
