"""
    _parse_chemical_formula(s::String)

Parse a chemical formula string like `"C6H12O6"` into a tuple expression of
`(element, count)` pairs. Elements are identified by an uppercase letter
optionally followed by lowercase letters, then an optional integer count
(defaults to 1).
"""
function _parse_chemical_formula(s::String)
    atoms = Expr(:tuple)
    matched_len = 0
    for m in eachmatch(r"([A-Z][a-z]*)(\d*)", s)
        elem = Symbol(m.captures[1]::SubString)
        cap2 = m.captures[2]::SubString
        count = isempty(cap2) ? 1 : parse(Int, cap2)
        push!(atoms.args, Expr(:tuple, QuoteNode(elem), count))
        matched_len += length(m.match)
    end
    matched_len == length(s) || error("Invalid chemical formula: \"$s\"")
    atoms
end

function _parse_species_tuple_expr(expr)
    if expr isa Symbol
        atoms = Expr(:tuple)
        return Expr(:tuple, QuoteNode(expr), atoms)
    elseif expr isa Expr && expr.head == :ref
        # S[C6H12O6] syntax: expr.args[1] is the name, expr.args[2] is the formula symbol
        name = expr.args[1]
        formula = string(expr.args[2])
        atoms = _parse_chemical_formula(formula)
        return Expr(:tuple, QuoteNode(name), atoms)
    else
        error("Cannot parse species definition: $expr")
    end
end

function _parse_label_species_tuple(expr)
    if expr isa Expr && expr.head == :call && expr.args[1] == :(:)
        label = expr.args[2]
        species = _parse_species_tuple_expr(expr.args[3])
        return label, species
    end
    error("Expected label: species, got $expr")
end

function _parse_labeled_block(block, valid_labels::Set{Symbol})
    result = Dict{Symbol, Expr}()
    for arg in block.args
        arg isa LineNumberNode && continue
        if arg isa Expr && arg.head == :tuple
            label, first_species = _parse_label_species_tuple(arg.args[1])
            label in valid_labels || error("Unknown label: $label")
            species_list = Expr(:tuple, first_species)
            for j in 2:length(arg.args)
                push!(species_list.args, _parse_species_tuple_expr(arg.args[j]))
            end
            result[label] = species_list
        elseif arg isa Expr && arg.head == :call && arg.args[1] == :(:)
            label = arg.args[2]
            label in valid_labels || error("Unknown label: $label")
            result[label] = Expr(:tuple, _parse_species_tuple_expr(arg.args[3]))
        end
    end
    result
end

"""
Extract the label symbol from a DSL expression, or `nothing` if not a labeled line.

Handles both standard `label: value` and the unparenthesized `max_sites: A => 2`
where Julia parses it as `(max_sites:A) => 2`.
"""
function _get_label(arg)
    # Standard: label: value
    if arg isa Expr && arg.head == :call && arg.args[1] == :(:)
        return arg.args[2]
    elseif arg isa Expr && arg.head == :tuple && length(arg.args) > 0
        first_arg = arg.args[1]
        # Standard multi: label: X, Y  →  tuple((: label X), Y)
        if first_arg isa Expr && first_arg.head == :call && first_arg.args[1] == :(:)
            return first_arg.args[2]
        end
        # Unparenthesized multi: max_sites: A => 2, B => 3
        # →  tuple((=> (: max_sites A) 2), (=> B 3))
        if first_arg isa Expr && first_arg.head == :call && first_arg.args[1] == :(=>)
            lhs = first_arg.args[2]
            if lhs isa Expr && lhs.head == :call && lhs.args[1] == :(:)
                return lhs.args[2]
            end
        end
    end
    # Unparenthesized single: max_sites: S => 2  →  (=> (: max_sites S) 2)
    if arg isa Expr && arg.head == :call && arg.args[1] == :(=>)
        lhs = arg.args[2]
        if lhs isa Expr && lhs.head == :call && lhs.args[1] == :(:)
            return lhs.args[2]
        end
    end
    nothing
end

"""Parse a `Name => Int` pair expression into the max_sites dict."""
function _parse_pair_expr!(result::Dict{Symbol,Int}, expr)
    if expr isa Expr && expr.head == :call && expr.args[1] == :(=>)
        name = expr.args[2]
        name isa Symbol || error("max_sites key must be a symbol, got $name")
        count = expr.args[3]
        count isa Integer || error("max_sites value must be an integer, got $count")
        count >= 1 || error("max_sites value must be ≥ 1, got $count for $name")
        result[name] = count
    else
        error("max_sites entry must be of the form Name => Int, got $expr")
    end
end

"""
Parse a max_sites block expression into the result dict.

Handles both parenthesized (`max_sites: (A => 2)`) and unparenthesized
(`max_sites: A => 2`) forms, including comma-separated multi-entry variants.
"""
function _parse_max_sites_block!(result::Dict{Symbol,Int}, arg)
    if arg isa Expr && arg.head == :call && arg.args[1] == :(:)
        # Parenthesized: max_sites: (S => 2) or max_sites: (S => 2, B => 3)
        value = arg.args[3]
        if value isa Expr && value.head == :call && value.args[1] == :(=>)
            _parse_pair_expr!(result, value)
        elseif value isa Expr && value.head == :tuple
            for entry in value.args
                _parse_pair_expr!(result, entry)
            end
        end
    elseif arg isa Expr && arg.head == :call && arg.args[1] == :(=>)
        # Unparenthesized single: max_sites: S => 2  →  (=> (: max_sites S) 2)
        lhs = arg.args[2]
        name = lhs.args[3]
        count = arg.args[3]
        name isa Symbol || error("max_sites key must be a symbol, got $name")
        count isa Integer || error("max_sites value must be an integer, got $count")
        count >= 1 || error("max_sites value must be ≥ 1, got $count for $name")
        result[name] = count
    elseif arg isa Expr && arg.head == :tuple
        # Unparenthesized multi: max_sites: A => 2, B => 3
        # →  tuple((=> (: max_sites A) 2), (=> B 3))
        first_arg = arg.args[1]
        if first_arg isa Expr && first_arg.head == :call && first_arg.args[1] == :(=>)
            lhs = first_arg.args[2]
            result[lhs.args[3]] = first_arg.args[3]
        end
        for j in 2:length(arg.args)
            _parse_pair_expr!(result, arg.args[j])
        end
    end
end

"""Add max_sites (3rd element) to each species tuple expression."""
function _apply_max_sites_expr(species_tuple_expr, max_sites_map)
    new_expr = Expr(:tuple)
    for arg in species_tuple_expr.args
        name = arg.args[1].value  # QuoteNode value
        ms = get(max_sites_map, name, 1)
        push!(new_expr.args, Expr(:tuple, arg.args..., ms))
    end
    new_expr
end

"""
    @enzyme_reaction begin
        substrates: S[C]
        products:   P[C]
        regulators: I[C5]
        max_sites:  S => 2
    end

Create an `EnzymeReaction` from a DSL block. Species atoms use chemical
formula bracket syntax: `S[C6H12O6]`. Bare symbols (no brackets) are allowed
when all metabolites omit atoms.

Multi-species lines use comma separation:
    substrates: S[C6H12O6], ATP[C10H16N5O13P3]

The `max_sites` label maps metabolite names to occupancy values. Any metabolite
not listed defaults to 1. Multiple entries use comma separation:
    max_sites: A => 2, B => 3
"""
macro enzyme_reaction(block)
    # Extract max_sites entries and build a cleaned block for species parsing
    max_sites_map = Dict{Symbol, Int}()
    filtered_args = Any[]
    for arg in block.args
        arg isa LineNumberNode && (push!(filtered_args, arg); continue)
        if _get_label(arg) == :max_sites
            _parse_max_sites_block!(max_sites_map, arg)
        else
            push!(filtered_args, arg)
        end
    end

    species_block = Expr(:block, filtered_args...)
    parsed = _parse_labeled_block(species_block, Set([:substrates, :products, :regulators]))

    haskey(parsed, :substrates) || error("substrates not specified")
    haskey(parsed, :products) || error("products not specified")
    subs = _apply_max_sites_expr(parsed[:substrates], max_sites_map)
    prods = _apply_max_sites_expr(parsed[:products], max_sites_map)
    regs = _apply_max_sites_expr(get(parsed, :regulators, Expr(:tuple)), max_sites_map)

    return esc(:(EnzymeReaction($subs, $prods, $regs)))
end

function _parse_species_block(block)
    parsed = _parse_labeled_block(block, Set([:substrates, :products, :regulators, :enzymes]))

    haskey(parsed, :substrates) || error("substrates not specified in species block")
    haskey(parsed, :products) || error("products not specified in species block")
    haskey(parsed, :enzymes) || error("enzymes not specified in species block")
    subs = parsed[:substrates]
    prods = parsed[:products]
    regs = get(parsed, :regulators, Expr(:tuple))
    enzs = parsed[:enzymes]

    return Expr(:tuple, subs, prods, regs, enzs)
end

function _parse_step_side_symbols(expr)
    if expr isa Expr && expr.head == :vect
        syms = Expr(:tuple)
        for a in expr.args
            a isa Symbol || error("Step sides must be symbols; define atoms in species block")
            push!(syms.args, QuoteNode(a))
        end
        return syms
    end
    error("Expected [...] on each side of <-->, got $expr")
end

"""
    @enzyme_mechanism begin
        species: begin
            substrates: S[C]
            products:   P[C]
            regulators: I[C]
            enzymes:    E, ES[C]
        end
        steps: begin
            [E, S] <--> [ES]
            [ES] <--> [E, P]
        end
    end

Create an `EnzymeMechanism` from explicit species and step definitions.
Species atoms use chemical formula bracket syntax: `S[C6H12O6]`. Bare symbols
(no brackets) are allowed when all metabolites omit atoms.
Steps use the `<-->` arrow.
"""
# Decompose a constraint RHS expression into (coeff::Int, factors_tuple_expr).
# Handles: symbols, `a * b`, `a / b`, integer literals, and combinations.
function _parse_constraint_rhs(expr)
    factors = Dict{Symbol,Int}()
    coeff = Ref(1)
    _walk_rhs!(expr, factors, coeff, 1)
    factors_expr = Expr(:tuple)
    for (sym, exp) in factors
        push!(factors_expr.args, Expr(:tuple, QuoteNode(sym), exp))
    end
    coeff[], factors_expr
end

function _walk_rhs!(expr, factors::Dict{Symbol,Int}, coeff::Ref{Int}, sign::Int)
    if expr isa Symbol
        factors[expr] = get(factors, expr, 0) + sign
    elseif expr isa Integer
        sign > 0 || error("Integer divisors not supported in constraints")
        coeff[] *= expr
    elseif expr isa Expr && expr.head == :call
        op = expr.args[1]
        if op == :*; for i in 2:length(expr.args); _walk_rhs!(expr.args[i], factors, coeff, sign); end
        elseif op == :/; _walk_rhs!(expr.args[2], factors, coeff, sign); _walk_rhs!(expr.args[3], factors, coeff, -sign)
        else error("Unsupported operator in constraint: $op"); end
    else error("Unsupported constraint expression: $expr"); end
end

macro enzyme_mechanism(block)
    species_block = nothing
    steps_block = nothing
    constraints_block = nothing

    for arg in block.args
        arg isa LineNumberNode && continue
        if arg isa Expr && arg.head == :call && arg.args[1] == :(:)
            label = arg.args[2]
            value = arg.args[3]
            if label == :species
                species_block = value
            elseif label == :steps
                steps_block = value
            elseif label == :constraints
                constraints_block = value
            else
                error("Unknown @mechanism block label: $label")
            end
        else
            error("Expected species: begin ... end, steps: begin ... end, and optional constraints: begin ... end blocks")
        end
    end

    species_block === nothing && error("species block not specified")
    steps_block === nothing && error("steps block not specified")

    species_tuple = _parse_species_block(species_block)
    reactions = Expr(:tuple)
    eq_steps = Expr(:tuple)
    for arg in steps_block.args
        arg isa LineNumberNode && continue
        is_re = false
        if arg isa Expr && arg.head == :call && arg.args[1] == :(<-->)
            is_re = false
        elseif arg isa Expr && arg.head == :call && arg.args[1] == :⇌
            is_re = true
        else
            error("Expected [lhs] <--> [rhs] or [lhs] ⇌ [rhs], got $arg")
        end
        lhs = _parse_step_side_symbols(arg.args[2])
        rhs = _parse_step_side_symbols(arg.args[3])
        push!(reactions.args, Expr(:tuple, lhs, rhs))
        push!(eq_steps.args, is_re)
    end

    # Parse constraints block if present
    if constraints_block !== nothing
        constraints = Expr(:tuple)
        for arg in constraints_block.args
            arg isa LineNumberNode && continue
            if !(arg isa Expr && arg.head == :(=))
                error("Each constraint must be an assignment: target = rhs_expr, got $arg")
            end
            target = arg.args[1]
            target isa Symbol || error("Constraint target must be a symbol, got $target")
            coeff, factors = _parse_constraint_rhs(arg.args[2])
            push!(constraints.args, Expr(:tuple, QuoteNode(target), coeff, factors))
        end
        return esc(:(EnzymeMechanism($species_tuple, $reactions, $eq_steps, $constraints)))
    end

    return esc(:(EnzymeMechanism($species_tuple, $reactions, $eq_steps)))
end
