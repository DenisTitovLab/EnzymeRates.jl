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
        # S[C6H12O6] or S[C6H12O6, 2] syntax
        name = expr.args[1]
        formula = string(expr.args[2])
        atoms = _parse_chemical_formula(formula)
        if length(expr.args) >= 3
            return Expr(:tuple, QuoteNode(name), atoms, expr.args[3])
        end
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
    @enzyme_reaction begin
        substrates: S[C]
        products:   P[C]
        regulators: I[C5]
    end

Create an `EnzymeReaction` from a DSL block. Species atoms use chemical
formula bracket syntax: `S[C6H12O6]`. Bare symbols (no brackets) are allowed
when all metabolites omit atoms.

Multi-species lines use comma separation:
    substrates: S[C6H12O6], ATP[C10H16N5O13P3]

Max binding sites use an optional integer in brackets: `S[C, 2]` (default 1).
"""
macro enzyme_reaction(block)
    parsed = _parse_labeled_block(block, Set([:substrates, :products, :regulators]))
    haskey(parsed, :substrates) || error("substrates not specified")
    haskey(parsed, :products) || error("products not specified")
    subs = parsed[:substrates]
    prods = parsed[:products]
    regs = get(parsed, :regulators, Expr(:tuple))
    return esc(:(EnzymeReaction($subs, $prods, $regs)))
end

function _parse_species_block(block)
    valid = Set([:substrates, :products, :regulators, :enzymes])
    parsed = _parse_labeled_block(block, valid)

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
            a isa Symbol || error(
                "Step sides must be symbols; " *
                "define atoms in species block"
            )
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
        if op == :*
            for i in 2:length(expr.args)
                _walk_rhs!(expr.args[i], factors, coeff, sign)
            end
        elseif op == :/
            _walk_rhs!(expr.args[2], factors, coeff, sign)
            _walk_rhs!(expr.args[3], factors, coeff, -sign)
        else
            error("Unsupported operator in constraint: $op")
        end
    else
        error("Unsupported constraint expression: $expr")
    end
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
            error(
                "Expected species: begin ... end, " *
                "steps: begin ... end, and optional " *
                "constraints: begin ... end blocks"
            )
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
