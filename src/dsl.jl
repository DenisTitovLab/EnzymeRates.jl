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
    @enzyme_reaction begin
        substrates: S[C]
        products:   P[C]
        regulators: I
    end

Create an `EnzymeReaction` from a DSL block. Species atoms use chemical
formula bracket syntax: `S[C6H12O6]`. Bare symbols (no brackets) are allowed
when all metabolites omit atoms. Regulators are plain symbol names.

Multi-species lines use comma separation:
    substrates: S[C6H12O6], ATP[C10H16N5O13P3]
    regulators: I, A
"""
macro enzyme_reaction(block)
    parsed = _parse_labeled_block(block,
        Set([:substrates, :products, :regulators,
             :dead_end_inhibitors, :allosteric_regulators]))
    haskey(parsed, :substrates) || error("substrates not specified")
    haskey(parsed, :products) || error("products not specified")
    subs = parsed[:substrates]
    prods = parsed[:products]
    regs = Expr(:tuple)
    for (label, role_sym) in [
        (:regulators, :unknown),
        (:dead_end_inhibitors, :dead_end),
        (:allosteric_regulators, :allosteric),
    ]
        if haskey(parsed, label)
            syms = _regulator_tuple_to_symbols(parsed[label])
            for s in syms.args
                push!(regs.args,
                    Expr(:tuple, s, QuoteNode(role_sym)))
            end
        end
    end
    return esc(:(EnzymeReaction($subs, $prods, $regs)))
end

"""Convert a parsed species tuple for regulators into a plain Symbol tuple.
Accepts both bare Symbols (:I) and bracket syntax (I[C5]) — extracts just the name."""
function _regulator_tuple_to_symbols(species_tuple::Expr)
    result = Expr(:tuple)
    for arg in species_tuple.args
        if arg isa Expr && arg.head == :tuple
            # (QuoteNode(:I), atoms...) → QuoteNode(:I)
            push!(result.args, arg.args[1])
        else
            push!(result.args, arg)
        end
    end
    result
end

function _parse_species_block(block)
    valid = Set([:substrates, :products, :regulators, :enzymes])
    parsed = _parse_labeled_block(block, valid)

    haskey(parsed, :substrates) || error("substrates not specified in species block")
    haskey(parsed, :products) || error("products not specified in species block")
    haskey(parsed, :enzymes) || error("enzymes not specified in species block")
    subs = parsed[:substrates]
    prods = parsed[:products]
    regs = get(parsed, :regulators, nothing)
    if regs === nothing
        regs = Expr(:tuple)
    else
        regs = _regulator_tuple_to_symbols(regs)
    end
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
            regulators: I
            enzymes:    E, ES[C]
        end
        steps: begin
            [E, S] <--> [ES]
            [ES] <--> [E, P]
        end
    end

Create an `EnzymeMechanism` from explicit species and step definitions.
Species atoms use chemical formula bracket syntax: `S[C6H12O6]`. Bare symbols
(no brackets) are allowed when all metabolites omit atoms. Regulators are
plain symbol names (no atoms).
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
    # Detect new OligomericEnzymeMechanism syntax (metabolites: or site(...):)
    for arg in block.args
        arg isa LineNumberNode && continue
        if _is_oligomeric_label(arg)
            return esc(_parse_oligomeric_mechanism(block))
        end
    end
    return esc(_parse_enzyme_mechanism(block))
end

"""Return true if an @enzyme_mechanism block arg is part of new oligomeric syntax."""
function _is_oligomeric_label(arg)
    # metabolites: ... (single or tuple form)
    if arg isa Expr && arg.head == :tuple
        inner = arg.args[1]
        return inner isa Expr && inner.head == :call && inner.args[1] == :(:) &&
               inner.args[2] == :metabolites
    end
    if arg isa Expr && arg.head == :call && arg.args[1] == :(:)
        label = arg.args[2]
        return label == :metabolites || label == :conformations ||
               (label isa Expr && label.head == :call && label.args[1] == :site)
    end
    false
end

"""Parse the original EnzymeMechanism DSL (species:/steps:/constraints: blocks)."""
function _parse_enzyme_mechanism(block)
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
    reactions, eq_steps = _parse_steps_block(steps_block)

    if constraints_block !== nothing
        constraints = _parse_constraints_block(constraints_block)
        return :(EnzymeMechanism($species_tuple, $reactions, $eq_steps, $constraints))
    end
    return :(EnzymeMechanism($species_tuple, $reactions, $eq_steps))
end

"""Parse steps: begin ... end into (reactions_expr, eq_steps_expr)."""
function _parse_steps_block(steps_block)
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
    reactions, eq_steps
end

"""Parse constraints: begin ... end into a constraints_expr tuple."""
function _parse_constraints_block(constraints_block)
    constraints = Expr(:tuple)
    for arg in constraints_block.args
        arg isa LineNumberNode && continue
        # Handle semicolon-separated constraints: K5=K3; K6=K1 (parsed as tuple)
        if arg isa Expr && arg.head == :tuple
            for a in arg.args
                a isa Expr && a.head == :(=) || continue
                _push_constraint!(constraints, a)
            end
        elseif arg isa Expr && arg.head == :(=)
            _push_constraint!(constraints, arg)
        else
            error("Each constraint must be an assignment: target = rhs_expr, got $arg")
        end
    end
    constraints
end

function _push_constraint!(constraints, arg)
    target = arg.args[1]
    target isa Symbol || error("Constraint target must be a symbol, got $target")
    coeff, factors = _parse_constraint_rhs(arg.args[2])
    push!(constraints.args, Expr(:tuple, QuoteNode(target), coeff, factors))
end

"""
Parse the OligomericEnzymeMechanism DSL syntax.
Handles: metabolites:, conformations:, site(:catalytic, N):, site(:regulatory, N):
"""
function _parse_oligomeric_mechanism(block)
    met_names = nothing   # vector of Symbol (metabolite names only, no atoms)
    nconf = 1
    catalytic_n = nothing
    catalytic_block = nothing
    reg_sites = Any[]     # vector of (ligand_syms, n_reg)

    for arg in block.args
        arg isa LineNumberNode && continue

        if arg isa Expr && arg.head == :tuple
            # metabolites: S[C], P[C] → Expr(:tuple, :(metabolites: S[C]), :(P[C]))
            inner = arg.args[1]
            if inner isa Expr && inner.head == :call && inner.args[1] == :(:) &&
               inner.args[2] == :metabolites
                met_names = [_met_sym(inner.args[3])]
                for i in 2:length(arg.args)
                    push!(met_names, _met_sym(arg.args[i]))
                end
            else
                error("Unexpected tuple in @enzyme_mechanism: $arg")
            end
        elseif arg isa Expr && arg.head == :call && arg.args[1] == :(:)
            label = arg.args[2]
            value = arg.args[3]

            if label == :metabolites
                met_names = [_met_sym(value)]
            elseif label == :conformations
                nconf = value  # integer literal
            elseif label isa Expr && label.head == :call && label.args[1] == :site
                site_kind = label.args[2]   # QuoteNode(:catalytic) or QuoteNode(:regulatory)
                site_n = label.args[3]      # integer literal
                if site_kind == QuoteNode(:catalytic)
                    catalytic_n = site_n
                    catalytic_block = value
                elseif site_kind == QuoteNode(:regulatory)
                    ligs = _parse_reg_ligands_block(value)
                    push!(reg_sites, (ligs, site_n))
                else
                    error("Unknown site kind: $site_kind")
                end
            else
                error("Unknown @enzyme_mechanism block label: $label")
            end
        else
            error("Unexpected expression in @enzyme_mechanism: $arg")
        end
    end

    met_names === nothing && error("metabolites: block not specified")
    catalytic_block === nothing && error("site(:catalytic, N): block not specified")

    # Build metabolites type parameter tuple (just Symbol names, no atoms)
    mets_tuple = Expr(:tuple, QuoteNode.(met_names)...)

    # Parse catalytic site block: species:, steps:, constraints:
    species_tuple, cat_steps_block, cat_constraints_block =
        _parse_catalytic_block(catalytic_block)

    # Parse reactions and eq_steps
    reactions, eq_steps = _parse_steps_block(cat_steps_block)

    # Build inner EnzymeMechanism expression
    if cat_constraints_block !== nothing
        constraints = _parse_constraints_block(cat_constraints_block)
        cm_expr = :(EnzymeMechanism($species_tuple, $reactions, $eq_steps, $constraints))
    else
        cm_expr = :(EnzymeMechanism($species_tuple, $reactions, $eq_steps))
    end

    # Build RegSites type parameter tuple: ((ligand_syms...,), n_reg) pairs
    reg_sites_elems = Any[]
    for (ligs, n_reg) in reg_sites
        ligs_tuple = Expr(:tuple, (QuoteNode(l) for l in ligs)...)
        push!(reg_sites_elems, Expr(:tuple, ligs_tuple, n_reg))
    end
    reg_sites_expr = Expr(:tuple, reg_sites_elems...)

    # Emit: let _cm = EnzymeMechanism(...)
    #           OligomericEnzymeMechanism{mets, typeof(_cm), CatN, RegSites, NConf}()
    #       end
    :(let _cm = $cm_expr
        OligomericEnzymeMechanism{$mets_tuple, typeof(_cm), $catalytic_n, $reg_sites_expr, $nconf}()
    end)
end

"""Parse regulatory site block: begin ligands: L1, L2 end → vector of ligand symbols."""
function _parse_reg_ligands_block(block)
    ligs = Symbol[]
    for arg in block.args
        arg isa LineNumberNode && continue
        if arg isa Expr && arg.head == :call && arg.args[1] == :(:) && arg.args[2] == :ligands
            push!(ligs, arg.args[3])
        elseif arg isa Expr && arg.head == :tuple
            inner = arg.args[1]
            if inner isa Expr && inner.head == :call && inner.args[1] == :(:) &&
               inner.args[2] == :ligands
                push!(ligs, inner.args[3])
                for i in 2:length(arg.args)
                    push!(ligs, arg.args[i])
                end
            else
                error("Expected ligands: L1, L2, ... in regulatory site block, got $arg")
            end
        else
            error("Expected ligands: in regulatory site block, got $arg")
        end
    end
    ligs
end

"""Extract the metabolite Symbol from a raw DSL expression (Symbol or S[C] ref)."""
function _met_sym(expr)
    expr isa Symbol && return expr
    expr isa Expr && expr.head == :ref && return expr.args[1]
    error("Cannot extract metabolite name from: $expr")
end

"""
Parse catalytic site block.
Returns (species_tuple_expr, steps_block, constraints_block_or_nothing).
The species_tuple_expr is produced by `_parse_species_block` and contains
(subs, prods, regs, enzymes) sub-expressions for the inner EnzymeMechanism.
"""
function _parse_catalytic_block(block)
    species_expr = nothing
    steps_block = nothing
    constraints_block = nothing

    for arg in block.args
        arg isa LineNumberNode && continue
        arg isa Expr && arg.head == :call && arg.args[1] == :(:) ||
            error("Unexpected expression in catalytic site block: $arg")
        label = arg.args[2]
        value = arg.args[3]
        if label == :species
            species_expr = _parse_species_block(value)
        elseif label == :steps
            steps_block = value
        elseif label == :constraints
            constraints_block = value
        else
            error("Unknown label in catalytic site block: $label")
        end
    end

    species_expr === nothing && error("species: not specified in catalytic site block")
    steps_block === nothing && error("steps: not specified in catalytic site block")

    species_expr, steps_block, constraints_block
end
