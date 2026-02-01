function _parse_species_tuple_expr(expr)
    if expr isa Symbol
        atoms = Expr(:tuple)
        return Expr(:tuple, QuoteNode(expr), atoms)
    elseif expr isa Expr && expr.head == :call
        name = expr.args[1]
        atoms = Expr(:tuple)
        for i in 2:length(expr.args)
            arg = expr.args[i]
            if arg isa Expr && arg.head == :kw
                push!(atoms.args, Expr(:tuple, QuoteNode(arg.args[1]), arg.args[2]))
            end
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

"""
    @enzyme_reaction begin
        substrates: S(C=1)
        products:   P(C=1)
        regulators: I(C=5)
    end

Create an `EnzymeReaction` from a DSL block.

Multi-species lines use comma separation:
    substrates: S(C=6, H=12), ATP(C=10)
"""
macro enzyme_reaction(block)
    subs = nothing
    prods = nothing
    regs = nothing

    items = Any[]
    for arg in block.args
        arg isa LineNumberNode && continue
        push!(items, arg)
    end

    i = 1
    while i <= length(items)
        item = items[i]
        if item isa Expr && item.head == :tuple
            first_elem = item.args[1]
            label, first_species = _parse_label_species_tuple(first_elem)
            species_list = Expr(:tuple, first_species)
            for j in 2:length(item.args)
                push!(species_list.args, _parse_species_tuple_expr(item.args[j]))
            end
            if label == :substrates
                subs = species_list
            elseif label == :products
                prods = species_list
            elseif label == :regulators
                regs = species_list
            end
        elseif item isa Expr && item.head == :call && item.args[1] == :(:)
            label = item.args[2]
            species_expr = item.args[3]
            species_list = Expr(:tuple, _parse_species_tuple_expr(species_expr))
            if label == :substrates
                subs = species_list
            elseif label == :products
                prods = species_list
            elseif label == :regulators
                regs = species_list
            end
        end
        i += 1
    end

    subs === nothing && error("substrates not specified")
    prods === nothing && error("products not specified")
    regs === nothing && (regs = Expr(:tuple))

    return esc(:(EnzymeReaction{$subs, $prods, $regs}()))
end

function _parse_species_block(block)
    subs = nothing
    prods = nothing
    regs = nothing
    enzs = nothing

    items = Any[]
    for arg in block.args
        arg isa LineNumberNode && continue
        push!(items, arg)
    end

    for item in items
        if item isa Expr && item.head == :tuple
            first_elem = item.args[1]
            label, first_species = _parse_label_species_tuple(first_elem)
            species_list = Expr(:tuple, first_species)
            for j in 2:length(item.args)
                push!(species_list.args, _parse_species_tuple_expr(item.args[j]))
            end
            if label == :substrates
                subs = species_list
            elseif label == :products
                prods = species_list
            elseif label == :regulators
                regs = species_list
            elseif label == :enzymes
                enzs = species_list
            end
        elseif item isa Expr && item.head == :call && item.args[1] == :(:)
            label = item.args[2]
            species_expr = item.args[3]
            species_list = Expr(:tuple, _parse_species_tuple_expr(species_expr))
            if label == :substrates
                subs = species_list
            elseif label == :products
                prods = species_list
            elseif label == :regulators
                regs = species_list
            elseif label == :enzymes
                enzs = species_list
            end
        end
    end

    subs === nothing && error("substrates not specified in species block")
    prods === nothing && error("products not specified in species block")
    enzs === nothing && error("enzymes not specified in species block")
    regs === nothing && (regs = Expr(:tuple))

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
    error("Expected [...] on each side of -->, got $expr")
end

"""
    @mechanism begin
        species: begin
            substrates: S(C=1)
            products:   P(C=1)
            regulators: I(C=1)
            enzymes:    E(), ES(C=1)
        end
        steps: begin
            [E, S] --> [ES]
            [ES] --> [E, P]
        end
    end

Create an `EnzymeMechanism` from explicit species and step definitions.
"""
macro mechanism(block)
    species_block = nothing
    steps_block = nothing

    for arg in block.args
        arg isa LineNumberNode && continue
        if arg isa Expr && arg.head == :call && arg.args[1] == :(:)
            label = arg.args[2]
            value = arg.args[3]
            if label == :species
                species_block = value
            elseif label == :steps
                steps_block = value
            else
                error("Unknown @mechanism block label: $label")
            end
        else
            error("Expected species: begin ... end and steps: begin ... end blocks")
        end
    end

    species_block === nothing && error("species block not specified")
    steps_block === nothing && error("steps block not specified")

    species_tuple = _parse_species_block(species_block)
    reactions = Expr(:tuple)
    for arg in steps_block.args
        arg isa LineNumberNode && continue
        if arg isa Expr && arg.head == :(-->)
            lhs = _parse_step_side_symbols(arg.args[1])
            rhs = _parse_step_side_symbols(arg.args[2])
            push!(reactions.args, Expr(:tuple, lhs, rhs))
        else
            error("Expected [lhs] --> [rhs], got $arg")
        end
    end

    return esc(:(EnzymeMechanism($species_tuple, $reactions)))
end
