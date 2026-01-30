"""
Parse a species expression like `S(C=6, H=12, O=6)` into a Species constructor Expr.
"""
function _parse_species_expr(expr)
    if expr isa Symbol
        return :(Species($(QuoteNode(expr)), metabolite, Dict{Symbol,Int}()))
    elseif expr isa Expr && expr.head == :call
        name = expr.args[1]
        atoms = Expr(:call, :(Dict{Symbol,Int}))
        for i in 2:length(expr.args)
            arg = expr.args[i]
            if arg isa Expr && arg.head == :kw
                push!(atoms.args, :($(Expr(:call, :(=>), QuoteNode(arg.args[1]), arg.args[2]))))
            end
        end
        return :(Species($(QuoteNode(name)), metabolite, $atoms))
    else
        error("Cannot parse species expression: $expr")
    end
end

"""
    @enzyme_reaction begin
        substrates: S(C=1)
        products:   P(C=1)
        regulators: I(C=5)
    end

Create a `ReactionSpec` from a DSL block.

Multi-species lines use comma separation:
    substrates: S(C=6, H=12), ATP(C=10)
which Julia parses as a tuple `(substrates: S(C=6,H=12), ATP(C=10))`.
"""
macro enzyme_reaction(block)
    substrates = nothing
    products = nothing
    regulators = nothing

    # Collect lines, handling both plain exprs and tuples
    items = Any[]
    for arg in block.args
        arg isa LineNumberNode && continue
        push!(items, arg)
    end

    # Process items: a tuple means multi-species line, a call means single-species line
    i = 1
    while i <= length(items)
        item = items[i]
        if item isa Expr && item.head == :tuple
            # First element is `label: Species(...)`, rest are additional species
            first_elem = item.args[1]
            label, first_species = _parse_label_species(first_elem)
            species_list = Expr(:vect, first_species)
            for j in 2:length(item.args)
                push!(species_list.args, _parse_species_expr(item.args[j]))
            end
            if label == :substrates
                substrates = species_list
            elseif label == :products
                products = species_list
            elseif label == :regulators
                regulators = species_list
            end
        elseif item isa Expr && item.head == :call && item.args[1] == :(:)
            label = item.args[2]
            species_expr = item.args[3]
            species_list = Expr(:vect, _parse_species_expr(species_expr))
            if label == :substrates
                substrates = species_list
            elseif label == :products
                products = species_list
            elseif label == :regulators
                regulators = species_list
            end
        end
        i += 1
    end

    substrates === nothing && error("substrates not specified")
    products === nothing && error("products not specified")
    regulators === nothing && (regulators = :(Species[]))

    return esc(:(ReactionSpec($substrates, $products, $regulators)))
end

"""
Parse a species in the @mechanism DSL.
Bare symbol → enzyme form; Name(K=V, ...) → metabolite with atoms.
"""
function _parse_mechanism_species(expr)
    if expr isa Symbol
        return :(Species($(QuoteNode(expr)), enzyme))
    elseif expr isa Expr && expr.head == :call
        name = expr.args[1]
        atoms = Expr(:call, :(Dict{Symbol,Int}))
        for i in 2:length(expr.args)
            arg = expr.args[i]
            if arg isa Expr && arg.head == :kw
                push!(atoms.args, Expr(:call, :(=>), QuoteNode(arg.args[1]), arg.args[2]))
            end
        end
        return :(Species($(QuoteNode(name)), metabolite, $atoms))
    else
        error("Cannot parse mechanism species: $expr")
    end
end

function _parse_side(expr)
    # expr is a :vect expression [A, B, ...]
    if expr isa Expr && expr.head == :vect
        return Expr(:vect, [_parse_mechanism_species(a) for a in expr.args]...)
    end
    error("Expected [...] on each side of -->, got $expr")
end

"""
    @mechanism begin
        [E, S(C=1)] --> [ES]
        [ES] --> [E, P(C=1)]
    end

Create an `EnzymeMechanism` from step definitions.

Bare symbols (e.g. `E`, `ES`) are enzyme forms.
Symbols with keyword arguments (e.g. `S(C=1)`) are metabolites with atomic composition.
"""
macro mechanism(block)
    steps = Expr[]
    for arg in block.args
        arg isa LineNumberNode && continue
        if arg isa Expr && arg.head == :(-->)
            lhs = _parse_side(arg.args[1])
            rhs = _parse_side(arg.args[2])
            push!(steps, :($lhs => $rhs))
        else
            error("Expected [lhs] --> [rhs], got $arg")
        end
    end
    return esc(:(EnzymeMechanism([$(steps...)])))
end

function _parse_label_species(expr)
    # expr is `label: Species(...)` which is Expr(:call, :(:), :label, species_expr)
    if expr isa Expr && expr.head == :call && expr.args[1] == :(:)
        label = expr.args[2]
        species = _parse_species_expr(expr.args[3])
        return label, species
    end
    error("Expected label: species, got $expr")
end
