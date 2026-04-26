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
        # Parse each ref arg individually and merge
        # (handles A[C,N] where Julia parses as
        # ref with args [:A, :C, :N]).
        # Each arg is parsed as a chemical formula
        # separately to avoid ambiguity (e.g.,
        # A[CO2,H] must not become "CO2H" → cobalt).
        atoms = Expr(:tuple)
        for arg in expr.args[2:end]
            parsed = _parse_chemical_formula(string(arg))
            for atom in parsed.args
                push!(atoms.args, atom)
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

function _parse_labeled_block(
    block, valid_labels::Set{Symbol},
    scalar_labels::Set{Symbol}=Set{Symbol}(),
)
    result = Dict{Symbol, Any}()
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
            if label in scalar_labels
                result[label] = arg.args[3]
            else
                result[label] = Expr(:tuple, _parse_species_tuple_expr(arg.args[3]))
            end
        end
    end
    result
end

"""
    @enzyme_reaction begin
        substrates: S[C]
        products:   P[C]
        regulators: I
        oligomeric_state: 4
    end

Create an `EnzymeReaction` from a DSL block. Species atoms use chemical
formula bracket syntax: `S[C6H12O6]`. Bare symbols (no brackets) are allowed
when all metabolites omit atoms. Regulators are plain symbol names.
`oligomeric_state` is an optional integer (defaults to 1).

Multi-species lines use comma separation:
    substrates: S[C6H12O6], ATP[C10H16N5O13P3]
    regulators: I, A
"""
macro enzyme_reaction(block)
    parsed = _parse_labeled_block(block,
        Set([:substrates, :products, :regulators,
             :dead_end_inhibitors, :allosteric_regulators,
             :oligomeric_state]),
        Set([:oligomeric_state]))
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
    oligo = haskey(parsed, :oligomeric_state) ? parsed[:oligomeric_state] : 1
    return esc(:(EnzymeReaction($subs, $prods, $regs; oligomeric_state=$oligo)))
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
        substrates: S
        products:   P
        regulators: I

        steps: begin
            ([E, S] ⇌ [ES], [EP, S] ⇌ [EPS])    # parenthesized → shared kinetics
            [ES, I] ⇌ [ESI]                      # dead-end
            [ES]   <--> [EP]
            [EP]   ⇌    [E, P]
        end
    end

Build a plain (non-allosteric) `EnzymeMechanism`.
- `substrates:`, `products:`, `regulators:` accept comma-separated bare symbols.
  Atom brackets (e.g. `S[C]`) are rejected.
- No `enzymes:` block (forms inferred from steps).
- No `constraints:` block — same-kinetics groups are expressed via parenthesized
  step-groups.
- Allosteric-only constructs (`site(...)` / `::Tag` / `allosteric_regulators:` /
  `catalytic_inhibitors:`) are rejected.
"""
macro enzyme_mechanism(block)
    _reject_allosteric_syntax!(block)
    mets_expr, rxns_expr = _parse_plain_mechanism_body(block)
    return esc(:(EnzymeMechanism($mets_expr, $rxns_expr)))
end

function _reject_allosteric_syntax!(block)
    for arg in block.args
        arg isa LineNumberNode && continue
        label_expr = _line_label_expr(arg)
        label_expr === nothing && continue
        if label_expr isa Expr && label_expr.head == :call && label_expr.args[1] == :site
            error("@enzyme_mechanism: `site(...)` belongs in @allosteric_mechanism")
        end
        label_expr in (:allosteric_regulators, :catalytic_inhibitors) &&
            error("@enzyme_mechanism: `$label_expr:` is allosteric-only; " *
                  "use @allosteric_mechanism instead")
    end
end

"""
Return the label of a line in the mechanism block, or `nothing` if not a
labeled line. Handles both Julia parse shapes:
  - `Expr(:call, :(:), label, value)` — single labeled value.
  - `Expr(:tuple, Expr(:call, :(:), label, first), rest...)` — multi-element labeled.
"""
function _line_label_expr(arg)
    if arg isa Expr && arg.head == :call && arg.args[1] == :(:)
        return arg.args[2]
    elseif arg isa Expr && arg.head == :tuple && !isempty(arg.args)
        first_arg = arg.args[1]
        if first_arg isa Expr && first_arg.head == :call && first_arg.args[1] == :(:)
            return first_arg.args[2]
        end
    end
    nothing
end

"""
Parse a labeled-line, returning `(label, values_vector)`. `values_vector` is the
list of args after the label, in source order. Each value is either a bare
Symbol or an `Expr(:(::), name, tag)` (for tagged lists, allosteric only).
"""
function _parse_labeled_line(arg)
    if arg isa Expr && arg.head == :call && arg.args[1] == :(:)
        label = arg.args[2]
        return label, [arg.args[3]]
    elseif arg isa Expr && arg.head == :tuple && !isempty(arg.args)
        first_arg = arg.args[1]
        if first_arg isa Expr && first_arg.head == :call && first_arg.args[1] == :(:)
            label = first_arg.args[2]
            values = Any[first_arg.args[3]]
            append!(values, arg.args[2:end])
            return label, values
        end
    end
    error("Expected `label: value` or `label: v1, v2, ...`; got $arg")
end

function _parse_plain_mechanism_body(block)
    subs_list, prods_list, regs_list = Symbol[], Symbol[], Symbol[]
    steps_block = nothing
    for arg in block.args
        arg isa LineNumberNode && continue
        if arg isa Expr && arg.head == :call && arg.args[1] == :(:) && arg.args[2] == :steps
            steps_block = arg.args[3]
            continue
        end
        label, values = _parse_labeled_line(arg)
        if label == :substrates
            append!(subs_list, _bare_symbols_from_values(values, label))
        elseif label == :products
            append!(prods_list, _bare_symbols_from_values(values, label))
        elseif label == :regulators
            append!(regs_list, _bare_symbols_from_values(values, label))
        else
            error("Unknown @enzyme_mechanism label: $label")
        end
    end
    isempty(subs_list) && error("substrates: not specified")
    isempty(prods_list) && error("products: not specified")
    steps_block === nothing && error("steps: not specified")

    rxns_expr = _parse_steps_block_with_groups(steps_block)
    mets_expr = Expr(:tuple,
        Expr(:tuple, QuoteNode.(subs_list)...),
        Expr(:tuple, QuoteNode.(prods_list)...),
        Expr(:tuple, QuoteNode.(regs_list)...),
    )
    mets_expr, rxns_expr
end

"""
Coerce labeled-line values to bare Symbols. Reject atom brackets and tag
annotations.
"""
function _bare_symbols_from_values(values, label)
    syms = Symbol[]
    for v in values
        if v isa Symbol
            push!(syms, v)
        elseif v isa Expr && v.head == :ref
            error("@enzyme_mechanism: atom bracket syntax `$v` is not allowed " *
                  "at the mechanism level; declare atoms in @enzyme_reaction.")
        elseif v isa Expr && v.head == :(::)
            error("@enzyme_mechanism: tag annotation `$v` is not allowed; " *
                  "tags are only valid in @allosteric_mechanism.")
        else
            error("@enzyme_mechanism `$label:` expects bare Symbol names; got $v")
        end
    end
    syms
end

"""
Parse the steps block. Each top-level expression is either:
  - `Expr(:(::), Expr(:tuple, step1, step2, ...), Tag)` — parenthesized group with tag
    (allosteric only).
  - `Expr(:tuple, step1, step2, ...)` — parenthesized group with no tag (plain mech).
  - `Expr(:call, ⇌|<-->, lhs, Expr(:(::), rhs, Tag))` — single tagged step (allosteric).
  - `Expr(:call, ⇌|<-->, lhs, rhs)` — single untagged step (plain).

With `allow_tag=false` (plain mechanism), reject any `::Tag` annotations and return
just the rxns tuple-Expr. With `allow_tag=true` (allosteric, future Task 2.5),
collect tags and return both.
"""
function _parse_steps_block_with_groups(steps_block; allow_tag::Bool=false)
    next_group = Ref(0)
    rxns = Expr(:tuple)
    tags = Pair{Int, Symbol}[]

    for arg in steps_block.args
        arg isa LineNumberNode && continue

        # Parenthesized-group-with-tag (allosteric)
        if arg isa Expr && arg.head == :(::) &&
           arg.args[1] isa Expr && arg.args[1].head == :tuple
            allow_tag ||
                error("@enzyme_mechanism: tag annotation `$arg` is not allowed")
            next_group[] += 1
            gnum = next_group[]
            tag = arg.args[2]
            tag isa Symbol || error("Step-group tag must be a Symbol; got $tag")
            push!(tags, gnum => tag)
            for step_expr in arg.args[1].args
                push!(rxns.args, _parse_single_step(step_expr, gnum))
            end
        # Parenthesized-group-without-tag (plain)
        elseif arg isa Expr && arg.head == :tuple
            next_group[] += 1
            gnum = next_group[]
            for step_expr in arg.args
                push!(rxns.args, _parse_single_step(step_expr, gnum))
            end
        # Single step (with or without tag)
        elseif arg isa Expr && arg.head == :call
            next_group[] += 1
            gnum = next_group[]
            original = string(arg)
            tag = _peel_step_tag!(arg)
            if tag !== nothing
                allow_tag ||
                    error("@enzyme_mechanism: tag annotation on `$original` " *
                          "is not allowed")
                push!(tags, gnum => tag)
            end
            push!(rxns.args, _parse_single_step(arg, gnum))
        else
            error("Expected step or step-group; got $arg")
        end
    end

    if allow_tag
        return rxns, tags
    else
        return rxns
    end
end

"""
If the step Expr has a `::Tag` attached to its RHS arg, remove the wrapper and
return the tag Symbol. Otherwise return `nothing`. Mutates `step_expr.args[3]`.

Single tagged step parses as:
  Expr(:call, op, Expr(:vect, lhs_syms...), Expr(:(::), Expr(:vect, rhs_syms...), Tag))
"""
function _peel_step_tag!(step_expr)
    rhs = step_expr.args[3]
    if rhs isa Expr && rhs.head == :(::)
        tag = rhs.args[2]
        tag isa Symbol || error("Step tag must be a Symbol; got $tag")
        step_expr.args[3] = rhs.args[1]
        return tag
    end
    nothing
end

"""
Parse a single (already-de-tagged) step `[lhs] ⇌ [rhs]` or `[lhs] <--> [rhs]`.
Returns the 4-tuple Expr `(lhs_syms, rhs_syms, is_eq, kinetic_group)`.
"""
function _parse_single_step(expr, gnum::Int)
    expr isa Expr && expr.head == :call ||
        error("Expected [lhs] ⇌ [rhs] or [lhs] <--> [rhs]; got $expr")
    op = expr.args[1]
    is_eq = op == :⇌
    is_eq || op == :(<-->) ||
        error("Expected ⇌ or <--> step operator; got $op")
    lhs = _parse_step_side_symbols(expr.args[2])
    rhs = _parse_step_side_symbols(expr.args[3])
    Expr(:tuple, lhs, rhs, is_eq, gnum)
end
