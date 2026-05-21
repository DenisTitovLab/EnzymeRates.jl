"""
    @enzyme_reaction begin
        substrates: S[C6H12O6], ATP[C10H16N5O13P3]
        products:   P[C6H13O9P]
        competitive_inhibitors: I            # bare OK; mults default to catalytic
        allosteric_regulators: A(1, 2, 4)    # per-reg mults required
        allowed_catalytic_multiplicities: (1, 2, 4)
        # or shorthand:
        # oligomeric_state: 2
    end

Emit an `EnzymeReaction`.

- `substrates:` / `products:` — comma-separated entries with required atom
  brackets (`S[C6H12O6]`). Multi-atom forms like `[C2,N]` are allowed.
- `competitive_inhibitors:` / `dead_end_inhibitors:` —
  `CompetitiveInhibitor` entries (catalytic-site binding). May be bare
  `I` (multiplicities default to `allowed_catalytic_multiplicities`) or
  `I(m1, m2, ...)` to override.
- `allosteric_regulators:` — `AllostericRegulator` entries. Each entry
  must declare per-regulator multiplicities as `A(m1, m2, ...)` or
  a single value `A(m)`.
- `allowed_catalytic_multiplicities:` — tuple of positive Ints (default `(1,)`).
- `oligomeric_state: N` — shorthand for `allowed_catalytic_multiplicities: (N,)`.
"""
macro enzyme_reaction(block)
    parsed = _parse_reaction_block(block)
    reactants_expr  = _build_reactants_expr(parsed.subs, parsed.prods)
    mults_expr      = _build_catalytic_mults_expr(parsed.mults)
    regulators_expr = _build_regulators_expr(parsed.regs, parsed.mults)
    return esc(:(EnzymeRates.EnzymeReaction(
        $reactants_expr, $regulators_expr, $mults_expr,
    )))
end

# Parse the @enzyme_reaction body. Returns a NamedTuple:
#   subs  ::Vector{Tuple{Symbol, Expr}}  — (name, atoms-tuple-Expr)
#   prods ::Vector{Tuple{Symbol, Expr}}
#   regs  ::Vector{Tuple{Symbol, Symbol, Union{Nothing, Vector{Int}}}}
#                                         — (name, kind, mults) where
#                                           kind ∈ (:competitive, :allosteric)
#                                           and mults is nothing if omitted.
#   mults ::Union{Nothing, Vector{Int}}   — allowed catalytic multiplicities;
#                                           nothing → default (1,).
const _VALID_REACTION_LABELS = Set([
    :substrates, :products,
    :dead_end_inhibitors, :competitive_inhibitors,
    :allosteric_regulators,
    :allowed_catalytic_multiplicities, :oligomeric_state,
])

function _parse_reaction_block(block)
    block isa Expr && block.head === :block ||
        error("@enzyme_reaction: expected a `begin ... end` block, got $block")
    subs  = Tuple{Symbol, Expr}[]
    prods = Tuple{Symbol, Expr}[]
    regs  = Tuple{Symbol, Symbol, Union{Nothing, Vector{Int}}}[]
    mults::Union{Nothing, Vector{Int}} = nothing

    for arg in block.args
        arg isa LineNumberNode && continue
        label, values = _parse_labeled_line(arg)
        label in _VALID_REACTION_LABELS ||
            error("@enzyme_reaction: unknown label `$label:`. Valid labels: " *
                  "$(sort(collect(_VALID_REACTION_LABELS))).")
        if label === :substrates
            append!(subs, _parse_atom_bracket_entries(values, label))
        elseif label === :products
            append!(prods, _parse_atom_bracket_entries(values, label))
        elseif label === :dead_end_inhibitors ||
               label === :competitive_inhibitors
            append!(regs, _parse_regulator_entries(values, :competitive))
        elseif label === :allosteric_regulators
            append!(regs, _parse_regulator_entries(values, :allosteric))
        elseif label === :allowed_catalytic_multiplicities
            mults = _parse_multiplicity_tuple(values, label)
        elseif label === :oligomeric_state
            mults === nothing ||
                error("@enzyme_reaction: cannot specify both `oligomeric_state:` " *
                      "and `allowed_catalytic_multiplicities:`.")
            length(values) == 1 ||
                error("@enzyme_reaction: `oligomeric_state:` takes a single Int.")
            v = values[1]
            v isa Integer && v >= 1 ||
                error("@enzyme_reaction: `oligomeric_state:` must be a positive " *
                      "Int, got $v.")
            mults = Int[v]
        end
    end

    isempty(subs)  && error("@enzyme_reaction: `substrates:` not specified.")
    isempty(prods) && error("@enzyme_reaction: `products:` not specified.")
    (; subs, prods, regs, mults)
end

# Parse `S[C6H12O6]` / `B[N, P]` entries. Returns Vector{Tuple{Symbol, Expr}}
# where the Expr is a tuple of `(elem, count)` pairs.
function _parse_atom_bracket_entries(values, label)
    out = Tuple{Symbol, Expr}[]
    for v in values
        v isa Expr && v.head === :ref && v.args[1] isa Symbol ||
            error("@enzyme_reaction `$label:`: expected `Sym[atoms]`; got $v.")
        atoms_expr = Expr(:tuple)
        for atom_arg in v.args[2:end]
            parsed = _parse_chemical_formula(string(atom_arg))
            for atom in parsed.args
                push!(atoms_expr.args, atom)
            end
        end
        push!(out, (v.args[1]::Symbol, atoms_expr))
    end
    out
end

"""
    _parse_chemical_formula(s::String)

Parse a chemical formula string like `"C6H12O6"` into a tuple expression
of `(element, count)` pairs. Elements are identified by an uppercase
letter optionally followed by lowercase letters, then an optional
integer count (defaults to 1).
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

# Parse `R`, `R(1, 2)`, or `R(4)` entries. Bare `R` produces `nothing` mults
# (filled by the macro from `allowed_catalytic_multiplicities` for competitive
# entries; rejected at emit time for allosteric entries).
function _parse_regulator_entries(values, kind::Symbol)
    out = Tuple{Symbol, Symbol, Union{Nothing, Vector{Int}}}[]
    for v in values
        if v isa Symbol
            push!(out, (v, kind, nothing))
        elseif v isa Expr && v.head === :call && length(v.args) >= 2 &&
               v.args[1] isa Symbol
            name = v.args[1]::Symbol
            ms = Int[]
            for a in v.args[2:end]
                a isa Integer && a >= 1 ||
                    error("@enzyme_reaction: regulator $name multiplicity " *
                          "must be a positive Int, got $a.")
                push!(ms, Int(a))
            end
            push!(out, (name, kind, ms))
        else
            error("@enzyme_reaction: cannot parse regulator entry $v. " *
                  "Use `R` or `R(1, 2)`.")
        end
    end
    out
end

# Parse `(1, 2, 4)` or `4` into Vector{Int}.
function _parse_multiplicity_tuple(values, label)
    length(values) == 1 ||
        error("@enzyme_reaction: `$label:` takes a single tuple, got $values.")
    v = values[1]
    if v isa Integer
        v >= 1 || error("@enzyme_reaction: `$label:` entry must be a positive " *
                        "Int, got $v.")
        return Int[v]
    elseif v isa Expr && v.head === :tuple
        ms = Int[]
        for a in v.args
            a isa Integer && a >= 1 ||
                error("@enzyme_reaction: `$label:` entry must be a positive " *
                      "Int, got $a.")
            push!(ms, Int(a))
        end
        return ms
    end
    error("@enzyme_reaction: `$label:` must be a tuple of positive Ints, got $v.")
end

# Build the reactants vector Expr: each entry is
# `ReactantAtoms(<Substrate|Product>(:Name), [:elem => count, ...])`.
function _build_reactants_expr(subs, prods)
    entries = Expr[]
    for (n, atoms) in subs
        push!(entries, :(EnzymeRates.ReactantAtoms(
            EnzymeRates.Substrate($(QuoteNode(n))),
            $(_atoms_pairs_expr(atoms)),
        )))
    end
    for (n, atoms) in prods
        push!(entries, :(EnzymeRates.ReactantAtoms(
            EnzymeRates.Product($(QuoteNode(n))),
            $(_atoms_pairs_expr(atoms)),
        )))
    end
    :([$(entries...)])
end

# Convert a `(elem, count)` tuple Expr into a `[Symbol => Int, ...]` Vector Expr.
function _atoms_pairs_expr(atoms::Expr)
    pairs = Expr[]
    for atom in atoms.args
        # atom is Expr(:tuple, QuoteNode(:C), 6)
        elem_q = atom.args[1]
        count  = atom.args[2]
        push!(pairs, :($(elem_q) => $(count)))
    end
    :(Pair{Symbol,Int}[$(pairs...)])
end

# Build the regulators vector Expr: each entry is
# `RegulatorMults(<CompetitiveInhibitor|AllostericRegulator>(:Name), [m1, m2, ...])`.
# Bare entries (mults === nothing) inherit `default_mults` (the parsed
# `allowed_catalytic_multiplicities` or its default of `[1]`).
function _build_regulators_expr(regs, default_mults)
    entries = Expr[]
    default_mults_resolved = default_mults === nothing ? Int[1] : default_mults
    for (name, kind, ms) in regs
        subtype_expr = if kind === :allosteric
            :(EnzymeRates.AllostericRegulator($(QuoteNode(name))))
        elseif kind === :competitive
            :(EnzymeRates.CompetitiveInhibitor($(QuoteNode(name))))
        else
            error("internal: unknown regulator kind $kind")
        end
        ms_resolved = ms === nothing ? default_mults_resolved : ms
        ms_expr = :(Int[$(ms_resolved...)])
        push!(entries, :(EnzymeRates.RegulatorMults(
            $(subtype_expr), $(ms_expr),
        )))
    end
    :(EnzymeRates.RegulatorMults[$(entries...)])
end

function _build_catalytic_mults_expr(mults)
    mults === nothing && return :(Int[1])
    :(Int[$(mults...)])
end

function _parse_step_side_symbols(expr, declared_mets::Set{Symbol})
    if expr isa Expr && expr.head == :call && expr.args[1] == :+
        # Multi-term side: `E + S` parses as Expr(:call, :+, …);
        # `E + S + ATP` parses as a single multi-arg `+` call. Each term is
        # one species, lowered to a Symbol via _step_side_term_to_symbol.
        syms = Expr(:tuple)
        for a in expr.args[2:end]
            push!(syms.args, QuoteNode(_step_side_term_to_symbol(a, declared_mets)))
        end
        return syms
    end
    Expr(:tuple, QuoteNode(_step_side_term_to_symbol(expr, declared_mets)))
end

# Resolve one term on a step side to a Symbol.
# - Bare Symbol → metabolite (if declared) or conformation-only species name.
# - Function call `E(S)` / `Estar(B; residual = A - P)` → synthesized species
#   name matching `name(::Species)` from src/types.jl.
function _step_side_term_to_symbol(expr, declared_mets::Set{Symbol})
    if expr isa Symbol
        return expr
    elseif expr isa Expr && expr.head == :call &&
           expr.args[1] isa Symbol
        return _synthesize_species_name(expr, declared_mets)
    end
    error("Expected metabolite Symbol or species expression on step side; got $expr")
end

# Build a Symbol matching `name(::Species)` from a `Conformation(...)` AST.
# Accepts:
#   - `E()` → conformation :E, no bound, no residual                → :E
#   - `E(S)` / `E(S, P)` → bound metabolites (sorted by name)        → :E_S / :E_P_S
#   - `Estar(; residual = A - P)` → residual only                    → :Estar_res_+A_-P
#   - `Estar(B; residual = A - P)` → bound + residual                → :Estar_B_res_+A_-P
function _synthesize_species_name(expr::Expr, declared_mets::Set{Symbol})
    conformation = expr.args[1]::Symbol
    conformation in declared_mets &&
        error("@enzyme_mechanism: conformation label `$conformation` collides " *
              "with declared metabolite `$conformation`; choose a different " *
              "conformation label.")
    bound_syms = Symbol[]
    added_syms = Symbol[]
    subtracted_syms = Symbol[]
    for arg in expr.args[2:end]
        if arg isa Symbol
            arg in declared_mets ||
                error("@enzyme_mechanism: bound metabolite `$arg` in species " *
                      "`$expr` is not declared. Declared: " *
                      "$(sort(collect(declared_mets))).")
            push!(bound_syms, arg)
        elseif arg isa Expr && arg.head === :parameters
            for kw in arg.args
                if kw isa Expr && kw.head === :kw && kw.args[1] === :residual
                    _walk_residual_expr(kw.args[2], true, added_syms,
                                        subtracted_syms, declared_mets)
                else
                    error("@enzyme_mechanism: unknown keyword in species " *
                          "`$expr`: $kw. Only `residual = ...` is allowed.")
                end
            end
        else
            error("@enzyme_mechanism: invalid entry in species `$expr`: $arg")
        end
    end
    parts = String[String(conformation)]
    for s in sort(bound_syms); push!(parts, String(s)); end
    if !(isempty(added_syms) && isempty(subtracted_syms))
        push!(parts, "res")
        for s in sort(added_syms);      push!(parts, "+" * String(s)); end
        for s in sort(subtracted_syms); push!(parts, "-" * String(s)); end
    end
    Symbol(join(parts, "_"))
end

# Walk a residual arithmetic expression (`A`, `A - P`, `S1 + S2 - P1 - P3`, etc.)
# and classify each metabolite Symbol as added (positive) or subtracted (negative).
function _walk_residual_expr(e, sign_positive, added, subtracted, declared_mets)
    if e isa Symbol
        e in declared_mets ||
            error("@enzyme_mechanism: residual entry `$e` is not a declared " *
                  "metabolite. Declared: $(sort(collect(declared_mets))).")
        sign_positive ? push!(added, e) : push!(subtracted, e)
    elseif e isa Expr && e.head === :call
        op = e.args[1]
        if op === :+ && length(e.args) >= 3
            for a in e.args[2:end]
                _walk_residual_expr(a, sign_positive, added, subtracted,
                                    declared_mets)
            end
        elseif op === :- && length(e.args) == 3
            _walk_residual_expr(e.args[2], sign_positive, added, subtracted,
                                declared_mets)
            _walk_residual_expr(e.args[3], !sign_positive, added, subtracted,
                                declared_mets)
        elseif op === :- && length(e.args) == 2
            _walk_residual_expr(e.args[2], !sign_positive, added, subtracted,
                                declared_mets)
        else
            error("@enzyme_mechanism: invalid residual expression: $e")
        end
    else
        error("@enzyme_mechanism: invalid residual expression: $e")
    end
end

"""
    @enzyme_mechanism begin
        substrates: S
        products:   P
        regulators: I

        steps: begin
            E + S ⇌ E(S)                   # function-call species notation
            (E(S) ⇌ E(P), E_alt(S) ⇌ E_alt(P))   # parenthesized → shared kinetics
            E(S) + I ⇌ E(S, I)             # dead-end
            E(P) ⇌ E + P
        end
    end

Build a plain (non-allosteric) `EnzymeMechanism`.

- `substrates:`, `products:`, `regulators:` accept comma-separated bare
  symbols. Atom brackets (`S[C]`) are rejected at the mechanism level.
- `regulators:` entries are treated as `CompetitiveInhibitor`s when later
  passed to `EnzymeReaction`.
- Same-kinetics groups are expressed via parenthesized step-groups; no
  `constraints:` block needed.
- Allosteric-only constructs (`site(...)` / `::Tag` /
  `allosteric_regulators:` / `catalytic_inhibitors:`) are rejected.

Species notation on step sides:

- Bare Symbol that matches a declared metabolite → that metabolite.
- Bare Symbol otherwise (e.g. `E`, `Estar`, `ES`) → conformation-only
  species named after the Symbol.
- `E(S)` / `E(S, P)` → species with conformation `:E` and bound
  metabolites; synthesized name is `:E_<bound...>` with bound names
  sorted alphabetically (matching `name(::Species)`).
- `Estar(; residual = A - P)` → species with empty bound and a residual
  recording `+A` / `−P`; synthesized name is `:Estar_res_+A_-P`.
- `Estar(B; residual = A - P)` → bound + residual; name
  `:Estar_B_res_+A_-P`.

Conformation labels cannot shadow declared metabolite names.
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

    declared_mets = Set{Symbol}(subs_list) ∪ Set{Symbol}(prods_list) ∪
                    Set{Symbol}(regs_list)
    rxns_expr = _parse_steps_block_with_groups(steps_block, declared_mets)
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
just the rxns tuple-Expr. With `allow_tag=true` (allosteric mechanism), collect
tags and return both.
"""
function _parse_steps_block_with_groups(steps_block, declared_mets::Set{Symbol};
                                        allow_tag::Bool=false)
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
                push!(rxns.args, _parse_single_step(step_expr, gnum, declared_mets))
            end
        # Parenthesized-group-without-tag (plain)
        elseif arg isa Expr && arg.head == :tuple
            allow_tag &&
                error("@allosteric_mechanism: parenthesized step group " *
                      "`$(arg)` is missing `:: <:OnlyR|:EqualRT|:NonequalRT>` " *
                      "annotation. Add `:: <state>` after the closing paren.")
            next_group[] += 1
            gnum = next_group[]
            for step_expr in arg.args
                push!(rxns.args, _parse_single_step(step_expr, gnum, declared_mets))
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
            elseif allow_tag
                error("@allosteric_mechanism: step `$(original)` is missing " *
                      "`:: <:OnlyR|:EqualRT|:NonequalRT>` annotation. Add " *
                      "`:: <state>` after the step expression.")
            end
            push!(rxns.args, _parse_single_step(arg, gnum, declared_mets))
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

Two RHS shapes carry a tag:
  - `Sym :: Tag` parses as `Expr(:(::), Sym, Tag)`.
  - `S1 + S2 + … + LastSym :: Tag` parses as
    `Expr(:call, :+, S1, …, Expr(:(::), LastSym, Tag))` because `::` binds
    tighter than `+`. We peel the inner `::` and put `LastSym` back in place.
"""
function _peel_step_tag!(step_expr)
    rhs = step_expr.args[3]
    if rhs isa Expr && rhs.head == :(::)
        tag = rhs.args[2]
        tag isa Symbol || error("Step tag must be a Symbol; got $tag")
        step_expr.args[3] = rhs.args[1]
        return tag
    elseif rhs isa Expr && rhs.head == :call && rhs.args[1] == :+
        last = rhs.args[end]
        if last isa Expr && last.head == :(::)
            tag = last.args[2]
            tag isa Symbol || error("Step tag must be a Symbol; got $tag")
            rhs.args[end] = last.args[1]
            return tag
        end
    end
    nothing
end

"""
Parse a single (already-de-tagged) step `lhs ⇌ rhs` or `lhs <--> rhs`, where
each side is either a bare Symbol or `Sym + Sym + …`. Returns the 4-tuple
Expr `(lhs_syms, rhs_syms, is_eq, kinetic_group)`.
"""
function _parse_single_step(expr, gnum::Int, declared_mets::Set{Symbol})
    expr isa Expr && expr.head == :call ||
        error("Expected lhs ⇌ rhs or lhs <--> rhs; got $expr")
    op = expr.args[1]
    is_eq = op == :⇌
    is_eq || op == :(<-->) ||
        error("Expected ⇌ or <--> step operator; got $op")
    lhs = _parse_step_side_symbols(expr.args[2], declared_mets)
    rhs = _parse_step_side_symbols(expr.args[3], declared_mets)
    Expr(:tuple, lhs, rhs, is_eq, gnum)
end

"""
    @allosteric_mechanism begin
        substrates: F6P
        products:   F16BP
        allosteric_regulators: I::OnlyT

        site(:catalytic, 2): begin
            steps: begin
                E + F6P ⇌ E_F6P       :: EqualRT
                E_F6P <--> E_F16BP    :: EqualRT
                E_F16BP ⇌ E + F16BP   :: EqualRT
            end
        end

        site(:regulatory, 2): begin
            ligands: A, I
        end
    end

Build an `AllostericEnzymeMechanism` (MWC, two conformations).
- `substrates:`, `products:`, `catalytic_inhibitors:` accept comma-separated
  bare symbols.
- `allosteric_regulators:` requires `name::Tag` per entry, where Tag is one of
  `OnlyR`, `OnlyT`, `EqualRT`, `NonequalRT`.
- `site(:catalytic, N): begin steps: ... end` is required (exactly once); each
  step or step-group must carry a `::Tag` from the same set.
- `site(:regulatory, N): begin ligands: A, I end` blocks group competing
  ligands into a single regulatory site (optional, multiple allowed). Ligands
  not listed in any `site(:regulatory, ...):` block default to a per-ligand
  site at multiplicity `N` (the catalytic multiplicity).
"""
macro allosteric_mechanism(block)
    return esc(_parse_allosteric_mechanism_body(block))
end

const _ALLOSTERIC_REG_STATES = Set([:OnlyR, :OnlyT, :EqualRT, :NonequalRT])

"""
Coerce labeled-line values to `(name, tag)` pairs. Each value must be
`Expr(:(::), name, tag)`; bare symbols are rejected. Used for
`allosteric_regulators:` and similar tagged lists.
"""
function _tagged_symbols_from_values(values, label, valid_tags)
    pairs = Pair{Symbol,Symbol}[]
    for v in values
        v isa Expr && v.head == :(::) ||
            error("@allosteric_mechanism `$label:` requires per-entry " *
                  "::Tag annotations (e.g., I::OnlyT); got $v")
        name, tag = v.args[1], v.args[2]
        name isa Symbol ||
            error("@allosteric_mechanism `$label:`: expected Symbol name " *
                  "in `name::Tag`; got $name")
        tag isa Symbol ||
            error("@allosteric_mechanism `$label:`: tag must be a Symbol; " *
                  "got $tag")
        tag in valid_tags ||
            error("@allosteric_mechanism `$label:`: tag :$tag not in " *
                  "($(_format_state_set(valid_tags)))")
        push!(pairs, name => tag)
    end
    pairs
end

"""Format a state set as a sorted, comma-joined list for error messages."""
_format_state_set(tags) = join((":$t" for t in sort(collect(tags))), ", ")

"""
Extract the inner `steps:` block from a `site(:catalytic, N):` body.
The body is `begin steps: <block> end`; reject any other label.
"""
function _extract_steps_block(value)
    value isa Expr && value.head == :block ||
        error("@allosteric_mechanism: site(:catalytic, ...) body must be " *
              "a `begin ... end` block; got $value")
    steps_block = nothing
    for arg in value.args
        arg isa LineNumberNode && continue
        if arg isa Expr && arg.head == :call && arg.args[1] == :(:) &&
           arg.args[2] == :steps
            steps_block === nothing ||
                error("@allosteric_mechanism: multiple `steps:` blocks " *
                      "inside site(:catalytic, ...)")
            steps_block = arg.args[3]
        else
            error("@allosteric_mechanism: site(:catalytic, ...) body " *
                  "expects `steps:`; got $arg")
        end
    end
    steps_block === nothing &&
        error("@allosteric_mechanism: site(:catalytic, ...) requires " *
              "a `steps:` block")
    steps_block
end

"""
Parse the body of a `site(:regulatory, N): begin ligands: A, B end` block.
Returns the ligand symbol vector.
"""
function _parse_regulatory_site_inner(value)
    value isa Expr && value.head == :block ||
        error("@allosteric_mechanism: site(:regulatory, ...) body must be " *
              "a `begin ... end` block; got $value")
    ligands = Symbol[]
    for arg in value.args
        arg isa LineNumberNode && continue
        label, values = _parse_labeled_line(arg)
        label == :ligands ||
            error("@allosteric_mechanism: site(:regulatory, ...) body " *
                  "expects `ligands:`; got `$label`")
        append!(ligands, _bare_symbols_from_values(values, label))
    end
    ligands
end

"""
Build the `RegSites` tuple expression. Each entry is
`(ligand_tuple, multiplicity, reg_allo_states)` where `reg_allo_states`
is a dense `Tuple{Symbol...}` parallel to `ligands`. Ligands not assigned
to any explicit `site(:regulatory, ...):` block become their own
single-ligand site at multiplicity `cat_n`.
"""
function _build_reg_sites_expr(allo_regs, reg_site_specs, cat_n)
    tag_of = Dict{Symbol,Symbol}(allo_regs)
    explicit = Set{Symbol}()
    for (_, ligs) in reg_site_specs
        for l in ligs
            l in explicit && error("@allosteric_mechanism: ligand $l " *
                                   "appears in multiple regulatory sites")
            haskey(tag_of, l) ||
                error("@allosteric_mechanism: ligand $l in " *
                      "site(:regulatory, ...) is not declared in " *
                      "`allosteric_regulators:`")
            push!(explicit, l)
        end
    end

    sites = Tuple{Any,Vector{Symbol}}[]
    for (mult, ligs) in reg_site_specs
        push!(sites, (mult, ligs))
    end
    for (name, _) in allo_regs
        name in explicit && continue
        push!(sites, (cat_n, [name]))
    end

    entries = Expr[]
    for (mult, ligs) in sites
        ligs_tuple = Expr(:tuple, QuoteNode.(ligs)...)
        states_tuple = Expr(:tuple, (QuoteNode(tag_of[l]) for l in ligs)...)
        entry = Expr(:tuple, ligs_tuple, mult, states_tuple)
        push!(entries, entry)
    end
    Expr(:tuple, entries...)
end

"""
Build the `CatSites` expression `(multiplicity, cat_allo_states)` for
the macro. `cat_allo_states` is a dense `Tuple{Symbol...}` with one
entry per kinetic group in source order.
"""
function _build_cat_sites_expr(cat_n, group_tags)
    for (_, tag) in group_tags
        tag in _ALLOSTERIC_REG_STATES ||
            error("@allosteric_mechanism: catalytic step tag :$tag not in " *
                  "($(_format_state_set(_ALLOSTERIC_REG_STATES)))")
    end
    tag_of = Dict{Int,Symbol}(group_tags)
    n_groups = isempty(group_tags) ? 0 : maximum(g for (g, _) in group_tags)
    states_tuple = Expr(:tuple,
        (QuoteNode(get(tag_of, g, :NonequalRT)) for g in 1:n_groups)...)
    Expr(:tuple, cat_n, states_tuple)
end

"""
Detect a `site(KIND, N): begin ... end` line. Returns `(kind::Symbol, n_expr,
body)` or `nothing` if `arg` isn't a site line.
"""
function _match_site_line(arg)
    arg isa Expr && arg.head == :call && arg.args[1] == :(:) || return nothing
    label, value = arg.args[2], arg.args[3]
    label isa Expr && label.head == :call && label.args[1] == :site ||
        return nothing
    site_kind = label.args[2]
    site_kind isa QuoteNode ||
        error("@allosteric_mechanism: site kind must be a Symbol literal; " *
              "got $site_kind")
    (site_kind.value, label.args[3], value)
end

function _parse_allosteric_mechanism_body(block)
    subs_list, prods_list, cat_inhibitors = Symbol[], Symbol[], Symbol[]
    allo_regs = Pair{Symbol,Symbol}[]
    cat_n = nothing
    cat_steps_block = nothing
    reg_site_specs = Tuple{Any,Vector{Symbol}}[]

    for arg in block.args
        arg isa LineNumberNode && continue
        site = _match_site_line(arg)
        if site !== nothing
            kind, n_expr, body = site
            if kind == :catalytic
                cat_n === nothing ||
                    error("@allosteric_mechanism: multiple " *
                          "site(:catalytic, ...) blocks")
                cat_n = n_expr
                cat_steps_block = _extract_steps_block(body)
            elseif kind == :regulatory
                push!(reg_site_specs,
                      (n_expr, _parse_regulatory_site_inner(body)))
            else
                error("@allosteric_mechanism: unknown site kind :$kind")
            end
            continue
        end
        label, values = _parse_labeled_line(arg)
        if label == :substrates
            append!(subs_list, _bare_symbols_from_values(values, label))
        elseif label == :products
            append!(prods_list, _bare_symbols_from_values(values, label))
        elseif label == :catalytic_inhibitors
            append!(cat_inhibitors,
                    _bare_symbols_from_values(values, label))
        elseif label == :allosteric_regulators
            append!(allo_regs,
                    _tagged_symbols_from_values(values, label,
                                                _ALLOSTERIC_REG_STATES))
        else
            error("@allosteric_mechanism: unknown label `$label:`")
        end
    end

    isempty(subs_list) &&
        error("@allosteric_mechanism: substrates: not specified")
    isempty(prods_list) &&
        error("@allosteric_mechanism: products: not specified")
    cat_n === nothing &&
        error("@allosteric_mechanism: site(:catalytic, N): is required")

    declared_mets = Set{Symbol}(subs_list) ∪ Set{Symbol}(prods_list) ∪
                    Set{Symbol}(cat_inhibitors) ∪
                    Set{Symbol}(name for (name, _) in allo_regs)
    rxns_expr, group_tags = _parse_steps_block_with_groups(
        cat_steps_block, declared_mets; allow_tag=true,
    )
    # Bare-step rejection now happens inside _parse_steps_block_with_groups.

    cm_mets_expr = Expr(:tuple,
        Expr(:tuple, QuoteNode.(subs_list)...),
        Expr(:tuple, QuoteNode.(prods_list)...),
        Expr(:tuple, QuoteNode.(cat_inhibitors)...),
    )
    cm_expr = :(EnzymeMechanism($cm_mets_expr, $rxns_expr))

    cat_sites_expr = _build_cat_sites_expr(cat_n, group_tags)
    reg_sites_expr = _build_reg_sites_expr(allo_regs, reg_site_specs, cat_n)

    :(AllostericEnzymeMechanism($cm_expr, $cat_sites_expr, $reg_sites_expr))
end
