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

# Parse one step side into a Vector{_StepSideTerm}. Mirrors
# `_parse_step_side_symbols` but preserves structural info (Call-form
# decomposition, declared-metabolite role) for the new-emission path.
function _parse_step_side_terms(expr, declared_mets::Set{Symbol})
    if expr isa Expr && expr.head == :call && expr.args[1] == :+
        return _StepSideTerm[
            _step_side_term_info(a, declared_mets) for a in expr.args[2:end]
        ]
    end
    _StepSideTerm[_step_side_term_info(expr, declared_mets)]
end

# Build a `_StepSideTerm` for a single term on a step side.
function _step_side_term_info(expr, declared_mets::Set{Symbol})
    if expr isa Symbol
        return expr in declared_mets ?
               _term_metabolite(expr) :
               _term_bare_enzyme(expr)
    elseif expr isa Expr && expr.head == :call &&
           expr.args[1] isa Symbol
        return _call_form_term_info(expr, declared_mets)
    end
    error("Expected metabolite Symbol or species expression on step side; got $expr")
end

# Build a `_StepSideTerm` for a Call-form species `E(S, ATP)` /
# `Estar(B; residual = A - P)`. The legacy synthesized Symbol is built
# alongside the structural fields so a single call to this function
# covers both emission paths.
function _call_form_term_info(expr::Expr, declared_mets::Set{Symbol})
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
    sort!(bound_syms)
    sort!(added_syms)
    sort!(subtracted_syms)
    parts = String[String(conformation)]
    for s in bound_syms; push!(parts, String(s)); end
    if !(isempty(added_syms) && isempty(subtracted_syms))
        push!(parts, "res")
        for s in added_syms;      push!(parts, "+" * String(s)); end
        for s in subtracted_syms; push!(parts, "-" * String(s)); end
    end
    sym = Symbol(join(parts, "_"))
    _StepSideTerm(sym, :call, conformation, bound_syms,
                  added_syms, subtracted_syms)
end

# Build a Symbol matching `name(::Species)` from a `Conformation(...)` AST.
# Accepts:
#   - `E()` → conformation :E, no bound, no residual                → :E
#   - `E(S)` / `E(S, P)` → bound metabolites (sorted by name)        → :E_S / :E_P_S
#   - `Estar(; residual = A - P)` → residual only                    → :Estar_res_+A_-P
#   - `Estar(B; residual = A - P)` → bound + residual                → :Estar_B_res_+A_-P
function _synthesize_species_name(expr::Expr, declared_mets::Set{Symbol})
    _call_form_term_info(expr, declared_mets).sym
end

# Structural side-term record collected during step parsing. Carries
# both the legacy Symbol form (the synthesized species name) and the
# decomposed-Species info needed to emit `Mechanism(...)` directly.
#
# `kind`:
#   :metabolite   — bare Symbol matching a declared metabolite (`S`).
#   :bare_enzyme  — bare Symbol enzyme-form name. Reclassified to
#                   :conformation or :opaque after all steps parsed,
#                   based on whether it appears as a Call-head elsewhere
#                   or matches the single-cap-then-lower conformation
#                   shape (`E`, `Estar`, `Eprime`).
#   :call         — call-form `E(S)` / `Estar(B; residual=A-P)`. Always
#                   decomposed-compatible; carries bound + residual data.
struct _StepSideTerm
    sym::Symbol                          # legacy synthesized Symbol
    kind::Symbol
    conformation::Symbol                 # for :call/bare-enzyme cases
    bound::Vector{Symbol}                # for :call (sorted by parser)
    residual_added::Vector{Symbol}       # for :call (sorted)
    residual_subtracted::Vector{Symbol}  # for :call (sorted)
end

_term_metabolite(sym::Symbol) = _StepSideTerm(
    sym, :metabolite, sym, Symbol[], Symbol[], Symbol[])
_term_bare_enzyme(sym::Symbol) = _StepSideTerm(
    sym, :bare_enzyme, sym, Symbol[], Symbol[], Symbol[])

# A bare Symbol is "conformation-shaped" iff it starts with a single
# capital letter followed by any mix of lowercase letters, digits, and
# underscore-separated lowercase/digit run: :E, :Estar, :Estar2, :E_c,
# :E_secondary. Multi-capital Symbols (:ES, :EAB) and underscore-then-
# uppercase Symbols (:E_S, :Estar_A_B) are opaque legacy enzyme-form
# names — rejected so the dual-grammar emission picks the legacy path.
_is_conformation_shape(sym::Symbol) =
    occursin(r"^[A-Z][a-z0-9]*(_[a-z0-9]+)*$", String(sym))

# True iff every bare-enzyme term in the parsed steps is decomposed-
# compatible. A bare-enzyme term `:X` is compatible iff `:X` is either
# (a) seen as a Call-form head in this steps block, or (b) matches the
# single-cap-then-lower conformation shape.
function _all_bare_terms_compatible(side_terms_per_step,
                                    call_heads::Set{Symbol})
    for (_, lhs, rhs) in side_terms_per_step
        for t in (lhs..., rhs...)
            t.kind === :bare_enzyme || continue
            (t.sym in call_heads || _is_conformation_shape(t.sym)) ||
                return false
        end
    end
    true
end

# True iff at least one Call-form term was used.
function _any_call_form(side_terms_per_step)
    for (_, lhs, rhs) in side_terms_per_step
        for t in (lhs..., rhs...)
            t.kind === :call && return true
        end
    end
    false
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
    return esc(_parse_plain_mechanism_body(block))
end

function _reject_allosteric_syntax!(block)
    for arg in block.args
        arg isa LineNumberNode && continue
        label_expr = _line_label_expr(arg)
        label_expr === nothing && continue
        if label_expr isa Expr && label_expr.head == :call &&
           label_expr.args[1] == :regulatory_site
            error("@enzyme_mechanism: `regulatory_site(...)` belongs in " *
                  "@allosteric_mechanism")
        end
        label_expr in (:allosteric_regulators, :catalytic_inhibitors,
                       :catalytic_multiplicity, :catalytic_steps) &&
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
    role_of = Dict{Symbol,Symbol}()
    for s in subs_list;  role_of[s] = :Substrate;            end
    for p in prods_list; role_of[p] = :Product;              end
    for r in regs_list;  role_of[r] = :CompetitiveInhibitor; end

    rxns_expr, side_terms_per_step =
        _parse_steps_block_with_groups(steps_block, declared_mets)

    if _should_emit_new_grammar(side_terms_per_step)
        return _build_mechanism_expr(subs_list, prods_list, regs_list,
                                     role_of, side_terms_per_step)
    end

    mets_expr = Expr(:tuple,
        Expr(:tuple, QuoteNode.(subs_list)...),
        Expr(:tuple, QuoteNode.(prods_list)...),
        Expr(:tuple, QuoteNode.(regs_list)...),
    )
    :(EnzymeMechanism($mets_expr, $rxns_expr))
end

# Decide whether to emit the new `EnzymeMechanism(Mechanism(...))` shape.
# Triggered by: at least one Call-form term AND every bare-enzyme term
# is decomposed-compatible (matches a Call-form head seen elsewhere or
# matches the single-cap-then-lower conformation shape).
function _should_emit_new_grammar(side_terms_per_step)
    _any_call_form(side_terms_per_step) || return false
    call_heads = Set{Symbol}()
    for (_, lhs, rhs) in side_terms_per_step
        for t in (lhs..., rhs...)
            t.kind === :call && push!(call_heads, t.conformation)
        end
    end
    _all_bare_terms_compatible(side_terms_per_step, call_heads)
end

# Build the `EnzymeMechanism(Mechanism(reaction, grouped_steps))` Expr
# from the structural per-step records collected during parsing.
#
# Atoms for each declared metabolite default to `[:C => 1]` — the same
# placeholder convention `_mechanism_from_legacy_sig` uses when lifting
# the legacy Sig shape. Real atom payloads live at the @enzyme_reaction
# level, not @enzyme_mechanism.
function _build_mechanism_expr(subs_list, prods_list, regs_list,
                               role_of::Dict{Symbol,Symbol},
                               side_terms_per_step)
    reactants_entries = Expr[]
    for s in subs_list
        push!(reactants_entries,
              :(EnzymeRates.ReactantAtoms(
                    EnzymeRates.Substrate($(QuoteNode(s))),
                    Pair{Symbol,Int}[:C => 1])))
    end
    for p in prods_list
        push!(reactants_entries,
              :(EnzymeRates.ReactantAtoms(
                    EnzymeRates.Product($(QuoteNode(p))),
                    Pair{Symbol,Int}[:C => 1])))
    end
    reactants_expr = :(EnzymeRates.ReactantAtoms[$(reactants_entries...)])

    regulator_entries = Expr[]
    for r in regs_list
        push!(regulator_entries,
              :(EnzymeRates.RegulatorMults(
                    EnzymeRates.CompetitiveInhibitor($(QuoteNode(r))),
                    Int[1])))
    end
    regulators_expr = :(EnzymeRates.RegulatorMults[$(regulator_entries...)])

    reaction_expr = :(EnzymeRates.EnzymeReaction(
        $reactants_expr, $regulators_expr, Int[1]))

    # Group structural step records by gnum (preserving source order).
    group_order = Int[]
    by_group = Dict{Int, Vector{Tuple{Vector{_StepSideTerm},
                                      Vector{_StepSideTerm}, Bool}}}()
    for (g, lhs, rhs, is_eq) in side_terms_per_step
        if !haskey(by_group, g)
            by_group[g] = Tuple{Vector{_StepSideTerm},
                                Vector{_StepSideTerm}, Bool}[]
            push!(group_order, g)
        end
        push!(by_group[g], (lhs, rhs, is_eq))
    end

    group_exprs = Expr[]
    for g in group_order
        step_exprs = Expr[]
        for (lhs, rhs, is_eq) in by_group[g]
            push!(step_exprs,
                  _build_step_expr(lhs, rhs, is_eq, role_of))
        end
        push!(group_exprs, :(EnzymeRates.Step[$(step_exprs...)]))
    end
    groups_expr = :(Vector{EnzymeRates.Step}[$(group_exprs...)])

    :(EnzymeRates.EnzymeMechanism(
        EnzymeRates.Mechanism($reaction_expr, $groups_expr)))
end

# Build a `Step(from_species, to_species, bound_metabolite, is_eq)` Expr
# from one step's LHS/RHS structural terms. Each side has exactly one
# enzyme-form term (bare conformation OR call-form) and zero or one
# metabolite terms.
function _build_step_expr(lhs::Vector{_StepSideTerm},
                          rhs::Vector{_StepSideTerm},
                          is_eq::Bool,
                          role_of::Dict{Symbol,Symbol})
    lhs_enzyme, lhs_met = _split_side(lhs)
    rhs_enzyme, rhs_met = _split_side(rhs)
    bound_met_term = lhs_met !== nothing ? lhs_met :
                     rhs_met !== nothing ? rhs_met : nothing
    from_expr = _species_expr_from_term(lhs_enzyme, role_of)
    to_expr   = _species_expr_from_term(rhs_enzyme, role_of)
    met_expr  = bound_met_term === nothing ? :nothing :
                _metabolite_expr(bound_met_term.sym, role_of)
    :(EnzymeRates.Step($from_expr, $to_expr, $met_expr, $is_eq))
end

# Split a step side into its (enzyme_term, optional_metabolite_term).
# Errors if there is not exactly one enzyme term or more than one
# metabolite term.
function _split_side(side::Vector{_StepSideTerm})
    enzyme_term = nothing
    met_term = nothing
    for t in side
        if t.kind === :metabolite
            met_term === nothing ||
                error("@enzyme_mechanism: step side has more than one " *
                      "metabolite term ($(met_term.sym), $(t.sym)); each " *
                      "elementary step binds at most one metabolite.")
            met_term = t
        else
            enzyme_term === nothing ||
                error("@enzyme_mechanism: step side has more than one " *
                      "enzyme-form term ($(enzyme_term.sym), $(t.sym)); " *
                      "each elementary step has exactly one enzyme form " *
                      "per side.")
            enzyme_term = t
        end
    end
    enzyme_term === nothing &&
        error("@enzyme_mechanism: step side has no enzyme-form term " *
              "(terms: $(Symbol[t.sym for t in side])).")
    enzyme_term, met_term
end

# Build a `Species(bound, conformation, residual)` Expr from an enzyme-
# form `_StepSideTerm` (either bare conformation or Call-form).
function _species_expr_from_term(t::_StepSideTerm,
                                 role_of::Dict{Symbol,Symbol})
    bound_entries = Expr[
        _metabolite_expr(b, role_of) for b in t.bound
    ]
    bound_expr = :(EnzymeRates.Metabolite[$(bound_entries...)])
    added_entries = Expr[
        _metabolite_expr(a, role_of) for a in t.residual_added
    ]
    sub_entries = Expr[
        _metabolite_expr(s, role_of) for s in t.residual_subtracted
    ]
    residual_expr = if isempty(added_entries) && isempty(sub_entries)
        :(EnzymeRates.Residual())
    else
        :(EnzymeRates.Residual(
            EnzymeRates.Substrate[$(added_entries...)],
            EnzymeRates.Product[$(sub_entries...)]))
    end
    :(EnzymeRates.Species($bound_expr, $(QuoteNode(t.conformation)),
                          $residual_expr))
end

# Build an Expr that constructs the appropriate `Metabolite` subtype for
# a declared name. The role is looked up from `role_of`.
function _metabolite_expr(name::Symbol, role_of::Dict{Symbol,Symbol})
    role = get(role_of, name, nothing)
    role === nothing &&
        error("@enzyme_mechanism: metabolite `$name` is not declared in " *
              "substrates:, products:, or regulators:.")
    if role === :Substrate
        :(EnzymeRates.Substrate($(QuoteNode(name))))
    elseif role === :Product
        :(EnzymeRates.Product($(QuoteNode(name))))
    elseif role === :CompetitiveInhibitor
        :(EnzymeRates.CompetitiveInhibitor($(QuoteNode(name))))
    elseif role === :AllostericRegulator
        :(EnzymeRates.AllostericRegulator($(QuoteNode(name))))
    else
        error("@enzyme_mechanism: unknown metabolite role $role for $name")
    end
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

Returns the rxns tuple-Expr (legacy emission shape) AND a Vector of
structural per-step records `(gnum, lhs_terms, rhs_terms, is_eq)` used
by the new-emission decision logic in `_parse_plain_mechanism_body` /
`_parse_allosteric_mechanism_body`. With `allow_tag=false` (plain
mechanism), reject any `::Tag` annotations. With `allow_tag=true`
(allosteric mechanism), also return collected `gnum => tag` pairs.
"""
function _parse_steps_block_with_groups(steps_block, declared_mets::Set{Symbol};
                                        allow_tag::Bool=false)
    next_group = Ref(0)
    rxns = Expr(:tuple)
    tags = Pair{Int, Symbol}[]
    side_terms_per_step = Tuple{Int, Vector{_StepSideTerm},
                                Vector{_StepSideTerm}, Bool}[]

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
                push!(rxns.args,
                      _parse_single_step(step_expr, gnum, declared_mets))
                push!(side_terms_per_step,
                      _step_struct_info(step_expr, gnum, declared_mets))
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
                push!(rxns.args,
                      _parse_single_step(step_expr, gnum, declared_mets))
                push!(side_terms_per_step,
                      _step_struct_info(step_expr, gnum, declared_mets))
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
            push!(side_terms_per_step,
                  _step_struct_info(arg, gnum, declared_mets))
        else
            error("Expected step or step-group; got $arg")
        end
    end

    if allow_tag
        return rxns, tags, side_terms_per_step
    else
        return rxns, side_terms_per_step
    end
end

"""
Return the structural per-step record `(gnum, lhs_terms, rhs_terms, is_eq)`
for a single (possibly already-de-tagged) step expression. Mirrors
`_parse_single_step` but produces `_StepSideTerm`s instead of QuoteNoded
Symbol tuple-Exprs.
"""
function _step_struct_info(expr, gnum::Int, declared_mets::Set{Symbol})
    expr isa Expr && expr.head == :call ||
        error("Expected lhs ⇌ rhs or lhs <--> rhs; got $expr")
    op = expr.args[1]
    is_eq = op == :⇌
    is_eq || op == :(<-->) ||
        error("Expected ⇌ or <--> step operator; got $op")
    lhs = _parse_step_side_terms(expr.args[2], declared_mets)
    rhs = _parse_step_side_terms(expr.args[3], declared_mets)
    (gnum, lhs, rhs, is_eq)
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
        catalytic_multiplicity: 2
        allosteric_regulators: A::OnlyR, I::OnlyT

        catalytic_steps: begin
            E + F6P ⇌ E(F6P)        :: EqualRT
            E(F6P) <--> E(F16BP)    :: EqualRT
            (E(F16BP) ⇌ E + F16BP)  :: EqualRT
        end

        regulatory_site(multiplicity = 4): begin
            ligands: A
        end
        regulatory_site(multiplicity = 4): begin
            ligands: I
        end
    end

Build an `AllostericEnzymeMechanism` (MWC, two conformations).

- `substrates:`, `products:`, `catalytic_inhibitors:` accept comma-separated
  bare symbols.
- `allosteric_regulators:` requires `name::Tag` per entry, where Tag is one of
  `OnlyR`, `OnlyT`, `EqualRT`, `NonequalRT`.
- `catalytic_multiplicity: N` is the subunit count for the catalytic site
  (default 1).
- `catalytic_steps: begin ... end` is required (exactly once); each step or
  parenthesized step-group must carry a `::Tag` from the same set. Function-
  call species notation (`E(F6P)`, `Estar(B; residual = A - P)`) is supported.
- `regulatory_site(multiplicity = N): begin ligands: L1, L2 end` declares one
  regulatory site per block with multiplicity `N` and the ligands listed
  inside. Each ligand must appear in `allosteric_regulators:`; a ligand may
  not appear in two sites. Ligands declared in `allosteric_regulators:` but
  not assigned to any `regulatory_site(...):` block default to a single-ligand
  site at the catalytic multiplicity.
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
Match a `regulatory_site(multiplicity = N): begin ligands: ... end` line.
Returns `(mult::Int, ligands::Vector{Symbol})` or `nothing` if the line is
not a regulatory-site declaration.
"""
function _match_regulatory_site_line(arg)
    arg isa Expr && arg.head == :call && length(arg.args) >= 3 &&
        arg.args[1] == :(:) || return nothing
    label, body = arg.args[2], arg.args[3]
    label isa Expr && label.head == :call &&
        label.args[1] == :regulatory_site || return nothing

    mult = nothing
    for kw in label.args[2:end]
        kw isa Expr && (kw.head == :kw || kw.head == :(=)) &&
            kw.args[1] == :multiplicity ||
            error("@allosteric_mechanism: `regulatory_site` only accepts " *
                  "`multiplicity = N` kwarg; got $kw")
        v = kw.args[2]
        v isa Integer && v >= 1 ||
            error("@allosteric_mechanism: `regulatory_site` multiplicity " *
                  "must be a positive Int, got $v")
        mult = Int(v)
    end
    mult === nothing &&
        error("@allosteric_mechanism: `regulatory_site` requires " *
              "`multiplicity = N`")

    body isa Expr && body.head == :block ||
        error("@allosteric_mechanism: `regulatory_site` body must be a " *
              "`begin ... end` block; got $body")
    ligands = Symbol[]
    for inner in body.args
        inner isa LineNumberNode && continue
        label, values = _parse_labeled_line(inner)
        label == :ligands ||
            error("@allosteric_mechanism: `regulatory_site` body expects " *
                  "`ligands:`; got `$label`")
        append!(ligands, _bare_symbols_from_values(values, label))
    end
    isempty(ligands) &&
        error("@allosteric_mechanism: `regulatory_site` has no `ligands:`")
    (mult, ligands)
end

"""
Build the `RegSites` tuple expression. Each entry is
`(ligand_tuple, multiplicity, reg_allo_states)` where `reg_allo_states`
is a dense `Tuple{Symbol...}` parallel to `ligands`. Ligands not assigned
to any explicit `regulatory_site(...):` block become their own single-ligand
site at multiplicity `cat_n`.
"""
function _build_reg_sites_expr(allo_regs, reg_site_specs, cat_n)
    tag_of = Dict{Symbol,Symbol}(allo_regs)
    explicit = Set{Symbol}()
    for (_, ligs) in reg_site_specs
        for l in ligs
            l in explicit && error("@allosteric_mechanism: ligand $l " *
                                   "appears in multiple regulatory sites")
            haskey(tag_of, l) ||
                error("@allosteric_mechanism: ligand $l on a " *
                      "`regulatory_site` is not declared in " *
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

function _parse_allosteric_mechanism_body(block)
    subs_list, prods_list, cat_inhibitors = Symbol[], Symbol[], Symbol[]
    allo_regs = Pair{Symbol,Symbol}[]
    cat_n::Int = 1
    cat_steps_block = nothing
    reg_site_specs = Tuple{Any,Vector{Symbol}}[]

    for arg in block.args
        arg isa LineNumberNode && continue
        reg_site = _match_regulatory_site_line(arg)
        if reg_site !== nothing
            push!(reg_site_specs, reg_site)
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
        elseif label == :catalytic_multiplicity
            length(values) == 1 ||
                error("@allosteric_mechanism: `catalytic_multiplicity:` " *
                      "takes a single Int.")
            v = values[1]
            v isa Integer && v >= 1 ||
                error("@allosteric_mechanism: `catalytic_multiplicity:` " *
                      "must be a positive Int, got $v.")
            cat_n = Int(v)
        elseif label == :catalytic_steps
            cat_steps_block === nothing ||
                error("@allosteric_mechanism: multiple " *
                      "`catalytic_steps:` blocks.")
            length(values) == 1 ||
                error("@allosteric_mechanism: `catalytic_steps:` takes a " *
                      "single `begin ... end` block.")
            cat_steps_block = values[1]
        else
            error("@allosteric_mechanism: unknown label `$label:`")
        end
    end

    isempty(subs_list) &&
        error("@allosteric_mechanism: substrates: not specified")
    isempty(prods_list) &&
        error("@allosteric_mechanism: products: not specified")
    cat_steps_block === nothing &&
        error("@allosteric_mechanism: `catalytic_steps:` block is required")

    declared_mets = Set{Symbol}(subs_list) ∪ Set{Symbol}(prods_list) ∪
                    Set{Symbol}(cat_inhibitors) ∪
                    Set{Symbol}(name for (name, _) in allo_regs)
    rxns_expr, group_tags, _ = _parse_steps_block_with_groups(
        cat_steps_block, declared_mets; allow_tag=true,
    )

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
