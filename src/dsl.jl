# ABOUTME: DSL macros @enzyme_reaction, @enzyme_mechanism, @allosteric_mechanism.
# ABOUTME: Parse reaction/mechanism declarations into concrete struct constructors.

"""
    @enzyme_reaction begin
        substrates: S[C6H12O6], ATP[C10H16N5O13P3]
        products:   P[C6H13O9P]
        competitive_inhibitors: I            # bare OK; mults default to catalytic
        allosteric_regulators: A(1, 2, 4)    # bare OK; mults default to catalytic
        allowed_catalytic_multiplicities: (1, 2, 4)
        # or shorthand:
        # oligomeric_state: 2
    end

Emit an `EnzymeReaction`.

- `substrates:` / `products:` — comma-separated entries with required atom
  brackets (`S[C6H12O6]`). Multi-atom forms like `[C2,N]` are allowed.
  An enzyme form's name joins its bound metabolites' names without a
  separator, so names whose joins coincide are rejected: with `Ac`, `CoA`
  and `AcCoA`, E with Ac and CoA bound and E with AcCoA bound are both
  `EAcCoA`. Rename one (`Acetyl`).
- `competitive_inhibitors:` / `dead_end_inhibitors:` —
  `CompetitiveInhibitor` entries (catalytic-site binding). May be bare
  `I` (multiplicities default to `allowed_catalytic_multiplicities`) or
  `I(m1, m2, ...)` to override.
- `allosteric_regulators:` — `AllostericRegulator` entries. May be bare
  `A` (multiplicities default to `allowed_catalytic_multiplicities`),
  `A(m1, m2, ...)` or a single value `A(m)`. An entry may carry a type
  tag, `A::Activator` or `A::Inhibitor`; untagged entries default to
  `:unspecified`. Type tags are only valid on `allosteric_regulators:`
  entries.
- `allowed_catalytic_multiplicities:` — tuple of positive Ints (default `(1,)`).
- `oligomeric_state: N` — shorthand for `allowed_catalytic_multiplicities: (N,)`.
- `shared_catalytic_site: (sub, prod), ...` — pairs of substrate/product
  names (in either order) that bind at the same catalytic site.
"""
macro enzyme_reaction(block)
    (; reactants, regs, mults, shared) = _parse_reaction_block(block)
    return esc(_reaction_expr(reactants, regs, mults, shared))
end

const _VALID_REACTION_LABELS = Set([
    :substrates, :products,
    :dead_end_inhibitors, :competitive_inhibitors,
    :allosteric_regulators,
    :allowed_catalytic_multiplicities, :oligomeric_state,
    :shared_catalytic_site,
])

"""
Parse the `@enzyme_reaction` body. Returns a `NamedTuple`:
- `reactants ::Vector{Tuple{Symbol, Symbol, Vector{Pair{Symbol,Int}}}}` —
  `(type, name, atoms)` with `type ∈ (:Substrate, :Product)` and `atoms` the
  `element => count` pairs.
- `regs  ::Vector{Tuple{Symbol, Symbol, Vector{Int}, Symbol}}` —
  `(type, name, mults, reg_type)` where `type ∈ (:CompetitiveInhibitor,
  :AllostericRegulator)`, a bare entry takes the catalytic multiplicities, and
  `reg_type ∈ (:activator, :inhibitor, :unspecified)`.
- `mults ::Vector{Int}` — allowed catalytic multiplicities (default `[1]`).
- `shared ::Vector{Tuple{Symbol,Symbol}}` — `shared_catalytic_site:` pairs.
"""
function _parse_reaction_block(block)
    block isa Expr && block.head === :block ||
        error("@enzyme_reaction: expected a `begin ... end` block, got $block")
    reactants = Tuple{Symbol, Symbol, Vector{Pair{Symbol,Int}}}[]
    regs  = Tuple{Symbol, Symbol, Union{Nothing, Vector{Int}}, Symbol}[]
    mults::Union{Nothing, Vector{Int}} = nothing
    shared = Tuple{Symbol,Symbol}[]

    for arg in block.args
        arg isa LineNumberNode && continue
        label, values = _parse_labeled_line(arg)
        label in _VALID_REACTION_LABELS ||
            error("@enzyme_reaction: unknown label `$label:`. Valid labels: " *
                  "$(sort(collect(_VALID_REACTION_LABELS))).")
        if label === :substrates || label === :products
            type = label === :substrates ? :Substrate : :Product
            for (name, atoms) in _parse_atom_bracket_entries(values, label)
                push!(reactants, (type, name, atoms))
            end
        elseif label === :dead_end_inhibitors ||
               label === :competitive_inhibitors
            append!(regs, _parse_regulator_entries(values, :CompetitiveInhibitor))
        elseif label === :allosteric_regulators
            append!(regs, _parse_regulator_entries(values, :AllostericRegulator))
        elseif label === :allowed_catalytic_multiplicities
            mults = _parse_multiplicity_tuple(values, label)
        elseif label === :oligomeric_state
            mults === nothing ||
                error("@enzyme_reaction: cannot specify both `oligomeric_state:` " *
                      "and `allowed_catalytic_multiplicities:`.")
            length(values) == 1 ||
                error("@enzyme_reaction: `oligomeric_state:` takes a single Int.")
            mults = Int[_positive_int(values[1], "@enzyme_reaction: `oligomeric_state:`")]
        elseif label === :shared_catalytic_site
            append!(shared, _parse_shared_site_pairs(values))
        end
    end

    any(r -> r[1] === :Substrate, reactants) ||
        error("@enzyme_reaction: `substrates:` not specified.")
    any(r -> r[1] === :Product, reactants) ||
        error("@enzyme_reaction: `products:` not specified.")
    catalytic_mults = something(mults, [1])
    regs = [(t, n, something(ms, catalytic_mults), rt) for (t, n, ms, rt) in regs]
    (; reactants, regs, mults = catalytic_mults, shared)
end

"""
Parse `S[C6H12O6]` / `B[N, P]` entries into `(name, atoms)` pairs, `atoms` a
`Vector{Pair{Symbol,Int}}` of `element => count`.
"""
function _parse_atom_bracket_entries(values, label)
    map(values) do v
        v isa Expr && v.head === :ref && v.args[1] isa Symbol ||
            error("@enzyme_reaction `$label:`: expected `Sym[atoms]`; got $v.")
        (v.args[1]::Symbol, Pair{Symbol,Int}[
            p for a in v.args[2:end] for p in _parse_chemical_formula(string(a))])
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

"""`v` as an `Int`; errors unless it is a positive integer literal. `what` names it."""
_positive_int(v, what) =
    v isa Integer && v >= 1 ? Int(v) : error("$what must be a positive Int, got $v.")

"""
Parse `R`, `R(1, 2)`, or `R(4)` entries, optionally tagged with a regulator
type: `R::Activator`, `R::Inhibitor`, or `R(1, 2)::Inhibitor`, into
`(type, name, mults, reg_type)`. Bare `R` produces `nothing` mults (filled from
`allowed_catalytic_multiplicities`, default `[1]`) and reg_type `:unspecified`.
A type tag is only valid on `type === :AllostericRegulator` entries.
"""
function _parse_regulator_entries(values, type::Symbol)
    map(values) do v
        reg_type, inner = :unspecified, v
        if v isa Expr && v.head === :(::)
            inner, tag = v.args
            reg_type = tag === :Activator ? :activator :
                   tag === :Inhibitor ? :inhibitor :
                   error("@enzyme_reaction: unknown regulator type tag ::$tag " *
                         "on $v; use ::Activator or ::Inhibitor.")
            type === :AllostericRegulator ||
                error("@enzyme_reaction: regulator type ::$tag on $v is only " *
                      "valid on allosteric_regulators: entries.")
        end
        inner isa Symbol && return (type, inner, nothing, reg_type)
        inner isa Expr && inner.head === :call && length(inner.args) >= 2 &&
            inner.args[1] isa Symbol ||
            error("@enzyme_reaction: cannot parse regulator entry $v. " *
                  "Use `R` or `R(1, 2)`.")
        name = inner.args[1]::Symbol
        what = "@enzyme_reaction: regulator $name multiplicity"
        (type, name, Int[_positive_int(a, what) for a in inner.args[2:end]], reg_type)
    end
end

"""Parse `(1, 2, 4)` or `4` into `Vector{Int}`."""
function _parse_multiplicity_tuple(values, label)
    length(values) == 1 ||
        error("@enzyme_reaction: `$label:` takes a single tuple, got $values.")
    v = values[1]
    Int[_positive_int(a, "@enzyme_reaction: `$label:` entry")
        for a in (v isa Expr && v.head === :tuple ? v.args : (v,))]
end

"""
Build the `EnzymeReaction` `Expr` from parsed declarations: `reactants` as
`(type, name, atoms)` and `regulators` as `(type, name, mults, reg_type)`, each
`type` naming the `Metabolite` subtype to construct.
"""
function _reaction_expr(reactants, regulators, mults, shared)
    reactant_exprs = [
        :(EnzymeRates.ReactantAtoms(EnzymeRates.$type($(QuoteNode(name))),
            Pair{Symbol,Int}[$((:($(QuoteNode(e)) => $c) for (e, c) in atoms)...)]))
        for (type, name, atoms) in reactants]
    regulator_exprs = [
        :(EnzymeRates.RegulatorMults(EnzymeRates.$type($(QuoteNode(name))),
                                     Int[$(ms...)], $(QuoteNode(reg_type))))
        for (type, name, ms, reg_type) in regulators]
    shared_exprs = [:(($(QuoteNode(a)), $(QuoteNode(b)))) for (a, b) in shared]
    :(EnzymeRates.EnzymeReaction(
        EnzymeRates.ReactantAtoms[$(reactant_exprs...)],
        EnzymeRates.RegulatorMults[$(regulator_exprs...)], Int[$(mults...)];
        shared_catalytic_site = Tuple{Symbol,Symbol}[$(shared_exprs...)]))
end

"""
Build the `EnzymeReaction` `Expr` of a mechanism macro, whose metabolites are
declared without atoms. Atoms are balanced placeholders: each substrate carries
`length(prods)` carbons and each product carries `length(subs)` carbons, so
total carbon is `length(subs) * length(prods)` on both sides and the
`EnzymeReaction` atom-balance check passes for any substrate/product counts.
Real atom payloads live at the `@enzyme_reaction` level. `regs` become
`CompetitiveInhibitor`s at multiplicity 1.
"""
_mechanism_reaction_expr(subs, prods, regs) = _reaction_expr(
    [[(:Substrate, s, [:C => length(prods)]) for s in subs];
     [(:Product, p, [:C => length(subs)]) for p in prods]],
    [(:CompetitiveInhibitor, r, [1], :unspecified) for r in regs], [1], ())

"""
Parse one step side into a `Vector{_StepSideTerm}`, preserving structural
info (Call-form decomposition, declared-metabolite role) for emission.
"""
function _parse_step_side_terms(expr, declared_mets::Set{Symbol})
    if expr isa Expr && expr.head == :call && expr.args[1] == :+
        return _StepSideTerm[
            _step_side_term_info(a, declared_mets) for a in expr.args[2:end]
        ]
    end
    _StepSideTerm[_step_side_term_info(expr, declared_mets)]
end

"""Build a `_StepSideTerm` for a single term on a step side."""
function _step_side_term_info(expr, declared_mets::Set{Symbol})
    if expr isa Symbol
        return expr in declared_mets ?
               _term_metabolite(expr) :
               _term_bare_enzyme(expr)
    elseif expr isa Expr && expr.head == :(::)
        name = expr.args[1]
        tag  = expr.args[2]
        (name isa Symbol) ||
            error("@enzyme_mechanism: expected `name::Inh`; got $expr")
        tag === :Inh ||
            error("@enzyme_mechanism: unknown role tag ::$tag on $name; " *
                  "only ::Inh is supported.")
        name in declared_mets ||
            error("@enzyme_mechanism: tagged metabolite `$name` in `$expr` is " *
                  "not declared. Declared: $(sort(collect(declared_mets))).")
        return _term_metabolite(name, :inh)
    elseif expr isa Expr && expr.head == :call &&
           expr.args[1] isa Symbol
        return _call_form_term_info(expr, declared_mets)
    end
    error("Expected metabolite Symbol or species expression on step side; got $expr")
end

"""
Build a `_StepSideTerm` for a Call-form species `E(S, ATP)` /
`Estar(B; residual = A - P)`, decomposing it into the conformation
label, bound metabolites (with roles), and residual atom deltas.
"""
function _call_form_term_info(expr::Expr, declared_mets::Set{Symbol})
    conformation = expr.args[1]::Symbol
    conformation in declared_mets &&
        error("@enzyme_mechanism: conformation label `$conformation` collides " *
              "with declared metabolite `$conformation`; choose a different " *
              "conformation label.")
    bound_syms = Symbol[]
    bound_roles = Symbol[]
    added_syms = Symbol[]
    subtracted_syms = Symbol[]
    for arg in expr.args[2:end]
        if arg isa Symbol
            arg in declared_mets ||
                error("@enzyme_mechanism: bound metabolite `$arg` in species " *
                      "`$expr` is not declared. Declared: " *
                      "$(sort(collect(declared_mets))).")
            push!(bound_syms, arg)
            push!(bound_roles, :default)
        elseif arg isa Expr && arg.head === :(::)
            name = arg.args[1]
            tag  = arg.args[2]
            (name isa Symbol) ||
                error("@enzyme_mechanism: expected `name::Inh` in species " *
                      "`$expr`; got $arg")
            tag === :Inh ||
                error("@enzyme_mechanism: unknown role tag ::$tag on $name; " *
                      "only ::Inh is supported.")
            name in declared_mets ||
                error("@enzyme_mechanism: bound metabolite `$name` in species " *
                      "`$expr` is not declared. Declared: " *
                      "$(sort(collect(declared_mets))).")
            push!(bound_syms, name)
            push!(bound_roles, :inh)
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
    perm = sortperm(collect(zip(bound_syms, bound_roles)))
    bound_syms = bound_syms[perm]
    bound_roles = bound_roles[perm]
    sort!(added_syms)
    sort!(subtracted_syms)
    _StepSideTerm(conformation, :call, conformation, bound_syms, bound_roles,
                  added_syms, subtracted_syms, :default)
end

"""
Structural side-term record collected during step parsing. Carries the
decomposed-Species info needed to emit `Mechanism(...)` directly.

`kind`:
- `:metabolite`  — bare `Symbol` matching a declared metabolite (`S`).
- `:bare_enzyme` — bare `Symbol` enzyme-form name. Reclassified to
  `:conformation` or `:opaque` after all steps parsed, based on whether
  it appears as a Call-head elsewhere or matches the single-cap-then-lower
  conformation shape (`E`, `Estar`, `Eprime`).
- `:call`        — call-form `E(S)` / `Estar(B; residual=A-P)`. Always
  decomposed-compatible; carries bound + residual data.
"""
struct _StepSideTerm
    sym::Symbol                          # metabolite/enzyme name, or :call conformation
    kind::Symbol
    conformation::Symbol                 # for :call/bare-enzyme cases
    bound::Vector{Symbol}                # for :call (sorted by parser)
    bound_roles::Vector{Symbol}          # parallel to `bound`: :default or :inh
    residual_added::Vector{Symbol}       # for :call (sorted)
    residual_subtracted::Vector{Symbol}  # for :call (sorted)
    role::Symbol                         # for :metabolite free term: :default or :inh
end

_term_metabolite(sym::Symbol, role::Symbol = :default) = _StepSideTerm(
    sym, :metabolite, sym, Symbol[], Symbol[], Symbol[], Symbol[], role)
_term_bare_enzyme(sym::Symbol) = _StepSideTerm(
    sym, :bare_enzyme, sym, Symbol[], Symbol[], Symbol[], Symbol[], :default)

"""
A bare `Symbol` is "conformation-shaped" iff it starts with a single
capital letter followed by any mix of lowercase letters, digits, and
underscore-separated lowercase/digit runs: `:E`, `:Estar`, `:Estar2`, `:E_c`,
`:E_secondary`. Multi-capital `Symbol`s (`:ES`, `:EAB`) and underscore-then-
uppercase `Symbol`s (`:E_S`, `:Estar_A_B`) are opaque bound-form names —
rejected in favor of decomposed call notation.
"""
_is_conformation_shape(sym::Symbol) =
    occursin(r"^[A-Z][a-z0-9]*(_[a-z0-9]+)*$", String(sym))


"""
Reject opaque bound-form bare-enzyme names. A bare-enzyme term `:X` is
acceptable iff `:X` is a call-form head seen in this steps block (`E` in
`E(S)`) or matches the conformation shape (`:E`, `:Estar`, `:E_c`).
Multi-capital (`:ES`) and underscore-then-uppercase (`:E_S`) names are
opaque and rejected in favor of decomposed call notation. `macro_name`
names the invoking macro so the error points at the right docs.
"""
function _reject_opaque_bound_forms(side_terms_per_step, macro_name::String)
    call_heads = Set{Symbol}()
    for (_, lhs, rhs, _) in side_terms_per_step
        for t in (lhs..., rhs...)
            t.kind === :call && push!(call_heads, t.conformation)
        end
    end
    for (_, lhs, rhs, _) in side_terms_per_step
        for t in (lhs..., rhs...)
            t.kind === :bare_enzyme || continue
            (t.sym in call_heads || _is_conformation_shape(t.sym)) && continue
            error("$macro_name: `$(t.sym)` looks like an opaque bound-form " *
                  "name; write it as decomposed call notation, e.g. `E(S)` " *
                  "or `E(A, B)`.")
        end
    end
end


"""
Walk a residual arithmetic expression (`A`, `A - P`, `S1 + S2 - P1 - P3`, etc.)
and classify each metabolite `Symbol` as added (positive) or subtracted (negative).
"""
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
- A name listed both as a substrate or product and in `regulators:` binds as
  the substrate or product in a bare `E(X)`; its inhibitor form is written
  `E(X::Inh)`.
- Same-kinetics groups are expressed via parenthesized step-groups; no
  `constraints:` block needed.
- Allosteric-only constructs (`regulatory_site(...)` / `::Tag` /
  `allosteric_regulators:` / `catalytic_inhibitors:`) are rejected.

Species notation on step sides:

- Bare Symbol that matches a declared metabolite → that metabolite.
- Bare Symbol otherwise (e.g. `E`, `Estar`, `E_c`) → conformation-only
  species named after the Symbol. A bare multi-capital name such as `ES`
  is rejected as an opaque bound-form name; write `E(S)`.
- `E(S)` / `E(S, P)` → species with conformation `:E` and bound
  metabolites; synthesized name is `:E<bound...>` (the conformation
  followed by the bound names, sorted alphabetically, with no separator:
  `:ES`, `:EPS`), matching `name(::Species)`.
- `Estar(; residual = A - P)` → species with empty bound and a residual
  recording `+A` / `−P`; synthesized name is `:Estar_res_+A_-P`.
- `Estar(B; residual = A - P)` → bound + residual; name
  `:EstarB_res_+A_-P`.

Conformation labels cannot shadow declared metabolite names.
"""
macro enzyme_mechanism(block)
    return esc(_parse_plain_mechanism_body(block)[1])
end

"""
Parse a labeled-line, returning `(label, values_vector)`. `values_vector` is the
list of args after the label, in source order. Each value is either a bare
Symbol or an `Expr(:(::), name, tag)` (for tagged lists, allosteric only).
Handles both Julia parse shapes:
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

function _parse_plain_mechanism_body(block)
    subs_list, prods_list, regs_list = Symbol[], Symbol[], Symbol[]
    steps_block = nothing
    for arg in block.args
        arg isa LineNumberNode && continue
        label, values = _parse_labeled_line(arg)
        if label == :substrates
            append!(subs_list, _bare_symbols_from_values(values, label))
        elseif label == :products
            append!(prods_list, _bare_symbols_from_values(values, label))
        elseif label == :regulators
            append!(regs_list, _bare_symbols_from_values(values, label))
        elseif label == :steps
            steps_block = only(values)
        else
            error("@enzyme_mechanism: unknown label `$label:`. Allosteric " *
                  "declarations (`allosteric_regulators:`, `catalytic_steps:`, " *
                  "`regulatory_site(...)`, ...) belong in @allosteric_mechanism.")
        end
    end
    isempty(subs_list) && error("substrates: not specified")
    isempty(prods_list) && error("products: not specified")
    steps_block === nothing && error("steps: not specified")

    declared_mets = Set{Symbol}(subs_list) ∪ Set{Symbol}(prods_list) ∪
                    Set{Symbol}(regs_list)
    # A metabolite that is both a substrate/product and its own competitive
    # inhibitor takes the substrate/product role for a bare `E(X)` binding; its
    # inhibitor form is written `E(X::Inh)`.
    role_of = Dict{Symbol,Symbol}()
    for r in regs_list;  role_of[r] = :CompetitiveInhibitor; end
    for s in subs_list;  role_of[s] = :Substrate;            end
    for p in prods_list; role_of[p] = :Product;              end

    side_terms_per_step =
        _parse_steps_block_with_groups(steps_block, declared_mets)

    _reject_opaque_bound_forms(side_terms_per_step, "@enzyme_mechanism")
    reaction_expr, groups_expr = _build_mechanism_expr(
        subs_list, prods_list, regs_list, role_of, side_terms_per_step)
    mech_expr = :(EnzymeRates.EnzymeMechanism(
        EnzymeRates.Mechanism($reaction_expr, $groups_expr)))
    # Second value is the SOURCE-order step groups (before the constructor
    # canonicalizes), for positional-oracle tests to bridge as-written step
    # indices to canonical stored order.
    (mech_expr, groups_expr)
end

"""
Build the `(reaction_expr, grouped_steps_expr)` pair from the structural
per-step records collected during parsing. `@enzyme_mechanism` wraps these
into `EnzymeMechanism(Mechanism(reaction, grouped_steps))`;
`@allosteric_mechanism` feeds the SAME source-order groups into
`AllostericMechanism`, which canonicalizes catalytic steps and allosteric
states together.
"""
function _build_mechanism_expr(subs_list, prods_list, regs_list,
                               role_of::Dict{Symbol,Symbol},
                               side_terms_per_step)
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

    (_mechanism_reaction_expr(subs_list, prods_list, regs_list), groups_expr)
end

"""
Build a `Step(from_species, to_species, consumed, released, is_eq)` `Expr`
from one step's LHS/RHS structural terms. Each side has exactly one
enzyme-form term (bare conformation OR call-form) and any number of
metabolite terms: the left-hand metabolites are consumed, the right-hand
ones released.
"""
function _build_step_expr(lhs::Vector{_StepSideTerm},
                          rhs::Vector{_StepSideTerm},
                          is_eq::Bool,
                          role_of::Dict{Symbol,Symbol})
    lhs_enzyme, lhs_mets = _split_side(lhs)
    rhs_enzyme, rhs_mets = _split_side(rhs)
    met_exprs(ts) = Expr[_metabolite_expr(t.sym, role_of, t.role) for t in ts]
    from_expr = _species_expr_from_term(lhs_enzyme, role_of)
    to_expr   = _species_expr_from_term(rhs_enzyme, role_of)
    :(EnzymeRates.Step($from_expr, $to_expr,
                       EnzymeRates.Metabolite[$(met_exprs(lhs_mets)...)],
                       EnzymeRates.Metabolite[$(met_exprs(rhs_mets)...)], $is_eq))
end

"""
Split a step side into its `(enzyme_term, metabolite_terms)`.
Errors if there is not exactly one enzyme term.
"""
function _split_side(side::Vector{_StepSideTerm})
    enzyme_term = nothing
    met_terms = _StepSideTerm[]
    for t in side
        if t.kind === :metabolite
            push!(met_terms, t)
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
    enzyme_term, met_terms
end

"""
Build a `Species(bound, conformation, residual)` `Expr` from an enzyme-form
`_StepSideTerm` (either bare conformation or Call-form).
"""
function _species_expr_from_term(t::_StepSideTerm,
                                 role_of::Dict{Symbol,Symbol})
    bound_entries = Expr[
        _metabolite_expr(b, role_of, r)
        for (b, r) in zip(t.bound, t.bound_roles)
    ]
    bound_expr = :(EnzymeRates.Metabolite[$(bound_entries...)])
    for (names, role) in ((t.residual_added, :Substrate),
                          (t.residual_subtracted, :Product)), n in names
        role_of[n] === role ||
            error("@enzyme_mechanism: a residual adds substrates and subtracts " *
                  "products; got `$n`.")
    end
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

"""
Build an `Expr` that constructs the appropriate `Metabolite` subtype for
a declared name. The role is looked up from `role_of`.
"""
function _metabolite_expr(name::Symbol, role_of::Dict{Symbol,Symbol},
                          override::Symbol = :default)
    if override === :inh
        return :(EnzymeRates.CompetitiveInhibitor($(QuoteNode(name))))
    end
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
    all(v -> v isa Symbol, values) ||
        error("`$label:` expects bare names; got $(join(values, ", ")); " *
              "atom brackets belong in @enzyme_reaction.")
    Symbol[values...]
end

"""
Parse the steps block. Each top-level expression is either:
  - `Expr(:(::), Expr(:tuple, step1, step2, ...), Tag)` — parenthesized group with tag
    (allosteric only).
  - `Expr(:tuple, step1, step2, ...)` — parenthesized group with no tag (plain mech).
  - `Expr(:call, ⇌|<-->, lhs, Expr(:(::), rhs, Tag))` — single tagged step (allosteric).
  - `Expr(:call, ⇌|<-->, lhs, rhs)` — single untagged step (plain).

Returns a Vector of structural per-step records
`(gnum, lhs_terms, rhs_terms, is_eq)` used by the emission decision logic
in `_parse_plain_mechanism_body` / `_parse_allosteric_mechanism_body`.
With `allow_tag=false` (plain mechanism), reject any `::Tag` annotations.
With `allow_tag=true` (allosteric mechanism), also return collected
`gnum => tag` pairs.
"""
function _parse_steps_block_with_groups(steps_block, declared_mets::Set{Symbol};
                                        allow_tag::Bool=false)
    next_group = Ref(0)
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
                push!(side_terms_per_step,
                      _step_struct_info(step_expr, gnum, declared_mets))
            end
        # Parenthesized-group-without-tag (plain)
        elseif arg isa Expr && arg.head == :tuple
            allow_tag &&
                error("@allosteric_mechanism: parenthesized step group " *
                      "`$(arg)` is missing `:: <:OnlyA|:EqualAI|:NonequalAI>` " *
                      "annotation. Add `:: <state>` after the closing paren.")
            next_group[] += 1
            gnum = next_group[]
            for step_expr in arg.args
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
                      "`:: <:OnlyA|:EqualAI|:NonequalAI>` annotation. Add " *
                      "`:: <state>` after the step expression.")
            end
            push!(side_terms_per_step,
                  _step_struct_info(arg, gnum, declared_mets))
        else
            error("Expected step or step-group; got $arg")
        end
    end

    if allow_tag
        return tags, side_terms_per_step
    else
        return side_terms_per_step
    end
end

"""
Return the structural per-step record `(gnum, lhs_terms, rhs_terms, is_eq)`
for a single (possibly already-de-tagged) step expression. Each side is
decomposed into a `Vector{_StepSideTerm}`.
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
    @allosteric_mechanism begin
        substrates: F6P
        products:   F16BP
        catalytic_multiplicity: 2
        allosteric_regulators: A::OnlyA, I::OnlyI

        catalytic_steps: begin
            E + F6P ⇌ E(F6P)        :: EqualAI
            E(F6P) <--> E(F16BP)    :: EqualAI
            (E(F16BP) ⇌ E + F16BP)  :: EqualAI
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
  `OnlyA`, `OnlyI`, `EqualAI`, `NonequalAI`.
- `catalytic_multiplicity: N` is the subunit count for the catalytic site
  (default 1).
- `catalytic_steps: begin ... end` is required (exactly once); each step or
  parenthesized step-group must carry a `::Tag`, one of `OnlyA`, `EqualAI` or
  `NonequalAI` (`OnlyI` is rejected on catalytic groups). Function-
  call species notation (`E(F6P)`, `Estar(B; residual = A - P)`) is supported.
- `regulatory_site(multiplicity = N): begin ligands: L1, L2 end` declares one
  regulatory site per block with multiplicity `N` and the ligands listed
  inside. Each ligand must appear in `allosteric_regulators:`; a ligand may
  not appear in two sites. Ligands declared in `allosteric_regulators:` but
  not assigned to any `regulatory_site(...):` block default to a single-ligand
  site at the catalytic multiplicity.
- A name with several roles binds in a bare catalytic-step `E(X)` as its
  substrate or product role first, then as a catalytic inhibitor; an
  allosteric regulator binds only at its regulatory site. A catalytic
  inhibitor's form is written `E(X::Inh)`.
"""
macro allosteric_mechanism(block)
    return esc(_parse_allosteric_mechanism_body(block)[1])
end

const _ALLOSTERIC_REG_STATES = Set([:OnlyA, :OnlyI, :EqualAI, :NonequalAI])

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
                  "::Tag annotations (e.g., I::OnlyI); got $v")
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
        mult = _positive_int(kw.args[2],
                             "@allosteric_mechanism: `regulatory_site` multiplicity")
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
Build the `Vector{RegulatorySite}` expression for the `AllostericMechanism`
constructor. Each site wraps its ligands (as `AllostericRegulator`s) with a
dense `Vector{Symbol}` of per-ligand allosteric states. Ligands not assigned
to any explicit `regulatory_site(...):` block become their own single-ligand
site at multiplicity `cat_n`.
"""
function _build_reg_sites_expr(allo_regs, reg_site_specs, cat_n)
    tag_of = Dict{Symbol,Symbol}(allo_regs)
    explicit = Set{Symbol}()
    for (_, ligs) in reg_site_specs, l in ligs
        l in explicit && error("@allosteric_mechanism: ligand $l " *
                               "appears in multiple regulatory sites")
        haskey(tag_of, l) ||
            error("@allosteric_mechanism: ligand $l on a " *
                  "`regulatory_site` is not declared in " *
                  "`allosteric_regulators:`")
        push!(explicit, l)
    end
    sites = [reg_site_specs; [(cat_n, [n]) for (n, _) in allo_regs if n ∉ explicit]]
    entries = map(sites) do (mult, ligs)
        :(EnzymeRates.RegulatorySite(
            EnzymeRates.AllostericRegulator[
                $((:(EnzymeRates.AllostericRegulator($(QuoteNode(l)))) for l in ligs)...)],
            $mult, Symbol[$((QuoteNode(tag_of[l]) for l in ligs)...)]))
    end
    :(EnzymeRates.RegulatorySite[$(entries...)])
end

"""
Build the `cat_allo_states` `Vector{Symbol}` expression for the
`AllostericMechanism` constructor: a dense vector with one entry per
catalytic kinetic group in source order (default `:NonequalAI`).
"""
function _build_cat_allo_states_expr(group_tags)
    for (_, tag) in group_tags
        tag in _ALLOSTERIC_REG_STATES ||
            error("@allosteric_mechanism: catalytic step tag :$tag not in " *
                  "($(_format_state_set(_ALLOSTERIC_REG_STATES)))")
    end
    tag_of = Dict{Int,Symbol}(group_tags)
    n_groups = isempty(group_tags) ? 0 : maximum(g for (g, _) in group_tags)
    :(Symbol[$((QuoteNode(get(tag_of, g, :NonequalAI)) for g in 1:n_groups)...)])
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
            cat_n = _positive_int(only(values),
                                  "@allosteric_mechanism: `catalytic_multiplicity:`")
        elseif label == :catalytic_steps
            cat_steps_block === nothing ||
                error("@allosteric_mechanism: multiple " *
                      "`catalytic_steps:` blocks.")
            cat_steps_block = only(values)
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

    # Order matters: a metabolite that is both a substrate/product and its own
    # competitive inhibitor (self-inhibition) takes the substrate/product role
    # for a bare `E(X)` binding; its inhibitor form is written `E(X::Inh)`. An
    # allosteric regulator binds only at its regulatory site, so every other
    # role of the same name wins in catalytic steps.
    role_of = Dict{Symbol,Symbol}()
    for (r, _) in allo_regs;  role_of[r] = :AllostericRegulator;  end
    for i in cat_inhibitors;  role_of[i] = :CompetitiveInhibitor; end
    for s in subs_list;       role_of[s] = :Substrate;            end
    for p in prods_list;      role_of[p] = :Product;              end

    group_tags, side_terms_per_step = _parse_steps_block_with_groups(
        cat_steps_block, declared_mets; allow_tag=true,
    )

    _reject_opaque_bound_forms(side_terms_per_step, "@allosteric_mechanism")
    reaction_expr, groups_expr = _build_mechanism_expr(
        subs_list, prods_list, cat_inhibitors, role_of, side_terms_per_step)

    cat_allo_states_expr = _build_cat_allo_states_expr(group_tags)
    reg_sites_expr = _build_reg_sites_expr(allo_regs, reg_site_specs, cat_n)

    # Route through AllostericMechanism so catalytic steps and their
    # allosteric-state tags canonicalize together; the lift back to
    # AllostericEnzymeMechanism keeps that alignment in the singleton.
    mech_expr = :(EnzymeRates.AllostericEnzymeMechanism(
        EnzymeRates.AllostericMechanism(
            $reaction_expr, $groups_expr, $cat_allo_states_expr,
            $cat_n, $reg_sites_expr)))
    # Extra values are SOURCE-order catalytic step groups and regulatory
    # sites (before the constructor canonicalizes), for positional-oracle
    # tests to bridge as-written indices to canonical stored order.
    (mech_expr, groups_expr, reg_sites_expr)
end
