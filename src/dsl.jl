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
  A name may appear under only one of the two labels.
  An enzyme form's name joins its bound metabolites' names without a
  separator, so two forms can render one name: with `Ac`, `CoA` and
  `AcCoA`, E with Ac and CoA bound and E with AcCoA bound are both
  `EAcCoA`. `@enzyme_reaction` accepts such names; the mechanism
  constructors reject a mechanism that holds two forms whose names
  coincide. For Ac + CoA ⇌ AcCoA every seed holds both E(Ac, CoA) and
  E(AcCoA), so `init_mechanisms` and `identify_rate_equation` stop with an
  error naming both forms. Rename one (`Acetyl`).
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
        elseif label === :allowed_catalytic_multiplicities || label === :oligomeric_state
            mults === nothing ||
                error("@enzyme_reaction: `$label:` sets the catalytic multiplicities, " *
                      "which an earlier `allowed_catalytic_multiplicities:` or " *
                      "`oligomeric_state:` line already set.")
            mults = _parse_multiplicity_tuple(values, label)
            label === :oligomeric_state && length(mults) != 1 &&
                error("@enzyme_reaction: `oligomeric_state:` takes a single Int.")
        elseif label === :shared_catalytic_site
            append!(shared, _parse_shared_site_pairs(values))
        else
            error("@enzyme_reaction: unknown label `$label:`. Valid labels: " *
                  "substrates, products, competitive_inhibitors, dead_end_inhibitors, " *
                  "allosteric_regulators, allowed_catalytic_multiplicities, " *
                  "oligomeric_state, shared_catalytic_site.")
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
        error("@enzyme_reaction: `$label:` takes a single " *
              "$(label === :oligomeric_state ? "Int" : "tuple"), " *
              "got $(join(values, ", ")).")
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
  `E(X::Inh)`. `X::Inh` requires `X` in `regulators:`.
- Same-kinetics groups are expressed via parenthesized step-groups; no
  `constraints:` block needed.
- Allosteric-only constructs (`regulatory_site(...)` / `::Tag` /
  `allosteric_regulators:` / `catalytic_inhibitors:`) are rejected.

Species notation on step sides:

- Bare Symbol that matches a declared metabolite → that metabolite.
- Bare Symbol otherwise (e.g. `E`, `Estar`, `E_c`) → conformation-only
  species named after the Symbol. A multi-capital conformation label such as `ES`,
  bare or as a call head (`ES(P)`), is rejected as an opaque bound-form name; write
  `E(S)`.
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
    return esc(_parse_mechanism_body(block, false)[1])
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
Parse the steps block into the `Vector{Vector{Step}}` `Expr` of its kinetic groups
and the groups' tags, both in source order. Each top-level expression is either:
  - `Expr(:(::), Expr(:tuple, step1, step2, ...), Tag)` — parenthesized group with tag
    (allosteric only); a parenthesized single step, `Expr(:(::), step, Tag)`, is a
    one-step group.
  - `Expr(:tuple, step1, step2, ...)` — parenthesized group with no tag (plain mech).
  - `Expr(:call, ⇌|<-->, lhs, Expr(:(::), rhs, Tag))` — single tagged step (allosteric).
  - `Expr(:call, ⇌|<-->, lhs, rhs)` — single untagged step (plain).

With `allow_tag=false` (plain mechanism), reject any `::Tag` annotations; the tags
come back empty. With `allow_tag=true` (allosteric mechanism), every group carries
one. `role_of` maps each declared metabolite to the `Metabolite` subtype it binds
as; `inhibitors` holds the `label` that declares competitive inhibitors and the
`names` it lists; `macro_name` names the invoking macro in error messages.
"""
function _parse_steps_block(steps_block, role_of, inhibitors, macro_name; allow_tag::Bool)
    groups, tags = Expr[], Symbol[]
    for arg in steps_block.args
        arg isa LineNumberNode && continue

        # Parenthesized-group-with-tag (allosteric); a lone parenthesized step is a group
        if arg isa Expr && arg.head == :(::) &&
           arg.args[1] isa Expr && arg.args[1].head in (:tuple, :call)
            allow_tag ||
                error("$macro_name: tag annotation `$arg` is not allowed")
            inner, tag = arg.args
            steps = inner.head == :tuple ? inner.args : Any[inner]
            tag isa Symbol ||
                error("$macro_name: step-group tag must be a Symbol; got $tag")
            push!(tags, tag)
        # Parenthesized-group-without-tag (plain)
        elseif arg isa Expr && arg.head == :tuple
            allow_tag &&
                error("$macro_name: parenthesized step group `$(arg)` is missing " *
                      "`:: <:OnlyA|:EqualAI|:NonequalAI>` annotation. Add " *
                      "`:: <state>` after the closing paren.")
            steps = arg.args
        # Single step (with or without tag): an operator call with two sides
        elseif arg isa Expr && arg.head == :call && length(arg.args) == 3
            original = string(arg)
            tag = _peel_step_tag!(arg)
            if tag !== nothing
                tag isa Symbol ||
                    error("$macro_name: step tag must be a Symbol; got $tag")
                allow_tag ||
                    error("$macro_name: tag annotation on `$original` is not allowed")
                push!(tags, tag)
            elseif allow_tag
                error("$macro_name: step `$(original)` is missing " *
                      "`:: <:OnlyA|:EqualAI|:NonequalAI>` annotation. Add " *
                      "`:: <state>` after the step expression.")
            end
            steps = Any[arg]
        else
            error("$macro_name: expected step or step-group; got $arg")
        end
        push!(groups, :(EnzymeRates.Step[
            $((_step_expr(s, role_of, inhibitors, macro_name) for s in steps)...)]))
    end
    :(Vector{EnzymeRates.Step}[$(groups...)]), tags
end

"""
If the step Expr has a `::Tag` attached to its RHS arg, remove the wrapper and
return the tag. Otherwise return `nothing`. Mutates the step Expr. `::Inh`
marks a metabolite's inhibitor role, so it is never taken as a step tag.

Two RHS shapes carry a tag:
  - `Sym :: Tag` parses as `Expr(:(::), Sym, Tag)`.
  - `S1 + S2 + … + LastSym :: Tag` parses as
    `Expr(:call, :+, S1, …, Expr(:(::), LastSym, Tag))` because `::` binds
    tighter than `+`. We peel the inner `::` and put `LastSym` back in place.
"""
function _peel_step_tag!(step_expr)
    rhs = step_expr.args[3]
    holder, i = rhs isa Expr && rhs.head == :call && rhs.args[1] == :+ ?
                (rhs, lastindex(rhs.args)) : (step_expr, 3)
    t = holder.args[i]
    t isa Expr && t.head == :(::) && t.args[2] !== :Inh || return nothing
    holder.args[i] = t.args[1]
    t.args[2]
end

"""
Build the `Step(from_species, to_species, consumed, released, is_eq)` `Expr` for one
`lhs ⇌ rhs` (rapid-equilibrium) or `lhs <--> rhs` (steady-state) step. Each side has
exactly one enzyme-form term and any number of metabolite terms (a declared name or
`X::Inh`): the left-hand metabolites are consumed, the right-hand ones released.
Every term of both sides is parsed before either side's enzyme forms are counted,
so a misspelled metabolite bound on either side is reported as undeclared rather
than as a second enzyme form or an opaque bound-form name. The enzyme form's
conformation label, bare or a call head, must be conformation-shaped
(`_is_conformation_shape`).
"""
function _step_expr(expr, role_of, inhibitors, macro_name)
    expr isa Expr && expr.head == :call ||
        error("$macro_name: expected lhs ⇌ rhs or lhs <--> rhs; got $expr")
    op = expr.args[1]
    op == :⇌ || op == :(<-->) ||
        error("$macro_name: expected ⇌ or <--> step operator; got $op")
    is_met(t) = t isa Symbol ? haskey(role_of, t) : t isa Expr && t.head == :(::)
    term_name(t) = t isa Expr ? t.args[1] : t
    sides = map(expr.args[2:3]) do side
        terms = side isa Expr && side.head == :call && side.args[1] == :+ ?
                side.args[2:end] : Any[side]
        terms, Expr[is_met(t) ? _metabolite_expr(t, role_of, inhibitors, macro_name) :
                    _species_expr(t, role_of, inhibitors, macro_name) for t in terms]
    end
    (from, consumed), (to, released) = map(sides) do (terms, exprs)
        i = findall(!is_met, terms)
        isempty(i) && error("$macro_name: step side has no enzyme-form term " *
                            "(terms: $(Symbol[term_name(t) for t in terms])).")
        length(i) == 1 ||
            error("$macro_name: step side has more than one enzyme-form term " *
                  "($(term_name(terms[i[1]])), $(term_name(terms[i[2]]))); each " *
                  "elementary step has exactly one enzyme form per side.")
        conformation = term_name(terms[only(i)])
        _is_conformation_shape(conformation) ||
            error("$macro_name: `$conformation` looks like an opaque bound-form name; " *
                  "write it as decomposed call notation, e.g. `E(S)` or `E(A, B)`.")
        exprs[only(i)], exprs[eachindex(exprs) .!= only(i)]
    end
    :(EnzymeRates.Step($from, $to, EnzymeRates.Metabolite[$(consumed...)],
                       EnzymeRates.Metabolite[$(released...)], $(op == :⇌)))
end

"""
Build the `Species(bound, conformation, residual)` `Expr` for an enzyme-form term:
a bare conformation `E`, or a call `E(S, X::Inh; residual = A - P)` whose head is
the conformation, whose positional arguments are the bound metabolites, and whose
`residual` lists the substrates added to and the products removed from the enzyme.
Conformation labels cannot shadow declared metabolite names. The `Species` and
`Residual` constructors sort what they hold.
"""
function _species_expr(t, role_of, inhibitors, macro_name)
    t isa Symbol || t isa Expr && t.head == :call && t.args[1] isa Symbol ||
        error("$macro_name: expected metabolite Symbol or species expression on " *
              "step side; got $t")
    conformation, args = t isa Symbol ? (t, Any[]) : (t.args[1], t.args[2:end])
    haskey(role_of, conformation) &&
        error("$macro_name: conformation label `$conformation` collides with " *
              "declared metabolite `$conformation`; choose a different " *
              "conformation label.")
    bound, added, subtracted = Expr[], Expr[], Expr[]
    for a in args
        if a isa Expr && a.head === :parameters
            for kw in a.args
                kw isa Expr && kw.head === :kw && kw.args[1] === :residual ||
                    error("$macro_name: unknown keyword in species `$t`: $kw. " *
                          "Only `residual = ...` is allowed.")
                _walk_residual_expr(kw.args[2], true, added, subtracted, role_of,
                                    macro_name)
            end
        else
            a isa Symbol || a isa Expr && a.head === :(::) ||
                error("$macro_name: invalid entry in species `$t`: $a")
            push!(bound, _metabolite_expr(a, role_of, inhibitors, macro_name, t))
        end
    end
    :(EnzymeRates.Species(EnzymeRates.Metabolite[$(bound...)],
                          $(QuoteNode(conformation)),
                          EnzymeRates.Residual(EnzymeRates.Substrate[$(added...)],
                                               EnzymeRates.Product[$(subtracted...)])))
end

"""
A `Symbol` is "conformation-shaped" iff it starts with a single
capital letter followed by any mix of lowercase letters, digits, and
underscore-separated lowercase/digit runs: `:E`, `:Estar`, `:Estar2`, `:E_c`,
`:E_secondary`. Multi-capital `Symbol`s (`:ES`, `:EAB`) and underscore-then-
uppercase `Symbol`s (`:E_S`, `:Estar_A_B`) are opaque bound-form names —
rejected in favor of decomposed call notation.
"""
_is_conformation_shape(sym::Symbol) =
    occursin(r"^[A-Z][a-z0-9]*(_[a-z0-9]+)*$", String(sym))

"""
Build the `Metabolite` `Expr` for a declared name `X` (the subtype `role_of[X]`) or
`X::Inh` (its `CompetitiveInhibitor` copy): a free term on a step side or, given
`species`, a metabolite bound in that species. `X::Inh` requires `X` among
`inhibitors.names`, the competitive inhibitors the `inhibitors.label` line declares. A
name declared only as an allosteric regulator is therefore rejected in either
spelling: it binds only at its regulatory site.
"""
function _metabolite_expr(t, role_of, inhibitors, macro_name, species = nothing)
    name, tag = t isa Expr ? t.args : (t, nothing)
    in_species = species === nothing ? "" : " in species `$species`"
    name isa Symbol ||
        error("$macro_name: expected `name::Inh`$in_species; got $t")
    tag === nothing || tag === :Inh ||
        error("$macro_name: unknown role tag ::$tag on $name; only ::Inh is " *
              "supported.")
    haskey(role_of, name) ||
        error("$macro_name: " * (species === nothing ?
                  "tagged metabolite `$name` in `$t`" :
                  "bound metabolite `$name` in species `$species`") *
              " is not declared. Declared: $(sort(collect(keys(role_of)))).")
    tag === :Inh && name ∉ inhibitors.names &&
        error("$macro_name: `$name::Inh` names `$name`, which `$(inhibitors.label):` " *
              "does not declare; add `$name` to `$(inhibitors.label):` to let it bind " *
              "the catalytic site as a competitive inhibitor.")
    type = tag === :Inh ? :CompetitiveInhibitor : role_of[name]
    type === :AllostericRegulator &&
        error("$macro_name: `$name` is an allosteric regulator; it binds only at " *
              "its regulatory site.")
    :(EnzymeRates.$type($(QuoteNode(name))))
end

"""
Walk a residual arithmetic expression (`A`, `A - P`, `S1 + S2 - P1 - P3`, etc.)
and push each metabolite's `Expr` onto `added` (positive) or `subtracted`
(negative). A residual adds substrates and subtracts products.
"""
function _walk_residual_expr(e, sign_positive, added, subtracted, role_of, macro_name)
    walk(x, s) = _walk_residual_expr(x, s, added, subtracted, role_of, macro_name)
    if e isa Symbol
        haskey(role_of, e) ||
            error("$macro_name: residual entry `$e` is not a declared metabolite. " *
                  "Declared: $(sort(collect(keys(role_of)))).")
        type = sign_positive ? :Substrate : :Product
        role_of[e] === type ||
            error("$macro_name: a residual adds substrates and subtracts " *
                  "products; got `$e`.")
        push!(sign_positive ? added : subtracted, :(EnzymeRates.$type($(QuoteNode(e)))))
    elseif e isa Expr && e.head === :call && e.args[1] === :+ && length(e.args) >= 3
        foreach(a -> walk(a, sign_positive), e.args[2:end])
    elseif e isa Expr && e.head === :call && e.args[1] === :- && length(e.args) in (2, 3)
        length(e.args) == 3 && walk(e.args[2], sign_positive)
        walk(e.args[end], !sign_positive)
    else
        error("$macro_name: invalid residual expression: $e")
    end
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
  inside, in any order. Each ligand must appear in `allosteric_regulators:`; a
  ligand may not appear in two sites. Ligands declared in
  `allosteric_regulators:` but not assigned to any `regulatory_site(...):` block
  default to a single-ligand site at the catalytic multiplicity.
- A name with several roles binds in a bare catalytic-step `E(X)` as its
  substrate or product role first, then as a catalytic inhibitor; an
  allosteric regulator binds only at its regulatory site. A catalytic
  inhibitor's form is written `E(X::Inh)`, and `X::Inh` requires `X` in
  `catalytic_inhibitors:`; list an allosteric regulator there as well to let it
  bind the catalytic site.
"""
macro allosteric_mechanism(block)
    return esc(_parse_mechanism_body(block, true)[1])
end

"""
Parse a `regulatory_site(multiplicity = N): begin ligands: ... end` line, given its
`regulatory_site(...)` label call and its body. Returns
`(mult::Int, ligands::Vector{Symbol})`.
"""
function _parse_regulatory_site(label, body)
    kws = label.args[2:end]
    length(kws) == 1 && kws[1] isa Expr && kws[1].head == :kw &&
        kws[1].args[1] == :multiplicity ||
        error("@allosteric_mechanism: `regulatory_site` takes exactly one " *
              "`multiplicity = N` argument; got `$label`")
    mult = _positive_int(kws[1].args[2],
                         "@allosteric_mechanism: `regulatory_site` multiplicity")

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
Parse the body of a mechanism macro into `(mech_expr, groups_expr, reg_sites_expr)`:
the `Expr` that builds the mechanism, its kinetic groups and, for an allosteric
mechanism, its regulatory sites (`nothing` for a plain one). The groups and sites are
in SOURCE order, before the mechanism constructor canonicalizes them, for
positional-oracle tests to bridge as-written indices to canonical stored order; each
site's ligands are in name order, which the `RegulatorySite` constructor sets.
`allosteric` selects the grammar of `@allosteric_mechanism` (`catalytic_inhibitors:`,
`catalytic_steps:`, `allosteric_regulators:`, `catalytic_multiplicity:`,
`regulatory_site(...)`) over that of `@enzyme_mechanism` (`regulators:`, `steps:`).
"""
function _parse_mechanism_body(block, allosteric::Bool)
    macro_name = allosteric ? "@allosteric_mechanism" : "@enzyme_mechanism"
    block isa Expr && block.head === :block ||
        error("$macro_name: expected a `begin ... end` block, got $block")
    inhibitors_label, steps_label =
        allosteric ? (:catalytic_inhibitors, :catalytic_steps) : (:regulators, :steps)
    subs_list, prods_list, inhibitors = Symbol[], Symbol[], Symbol[]
    allo_regs = Pair{Symbol,Symbol}[]
    cat_n = nothing
    steps_block = nothing
    reg_site_specs = Tuple{Int,Vector{Symbol}}[]

    for arg in block.args
        arg isa LineNumberNode && continue
        label, values = _parse_labeled_line(arg)
        if label == :substrates
            append!(subs_list, _bare_symbols_from_values(values, label))
        elseif label == :products
            append!(prods_list, _bare_symbols_from_values(values, label))
        elseif label == inhibitors_label
            append!(inhibitors, _bare_symbols_from_values(values, label))
        elseif label == steps_label
            steps_block === nothing ||
                error("$macro_name: `$label:` given more than once.")
            steps_block = only(values)
        elseif allosteric && label == :allosteric_regulators
            for v in values
                v isa Expr && v.head == :(::) && all(a -> a isa Symbol, v.args) ||
                    error("@allosteric_mechanism `$label:` requires per-entry " *
                          "`name::Tag` annotations (e.g., I::OnlyI); got $v")
                push!(allo_regs, v.args[1] => v.args[2])
            end
        elseif allosteric && label isa Expr && label.head == :call &&
               label.args[1] == :regulatory_site
            push!(reg_site_specs, _parse_regulatory_site(label, only(values)))
        elseif allosteric && label == :catalytic_multiplicity
            cat_n === nothing ||
                error("@allosteric_mechanism: `catalytic_multiplicity:` given more " *
                      "than once.")
            cat_n = _positive_int(only(values),
                                  "@allosteric_mechanism: `catalytic_multiplicity:`")
        else
            error("$macro_name: unknown label `$label:`." * (allosteric ? "" :
                  " Allosteric declarations (`allosteric_regulators:`, " *
                  "`catalytic_steps:`, `regulatory_site(...)`, ...) belong in " *
                  "@allosteric_mechanism."))
        end
    end

    isempty(subs_list) && error("$macro_name: `substrates:` not specified.")
    isempty(prods_list) && error("$macro_name: `products:` not specified.")
    steps_block === nothing && error("$macro_name: `$steps_label:` not specified.")

    # Order matters: a metabolite that is both a substrate/product and its own
    # competitive inhibitor (self-inhibition) takes the substrate/product role
    # for a bare `E(X)` binding; its inhibitor form is written `E(X::Inh)`. An
    # allosteric regulator binds only at its regulatory site, so every other
    # role of the same name wins in catalytic steps.
    role_of = Dict{Symbol,Symbol}([first.(allo_regs) .=> :AllostericRegulator;
                                   inhibitors .=> :CompetitiveInhibitor;
                                   subs_list .=> :Substrate; prods_list .=> :Product])
    groups_expr, group_tags = _parse_steps_block(
        steps_block, role_of, (label = inhibitors_label, names = inhibitors), macro_name;
        allow_tag = allosteric)
    reaction_expr = _mechanism_reaction_expr(subs_list, prods_list, inhibitors)
    if allosteric
        allo_names = first.(allo_regs)
        repeated = unique(n for (i, n) in enumerate(allo_names)
                          if n in view(allo_names, 1:i-1))
        isempty(repeated) ||
            error("@allosteric_mechanism: `allosteric_regulators:` lists " *
                  "$(join(("`$n`" for n in repeated), ", ")) more than once.")
        cat_allo_states_expr = :(Symbol[$(QuoteNode.(group_tags)...)])
        cat_n = something(cat_n, 1)
        reg_sites_expr = _build_reg_sites_expr(allo_regs, reg_site_specs, cat_n)
        # Route through AllostericMechanism so catalytic steps and their
        # allosteric-state tags canonicalize together; the lift back to
        # AllostericEnzymeMechanism keeps that alignment in the singleton.
        mech_expr = :(EnzymeRates.AllostericEnzymeMechanism(
            EnzymeRates.AllostericMechanism(
                $reaction_expr, $groups_expr, $cat_allo_states_expr, $cat_n,
                $reg_sites_expr)))
    else
        reg_sites_expr = nothing
        mech_expr = :(EnzymeRates.EnzymeMechanism(
            EnzymeRates.Mechanism($reaction_expr, $groups_expr)))
    end
    (mech_expr, groups_expr, reg_sites_expr)
end
