# ABOUTME: @enzyme_reaction, which builds an EnzymeReaction, and the reader of its block:
# ABOUTME: reactants with their atom brackets, regulators, oligomeric states, shared sites.

# ─── @enzyme_reaction ───────────────────────────────────────────────────

"""
    @enzyme_reaction begin
        substrates: S[C6H12O6], ...
        products: P[C6H13O9P], ...
        competitive_inhibitors: I, ...
        allosteric_regulators: A, ...
        oligomeric_state: 2
        shared_catalytic_site: (S, P), ...
    end

Build the [`EnzymeReaction`](@ref EnzymeRates.Mechanisms.EnzymeReaction) that the
mechanism search starts from: its reactants with their atoms, and the regulators and
oligomeric states the search may use.

- `substrates:` and `products:` (required) list each reactant with its atoms, one
  chemical formula in brackets: `S[C6H12O6]`. The atoms of the two sides balance, and a
  name is a substrate or a product, not both.
- `competitive_inhibitors:` lists the competitive inhibitors, names that bind the
  catalytic site without reacting. A substrate or product listed here may also bind
  where it does not react, as an inhibitor.
- `allosteric_regulators:` lists the allosteric regulators, names that bind regulatory
  sites. A type tag, `A::Activator` or `A::Inhibitor`, pins the regulator to the
  allosteric state that activates or inhibits; an untagged one (`:unspecified`) may take
  either. Only these entries take a type tag. A name may be both a competitive
  inhibitor and an allosteric regulator.
- `oligomeric_state:` gives the oligomeric states the search may assign: one
  positive Int, or a tuple of distinct ones such as `(1, 2, 4)`. It defaults to 1 and is
  given at most once.
- `shared_catalytic_site:` lists `(substrate, product)` pairs, in either order, that bind
  one catalytic site, so never both at once.
- The list labels may repeat; their entries accumulate.

No metabolite is named `Keq`, `L`, `N`, `params`, `concs` or `v`, or holds `_`: a
rate equation keeps those names, and its constant names join names with `_`. Nor
is one named `E`, the name of the free enzyme form; a name is an identifier that is no
keyword. The mechanism search names an enzyme form by joining its bound metabolites'
names without a separator, so the names may not spell one string two ways: with `Ac`,
`CoA` and `AcCoA`, E with Ac and CoA bound and E with AcCoA bound would both be
`EAcCoA`. The search rejects such a reaction when it starts and names both spellings
(`Ac + CoA and AcCoA both spell AcCoA`); rename one (`Acetyl`). A reactant bound where it
does not react is written `Xmod` in a form's name, which counts as a name too.

A reaction prints as the block that builds it, its names sorted:

```jldoctest
julia> using EnzymeRates

julia> @enzyme_reaction begin
           substrates: S[C6H12O6], ATP[C10H16N5O13P3]
           products: P[C6H13O9P], ADP[C10H15N5O10P2]
           competitive_inhibitors: P
           allosteric_regulators: A::Activator, I::Inhibitor
           oligomeric_state: (1, 2)
       end
@enzyme_reaction begin
    substrates: ATP[C10H16N5O13P3], S[C6H12O6]
    products: ADP[C10H15N5O10P2], P[C6H13O9P]
    competitive_inhibitors: P
    allosteric_regulators: A::Activator, I::Inhibitor
    oligomeric_state: (1, 2)
end
```
"""
macro enzyme_reaction(block)
    return reaction_expr(block)
end

# ─── Reaction macro blocks as written ───────────────────────────────────

"""
Parse the `@enzyme_reaction` body. Returns a `NamedTuple`:
- `reactants ::Vector{Tuple{Symbol, Symbol, Vector{Pair{Symbol,Int}}}}` —
  `(role, name, atoms)` with `role ∈ (:substrate, :product)` and `atoms` the
  `element => count` pairs.
- `regs  ::Vector{Tuple{Symbol, Symbol, Symbol}}` —
  `(kind, name, effect)` where `kind ∈ (:competitive, :allosteric)` and
  `effect ∈ (:activator, :inhibitor, :unspecified)`.
- `mults ::Vector{Int}` — allowed oligomeric states (default `[1]`).
- `shared ::Vector{Tuple{Symbol,Symbol}}` — `shared_catalytic_site:` pairs.
"""
function _parse_reaction_block(block)
    block isa Expr && block.head === :block ||
        error("@enzyme_reaction: expected a `begin ... end` block, got $block")
    reactants = Tuple{Symbol, Symbol, Vector{Pair{Symbol,Int}}}[]
    regs  = Tuple{Symbol, Symbol, Symbol}[]
    mults::Union{Nothing, Vector{Int}} = nothing
    shared = Tuple{Symbol,Symbol}[]

    for arg in block.args
        arg isa LineNumberNode && continue
        label, values = _parse_labeled_line(arg)
        if label === :substrates || label === :products
            role = label === :substrates ? :substrate : :product
            for (name, atoms) in _parse_atom_bracket_entries(values, label)
                push!(reactants, (role, name, atoms))
            end
        elseif label === :competitive_inhibitors
            append!(regs, _parse_regulator_entries(values, :competitive))
        elseif label === :allosteric_regulators
            append!(regs, _parse_regulator_entries(values, :allosteric))
        elseif label === :oligomeric_state
            mults === nothing ||
                error("@enzyme_reaction: `oligomeric_state:` given more than once.")
            mults = _parse_multiplicity_tuple(values, label)
        elseif label === :shared_catalytic_site
            append!(shared, _parse_shared_site_pairs(values))
        else
            error("@enzyme_reaction: unknown label `$label:`. Valid labels: " *
                  "substrates, products, competitive_inhibitors, " *
                  "allosteric_regulators, oligomeric_state, shared_catalytic_site.")
        end
    end
    (; reactants, regs, mults = something(mults, [1]), shared)
end

"""
Parse `S[C6H12O6]` entries, each a name and one chemical formula in brackets, into
`(name, atoms)` pairs, `atoms` a `Vector{Pair{Symbol,Int}}` of `element => count`.
"""
function _parse_atom_bracket_entries(values, label)
    map(values) do v
        v isa Expr && v.head === :ref && length(v.args) == 2 && v.args[1] isa Symbol ||
            error("@enzyme_reaction: `$label:` lists each reactant with one chemical " *
                  "formula in brackets, such as `S[C6H12O6]`; got `$v`.")
        (v.args[1]::Symbol, _parse_chemical_formula(string(v.args[2])))
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

"""
Parse `R` entries, optionally tagged with a regulator type, `R::Activator` or
`R::Inhibitor`, into `(kind, name, effect)`; an untagged `R` has effect
`:unspecified`. A multiplicity list such as `R(1, 2)` is an error: a regulator takes the
oligomeric states that `oligomeric_state:` lists. `Regulator` rejects a type tag on a
competitive inhibitor.
"""
function _parse_regulator_entries(values, kind::Symbol)
    map(values) do v
        effect, inner = :unspecified, v
        if v isa Expr && v.head === :(::)
            inner, tag = v.args
            effect = tag === :Activator ? :activator :
                   tag === :Inhibitor ? :inhibitor :
                   error("@enzyme_reaction: unknown regulator type tag ::$tag " *
                         "on $v; use ::Activator or ::Inhibitor.")
        end
        inner isa Symbol && return (kind, inner, effect)
        inner isa Expr && inner.head === :call && length(inner.args) >= 2 &&
            inner.args[1] isa Symbol && all(a -> a isa Integer, inner.args[2:end]) &&
            error("@enzyme_reaction: regulator $(inner.args[1]) lists multiplicities, " *
                  "got `$v`; a regulator takes the oligomeric states that " *
                  "`oligomeric_state:` lists, so list $(inner.args[1]) bare.")
        error("@enzyme_reaction: cannot parse regulator entry $v. Use `R`.")
    end
end

"""Parse `(1, 2, 4)` or `4` into `Vector{Int}`."""
function _parse_multiplicity_tuple(values, label)
    v = values[1]
    length(values) == 1 && v != Expr(:tuple) ||
        error("@enzyme_reaction: `$label:` takes one Int or a tuple of one or more " *
              "Ints, such as `(1, 2, 4)`; got $(join(values, ", ")).")
    Int[_positive_int(a, "@enzyme_reaction: `$label:` entry")
        for a in (v isa Expr && v.head === :tuple ? v.args : (v,))]
end

"""
    reaction_expr(block) -> Expr

The `EnzymeReaction` `Expr` that `@enzyme_reaction`'s `block` writes, from its parsed
declarations (`_parse_reaction_block`): `reactants` as `(role, name, atoms)` and
`regs` as `(kind, name, effect)`. The macro returns it unescaped, so hygiene resolves its
names in `DSL` itself and the macro works in a module that imports the package under
another name or imports only the macro.
"""
function reaction_expr(block)
    (; reactants, regs, mults, shared) = _parse_reaction_block(block)
    reactant_exprs(role) = [
        :(Mechanisms.Reactant($(QuoteNode(name)),
            Pair{Symbol,Int}[$((:($(QuoteNode(e)) => $c) for (e, c) in atoms)...)]))
        for (r, name, atoms) in reactants if r === role]
    regulator_exprs = [
        :(Mechanisms.Regulator(
            $(QuoteNode(name)), $(QuoteNode(kind)), $(QuoteNode(effect))))
        for (kind, name, effect) in regs]
    shared_exprs = [:(($(QuoteNode(a)), $(QuoteNode(b)))) for (a, b) in shared]
    :(Mechanisms.EnzymeReaction(
        Mechanisms.Reactant[$(reactant_exprs(:substrate)...)],
        Mechanisms.Reactant[$(reactant_exprs(:product)...)],
        Mechanisms.Regulator[$(regulator_exprs...)], Int[$(mults...)];
        shared_catalytic_site = Tuple{Symbol,Symbol}[$(shared_exprs...)]))
end
