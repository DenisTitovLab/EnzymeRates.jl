# ABOUTME: DSL macros @enzyme_reaction, @enzyme_mechanism, @allosteric_mechanism.
# ABOUTME: Parse reaction/mechanism declarations into concrete struct constructors.

"""
    @enzyme_reaction begin
        substrates: S[C6H12O6], ATP[C10H16N5O13P3]
        products:   P[C6H13O9P], ADP[C10H15N5O10P2]
        competitive_inhibitors: I            # bare OK; mults default to catalytic
        allosteric_regulators: A(1, 2, 4)    # bare OK; mults default to catalytic
        oligomeric_state: (1, 2, 4)          # or one Int: oligomeric_state: 2
    end

Emit an `EnzymeReaction`.

- `substrates:` / `products:` — comma-separated entries with required atom
  brackets (`S[C6H12O6]`). Multi-atom forms like `[C2,N]` are allowed.
  A name may appear under only one of the two labels.
  An enzyme form's name joins its bound metabolites' names without a
  separator, so the names may not spell one string two ways: with `Ac`,
  `CoA` and `AcCoA`, E with Ac and CoA bound and E with AcCoA bound would
  both be `EAcCoA`. `@enzyme_reaction` rejects such a reaction and names
  both spellings (`Ac + CoA and AcCoA both spell AcCoA`); rename one
  (`Acetyl`). A competitive inhibitor that is also a substrate or product
  binds as its modifier copy, named `Xmod`, which counts as a name too.
- No metabolite is named `E_total`, `Keq`, `L`, `N`, `params`, `concs` or `v`,
  or holds `_`: a rate equation keeps those names, and its constant names join
  names with `_`.
- `competitive_inhibitors:` —
  `CompetitiveInhibitor` entries (catalytic-site binding). May be bare
  `I` (multiplicities default to the `oligomeric_state:` ones) or
  `I(m1, m2, ...)` to override. Enumeration reads only the mechanism's
  catalytic multiplicity; the per-inhibitor list is stored but not read.
- `allosteric_regulators:` — `AllostericRegulator` entries. May be bare
  `A` (multiplicities default to the `oligomeric_state:` ones),
  `A(m1, m2, ...)` or a single value `A(m)`. Enumeration gives every
  regulatory site the mechanism's catalytic multiplicity; the per-regulator
  list is stored but not read. An entry may carry a type
  tag, `A::Activator` or `A::Inhibitor`; untagged entries default to
  `:unspecified`. Type tags are only valid on `allosteric_regulators:`
  entries.
- `oligomeric_state:` — the catalytic multiplicities the search may assign: one positive
  Int, or a tuple of them such as `(1, 2, 4)` (default 1).
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
        elseif label === :competitive_inhibitors
            append!(regs, _parse_regulator_entries(values, :CompetitiveInhibitor))
        elseif label === :allosteric_regulators
            append!(regs, _parse_regulator_entries(values, :AllostericRegulator))
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
`(type, name, mults, reg_type)`. Bare `R` produces `nothing` mults (filled from the
`oligomeric_state:` multiplicities, default `[1]`) and reg_type `:unspecified`.
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
        error("@enzyme_reaction: `$label:` takes one Int or one tuple such as " *
              "`(1, 2, 4)`; got $(join(values, ", ")).")
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
    @enzyme_mechanism begin
        substrates: S
        products:   P
        modifiers: I

        steps: begin
            E + S ⇌ E(S)                   # function-call species notation
            E(S) <--> E(P)
            E(P) ⇌ E + P
            # parenthesized → shared kinetics: dead-end I binds E and E(S) alike
            (E + I ⇌ E(I), E(S) + I ⇌ E(S, I))
        end
    end

Build a plain (non-allosteric) `EnzymeMechanism`.

- `substrates:`, `products:`, `modifiers:` accept comma-separated bare
  symbols. Atom brackets (`S[C]`) are rejected at the mechanism level.
- `modifiers:` lists names that bind the catalytic graph without reacting; they
  are treated as `CompetitiveInhibitor`s when later passed to `EnzymeReaction`.
- A name listed both as a substrate or product and in `modifiers:` binds as
  the substrate or product in a bare `E(X)`; its modifier form is written
  `E(X::Mod)`. `X::Mod` requires `X` in `modifiers:`; a name that is only a
  modifier is written bare.
- Same-kinetics groups are expressed via parenthesized step-groups; no
  `constraints:` block needed.
- Allosteric-only constructs (`regulatory_site(n):` / `::Tag` /
  `allosteric_regulators:` / `catalytic_multiplicity:`) are rejected.

Species notation on step sides:

- Bare Symbol that matches a declared metabolite → that metabolite.
- Bare Symbol otherwise (e.g. `E`, `Estar`, `F2`) → conformation-only species
  named after the Symbol; a conformation label is a capital letter followed by
  lowercase letters and digits. A multi-capital conformation label such as `ES`,
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
    return :($_compiled_mechanism($(_raw_mechanism(block, false))))
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
    @allosteric_mechanism begin
        substrates: F6P
        products:   F16BP
        catalytic_multiplicity: 2
        allosteric_regulators: A, I

        steps: begin
            E + F6P ⇌ E(F6P)        :: EqualAI
            E(F6P) <--> E(F16BP)    :: EqualAI
            E(F16BP) ⇌ E + F16BP    :: EqualAI
        end

        regulatory_site(4): A::OnlyA
        regulatory_site(4): I::OnlyI
    end

Build an `AllostericEnzymeMechanism` (MWC, two conformations).

- `substrates:`, `products:`, `modifiers:` accept comma-separated
  bare symbols.
- `allosteric_regulators:` lists the allosteric regulators by name; each sits on
  one `regulatory_site(n):` line, which carries its tag.
- `catalytic_multiplicity: N` is the subunit count for the catalytic site
  (default 1).
- `steps: begin ... end` is required (exactly once); each step or
  parenthesized step-group must carry a `::Tag`, one of `OnlyA`, `EqualAI` or
  `NonequalAI` (`OnlyI` is rejected on catalytic groups). A single step carries
  its tag after its right-hand side, a group after its closing parenthesis.
  Function-call species notation (`E(F6P)`, `Estar(B; residual = A - P)`) is
  supported.
- `regulatory_site(N): L1::Tag, L2::Tag` declares one regulatory site with
  multiplicity `N` and its ligands, in any order, each with its tag: one of
  `OnlyA`, `OnlyI`, `EqualAI`, `NonequalAI`. Each ligand must appear in
  `allosteric_regulators:`, and each allosteric regulator sits at exactly one
  site.
- A name with several roles binds in a bare catalytic-step `E(X)` as its
  substrate or product role first, then as a modifier; an allosteric regulator
  binds only at its regulatory site. A substrate or product bound as a modifier
  is written `E(X::Mod)`, and `X::Mod` requires `X` in `modifiers:`; list an
  allosteric regulator under `modifiers:` as well to let it bind the catalytic
  site.
"""
macro allosteric_mechanism(block)
    return :($_compiled_mechanism($(_raw_mechanism(block, true))))
end

# ─── Mechanism macro blocks as written ─────────────────────────────────

"""
    _raw_mechanism(block, allosteric::Bool) -> RawMechanism or RawAllostericMechanism

The mechanism a mechanism macro's `begin ... end` block writes, as raw input:
`@allosteric_mechanism`'s grammar when `allosteric`, else `@enzyme_mechanism`'s. It reads
the block's expression without evaluating it or building code, and throws the DSL's
syntax errors; the constructors check the mechanism (`Mechanism(raw)`,
`AllostericMechanism(raw)`).

`substrates:`, `products:` and `modifiers:` list bare names and may repeat; `steps:`
holds one `begin ... end` block (`_raw_steps`). `@allosteric_mechanism` adds
`allosteric_regulators:` (bare names), `catalytic_multiplicity:` (one positive Int, 1
when omitted) and `regulatory_site(n):` lines, each listing its ligands as `X::Tag`
(`_raw_site`); every allosteric regulator sits on exactly one of them.
"""
function _raw_mechanism(block, allosteric::Bool)
    macro_name = allosteric ? "@allosteric_mechanism" : "@enzyme_mechanism"
    block isa Expr && block.head === :block ||
        error("$macro_name: expected a `begin ... end` block, got $block")
    substrates, products, modifiers, regulators = Symbol[], Symbol[], Symbol[], Symbol[]
    steps_block, multiplicity, sites = nothing, nothing, RawSite[]
    for arg in block.args
        arg isa LineNumberNode && continue
        label, values = _parse_labeled_line(arg)
        if label === :substrates
            append!(substrates, _bare_symbols_from_values(values, label))
        elseif label === :products
            append!(products, _bare_symbols_from_values(values, label))
        elseif label === :modifiers
            append!(modifiers, _bare_symbols_from_values(values, label))
        elseif label === :steps
            steps_block === nothing ||
                error("$macro_name: `steps:` given more than once.")
            length(values) == 1 || error("$macro_name: `steps:` takes one " *
                                         "`begin ... end` block; got $(length(values)) " *
                                         "values.")
            steps_block = only(values)
            steps_block isa Expr && steps_block.head === :block ||
                error("$macro_name: `steps:` takes a `begin ... end` block of steps; " *
                      "got $steps_block.")
        elseif allosteric && label === :allosteric_regulators
            tagged = findfirst(v -> v isa Expr && v.head === :(::), values)
            tagged === nothing ||
                error("$macro_name: `allosteric_regulators:` lists names; got " *
                      "`$(values[tagged])`. A ligand's tag goes on its " *
                      "`regulatory_site(n):` line, as in " *
                      "`regulatory_site(2): $(values[tagged])`.")
            append!(regulators, _bare_symbols_from_values(values, label))
        elseif allosteric && label === :catalytic_multiplicity
            multiplicity === nothing ||
                error("$macro_name: `catalytic_multiplicity:` given more than once.")
            length(values) == 1 ||
                error("$macro_name: `catalytic_multiplicity:` takes one positive Int; " *
                      "got $(join(values, ", ")).")
            multiplicity = _positive_int(only(values),
                                         "$macro_name: `catalytic_multiplicity:`")
        elseif allosteric && label isa Expr && label.head === :call &&
               label.args[1] === :regulatory_site
            push!(sites, _raw_site(label, values))
        else
            error("$macro_name: unknown label `$label:`. Valid labels: " * (allosteric ?
                  "substrates, products, modifiers, allosteric_regulators, " *
                  "catalytic_multiplicity, steps, regulatory_site(n)." :
                  "substrates, products, modifiers, steps. Allosteric declarations " *
                  "(`allosteric_regulators:`, `catalytic_multiplicity:`, " *
                  "`regulatory_site(n):`) belong in @allosteric_mechanism."))
        end
    end
    isempty(substrates) && error("$macro_name: `substrates:` not specified.")
    isempty(products) && error("$macro_name: `products:` not specified.")
    steps_block === nothing && error("$macro_name: `steps:` not specified.")
    names = (; substrates, products, modifiers, regulators)
    steps, group, tags = _raw_steps(steps_block, names, macro_name, allosteric)
    raw = RawMechanism(substrates, products, modifiers, steps, group)
    allosteric || return raw
    _check_site_ligands(regulators, sites)
    RawAllostericMechanism(raw, tags, something(multiplicity, 1), sites)
end

"""
A `regulatory_site(n): X::Tag, ...` line as written: `label` is its `regulatory_site(n)`
call, whose one argument is the site's multiplicity, and `values` its ligand entries,
each a name with its tag.
"""
function _raw_site(label, values)
    length(label.args) == 2 && !(label.args[2] isa Expr) ||
        error("@allosteric_mechanism: `regulatory_site(n):` takes one positive Int " *
              "multiplicity, as in `regulatory_site(2): A::OnlyA, I::OnlyI`; got " *
              "`$label`.")
    multiplicity = _positive_int(label.args[2],
                                 "@allosteric_mechanism: `regulatory_site` multiplicity")
    ligands, tags = Symbol[], Symbol[]
    for v in values
        v isa Expr && v.head === :(::) && length(v.args) == 2 &&
            v.args[1] isa Symbol && v.args[2] isa Symbol ||
            error("@allosteric_mechanism: `regulatory_site(n):` lists each ligand with " *
                  "its tag, as in `regulatory_site(2): A::OnlyA, I::OnlyI`; got " *
                  (v isa Expr && v.head === :block ? "a `begin ... end` block." : "`$v`."))
        push!(ligands, v.args[1])
        push!(tags, v.args[2])
    end
    RawSite(ligands, multiplicity, tags)
end

"""
Error unless the `regulatory_site(n):` lines `sites` place each allosteric regulator of
`regulators`, the names `allosteric_regulators:` lists once each, exactly once, and list
no other ligand.
"""
function _check_site_ligands(regulators::Vector{Symbol}, sites::Vector{RawSite})
    repeated = unique(n for (i, n) in enumerate(regulators)
                      if n in view(regulators, 1:i-1))
    isempty(repeated) ||
        error("@allosteric_mechanism: `allosteric_regulators:` lists " *
              "$(join(("`$n`" for n in repeated), ", ")) more than once.")
    for site in sites, l in site.ligands
        l in regulators ||
            error("@allosteric_mechanism: ligand `$l` on a `regulatory_site(n):` line is " *
                  "not declared in `allosteric_regulators:`.")
    end
    for r in regulators
        n = sum(count(==(r), site.ligands) for site in sites; init = 0)
        n == 0 && error("@allosteric_mechanism: allosteric regulator `$r` sits at no " *
                        "regulatory site; list it on one `regulatory_site(n):` line, as " *
                        "in `regulatory_site(2): $r::OnlyA`.")
        n == 1 || error("@allosteric_mechanism: allosteric regulator `$r` is listed $n " *
                        "times on `regulatory_site(n):` lines; an allosteric regulator " *
                        "sits at exactly one regulatory site.")
    end
end

"""
The steps of a `steps:` block as written, each statement's group label per step and,
for an allosteric block, each statement's tag: `(steps, group, tags)`. A statement is
one kinetic group: a step (`lhs ⇌ rhs` rapid equilibrium, `lhs <--> rhs` steady state)
or a parenthesized tuple of steps (`(step,)` for one). An allosteric statement ends in
its tag (`step :: Tag`, `(step, step) :: Tag`); a plain one carries none. `names` holds
the block's declared names.
"""
function _raw_steps(steps_block, names, macro_name, allosteric::Bool)
    steps, group, tags = RawStep[], Int[], Symbol[]
    missing_tag = "`:: <:OnlyA|:EqualAI|:NonequalAI>` annotation. Add `:: <state>` after"
    n = 0
    for arg in steps_block.args
        arg isa LineNumberNode && continue
        n += 1
        if arg isa Expr && arg.head === :(::) && length(arg.args) == 2 &&
           arg.args[1] isa Expr && arg.args[1].head in (:tuple, :call)
            inner, tag = arg.args
            allosteric || error("$macro_name: tag annotation `$arg` is not allowed")
            inner isa Expr && inner.head === :tuple ||
                error("$macro_name: a single step carries its tag after its right-hand " *
                      "side, as in `E + S ⇌ E(S) :: EqualAI`. Got `$arg`.")
            tag isa Symbol ||
                error("$macro_name: step-group tag must be a Symbol; got $tag")
            push!(tags, tag)
            members = inner.args
        elseif arg isa Expr && arg.head === :tuple
            allosteric && error("$macro_name: parenthesized step group `$arg` is missing " *
                                "$missing_tag the closing paren.")
            members = arg.args
        elseif arg isa Expr && arg.head === :call && length(arg.args) == 3
            op, lhs, rhs = _step_parts(arg, macro_name)
            tag = _peel_tag!(rhs)
            if tag !== nothing
                tag isa Symbol || error("$macro_name: step tag must be a Symbol; got $tag")
                allosteric || error("$macro_name: tag annotation on `$arg` is not allowed")
                push!(tags, tag)
            elseif allosteric
                error("$macro_name: step `$arg` is missing $missing_tag the step " *
                      "expression.")
            end
            push!(steps, _raw_step(op, lhs, rhs, names, macro_name))
            push!(group, n)
            continue
        else
            error("$macro_name: expected step or step-group; got $arg")
        end
        isempty(members) && error("$macro_name: a parenthesized group holds no step; a " *
                                  "kinetic group holds at least one step.")
        for s in members
            push!(steps, _raw_step(_step_parts(s, macro_name)..., names, macro_name))
            push!(group, n)
        end
    end
    steps, group, tags
end

"""A step expression's operator, `⇌` or `<-->`, and the terms of its two sides; a side
is a `+` chain of terms or a single term."""
function _step_parts(expr, macro_name)
    expr isa Expr && expr.head === :call ||
        error("$macro_name: expected lhs ⇌ rhs or lhs <--> rhs; got $expr")
    op = expr.args[1]
    op === :⇌ || op === :(<-->) ||
        error("$macro_name: expected ⇌ or <--> step operator; got $op")
    length(expr.args) == 3 ||
        error("$macro_name: expected lhs ⇌ rhs or lhs <--> rhs; got $expr")
    op, _side_terms(expr.args[2]), _side_terms(expr.args[3])
end

"""The terms of a step side: the operands of a `+` chain, or the side itself."""
_side_terms(side) = side isa Expr && side.head === :call && side.args[1] === :+ ?
                    side.args[2:end] : Any[side]

"""
The tag a step's right-hand terms `rhs` end in, taken off the last term, or `nothing`.
`::` binds tighter than `+`, so `E(P) ⇌ E + P :: EqualAI` reads `P :: EqualAI` as its
last term, while `E(P) ⇌ (E + P) :: EqualAI` tags the whole side, whose terms then
replace it. A role tag (`::Mod`, `::Inh`) marks a metabolite, never a step.
"""
function _peel_tag!(rhs::Vector{Any})
    t = rhs[end]
    t isa Expr && t.head === :(::) && length(t.args) == 2 &&
        t.args[2] ∉ (:Mod, :Inh) || return nothing
    rhs[end] = t.args[1]
    length(rhs) == 1 && splice!(rhs, 1, _side_terms(t.args[1]))
    t.args[2]
end

"""
One step as written: `op` is `⇌` (rapid equilibrium) or `<-->` (steady state), `lhs`
and `rhs` the terms of its two sides. Each side holds exactly one enzyme form and any
number of metabolites (`_raw_met`): the left-hand metabolites are consumed, the
right-hand ones released. Every term of both sides is read before either side's enzyme
forms are counted, so a misspelled bound metabolite is reported as undeclared rather
than as a second enzyme form. Each form's conformation label is checked
(`_check_conformation_label`).
"""
function _raw_step(op, lhs, rhs, names, macro_name)
    is_met(t) = t isa Symbol ? _is_declared(t, names) : t isa Expr && t.head === :(::)
    read_term(t) = is_met(t) ? _raw_met(t, names, macro_name) :
                   _raw_form(t, names, macro_name)
    term_name(t) = t isa Expr && !isempty(t.args) ? t.args[1] : t
    parsed = (Any[read_term(t) for t in lhs], Any[read_term(t) for t in rhs])
    (from, consumed), (to, released) = map((lhs, rhs), parsed) do terms, side
        i = findall(x -> x isa RawForm, side)
        isempty(i) && error("$macro_name: step side has no enzyme-form term (terms: " *
                            "$(join(term_name.(terms), ", "))).")
        length(i) == 1 ||
            error("$macro_name: step side has more than one enzyme-form term " *
                  "($(term_name(terms[i[1]])), $(term_name(terms[i[2]]))); each " *
                  "elementary step has exactly one enzyme form per side.")
        _check_conformation_label(side[only(i)].conformation, macro_name)
        side[only(i)], Met[x for x in side if x isa Met]
    end
    RawStep(from, to, consumed, released, op === :⇌)
end

"""Whether `x` is a name the block declares (`names`: substrates, products, modifiers
and allosteric regulators)."""
_is_declared(x, names) = x in names.substrates || x in names.products ||
                         x in names.modifiers || x in names.regulators

"""The block's declared names (`names`), each once, in name order, for error messages."""
_declared_names(names) =
    sort!(unique!([names.substrates; names.products; names.modifiers; names.regulators]))

"""
A metabolite term as written: a declared name `X`, or `X::Mod`. A substrate or product
written bare binds as itself, a name under `modifiers:` that is no reactant as a modifier,
and `X::Mod` as the modifier copy of a reactant `X` that `modifiers:` also lists. An
allosteric regulator binds only at its regulatory site. `species` is the enzyme-form term
that binds the metabolite, or `nothing` for a free metabolite of a step.
"""
function _raw_met(t, names, macro_name, species = nothing)
    name, tag = t isa Expr && length(t.args) == 2 ? (t.args[1], t.args[2]) : (t, nothing)
    in_species = species === nothing ? "" : " in species `$species`"
    name isa Symbol || error("$macro_name: expected `name::Mod`$in_species; got $t")
    tag === nothing || tag === :Mod ||
        error("$macro_name: unknown role tag ::$tag on $name; the one role tag is " *
              "::Mod, which marks a substrate or product bound as a modifier.")
    _is_declared(name, names) ||
        error("$macro_name: " * (species === nothing ? "tagged metabolite `$name` in `$t`" :
                                 "bound metabolite `$name` in species `$species`") *
              " is not declared. Declared: $(_declared_names(names)).")
    reactant = name in names.substrates || name in names.products
    if tag === :Mod
        reactant ||
            error("$macro_name: `$name::Mod` marks a substrate or product bound as a " *
                  "modifier, and `$name` is neither; write a name under `modifiers:` " *
                  "bare.")
        name in names.modifiers ||
            error("$macro_name: `$name::Mod` names `$name`, which `modifiers:` does not " *
                  "declare; add `$name` to `modifiers:` to let it bind the catalytic " *
                  "site as a modifier.")
        return Met(name, true)
    end
    reactant && return Met(name, false)
    name in names.modifiers && return Met(name, true)
    error("$macro_name: `$name` is an allosteric regulator; it binds only at its " *
          "regulatory site.")
end

"""
An enzyme-form term as written: a bare conformation label `E`, or a call
`E(S, X::Mod; residual = A - P)` whose head is the conformation label, whose positional
arguments are the bound metabolites (`_raw_met`) and whose `residual` lists the
substrates added to and the products subtracted from the enzyme (`_read_residual!`). A
conformation label may not be a declared name.
"""
function _raw_form(t, names, macro_name)
    t isa Symbol || t isa Expr && t.head === :call && t.args[1] isa Symbol ||
        error("$macro_name: expected metabolite Symbol or species expression on step " *
              "side; got $t")
    conformation, args = t isa Symbol ? (t, Any[]) : (t.args[1], t.args[2:end])
    _is_declared(conformation, names) &&
        error("$macro_name: conformation label `$conformation` collides with declared " *
              "metabolite `$conformation`; choose a different conformation label.")
    bound, added, subtracted = Met[], Symbol[], Symbol[]
    for a in args
        if a isa Expr && a.head === :parameters
            for kw in a.args
                kw isa Expr && kw.head === :kw && kw.args[1] === :residual ||
                    error("$macro_name: unknown keyword in species `$t`: $kw. Only " *
                          "`residual = ...` is allowed.")
                _read_residual!(added, subtracted, kw.args[2], true, names, macro_name)
            end
        else
            a isa Symbol || a isa Expr && a.head === :(::) ||
                error("$macro_name: invalid entry in species `$t`: $a")
            push!(bound, _raw_met(a, names, macro_name, t))
        end
    end
    RawForm(conformation, bound, added, subtracted)
end

"""
Read a residual expression (`A`, `A - P`, `S1 + S2 - P1 - P3`, ...) into `added` (its
positive names) and `subtracted` (its negative ones); `positive` is the sign of `e`. A
residual adds substrates and subtracts products.
"""
function _read_residual!(added, subtracted, e, positive::Bool, names, macro_name)
    if e isa Symbol
        _is_declared(e, names) ||
            error("$macro_name: residual entry `$e` is not a declared metabolite. " *
                  "Declared: $(_declared_names(names)).")
        e in (positive ? names.substrates : names.products) ||
            error("$macro_name: a residual adds substrates and subtracts products; got " *
                  "`$e`.")
        push!(positive ? added : subtracted, e)
    elseif e isa Expr && e.head === :call && e.args[1] === :+ && length(e.args) >= 3
        for a in e.args[2:end]
            _read_residual!(added, subtracted, a, positive, names, macro_name)
        end
    elseif e isa Expr && e.head === :call && e.args[1] === :- && length(e.args) in (2, 3)
        length(e.args) == 3 &&
            _read_residual!(added, subtracted, e.args[2], positive, names, macro_name)
        _read_residual!(added, subtracted, e.args[end], !positive, names, macro_name)
    else
        error("$macro_name: invalid residual expression: $e")
    end
end

"""
Error unless the conformation label `c` is a capital letter followed by lowercase
letters and digits (`E`, `Estar`, `F2`). A label with a second capital reads as an
opaque bound-form name (`ES`), which is written in call notation (`E(S)`). Constant
names join form names with `_`, so a label holds none.
"""
function _check_conformation_label(c::Symbol, macro_name)
    occursin(r"^[A-Z][a-z0-9]*$", String(c)) && return nothing
    occursin(r"^[A-Z].*[A-Z]", String(c)) &&
        error("$macro_name: `$c` looks like an opaque bound-form name; write it as " *
              "decomposed call notation, e.g. `E(S)` or `E(A, B)`.")
    error("$macro_name: `$c` is not a conformation label: a conformation label is a " *
          "capital letter followed by lowercase letters and digits, such as `E`, `Estar` " *
          "or `F2`.")
end

"""
    _compiled_mechanism(raw) -> EnzymeMechanism or AllostericEnzymeMechanism

The compiled mechanism of `raw`, the walker's value for a mechanism macro's block
(`_raw_mechanism`): the macros' expansion calls it, so the constructors check the
mechanism when the expanded code runs.
"""
_compiled_mechanism(raw::RawMechanism) = EnzymeMechanism(Mechanism(raw))
_compiled_mechanism(raw::RawAllostericMechanism) =
    AllostericEnzymeMechanism(AllostericMechanism(raw))
