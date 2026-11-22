# ABOUTME: EnzymeReaction, a reaction as a value: its reactants with their atoms, its
# ABOUTME: regulators, the checks of metabolite names, and its text as the macro block.

"""
Reactant: a substrate or product of a reaction, with its atoms as `element => count`
pairs. The constructor sorts the atoms by element so two equivalent constructions compare
equal under `==` / `hash`; an element listed twice (`S[CHC]`) counts once with the summed
count, as `S[C2H]` writes it.
"""
struct Reactant
    name::Symbol
    atoms::Vector{Pair{Symbol,Int}}
    function Reactant(name::Symbol, atoms::Vector{<:Pair{Symbol,<:Integer}})
        isempty(atoms) && error(
            "Reactant: $name has no declared atoms; atoms " *
            "are mandatory (use `[C…]` bracket syntax in @enzyme_reaction).")
        bad = findfirst(a -> a.second isa Bool || a.second ≤ 0, atoms)
        bad === nothing ||
            error("Reactant: $name atom count for element " *
                  "$(atoms[bad].first) must be a positive integer, got " *
                  "$(atoms[bad].second).")
        bad = findfirst(a -> !occursin(r"^[A-Z][a-z]*$", String(a.first)), atoms)
        bad === nothing ||
            error("Reactant: $name has the element " *
                  "$(atoms[bad].first); an element is a capital letter and any lowercase " *
                  "letters, such as C or Cl")
        merged = Pair{Symbol,Int}[]
        for (e, c) in sort(Vector{Pair{Symbol,Int}}(atoms); by=first)
            if !isempty(merged) && merged[end].first === e
                merged[end] = e => merged[end].second + c
            else
                push!(merged, e => c)
            end
        end
        new(name, merged)
    end
end

"""
Regulator: a metabolite a reaction declares as a regulator. `kind` is `:allosteric` (it
binds a regulatory site) or `:competitive` (a competitive inhibitor of the catalytic
site). `effect` is `:activator` or `:inhibitor` when the declaration pins the regulator
to the conformation it favors, else `:unspecified`; a competitive inhibitor takes no
effect.
"""
struct Regulator
    name::Symbol
    kind::Symbol
    effect::Symbol
    function Regulator(name::Symbol, kind::Symbol, effect::Symbol = :unspecified)
        kind in (:allosteric, :competitive) ||
            error("Regulator: kind must be :allosteric or :competitive, got $kind")
        effect in (:activator, :inhibitor, :unspecified) ||
            error("Regulator: effect must be :activator, :inhibitor, or " *
                  ":unspecified, got $effect")
        kind === :competitive && effect !== :unspecified &&
            error("Regulator: $name is a competitive inhibitor, " *
                  "which takes no regulator type, got $(repr(effect)); only an " *
                  "allosteric regulator is an activator (::Activator) or an inhibitor " *
                  "(::Inhibitor)")
        new(name, kind, effect)
    end
end

"""
The names no metabolite takes: `Keq`, `L` and `N`, which name or are kept for
rate-equation constants, and `params`, `concs` and `v`, the local names of a rate
equation's code (`rate_equation_string`). A rate equation binds its constants, its
locals and the concentrations in one scope, so a metabolite of one of these names would
overwrite another or be overwritten.
"""
const _RESERVED_NAMES = (:Keq, :L, :N, :params, :concs, :v)

"""
Julia's reserved words and the infix keywords `in`, `isa` and `where`. `Base.isidentifier`
accepts each of them, but a DSL text cannot write one as a metabolite's name.
"""
const _KEYWORDS = (:baremodule, :begin, :break, :catch, :const, :continue, :do, :else,
                   :elseif, :end, :export, :finally, :for, :function, :global, :if,
                   :import, :in, :isa, :let, :local, :macro, :module, :quote, :return,
                   :struct, :try, :using, :where, :while)

"""
Error unless each name of `names`, the names of a `kind` of thing (`"metabolite"`), is one
a DSL text can write: an identifier that is no keyword, in the form Julia reads it. Julia's
parser normalizes an identifier (NFC, and `Unicode.julia_chartransform`: a micro sign reads
as `μ`), so a name in another form would print as text that reads as another name. `who`
names the constructor that checks.
"""
function _check_identifiers(who::String, kind::String, names::Vector{Symbol})
    article = first(kind) in "aeiou" ? "an" : "a"
    for n in names
        s = String(n)
        Base.isidentifier(s) && !(n in _KEYWORDS) ||
            error("$who: \"$s\" cannot name $article $kind; $article $kind name is an " *
                  "identifier that is no keyword")
        read = isascii(s) ? s :
               Unicode.normalize(s; stable = true, compose = true,
                                 chartransform = Unicode.julia_chartransform)
        read == s || error("$who: Julia reads the name \"$s\" as \"$read\", so no DSL " *
                           "text can write it; use \"$read\"")
    end
end

"""
Error when a metabolite name of `names` is no name a DSL text can write
(`_check_identifiers`), is `E`, the name of the free enzyme, is reserved
(`_RESERVED_NAMES`) or holds `_`: constant names join form and metabolite names with `_`,
so a metabolite name holds none. `who` names the constructor that checks.
"""
function _check_metabolite_names(who::String, names::Vector{Symbol})
    _check_identifiers(who, "metabolite", names)
    for n in names
        n === :E &&
            error("$who: E cannot name a metabolite; it labels the free enzyme form")
        n in _RESERVED_NAMES &&
            error("$who: $n is a reserved name; $(join(_RESERVED_NAMES, ", ", " and ")) " *
                  "are kept for a rate equation's constants and local names")
        occursin('_', String(n)) &&
            error("$who: the metabolite name $n holds `_`; constant names join " *
                  "metabolite names with `_`, so a metabolite name holds none")
    end
end

"""
    EnzymeReaction

A reaction: its substrates, products, optional regulators, and the oligomeric states the
mechanism enumerator may consider. It is the entry point to the package —
pair it with rate data in an
[`IdentifyRateEquationProblem`](@ref EnzymeRates.IdentifyRateEquationProblem), or pass
it to `EnzymeRates.Enumeration.candidate_mechanisms` to enumerate candidate mechanisms.

# Fields
- `substrates::Vector{Reactant}` and `products::Vector{Reactant}` — each reactant
  carrying its per-atom inventory.
- `regulators::Vector{Regulator}` — competitive inhibitors and
  allosteric regulators, each with its `kind` and `effect`.
- `allowed_oligomeric_states::Vector{Int}` — the oligomeric states, one
  or more, that the allosteric enumerator may assign to the catalytic core.
- `shared_catalytic_site::Vector{Tuple{Symbol,Symbol}}` — `(substrate,
  product)` pairs that bind the same catalytic site, and so are never both
  bound to it at once.

Construct one with the [`@enzyme_reaction`](@ref EnzymeRates.DSL.@enzyme_reaction) DSL. Each
reactant takes an atom bracket (`S[C]`, `A[C1H1]`); the brackets are load-bearing for
ping-pong and multi-substrate reactions, whose residues and atom balance the search reads
from them. The constructor sorts substrates, products and regulators by name, a name's
competitive inhibitor before its allosteric regulator, so two equivalent declarations
compare equal under `==`/`hash`. It rejects unbalanced atoms, a name listed as both a
substrate and a product, the names that `@enzyme_reaction` says no metabolite takes, and a
name that reads in Julia code as another (a micro sign as `μ`), which no DSL text can
write. A reaction prints as the `@enzyme_reaction` block that builds it.

```jldoctest
julia> using EnzymeRates

julia> rxn = @enzyme_reaction begin
           substrates: S[C]
           products:   P[C]
       end;

julia> rxn isa EnzymeReaction
true

julia> rxn
@enzyme_reaction begin
    substrates: S[C]
    products: P[C]
end
```
"""
struct EnzymeReaction
    substrates::Vector{Reactant}
    products::Vector{Reactant}
    regulators::Vector{Regulator}
    allowed_oligomeric_states::Vector{Int}
    shared_catalytic_site::Vector{Tuple{Symbol,Symbol}}

    function EnzymeReaction(substrates::Vector{Reactant}, products::Vector{Reactant},
                            regulators::Vector{Regulator},
                            allowed_oligomeric_states::Vector{Int};
                            shared_catalytic_site = Tuple{Symbol,Symbol}[])
        sorted_substrates = sort(substrates; by = r -> r.name)
        sorted_products = sort(products; by = r -> r.name)
        sorted_regulators = sort(regulators; by = r -> (r.name, r.kind === :allosteric))
        # A nonempty list of distinct Ints ≥ 1, as a DSL text writes one.
        ms = allowed_oligomeric_states
        !isempty(ms) && all(≥(1), ms) ||
            error("EnzymeReaction: the allowed oligomeric states must be a nonempty " *
                  "list of Ints ≥ 1, got $ms")
        allunique(ms) || error("EnzymeReaction: the allowed oligomeric states list a " *
                               "multiplicity more than once; got $ms.")
        sorted_mults = sort(allowed_oligomeric_states)

        sub_names  = Symbol[r.name for r in sorted_substrates]
        prod_names = Symbol[r.name for r in sorted_products]
        isempty(sub_names) && error("EnzymeReaction: no substrate is given; a reaction " *
                                    "has at least one substrate.")
        isempty(prod_names) && error("EnzymeReaction: no product is given; a reaction " *
                                     "has at least one product.")

        reg_keys   = [(r.name, r.kind) for r in sorted_regulators]
        allunique(sub_names)  || error("EnzymeReaction: duplicate substrate names")
        allunique(prod_names) || error("EnzymeReaction: duplicate product names")
        allunique(reg_keys)   ||
            error("EnzymeReaction: duplicate regulator of the same kind")
        both_roles = intersect(sub_names, prod_names)
        isempty(both_roles) || error(
            "EnzymeReaction: $(join(both_roles, ", ")) named as both a substrate and " *
            "a product; concentrations and constants are keyed by name")
        _check_metabolite_names("EnzymeReaction",
                                Symbol[sub_names; prod_names; first.(reg_keys)])

        sub_set  = Set(sub_names)
        prod_set = Set(prod_names)
        normalized_shared = Tuple{Symbol,Symbol}[]
        for pair in shared_catalytic_site
            length(pair) == 2 || error(
                "EnzymeReaction: each shared_catalytic_site entry must name " *
                "exactly two metabolites, got $pair")
            a, b = pair[1], pair[2]
            if a in sub_set && b in prod_set
                push!(normalized_shared, (a, b))
            elseif b in sub_set && a in prod_set
                push!(normalized_shared, (b, a))
            else
                error("EnzymeReaction: shared_catalytic_site pair $pair must " *
                      "name one declared substrate and one declared product")
            end
        end
        sorted_shared = sort(unique(normalized_shared))
        length(sorted_shared) == length(normalized_shared) || error(
            "EnzymeReaction: duplicate shared_catalytic_site pair")

        # A name may be declared once per regulator kind — the same metabolite
        # may bind BOTH as an allosteric regulator and as a competitive inhibitor
        # (the two roles render to distinct parameter names). A regulator may
        # also share a substrate/product name (so `concs.X` drives every role);
        # cross-category names need not be unique.

        # Atom mass-balance: per-element sum over substrates vs products.
        sub_atoms  = Dict{Symbol,Int}()
        prod_atoms = Dict{Symbol,Int}()
        for (rs, tgt) in ((sorted_substrates, sub_atoms), (sorted_products, prod_atoms)),
            r in rs, (elem, c) in r.atoms
            tgt[elem] = get(tgt, elem, 0) + c
        end
        for elem in union(keys(sub_atoms), keys(prod_atoms))
            s_c = get(sub_atoms, elem, 0)
            p_c = get(prod_atoms, elem, 0)
            s_c == p_c || error(
                "EnzymeReaction: atom imbalance — element $elem appears " *
                "$s_c time(s) on substrate side and $p_c on product side. " *
                "Declared atoms must balance.")
        end

        new(sorted_substrates, sorted_products, sorted_regulators, sorted_mults,
            sorted_shared)
    end
end

"""
Print `r` as the `@enzyme_reaction` text that builds it: each reactant with its atoms as
one formula, the competitive inhibitors and allosteric regulators, the oligomeric states
unless they are the default 1, and the shared catalytic sites.
"""
function Base.show(io::IO, r::EnzymeReaction)
    mults = r.allowed_oligomeric_states
    formula(re) = join("$e" * (c == 1 ? "" : "$c") for (e, c) in re.atoms)
    tag(reg) = reg.effect === :activator ? "::Activator" :
               reg.effect === :inhibitor ? "::Inhibitor" : ""
    entry(reg) = "$(reg.name)" * tag(reg)
    line(label, entries) = isempty(entries) ||
                           print(io, "\n    ", label, ": ", join(entries, ", "))
    print(io, "@enzyme_reaction begin")
    line("substrates", ["$(re.name)[$(formula(re))]" for re in r.substrates])
    line("products", ["$(re.name)[$(formula(re))]" for re in r.products])
    for (label, kind) in (("competitive_inhibitors", :competitive),
                          ("allosteric_regulators", :allosteric))
        line(label, [entry(reg) for reg in r.regulators if reg.kind === kind])
    end
    mults == [1] || print(io, "\n    oligomeric_state: ",
                          length(mults) == 1 ? only(mults) : "(" * join(mults, ", ") * ")")
    line("shared_catalytic_site", ["($s, $p)" for (s, p) in r.shared_catalytic_site])
    print(io, "\nend")
end
