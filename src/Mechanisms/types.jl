# ABOUTME: Value types: EnzymeReaction, and the mechanism as a reaction graph in canonical
# ABOUTME: form, built from nodes and steps, with its checks, identity and Sig.

# ─── Concrete types ───────────────────────────────────────────────────

"""
RegulatorySite: a binding site (possibly multimeric) for one or
more allosteric ligands, each given by its name, once. `ligands[i]` and `mwc_tags[i]`
are parallel. The constructor sorts the ligands by name and carries each ligand's state
along, so the order the ligands are given in never changes `==`, `hash` or
the compiled type.
"""
struct RegulatorySite
    ligands::Vector{Symbol}
    multiplicity::Int
    mwc_tags::Vector{Symbol}
    function RegulatorySite(ligands::Vector{Symbol},
                            multiplicity::Int, mwc_tags::Vector{Symbol})
        isempty(ligands) && error("RegulatorySite: needs at least one ligand")
        length(ligands) == length(mwc_tags) ||
            error("RegulatorySite: length(ligands)=$(length(ligands)) " *
                  "must equal length(mwc_tags)=$(length(mwc_tags))")
        multiplicity ≥ 1 ||
            error("RegulatorySite: multiplicity must be ≥ 1, got $multiplicity")
        bad = findfirst(∉(_TAG_CODES), mwc_tags)
        bad === nothing ||
            error("RegulatorySite: MWC tag $(mwc_tags[bad]) must be one of " *
                  "$_TAG_CODES")
        dup = findfirst(i -> ligands[i] in view(ligands, 1:(i - 1)), eachindex(ligands))
        dup === nothing || error("RegulatorySite: ligand $(ligands[dup]) appears twice " *
                                 "in one regulatory site")
        p = sortperm(ligands)
        new(ligands[p], multiplicity, mwc_tags[p])
    end
end

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
The names no metabolite takes: `E_total`, `Keq`, `L` and `N`, which name or are kept for
rate-equation constants, and `params`, `concs` and `v`, the local names of a rate
equation's code (`rate_equation_string`). A rate equation binds its constants, its
locals and the concentrations in one scope, so a metabolite of one of these names would
overwrite another or be overwritten.
"""
const _RESERVED_NAMES = (:E_total, :Keq, :L, :N, :params, :concs, :v)

"""
Julia's reserved words and the infix keywords `in`, `isa` and `where`. `Base.isidentifier`
accepts each of them, but a DSL text cannot write one as a metabolite's name.
"""
const _KEYWORDS = (:baremodule, :begin, :break, :catch, :const, :continue, :do, :else,
                   :elseif, :end, :export, :finally, :for, :function, :global, :if,
                   :import, :in, :isa, :let, :local, :macro, :module, :quote, :return,
                   :struct, :try, :using, :where, :while)

"""
Error when a metabolite name of `names` is no name a DSL text can write (an identifier
that is no keyword, in the form Julia reads it), is `E`, the label of the free enzyme
form, is reserved (`_RESERVED_NAMES`) or holds `_`: constant names join form and
metabolite names with `_`, so a metabolite name holds none.
Julia's parser normalizes an identifier (NFC, and
`Unicode.julia_chartransform`: a micro sign reads as `μ`), so a name in another form would
print as text that reads as another name. `who` names the constructor that checks.
"""
function _check_metabolite_names(who::String, names::Vector{Symbol})
    for n in names
        s = String(n)
        Base.isidentifier(s) && !(n in _KEYWORDS) ||
            error("$who: \"$s\" cannot name a metabolite; a metabolite name is an " *
                  "identifier that is no keyword")
        read = isascii(s) ? s :
               Unicode.normalize(s; stable = true, compose = true,
                                 chartransform = Unicode.julia_chartransform)
        read == s || error("$who: Julia reads the name \"$s\" as \"$read\", so no DSL " *
                           "text can write it; use \"$read\"")
        n === :E &&
            error("$who: E cannot name a metabolite; it labels the free enzyme form")
        n in _RESERVED_NAMES &&
            error("$who: $n is a reserved name; $(join(_RESERVED_NAMES, ", ", " and ")) " *
                  "are kept for a rate equation's constants and local names")
        occursin('_', s) &&
            error("$who: the metabolite name $n holds `_`; constant names join " *
                  "metabolite names with `_`, so a metabolite name holds none")
    end
end

"""
Error unless the conformation label `c` is a capital letter followed by lowercase
letters and digits (`E`, `Estar`, `F2`), other than `A` and `I`. A label with a second
capital reads as an opaque bound-form name (`ES`), which is written in call notation
(`E(S)`). Constant names join form names with `_`, so a label holds none, and they tag a
conformation with `A_` or `I_`, so a step side that starts with a form labelled `A` or
`I` could give an `:EqualAI` step's constant the tagged name of another step's
(`A(Eb) <--> A + Eb` and the active constant of `E(b) <--> A(Eb)` both render
`k_A_Eb_to_AEb`). `who` names the macro or the constructor that checks.
"""
function check_conformation_label(c::Symbol, who::String)
    c in (:A, :I) &&
        error("$who: `$c` cannot be a conformation label: constant names tag a " *
              "conformation with A_ or I_ (`k_A_ES_to_EP`), so the labels A and I " *
              "could give two constants one name; choose another label.")
    occursin(r"^[A-Z][a-z0-9]*$", String(c)) && return nothing
    occursin(r"^[A-Z].*[A-Z]", String(c)) &&
        error("$who: `$c` looks like an opaque bound-form name; write it as " *
              "decomposed call notation, e.g. `E(S)` or `E(A, B)`.")
    error("$who: `$c` is not a conformation label: a conformation label is a " *
          "capital letter followed by lowercase letters and digits, such as `E`, `Estar` " *
          "or `F2`.")
end

"""Error that the conformation label `c` also names a metabolite, which `kind` calls a
`"declared metabolite"` or a `"ligand"`: a form's text starts with its label, so a DSL
text would read the label as the metabolite. `who` names the constructor that checks."""
_label_collision_error(c::Symbol, who::String, kind::String) =
    error("$who: conformation label `$c` collides with $kind `$c`; choose a different " *
          "conformation label.")

"""
Error unless the labels of `reactants` (substrate and product names) and `modifiers`
(names that bind as modifiers) are uniquely decodable: no string splits into labels two
ways. Enzyme forms and their constants join the labels of the metabolites they bind
without a separator (`_form_name`), so with `Ac`, `CoA` and `AcCoA` the forms
E(Ac, CoA) and E(AcCoA) would share a name. Each name is its own label, but a reactant
that also binds as a modifier does so as its copy `Xmod`, written `X::Mod`; an allosteric
regulator binds no form, so its name takes no part. The Sardinas–Patterson test decides
it, and the error, which `who` begins, spells both splits.
"""
function check_decodable_names(reactants::Vector{Symbol}, modifiers::Vector{Symbol},
                                   who::String)
    spelled = Dict{String, String}()
    modifier_label(n) = n in reactants ? "$(n)mod" => "$n::Mod" : String(n) => String(n)
    for (label, written) in [[String(n) => String(n) for n in reactants];
                             modifier_label.(modifiers)]
        other = get!(spelled, label, written)
        other == written || _spelled_twice_error([other], [written], who)
    end
    words = sort!(collect(keys(spelled)))
    tail(w, prefix) = w[ncodeunits(prefix) + 1:end]
    # Each entry holds a dangling suffix `d` and two splits: join(long) == join(short) * d.
    queue = [(tail(v, u), [v], [u]) for u in words for v in words
             if u != v && startswith(v, u)]
    seen = Set{String}()
    while !isempty(queue)
        d, long, short = popfirst!(queue)
        d in seen && continue
        push!(seen, d)
        haskey(spelled, d) && _spelled_twice_error(getindex.(Ref(spelled), [short; d]),
                                                   getindex.(Ref(spelled), long), who)
        for w in words
            startswith(w, d) && push!(queue, (tail(w, d), [short; w], long))
            startswith(d, w) && push!(queue, (tail(d, w), long, [short; w]))
        end
    end
end

"""Error, which `who` begins, that the label splits `a` and `b`, each a list of DSL
spellings, spell one string."""
function _spelled_twice_error(a::Vector{String}, b::Vector{String}, who::String)
    spelled = join(replace.(a, "::Mod" => "mod"))
    error("$who: $(join(a, " + ")) and $(join(b, " + ")) both spell $spelled: " *
          "enzyme forms and their constants join the names of the metabolites they bind " *
          "without a separator, so two forms would share a name; rename a metabolite")
end

"""
    EnzymeReaction

The public concrete reaction descriptor: substrates, products, optional
regulators, and the set of catalytic multiplicities the mechanism
enumerator is allowed to consider. It is the entry point to the package —
pair it with rate data in an
[`IdentifyRateEquationProblem`](@ref EnzymeRates.IdentifyRateEquationProblem), or pass
it to `EnzymeRates.Enumeration.init_mechanisms` to enumerate candidate mechanisms.

# Fields
- `substrates::Vector{Reactant}` and `products::Vector{Reactant}` — each reactant
  carrying its per-atom inventory (used for ping-pong residual bookkeeping
  and atom conservation across steps).
- `regulators::Vector{Regulator}` — competitive inhibitors and
  allosteric regulators, each with its `kind` and `effect`.
- `allowed_catalytic_multiplicities::Vector{Int}` — the oligomeric states, one
  or more, that the allosteric enumerator may assign to the catalytic core.
- `shared_catalytic_site::Vector{Tuple{Symbol,Symbol}}` — `(substrate,
  product)` pairs that bind the same catalytic site, and so are never both
  bound to it at once.

Construct one with the [`@enzyme_reaction`](@ref EnzymeRates.DSL.@enzyme_reaction) DSL. Each
reactant takes an atom bracket (`S[C]`, `A[C1H1]`); the brackets are load-bearing for
ping-pong and multi-substrate reactions. Substrates, products and regulators are
sorted by name in the constructor, a name's competitive inhibitor before its
allosteric regulator, so two equivalent declarations compare equal under
`==`/`hash`. The constructor rejects a name listed as both a substrate
and a product: concentrations and constants are keyed by name. For the same
reason it rejects a reactant or regulator named `E_total`, `Keq`, `L`, `N`,
`params`, `concs` or `v`, names a rate equation keeps for its constants and the
local names of its code, a name holding `_`, since constant names join
names with `_`, the name `E`, which labels the free enzyme form, and a name that is
no identifier, is a keyword or reads in Julia code as another name (a micro sign as `μ`),
which no DSL text can write. It also rejects names
that spell one string two ways, such as
`Ac`, `CoA` and `AcCoA`: an enzyme form's name joins the names of the
metabolites it binds without a separator, so E with `Ac` and `CoA` bound and E
with `AcCoA` bound would share a name. A reaction prints as the `@enzyme_reaction`
block that builds it.

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
    allowed_catalytic_multiplicities::Vector{Int}
    shared_catalytic_site::Vector{Tuple{Symbol,Symbol}}

    function EnzymeReaction(substrates::Vector{Reactant}, products::Vector{Reactant},
                            regulators::Vector{Regulator},
                            allowed_catalytic_multiplicities::Vector{Int};
                            shared_catalytic_site = Tuple{Symbol,Symbol}[])
        sorted_substrates = sort(substrates; by = r -> r.name)
        sorted_products = sort(products; by = r -> r.name)
        sorted_regulators = sort(regulators; by = r -> (r.name, r.kind === :allosteric))
        # A nonempty list of distinct Ints ≥ 1, as a DSL text writes one.
        ms = allowed_catalytic_multiplicities
        !isempty(ms) && all(≥(1), ms) && allunique(ms) ||
            error("EnzymeReaction: allowed_catalytic_multiplicities must be a nonempty " *
                  "list of distinct Ints ≥ 1, got $ms")
        sorted_mults = sort(allowed_catalytic_multiplicities)

        sub_names  = Symbol[r.name for r in sorted_substrates]
        prod_names = Symbol[r.name for r in sorted_products]
        isempty(sub_names)  && error("EnzymeReaction: substrates must not be empty")
        isempty(prod_names) && error("EnzymeReaction: products must not be empty")

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
        check_decodable_names([sub_names; prod_names],
            Symbol[n for (n, kind) in reg_keys if kind === :competitive],
            "EnzymeReaction")

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
one formula, the competitive inhibitors and allosteric regulators, the catalytic
multiplicities unless they are the default 1, and the shared catalytic sites.
"""
function Base.show(io::IO, r::EnzymeReaction)
    mults = r.allowed_catalytic_multiplicities
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

# ─── The reaction graph ───────────────────────────────────────────────

"""
A metabolite as a form holds it or a step exchanges it: its name, and whether it binds
as a modifier (a name under `modifiers:` written bare, or a reactant's copy `A::Mod`). A
reactant and its modifier copy are two different `Met`s.
"""
struct Met
    name::Symbol
    modifier::Bool
end
Base.isless(x::Met, y::Met) = isless((x.name, x.modifier), (y.name, y.modifier))
Base.hash(x::Met, h::UInt) = hash(x.modifier, hash(x.name, hash(:Met, h)))

"""
    Node

An enzyme form: its conformation label, the metabolites it binds, and its residual: the
substrates it adds and the products it subtracts. Each list is sorted and holds a name as
often as the form does.
"""
struct Node
    conformation::Symbol
    bound::Vector{Met}
    added::Vector{Symbol}
    subtracted::Vector{Symbol}
end

"""
    Step

A step between nodes `from` and `to` of its mechanism, stored in its physical-forward
direction: going `from → to` it takes up `taken_up` and gives off `given_off`, each list
sorted. `rapid_equilibrium` marks a rapid-equilibrium step, and `group` is the step's
kinetic group.
"""
struct Step
    from::Int
    to::Int
    taken_up::Vector{Met}
    given_off::Vector{Met}
    rapid_equilibrium::Bool
    group::Int
end

"""
    Mechanism

A mechanism as its reaction graph, in canonical form: two writings of one mechanism give
equal fields. `substrates` and `products` are sorted; `nodes` hold free E first, then the
other nodes in content order (`isless(::Node)`); `steps` run physical-forward, in content
order (`isless(::Step)`), the kinetic groups numbered by their first step.

It also stores facts its canonical fields give, computed once when it is built:
`modifiers`, every name bound as a modifier, sorted; `groups`, each kinetic group's steps
in step order; and per step, `step_kinds` (`_step_kind`), `toward_products` (whether it
runs toward products as stored, `_runs_toward_products`) and `step_sides` (the names of
its two sides, which name its constants, `_sides`). They take no part in `==` and `hash`.

`Mechanism(substrates, products, nodes, steps)` builds one from nodes and steps written
any way and checks it; `_unchecked_mechanism` takes canonical fields and checks nothing,
for canonical form and the Sig decoder. Do not mutate the vectors of a built mechanism in
place: its canonical form, stored facts, hash, compiled Sig and memo entries depend on
them; build a new one with the constructor. Lift it to its compiled type with
`EnzymeRates.compile_mechanism(m)` or `EnzymeMechanism(m)`.
"""
struct Mechanism
    substrates::Vector{Symbol}
    products::Vector{Symbol}
    nodes::Vector{Node}
    steps::Vector{Step}
    modifiers::Vector{Symbol}
    groups::Vector{Vector{Int}}
    step_kinds::Vector{Symbol}
    toward_products::BitVector
    step_sides::Vector{Tuple{String, String}}
    # `graph` is partially initialised: it holds the canonical fields alone. The fact
    # helpers may read only `substrates`, `products`, `nodes` and `steps` of it.
    global function _unchecked_mechanism(substrates::Vector{Symbol},
                                         products::Vector{Symbol}, nodes::Vector{Node},
                                         steps::Vector{Step})
        graph = new(substrates, products, nodes, steps)
        new(substrates, products, nodes, steps, copy(_modifiers(nodes, steps)),
            _groups(steps), Symbol[_step_kind(graph, e) for e in steps],
            _runs_toward_products.(Ref(graph), steps), _step_sides(graph))
    end
end

"""
    AllostericMechanism

An MWC enzyme: its catalytic `Mechanism`, one allosteric tag per kinetic group of it
(`mwc_tags[g]`), the catalytic multiplicity and the regulatory sites, sorted. Its
catalytic graph is a plain `Mechanism`, so a model with another derivation can wrap the
same graph in a type of its own.

It also stores, one flag per catalytic step, the steps each conformation excludes
(`excluded_active`, the steps of `:OnlyI` groups; `excluded_inactive`, those of `:OnlyA`
groups) and the steps of each conformation (`active_steps`, `inactive_steps`:
`_conformation_steps`). They take no part in `==` and `hash`.
"""
struct AllostericMechanism
    catalytic::Mechanism
    mwc_tags::Vector{Symbol}
    catalytic_multiplicity::Int
    regulatory_sites::Vector{RegulatorySite}
    active_steps::BitVector
    inactive_steps::BitVector
    excluded_active::BitVector
    excluded_inactive::BitVector
    function AllostericMechanism(catalytic::Mechanism, mwc_tags::Vector{Symbol},
                                 catalytic_multiplicity::Int,
                                 regulatory_sites::Vector{RegulatorySite})
        _check_tags(mwc_tags, length(catalytic.groups), catalytic_multiplicity)
        _check_onlyi_groups(catalytic, mwc_tags)
        sites = sort(regulatory_sites; by = _site_key)
        _check_ligands(sites, catalytic)
        excluded_active = (e -> excludes(mwc_tags[e.group], :active)).(catalytic.steps)
        excluded_inactive = (e -> excludes(mwc_tags[e.group], :inactive)).(catalytic.steps)
        new(catalytic, copy(mwc_tags), catalytic_multiplicity, sites,
            _conformation_steps(catalytic, excluded_active),
            _conformation_steps(catalytic, excluded_inactive), excluded_active,
            excluded_inactive)
    end
end

"""
    excludes(tag::Symbol, c::Symbol) -> Bool

Whether conformation `c` (`:active` or `:inactive`) excludes a step or a regulator ligand
tagged `tag`: `:OnlyI` excludes the active conformation, `:OnlyA` the inactive one (spec
section 7). The one rule by which the constructor's masks, the regulator constants and the
moves' ligand states read the tags.
"""
excludes(tag::Symbol, c::Symbol) =
    c === :active ? tag === :OnlyI : c === :inactive ? tag === :OnlyA :
    throw(ArgumentError("a conformation is :active or :inactive, not :$c"))

"""A site's place in the order of a mechanism's sites: its ligands, its multiplicity,
then its ligands' tags."""
_site_key(s::RegulatorySite) = (s.ligands, s.multiplicity, s.mwc_tags)

abstract type AbstractEnzymeMechanism end

"""
    EnzymeMechanism{Sig}

A plain enzyme mechanism as a singleton type, the form `rate_equation` and the fitter
take. `Sig`, its only parameter, is the mechanism encoded as one tuple of integers
(`_sig_of`), so equal mechanisms share one type and the derivation compiles once per
mechanism. `@enzyme_mechanism` returns one; `EnzymeMechanism(m)` or
`EnzymeRates.compile_mechanism(m)` lifts a `Mechanism` to it and
`EnzymeRates.Mechanism(em)` reads the mechanism back. The layout of `Sig` is internal.
"""
struct EnzymeMechanism{Sig} <: AbstractEnzymeMechanism end

"""
    AllostericEnzymeMechanism{Sig}

An allosteric (MWC) enzyme mechanism as a singleton type: `Sig`, its only parameter,
is the mechanism encoded as one tuple of integers (`_sig_of`), its catalytic graph
first. `@allosteric_mechanism` returns one; `AllostericEnzymeMechanism(am)` or
`EnzymeRates.compile_mechanism(am)` lifts an `AllostericMechanism` to it and
`EnzymeRates.AllostericMechanism(aem)` reads the mechanism back.
"""
struct AllostericEnzymeMechanism{Sig} <: AbstractEnzymeMechanism end

# ─── Content order ──────────────────────────────────────────────────────

"""Nodes by label, then bound metabolites, residual additions and residual subtractions,
each list in lexicographic order, a list coming before every longer list it begins."""
function Base.isless(f::Node, g::Node)
    f.conformation === g.conformation || return isless(f.conformation, g.conformation)
    c = cmp(f.bound, g.bound)
    c == 0 && (c = cmp(f.added, g.added))
    c == 0 && (c = cmp(f.subtracted, g.subtracted))
    c < 0
end

"""Steps by node indices, then taken-up and given-off metabolites as lists (as
`isless(::Node)` orders them), then the flag, steady state first; the group takes no part.
Tuples and `argmin` over steps compare by it."""
function Base.isless(e::Step, f::Step)
    e.from != f.from && return e.from < f.from
    e.to != f.to && return e.to < f.to
    c = cmp(e.taken_up, f.taken_up)
    c == 0 && (c = cmp(e.given_off, f.given_off))
    c == 0 ? e.rapid_equilibrium < f.rapid_equilibrium : c < 0
end

"""`e` stored in the opposite direction: nodes swapped, taken-up and given-off swapped."""
_reversed(e::Step) =
    Step(e.to, e.from, e.given_off, e.taken_up, e.rapid_equilibrium, e.group)

# ─── Canonical form ─────────────────────────────────────────────────────

"""
    _canonicalize(substrates, products, nodes, steps)
        -> (; mechanism, edge_of, flipped, group_of, node_of, misaligned)

The canonical `Mechanism` that `nodes` and `steps` write, with no representation check:
the nodes in any order, a node any number of times, each step between two of them in
either direction, its `group` a written kinetic group label (1, 2, …, each used), and
every list of a node or a step sorted. `edge_of[j]` is the step written step `j` became,
`flipped[j]` whether that step runs against step `j` as written, `group_of[l]` the
kinetic group written label `l` became, `node_of[i]` the node written node `i` became,
and `misaligned` the lead and member of a tied kinetic group that no orientation aligns,
or no step. Alignment stops at that group, so a result with a misaligned pair is not
canonical. The checked constructors build on it. The nodes and steps share the lists of
`nodes` and `steps`.
"""
function _canonicalize(substrates::Vector{Symbol}, products::Vector{Symbol},
                       nodes::Vector{Node}, steps::Vector{Step})
    distinct, node_of = _sorted_nodes(nodes)
    written = [Step(node_of[e.from], node_of[e.to], e.taken_up, e.given_off,
                    e.rapid_equilibrium, e.group) for e in steps]
    substrates = sort(substrates)
    fwd, misaligned = _orientation(distinct, written, substrates)
    stored = [fwd[j] ? written[j] : _reversed(written[j]) for j in eachindex(written)]
    p = sortperm(stored)
    group, group_of = _number_groups([stored[i].group for i in p],
                                     maximum(e -> e.group, steps; init = 0))
    canonical = [Step(e.from, e.to, e.taken_up, e.given_off, e.rapid_equilibrium, group[k])
                 for (k, e) in enumerate(view(stored, p))]
    m = _unchecked_mechanism(substrates, sort(products), distinct, canonical)
    (; mechanism = m, edge_of = invperm(p), flipped = .!fwd, group_of, node_of,
     misaligned)
end

"""The distinct nodes of `written` in canonical order, free E (label `E`, nothing bound,
no residual) first, and the index each written node became. A search keeps every
mechanism it emits, so the nodes, built in a vector grown past its length, are stored as
a copy at its length."""
function _sorted_nodes(written::Vector{Node})
    nodes, node_of = Node[], Vector{Int}(undef, length(written))
    for t in sortperm(written; by = n -> (!_is_free_E(n), n))
        (isempty(nodes) || nodes[end] != written[t]) && push!(nodes, written[t])
        node_of[t] = length(nodes)
    end
    copy(nodes), node_of
end

"""Canonical form's tier-1 side order: substrates, bound plus free, then fewer products.
`free` lists the metabolites a step takes up or gives off at node `n`; a modifier counts on
neither side, and every other name that is no substrate counts as a product."""
function _side_score(n::Node, free::Vector{Met}, substrates::Vector{Symbol})
    s = p = 0
    for list in (n.bound, free), x in list
        x.modifier && continue
        x.name in substrates ? (s += 1) : (p += 1)
    end
    (s, -p)
end

"""
    _orientation(nodes, written, substrates) -> (fwd::BitVector, misaligned::Vector{Step})

Whether each written step runs in its physical-forward direction. A step with one
metabolite in all takes it up (tier 0). Otherwise the side holding more substrates,
bound plus free, then fewer products, comes first (tier 1); otherwise the step runs from
a node where only products enter or leave to one where only substrates do (tier 2);
otherwise from the earlier node (tier 3). A modifier counts on neither side. Entry
sides are read over every step as written; reversing a step swaps its ends and lists
together, so they do not depend on any written direction. Tied groups then turn as one
(`_align_tied_groups!`).
"""
function _orientation(nodes::Vector{Node}, written::Vector{Step},
                      substrates::Vector{Symbol})
    substrate(x::Met) = !x.modifier && x.name in substrates
    product(x::Met) = !x.modifier && !(x.name in substrates)
    side(list::Vector{Met}) = UInt8(any(substrate, list) + 2 * any(product, list))
    enters = zeros(UInt8, length(nodes))
    for e in written
        enters[e.from] |= side(e.taken_up)
        enters[e.to] |= side(e.given_off)
    end
    n = length(written)
    fwd, tier = trues(n), zeros(Int, n)
    for (k, e) in enumerate(written)
        a = _side_score(nodes[e.from], e.taken_up, substrates)
        b = _side_score(nodes[e.to], e.given_off, substrates)
        if length(e.taken_up) + length(e.given_off) == 1
            fwd[k] = !isempty(e.taken_up)
        elseif a != b
            fwd[k], tier[k] = a > b, 1
        elseif (enters[e.from], enters[e.to]) in ((0x02, 0x01), (0x01, 0x02))
            fwd[k], tier[k] = enters[e.from] == 0x02, 2
        else
            fwd[k], tier[k] = e.from < e.to, 3
        end
    end
    misaligned = _align_tied_groups!(fwd, tier, nodes, written)
    fwd, misaligned
end

"""
Turn the steps tiers 2 and 3 decide in each kinetic group (each written step's `group`
label) with the group's lead: the tied step that changes conformation first, then by
tier, then by its oriented step. A group's steps share their constants, so a
conformational change and its mirror at other bound nodes must run the same way. Every
other tied step that takes up and gives off the lead's metabolites takes the orientation
that does, and, if it changes conformation, makes the lead's change; a step that fits both
ways keeps its own. Returns the lead and the first step that fits the lead's metabolites
but not its change, or no step.
"""
function _align_tied_groups!(fwd::BitVector, tier::Vector{Int}, nodes::Vector{Node},
                             written::Vector{Step})
    label_pair(e::Step) = (nodes[e.from].conformation, nodes[e.to].conformation)
    turns(e::Step) = label_pair(e)[1] !== label_pair(e)[2]
    oriented(k::Int) = fwd[k] ? written[k] : _reversed(written[k])
    for l in unique(e.group for e in written)
        tied = [k for k in eachindex(written) if written[k].group == l && tier[k] > 1]
        length(tied) < 2 && continue
        lead = oriented(argmin(k -> (!turns(oriented(k)), tier[k], oriented(k)), tied))
        for k in tied
            o = oriented(k)
            same = [x for x in (o, _reversed(o))
                    if (x.taken_up, x.given_off) == (lead.taken_up, lead.given_off)]
            isempty(same) && continue
            i = findfirst(x -> !turns(x) || label_pair(x) == label_pair(lead), same)
            i === nothing && return [lead, same[1]]
            fwd[k] = same[i] == written[k]
        end
    end
    Step[]
end

"""Group labels in step order renumbered by first occurrence, and the number each
written label took."""
function _number_groups(labels::Vector{Int}, G::Int)
    group_of, group = zeros(Int, G), similar(labels)
    next = 0
    for (k, l) in enumerate(labels)
        if group_of[l] == 0
            next += 1
            group_of[l] = next
        end
        group[k] = group_of[l]
    end
    group, group_of
end

# ─── Constructors and their checks ──────────────────────────────────────

"""
    Mechanism(substrates, products, nodes::Vector{Node}, steps::Vector{Step})

The mechanism that `nodes` and `steps` write, in canonical form, checked. The nodes come
in any order, a node any number of times, and each list of a node or a step in any
order; a step joins two of the nodes (`from` and `to` index `nodes`) in either direction,
and its `group` is a kinetic group label, the labels running 1, 2, … with each used.
Every name a node or a step binds as a modifier is a modifier. It checks the
declarations (`_check_declarations`), then the representation (`_check_representation`).
Whether its rate law is defined is the boundary check's question, not the constructor's.
"""
Mechanism(substrates::Vector{Symbol}, products::Vector{Symbol}, nodes::Vector{Node},
          steps::Vector{Step}) = _checked(substrates, products, nodes, steps).mechanism

"""
    AllostericMechanism(substrates, products, nodes, steps, mwc_tags,
                        catalytic_multiplicity, regulatory_sites)

The allosteric mechanism whose catalytic part `nodes` and `steps` write, as
`Mechanism(substrates, products, nodes, steps)` builds and checks it, with
`mwc_tags[l]` the tag of written group label `l`: it checks the tags, then the catalytic
part, then moves each tag to the canonical group its label became and composes.
"""
function AllostericMechanism(substrates::Vector{Symbol}, products::Vector{Symbol},
                             nodes::Vector{Node}, steps::Vector{Step},
                             mwc_tags::Vector{Symbol}, catalytic_multiplicity::Int,
                             regulatory_sites::Vector{RegulatorySite})
    _check_tags(mwc_tags, maximum(e -> e.group, steps; init = 0), catalytic_multiplicity)
    c = _checked(substrates, products, nodes, steps)
    AllostericMechanism(c.mechanism, _tags_by_group(c.group_of, mwc_tags),
                        catalytic_multiplicity, regulatory_sites)
end

"""The canonical form (`_canonicalize`) of the mechanism that `nodes` and `steps` write,
checked: each step's ends index `nodes`, each list is sorted into a fresh vector, the
names bound as modifiers are taken as its modifiers, then its declarations and its
representation are checked."""
function _checked(substrates::Vector{Symbol}, products::Vector{Symbol},
                  nodes::Vector{Node}, steps::Vector{Step})
    for (k, e) in enumerate(steps)
        e.from in eachindex(nodes) && e.to in eachindex(nodes) && continue
        error("Step $k: its ends are nodes $(e.from) and $(e.to), but there are " *
              "$(length(nodes)) nodes; a step's ends index the nodes")
    end
    nodes = [Node(n.conformation, sort(n.bound), sort(n.added), sort(n.subtracted))
             for n in nodes]
    steps = [Step(e.from, e.to, sort(e.taken_up), sort(e.given_off), e.rapid_equilibrium,
                  e.group) for e in steps]
    modifiers = _modifiers(nodes, steps)
    _check_declarations(substrates, products, modifiers, Int[e.group for e in steps])
    c = _canonicalize(substrates, products, nodes, steps)
    _check_representation(c.mechanism, c.misaligned, modifiers)
    c
end

"""The tags by canonical group: `mwc_tags[l]`, the tag of written label `l`, moved to
group `group_of[l]`."""
function _tags_by_group(group_of::Vector{Int}, mwc_tags::Vector{Symbol})
    tags = similar(mwc_tags)
    for (l, g) in enumerate(group_of)
        tags[g] = mwc_tags[l]
    end
    tags
end

# ─── Checks ─────────────────────────────────────────────────────────────

"""Error unless `mwc_tags` holds one valid tag (`_TAG_CODES`) for each of `G`
kinetic groups and the catalytic multiplicity is at least 1."""
function _check_tags(mwc_tags::Vector{Symbol}, G::Int, catalytic_multiplicity::Int)
    length(mwc_tags) == G ||
        error("AllostericMechanism: mwc_tags length $(length(mwc_tags)) " *
              "must match the number of kinetic groups $G")
    catalytic_multiplicity ≥ 1 ||
        error("AllostericMechanism: catalytic_multiplicity must be ≥ 1, " *
              "got $catalytic_multiplicity")
    bad = findfirst(∉(_TAG_CODES), mwc_tags)
    bad === nothing ||
        error("AllostericMechanism: catalytic group $bad has invalid " *
              "MWC tag $(mwc_tags[bad]) (must be one of $_TAG_CODES)")
    nothing
end

"""Error when a kinetic group of `catalytic` tagged `:OnlyI` (`mwc_tags`) holds a
step that does more than take up or give off modifiers, or whose two nodes differ in
anything but their modifiers and label: only the active conformation catalyzes, so a
group it lacks binds modifiers alone (spec section 7)."""
function _check_onlyi_groups(catalytic::Mechanism, mwc_tags::Vector{Symbol})
    reactants(list::Vector{Met}) = filter(x -> !x.modifier, list)
    for (k, e) in enumerate(catalytic.steps)
        mwc_tags[e.group] === :OnlyI || continue
        f, t = catalytic.nodes[e.from], catalytic.nodes[e.to]
        !_is_iso(e) && all(x -> x.modifier, e.taken_up) &&
            all(x -> x.modifier, e.given_off) &&
            (reactants(f.bound), f.added, f.subtracted) ==
            (reactants(t.bound), t.added, t.subtracted) && continue
        error("AllostericMechanism: $(step_text(catalytic, k)) is tagged :OnlyI but " *
              "does more than bind a modifier. Only the active conformation " *
              "catalyzes, so :OnlyI fits only a group whose steps take up or give off " *
              "modifiers alone. Tag its group :NonequalAI.")
    end
    nothing
end

"""Error when a site ligand cannot name a metabolite (`_check_metabolite_names`) or is a
conformation label of `catalytic`, or when a ligand sits on two sites: each ligand's
constant is named after the ligand alone."""
function _check_ligands(sites::Vector{RegulatorySite}, catalytic::Mechanism)
    ligs = Symbol[l for s in sites for l in s.ligands]
    _check_metabolite_names("AllostericMechanism", ligs)
    for l in ligs
        any(n -> n.conformation === l, catalytic.nodes) &&
            _label_collision_error(l, "AllostericMechanism", "ligand")
    end
    dup = findfirst(i -> ligs[i] in view(ligs, 1:(i - 1)), eachindex(ligs))
    dup === nothing ||
        error("AllostericMechanism: ligand $(ligs[dup]) appears in two " *
              "distinct regulatory sites; rendered Kreg names would " *
              "collide. Same-ligand-two-sites is not enumerated.")
    nothing
end

"""Error unless `substrates` and `products` are not empty, the declared `substrates` and
`products` name each name once, no name is both a substrate and a product, `labels`, one
per step, are kinetic group labels from 1 to the number of groups, every label used, and
every declared name and every name of `modifiers` is one a metabolite can take
(`_check_metabolite_names`)."""
function _check_declarations(substrates::Vector{Symbol}, products::Vector{Symbol},
                             modifiers::Vector{Symbol}, labels::Vector{Int})
    isempty(substrates) && error("Mechanism: substrates must not be empty")
    isempty(products) && error("Mechanism: products must not be empty")
    for (names, kind) in ((substrates, "substrate"), (products, "product"))
        allunique(names) || error("Mechanism: duplicate $kind names")
    end
    both = findfirst(in(products), substrates)
    both === nothing ||
        error("Mechanism: $(substrates[both]) named as both a substrate and a " *
              "product; concentrations and constants are keyed by name")
    all(≥(1), labels) ||
        error("Mechanism: each step takes one kinetic group label from 1 to the number " *
              "of groups")
    empty = findfirst(∉(labels), 1:maximum(labels; init = 0))
    empty === nothing ||
        error("Mechanism: kinetic group $empty is empty; a kinetic group holds at least " *
              "one step")
    for names in (substrates, products, modifiers)
        _check_metabolite_names("Mechanism", names)
    end
end

"""
Error unless the canonical mechanism `m` represents one: conformation labels shaped as
labels (`check_conformation_label`), none a metabolite's name; every name declared (in
`modifiers`, the names bound as modifiers, or among its substrates and products), a plain
entry only for a reactant, residuals that add substrates and subtract products; no step
joining a node to itself, taking up and giving off one metabolite, or binding a
metabolite into a node that holds it in another role; tied steps that make one
conformational change (`misaligned`, from `_canonicalize`); no two nodes that render one
name; each kinetic group one kind of step with one flag; each reaction once; free E; one
connected graph; and a step that exchanges each declared name and each name of
`modifiers`, a modifier as a modifier.
"""
function _check_representation(m::Mechanism, misaligned::Vector{Step},
                               modifiers::Vector{Symbol})
    table = _table(m)
    for c in sort!(unique!([n.conformation for n in m.nodes]))
        check_conformation_label(c, "Mechanism")
        any(x -> x.name === c, table) &&
            _label_collision_error(c, "Mechanism", "declared metabolite")
    end
    for x in table
        _is_reactant(m, x.name) || x.name in modifiers ||
            error("Mechanism: $(x.name) is not declared; declare it under substrates:, " *
                  "products: or modifiers:")
        x.modifier || _is_reactant(m, x.name) ||
            error("Mechanism: $(x.name) is a modifier, so it binds only as a modifier")
    end
    for n in m.nodes, (list, verb, wrong) in ((n.added, "adds", ∉(m.substrates)),
                                              (n.subtracted, "subtracts", ∈(m.substrates)))
        i = findfirst(wrong, list)
        i === nothing ||
            error("Mechanism: a residual adds substrates and subtracts products; one " *
                  "$verb $(list[i])")
    end
    for e in m.steps
        e.from == e.to &&
            error("Step: both ends are the form $(_form_name(m, e.from)); a step joins " *
                  "two forms")
        shared = unique!([x for x in e.taken_up if x in e.given_off])
        isempty(shared) ||
            error("Step $(_form_name(m, e.from)) → $(_form_name(m, e.to)): " *
                  "$(join([x.name for x in shared], ", ")) both " *
                  "consumed and released")
        _binds_by_name(m, e) && !_is_plain_binding(m, e) &&
            error("Step $(_form_name(m, e.from)) → $(_form_name(m, e.to)) changes a " *
                  "metabolite's role: $(_form_name(m, e.to)) holds " *
                  "$(_form_name(m, e.from))'s metabolites plus " *
                  "$(only(e.taken_up).name) by name, not by role; write a " *
                  "modifier copy as X::Mod on both sides")
    end
    _check_aligned(m, misaligned)
    _check_distinct_form_names(m)
    _check_uniform_groups(m, m.groups)
    _check_each_reaction_once(m, m.groups)
    _check_free_E(m)
    _check_connected(m)
    _check_declared_used(m, modifiers)
end

"""Error unless node 1 of `m` is free E: label `E`, nothing bound, no residual. Canonical
form puts free E first wherever a mechanism holds it."""
function _check_free_E(m::Mechanism)
    !isempty(m.nodes) && _is_free_E(m.nodes[1]) ||
        error("Mechanism: no free enzyme form E; every mechanism holds the free form E, " *
              "with nothing bound and no residual")
end

"""Error unless the steps of `m` join every node to free E, node 1: nodes that never
meet free E split the enzyme between parts with no single steady state."""
function _check_connected(m::Mechanism)
    f = findfirst(!, _joined_to_free_E(m, falses(length(m.steps))))
    f === nothing ||
        error("Mechanism: the steps do not form one connected graph; no steps lead from " *
              "free E to $(_form_name(m, f))")
end

"""Error unless some step of `m` takes up or gives off each declared substrate and
product, and some step takes up or gives off each modifier `modifiers` declares as a
modifier: a declared name that no step exchanges takes no part in the rate law."""
function _check_declared_used(m::Mechanism, modifiers::Vector{Symbol})
    used(x::Met) = any(e -> x in e.taken_up || x in e.given_off, m.steps)
    for (names, kind) in ((m.substrates, "substrate"), (m.products, "product"))
        k = findfirst(n -> !used(Met(n, false)), names)
        k === nothing ||
            error("Mechanism: $kind $(names[k]) is declared, but no step uses it")
    end
    k = findfirst(n -> !used(Met(n, true)), modifiers)
    k === nothing ||
        error("Mechanism: modifier $(modifiers[k]) is bound at a node, " *
              "but no step binds it")
end

"""Error unless `misaligned`, from `_canonicalize`, is empty: it holds the lead and a member
of a tied kinetic group that no orientation aligns."""
function _check_aligned(m::Mechanism, misaligned::Vector{Step})
    isempty(misaligned) ||
        error("Mechanism: a kinetic group holds $(join(_sides(m, misaligned[1]), " → ")) " *
              "and $(join(_sides(m, misaligned[2]), " → ")), which make different " *
              "conformational changes; the steps of a kinetic group share their " *
              "constants, so they must make the same conformational change")
    nothing
end

"""Error when two nodes of `m` render one name: every constant is named after the nodes
its step joins (`_form_name`)."""
function _check_distinct_form_names(m::Mechanism)
    seen = Dict{String, Int}()
    for f in eachindex(m.nodes)
        other = get!(seen, _form_name(m, f), f)
        other == f ||
            error("Mechanism: the enzyme forms $(form_text(m, other)) and " *
                  "$(form_text(m, f)) both render the name $(_form_name(m, f)); every " *
                  "constant is named after the forms its step joins, so rename a " *
                  "metabolite to tell the two forms apart")
    end
end

"""Error unless the steps of each kinetic group (`groups`) take up and give off the same
metabolites and carry one flag: they share their constants."""
function _check_uniform_groups(m::Mechanism, groups::Vector{Vector{Int}})
    label(e::Step) = join(_sides(m, e), " → ") * (e.rapid_equilibrium ? " (RE)" : " (SS)")
    for members in groups
        e1 = m.steps[first(members)]
        for k in members
            e = m.steps[k]
            (e.taken_up, e.given_off, e.rapid_equilibrium) ==
                (e1.taken_up, e1.given_off, e1.rapid_equilibrium) && continue
            error("Mechanism: a kinetic group holds $(label(e1)) and $(label(e)), " *
                  "which differ in kind or in RE/SS; the steps of a kinetic group " *
                  "share their constants, so they must be the same kind of step " *
                  "with the same flag")
        end
    end
end

"""Error when one reaction, the unordered pair of a step's two sides, appears in two
steps: a reaction belongs to one kinetic group, whose constants are named after it, and
one group holding it twice would count its step twice. Reads the groups (`groups`) in
order, each in step order, and names the reaction as the step seen first stores it."""
function _check_each_reaction_once(m::Mechanism, groups::Vector{Vector{Int}})
    seen = Dict{NTuple{2, Tuple{Int, Vector{Met}}}, Tuple{Int, Int}}()
    flag(e::Step) = e.rapid_equilibrium ? "RE" : "SS"
    for (g, members) in enumerate(groups), k in members
        e = m.steps[k]
        key = minmax((e.from, e.taken_up), (e.to, e.given_off))
        if haskey(seen, key)
            g0, k0 = seen[key]
            text = join(m.step_sides[k0], " ⇌ ")
            g0 == g && error("Mechanism: kinetic group $g holds the reaction $text " *
                             "twice; a mechanism holds each reaction once")
            e0 = m.steps[k0]
            flags = flag(e0) == flag(e) ? "" :
                " ($(flag(e0)) in group $g0, $(flag(e)) in group $g)"
            error("Mechanism: kinetic groups $g0 and $g both hold the reaction " *
                  "$text$flags; a reaction belongs to one kinetic group, whose " *
                  "constants are named after it")
        end
        seen[key] = (g, k)
    end
end

# ─── Identity ───────────────────────────────────────────────────────────

# `==` and `hash` by content for every value type whose fields hold vectors or structs,
# whose default `==` would compare identity: each field in declaration order, `hash`
# seeded with the type's name. A mechanism's stored facts take no part: they are
# functions of its other fields. A struct's default hash mixes in its type's identity,
# which changes between precompile sessions; a value's hash reads only symbols, integers
# and flags, so every process agrees on it. `Met` keeps Julia's `==`,
# which compares its two fields, and has its own `hash`.
for T in (RegulatorySite, Reactant, Regulator, EnzymeReaction, Node, Step, Mechanism,
          AllostericMechanism)
    fs = T === Mechanism ? (:substrates, :products, :nodes, :steps) :
         T === AllostericMechanism ?
         (:catalytic, :mwc_tags, :catalytic_multiplicity, :regulatory_sites) : fieldnames(T)
    @eval Base.:(==)(a::$T, b::$T) =
        $(reduce((l, r) -> :($l && $r), [:(a.$f == b.$f) for f in fs]))
    @eval Base.hash(x::$T, h::UInt) =
        $(foldl((acc, f) -> :(hash(x.$f, $acc)), fs;
                init = :(hash($(QuoteNode(nameof(T))), h))))
end

# ─── Sig ────────────────────────────────────────────────────────────────
#
# A mechanism's Sig, every entry an `Int`:
#   plain       its names: their number, then each name, the sorted distinct names of its
#               reactants, metabolites and conformation labels; the substrates and the
#               products; the number of nodes, then each node (label, bound, added,
#               subtracted); the number of steps, then each step (from, to, taken up,
#               given off, flag, group). After the names, a name is its index among them,
#               a list is its length followed by its entries, and a metabolite is its name
#               followed by its modifier flag.
#   name        its number of characters, then each character's code point.
#   allosteric  the plain Sig of its catalytic mechanism; one tag code per group; the
#               catalytic multiplicity; the number of sites, then per site its
#               multiplicity, its number of ligands and per ligand its name and tag code.

"""The allosteric tags, coded by position: catalytic groups and site ligands share them,
and the constructors accept no other (`_check_tags`, `RegulatorySite`)."""
const _TAG_CODES = (:OnlyA, :OnlyI, :EqualAI, :NonequalAI)
_tag_code(t::Symbol) = findfirst(==(t), _TAG_CODES)::Int

function _sig_name!(out::Vector{Int}, s::Symbol)
    str = String(s)
    push!(out, length(str))
    for c in str
        push!(out, Int(c))
    end
    out
end

"""The sorted distinct names of `m`: its reactants, its metabolite names and its
conformation labels."""
_sig_names(m::Mechanism) =
    sort!(unique!([m.substrates; m.products; [x.name for x in _table(m)];
                   [n.conformation for n in m.nodes]]))

function _sig_mechanism!(out::Vector{Int}, m::Mechanism)
    names = _sig_names(m)
    push!(out, length(names))
    foreach(s -> _sig_name!(out, s), names)
    index(s::Symbol) = searchsortedfirst(names, s)
    function symbols!(list::Vector{Symbol})
        push!(out, length(list))
        foreach(s -> push!(out, index(s)), list)
    end
    function mets!(list::Vector{Met})
        push!(out, length(list))
        foreach(x -> push!(out, index(x.name), x.modifier), list)
    end
    symbols!(m.substrates)
    symbols!(m.products)
    push!(out, length(m.nodes))
    for n in m.nodes
        push!(out, index(n.conformation))
        mets!(n.bound)
        symbols!(n.added)
        symbols!(n.subtracted)
    end
    push!(out, length(m.steps))
    for e in m.steps
        push!(out, e.from, e.to)
        mets!(e.taken_up)
        mets!(e.given_off)
        push!(out, e.rapid_equilibrium, e.group)
    end
    out
end

"""`m` as one homogeneous integer tuple, the only type parameter of its compiled type."""
_sig_of(m::Mechanism) = Tuple(_sig_mechanism!(Int[], m))
function _sig_of(am::AllostericMechanism)
    out = _sig_mechanism!(Int[], am.catalytic)
    foreach(t -> push!(out, _tag_code(t)), am.mwc_tags)
    push!(out, am.catalytic_multiplicity, length(am.regulatory_sites))
    for s in am.regulatory_sites
        push!(out, s.multiplicity, length(s.ligands))
        for (l, t) in zip(s.ligands, s.mwc_tags)
            _sig_name!(out, l)
            push!(out, _tag_code(t))
        end
    end
    Tuple(out)
end

"""A Sig read from its start: a struct, not captured variables, so nothing is boxed."""
mutable struct _SigReader
    sig::Vector{Int}
    at::Int
end
_next!(r::_SigReader) = r.sig[r.at += 1]
function _next_name!(r::_SigReader)
    n = _next!(r)
    Symbol(String([Char(_next!(r)) for _ in 1:n]))
end

function _read_mechanism!(r::_SigReader)
    names = [_next_name!(r) for _ in 1:_next!(r)]
    symbols() = Symbol[names[_next!(r)] for _ in 1:_next!(r)]
    mets() = Met[Met(names[_next!(r)], _next!(r) == 1) for _ in 1:_next!(r)]
    substrates = symbols()
    products = symbols()
    nodes = [Node(names[_next!(r)], mets(), symbols(), symbols()) for _ in 1:_next!(r)]
    steps = [Step(_next!(r), _next!(r), mets(), mets(), _next!(r) == 1, _next!(r))
             for _ in 1:_next!(r)]
    _unchecked_mechanism(substrates, products, nodes, steps)
end

function _read_site!(r::_SigReader)
    multiplicity = _next!(r)
    n = _next!(r)
    ligands, states = Symbol[], Symbol[]
    for _ in 1:n
        push!(ligands, _next_name!(r))
        push!(states, _TAG_CODES[_next!(r)])
    end
    RegulatorySite(ligands, multiplicity, states)
end

"""The mechanism whose Sig is `sig`, splatted. Decoding trusts the Sig, which only
`_sig_of` writes, so it canonicalizes and checks nothing. The varargs method compiles
once, whatever the Sig length, as long as callers splat a Sig whose length their code
leaves open."""
_mechanism_from_sig(sig::Int...) = _read_mechanism!(_SigReader(collect(sig), 0))
"""The allosteric mechanism whose Sig is `sig`, splatted; only the composition
constructor's tag and site checks run."""
function _allosteric_mechanism_from_sig(sig::Int...)
    r = _SigReader(collect(sig), 0)
    m = _read_mechanism!(r)
    tags = [_TAG_CODES[_next!(r)] for _ in 1:length(m.groups)]
    multiplicity = _next!(r)
    sites = [_read_site!(r) for _ in 1:_next!(r)]
    AllostericMechanism(m, tags, multiplicity, sites)
end

"""Lift a `Mechanism` to its singleton `EnzymeMechanism` type (`_sig_of`)."""
EnzymeMechanism(m::Mechanism) = EnzymeMechanism{_sig_of(m)}()
"""Lift an `AllostericMechanism` to its singleton `AllostericEnzymeMechanism` type."""
AllostericEnzymeMechanism(am::AllostericMechanism) =
    AllostericEnzymeMechanism{_sig_of(am)}()

# The lifts back read the Sig at run time, so they compile once, whatever the mechanism.
# Each gives back the mechanism its Sig encodes.
"""The `Mechanism` the compiled mechanism `em` encodes."""
Mechanism(@nospecialize(em::EnzymeMechanism)) =
    _mechanism_from_sig(typeof(em).parameters[1]...)
"""The `AllostericMechanism` the compiled mechanism `aem` encodes."""
AllostericMechanism(@nospecialize(aem::AllostericEnzymeMechanism)) =
    _allosteric_mechanism_from_sig(typeof(aem).parameters[1]...)
"""Lift a singleton mechanism to its concrete `Mechanism` or `AllostericMechanism`."""
_concrete(@nospecialize(em::EnzymeMechanism)) = Mechanism(em)
_concrete(@nospecialize(aem::AllostericEnzymeMechanism)) = AllostericMechanism(aem)

"""
    compile_mechanism(m::Mechanism)
    compile_mechanism(am::AllostericMechanism)

Convert a `Mechanism` to an `EnzymeMechanism`, or an
`AllostericMechanism` to an `AllostericEnzymeMechanism`.
"""
compile_mechanism(m::Mechanism) = EnzymeMechanism(m)
compile_mechanism(am::AllostericMechanism) = AllostericEnzymeMechanism(am)

# ─── Reads ──────────────────────────────────────────────────────────────

"""
Every metabolite entry `m` holds, sorted by (name, modifier) and once each: the
metabolites its nodes bind and its steps take up and give off, and its residual names as
non-modifier entries. A reactant and its modifier copy are two entries.
"""
function _table(m::Mechanism)
    table = Met[]
    for n in m.nodes
        append!(table, n.bound)
        foreach(x -> push!(table, Met(x, false)), n.added)
        foreach(x -> push!(table, Met(x, false)), n.subtracted)
    end
    for e in m.steps
        append!(table, e.taken_up)
        append!(table, e.given_off)
    end
    unique!(sort!(table))
end
"""Every name `nodes` and `steps` bind as a modifier, sorted, each once."""
function _modifiers(nodes::Vector{Node}, steps::Vector{Step})
    names = Symbol[]
    add(list::Vector{Met}) =
        foreach(x -> x.modifier && !(x.name in names) && push!(names, x.name), list)
    foreach(n -> add(n.bound), nodes)
    foreach(e -> (add(e.taken_up); add(e.given_off)), steps)
    sort!(names)
end
# The concentration names in the order the rate law reads them: substrates, products,
# then the modifier-only names, each by name; an allosteric mechanism appends each site's
# ligands in site order, once each.
metabolites(m::Mechanism) = unique!([m.substrates; m.products; m.modifiers])
metabolites(am::AllostericMechanism) =
    unique!([metabolites(am.catalytic);
             [l for s in am.regulatory_sites for l in s.ligands]])
"""Whether `x` names a substrate or a product of `m`."""
_is_reactant(m::Mechanism, x::Symbol) = insorted(x, m.substrates) || insorted(x, m.products)

"""The ligands of `am`'s regulatory sites, site by site, each once."""
allosteric_regulators(am::AllostericMechanism) =
    unique(Symbol[l for s in am.regulatory_sites for l in s.ligands])

catalytic(m::Mechanism) = m
catalytic(am::AllostericMechanism) = am.catalytic
_has_residual(n::Node) = !isempty(n.added) || !isempty(n.subtracted)
"""Whether step `e` takes up one metabolite and gives off none."""
_is_binding(e::Step) = isempty(e.given_off) && length(e.taken_up) == 1
"""Whether step `e` takes up and gives off nothing: an isomerization."""
_is_iso(e::Step) = isempty(e.taken_up) && isempty(e.given_off)
"""Whether binding `e` joins two nodes of `m` that carry one residual."""
function _keeps_residual(m::Mechanism, e::Step)
    f, t = m.nodes[e.from], m.nodes[e.to]
    f.added == t.added && f.subtracted == t.subtracted
end
"""Whether the sorted list `t` is the sorted list `f` with one more `x`."""
function _holds_one_more(t::Vector{Met}, f::Vector{Met}, x::Met)
    length(t) == length(f) + 1 || return false
    j, skipped = 1, false
    for y in t
        if j <= length(f) && y == f[j]
            j += 1
        elseif !skipped && y == x
            skipped = true
        else
            return false
        end
    end
    skipped
end
"""Whether binding `e` takes up its metabolite into a node that holds its source's
metabolites plus that one, with the same residual; the label may change."""
_is_plain_binding(m::Mechanism, e::Step) =
    _is_binding(e) && _keeps_residual(m, e) &&
    _holds_one_more(m.nodes[e.to].bound, m.nodes[e.from].bound, only(e.taken_up))
"""Whether binding `e` takes up its metabolite into a node that holds its source's
metabolites plus that one by name, with the same residual: a plain binding, or one that
binds a reactant where its modifier copy sits, or the reverse."""
function _binds_by_name(m::Mechanism, e::Step)
    _is_binding(e) && _keeps_residual(m, e) || return false
    names(n::Int) = Symbol[x.name for x in m.nodes[n].bound]
    names(e.to) == sort!(push!(names(e.from), only(e.taken_up).name))
end
"""
Whether step `e` of `m` runs, as stored, toward products: from the side that holds more
substrates, bound or free, then fewer products, the order canonical form gives a step that
exchanges several metabolites (`_orientation`). A step that takes up or gives off one
metabolite is stored as its uptake, so a half-reaction that releases a product,
`E(A) <--> E(; residual = A - P) + P`, runs toward products against its stored direction.
A step whose two sides tie, such as an isomerization of free forms, runs as stored.
"""
_runs_toward_products(m::Mechanism, e::Step) =
    _side_score(m.nodes[e.from], e.taken_up, m.substrates) >=
    _side_score(m.nodes[e.to], e.given_off, m.substrates)
"""Whether `n` is a free node: nothing bound and no residual, whatever its label."""
_is_free(n::Node) = isempty(n.bound) && !_has_residual(n)
"""Whether `n` is free E: label `E`, nothing bound and no residual."""
_is_free_E(n::Node) = n.conformation === :E && _is_free(n)

# ─── Accessors ─────────────────────────────────────────────────
#
# A compiled mechanism is read by lifting it to its concrete `Mechanism` or
# `AllostericMechanism` (`Mechanism(em)`, `AllostericMechanism(aem)`). Only
# `metabolites` is read on the compiled type itself.

"""
    metabolites(m::EnzymeMechanism) → Tuple{Symbol,...}

Return distinct metabolite names as a tuple of `Symbol`s: substrates, then products,
then the names the steps bind as modifiers, each sorted by name, then the ligands of an
allosteric mechanism's regulatory sites, site by site, without duplicates. The fitter
uses this as the key set for the per-datapoint concentration `NamedTuple` passed to
[`rate_equation`](@ref EnzymeRates.Derivation.rate_equation).

```jldoctest
julia> using EnzymeRates

julia> m = @enzyme_mechanism begin
           substrates: S
           products: P
           steps: begin
               E + S <--> E(S)
               E(S) <--> E + P
           end
       end;

julia> metabolites(m)
(:S, :P)
```
"""
@generated function metabolites(::M) where {M <: AbstractEnzymeMechanism}
    # Kept `@generated`: `loss!` uses `metabolites(m)` as a compile-time-constant
    # tuple to build the per-datapoint `NamedTuple{MetNames}` concs on the fitting
    # hot path. A runtime body would make that NamedTuple type-unstable and allocate.
    Tuple(metabolites(_concrete(M())))
end
