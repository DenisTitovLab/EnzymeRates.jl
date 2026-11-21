# ABOUTME: Value types: EnzymeReaction, and the mechanism as its catalytic scheme in
# ABOUTME: canonical form and its allostery, with their checks and identity.

# ─── Concrete types ───────────────────────────────────────────────────

"""The allosteric states a kinetic group or a regulatory site's ligand takes; the
constructors accept no other (`_check_tags`, `RegulatorySite`)."""
const _ALLOSTERIC_STATES = (:OnlyA, :OnlyI, :EqualAI, :NonequalAI)

"""
RegulatorySite: a binding site (possibly multimeric) for one or
more allosteric ligands, each given by its name, once. `ligands[i]` and
`allosteric_states[i]` are parallel. The constructor sorts the ligands by name and
carries each ligand's state along, so the order the ligands are given in never changes
`==`, `hash` or the compiled rate law.
"""
struct RegulatorySite
    ligands::Vector{Symbol}
    multiplicity::Int
    allosteric_states::Vector{Symbol}
    function RegulatorySite(ligands::Vector{Symbol},
                            multiplicity::Int, allosteric_states::Vector{Symbol})
        isempty(ligands) && error("RegulatorySite: needs at least one ligand")
        length(ligands) == length(allosteric_states) ||
            error("RegulatorySite: length(ligands)=$(length(ligands)) " *
                  "must equal length(allosteric_states)=$(length(allosteric_states))")
        multiplicity ≥ 1 ||
            error("RegulatorySite: multiplicity must be ≥ 1, got $multiplicity")
        bad = findfirst(∉(_ALLOSTERIC_STATES), allosteric_states)
        bad === nothing ||
            error("RegulatorySite: allosteric state $(allosteric_states[bad]) " *
                  "must be one of $_ALLOSTERIC_STATES")
        dup = findfirst(i -> ligands[i] in view(ligands, 1:(i - 1)), eachindex(ligands))
        dup === nothing || error("RegulatorySite: ligand $(ligands[dup]) appears twice " *
                                 "in one regulatory site")
        p = sortperm(ligands)
        new(ligands[p], multiplicity, allosteric_states[p])
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
    EnzymeReaction

The public concrete reaction descriptor: substrates, products, optional
regulators, and the set of oligomeric states the mechanism
enumerator is allowed to consider. It is the entry point to the package —
pair it with rate data in an
[`IdentifyRateEquationProblem`](@ref EnzymeRates.IdentifyRateEquationProblem), or pass
it to `EnzymeRates.Enumeration.candidate_mechanisms` to enumerate candidate mechanisms.

# Fields
- `substrates::Vector{Reactant}` and `products::Vector{Reactant}` — each reactant
  carrying its per-atom inventory (used for ping-pong residual bookkeeping
  and atom conservation across steps).
- `regulators::Vector{Regulator}` — competitive inhibitors and
  allosteric regulators, each with its `kind` and `effect`.
- `allowed_oligomeric_states::Vector{Int}` — the oligomeric states, one
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
reason it rejects a reactant or regulator named `Keq`, `L`, `N`, `params`,
`concs` or `v`, names a rate equation keeps for its constants and the
local names of its code, a name holding `_`, since constant names join
names with `_`, the name `E`, which labels the free enzyme form, and a name that is
no identifier, is a keyword or reads in Julia code as another name (a micro sign as `μ`),
which no DSL text can write. A reaction prints as the `@enzyme_reaction` block that
builds it.

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
        !isempty(ms) && all(≥(1), ms) && allunique(ms) ||
            error("EnzymeReaction: allowed_oligomeric_states must be a nonempty " *
                  "list of distinct Ints ≥ 1, got $ms")
        sorted_mults = sort(allowed_oligomeric_states)

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

# ─── The reaction graph ───────────────────────────────────────────────

"""
A metabolite as a form holds it or a step exchanges it: its name, and whether it binds
as a modifier (a name under `modifiers:` written bare, or a reactant's copy `A::Mod`). A
reactant and its modifier copy are two different `Metabolite`s.
"""
struct Metabolite
    name::Symbol
    modifier::Bool
end
Base.isless(x::Metabolite, y::Metabolite) =
    isless((x.name, x.modifier), (y.name, y.modifier))
Base.hash(x::Metabolite, h::UInt) = hash(x.modifier, hash(x.name, hash(:Metabolite, h)))

"""
    EnzymeForm

An enzyme form: its conformation label, the metabolites it binds, and its residual: the
substrates it adds and the products it subtracts. Each list is sorted and holds a name as
often as the form does.
"""
struct EnzymeForm
    conformation::Symbol
    bound::Vector{Metabolite}
    added::Vector{Symbol}
    subtracted::Vector{Symbol}
end

"""
    Step

A step between forms `from` and `to` of its catalytic scheme: going `from → to` it takes
up the metabolites `taken_up` and gives off `given_off`, each list sorted.
`rapid_equilibrium` marks a rapid-equilibrium step, and `group` is the step's kinetic
group, whose steps share their constants.
"""
struct Step
    from::Int
    to::Int
    taken_up::Vector{Metabolite}
    given_off::Vector{Metabolite}
    rapid_equilibrium::Bool
    group::Int
end

"""
    CatalyticScheme(substrates, products, nodes, steps)

The catalytic scheme of a mechanism: its reactants, its enzyme forms and the steps that
join them. `substrates` and `products` name the reactants; `nodes` holds the enzyme
forms; each step's `from` and `to` index `nodes`, and its `group` labels its kinetic group,
the labels running 1, 2, … with each used. The nodes may come in any order and repeat,
and a step may be written in either direction.

The constructor checks the scheme and stores it in canonical order: the substrates and
products sorted, free E (`E`, nothing bound, no residual) first and the other forms in
content order, and each step in the direction the reaction gives it, in content order.
Each step keeps the kinetic group label it is written with; [`Mechanism`](@ref) reads the
labels to give each group its allosteric state, then numbers the groups by their first
step. It throws an `ErrorException` when a step's ends do not index `nodes`, a name is
declared twice, used undeclared or cannot name a metabolite, a form is no free E's, a
form's name is shared, the steps leave a form unjoined to free E or hold one reaction
twice, a kinetic group mixes kinds of step, or a declared name takes part in no step.

Besides its fields `substrates`, `products`, `nodes` and `steps`, it stores what they give,
for reading: `enzyme_forms`, each node's name; `modifiers`, the names bound as modifiers,
sorted; `groups`, each kinetic group's steps; and per step `step_kinds`, `toward_products`
and `step_sides`. Two schemes are `==` when their four fields are.
"""
struct CatalyticScheme
    substrates::Vector{Symbol}
    products::Vector{Symbol}
    nodes::Vector{EnzymeForm}
    steps::Vector{Step}
    enzyme_forms::Vector{Symbol}
    modifiers::Vector{Symbol}
    groups::Vector{Vector{Int}}
    step_kinds::Vector{Symbol}
    toward_products::BitVector
    step_sides::Vector{Tuple{String, String}}
    # The scheme with these four fields and the facts they give, unchecked; `forms` holds
    # each node's name, passed in by a caller that has computed them already. `graph` is
    # partially initialised, so the fact helpers read only its first four fields.
    global function _unchecked_scheme(substrates::Vector{Symbol},
                                      products::Vector{Symbol},
                                      nodes::Vector{EnzymeForm},
                                      steps::Vector{Step},
                                      forms::Vector{String} =
                                          [_form_name(n, substrates, products)
                                           for n in nodes])
        graph = new(substrates, products, nodes, steps)
        new(substrates, products, nodes, steps, Symbol.(forms),
            copy(_modifiers(nodes, steps)), _groups(steps),
            Symbol[_step_kind(graph, e) for e in steps],
            _runs_toward_products.(Ref(graph), steps), _step_sides(graph, forms))
    end
    # `c` with each step's kinetic group label `l` changed to `label[l]`, a permutation of
    # the labels. No fact `c` stores but its groups depends on the labels.
    global function _relabelled(c::CatalyticScheme, label::Vector{Int})
        groups = similar(c.groups)
        for (l, members) in enumerate(c.groups)
            groups[label[l]] = members
        end
        new(c.substrates, c.products, c.nodes,
            [Step(e.from, e.to, e.taken_up, e.given_off, e.rapid_equilibrium,
                  label[e.group]) for e in c.steps],
            c.enzyme_forms, c.modifiers, groups, c.step_kinds, c.toward_products,
            c.step_sides)
    end
end

"""
    Allostery(allosteric_states, oligomeric_state, regulatory_sites)

What a mechanism adds to its catalytic scheme as an MWC enzyme, whose subunits switch
together between an active and an inactive conformation: `allosteric_states[g]`, the
allosteric state of kinetic group `g` (`:OnlyA`, `:OnlyI`, `:EqualAI` or `:NonequalAI`);
the `oligomeric_state`, the number of subunits that switch together; and the
`regulatory_sites`. A plain mechanism has no allosteric state, oligomeric state 1 and no
regulatory site: `Allostery(Symbol[], 1, RegulatorySite[])`. [`Mechanism`](@ref) checks
it against its catalytic scheme.
"""
struct Allostery
    allosteric_states::Vector{Symbol}
    oligomeric_state::Int
    regulatory_sites::Vector{RegulatorySite}
    # Typed arguments, so that no `Bool` converts to an oligomeric state.
    Allostery(allosteric_states::Vector{Symbol}, oligomeric_state::Int,
              regulatory_sites::Vector{RegulatorySite}) =
        new(allosteric_states, oligomeric_state, regulatory_sites)
end

"""
    Mechanism(catalytic::CatalyticScheme, allosteric::Allostery)

An enzyme mechanism: its catalytic scheme (`catalytic`) and its allostery (`allosteric`).
It is plain when `allosteric.allosteric_states` is empty, and allosteric (MWC) otherwise.
`@enzyme_mechanism` builds one from its text, and `show` prints that text.

`allosteric.allosteric_states[l]` is the allosteric state of the steps of `catalytic`
whose kinetic group label is `l`. The constructor stores the mechanism in canonical form,
so two writings of one mechanism are `==` and hash alike: the kinetic groups numbered by
their first step, the states moved with them, and the regulatory sites sorted. A plain
mechanism and an allosteric one with the same catalytic scheme differ. It throws an
`ErrorException` when an allosteric mechanism does not give each kinetic group one known
state or has an oligomeric state below 1, when a plain one has an oligomeric state other
than 1 or a regulatory site, or when a ligand cannot name a metabolite, is a conformation
label or sits at two sites. Whether its rate law is defined is
`check_valid_rate_equation`'s question, not the constructor's.
"""
struct Mechanism
    catalytic::CatalyticScheme
    allosteric::Allostery
    function Mechanism(catalytic::CatalyticScheme, allosteric::Allostery)
        states = allosteric.allosteric_states
        if isempty(states)
            n = length(allosteric.regulatory_sites)
            allosteric.oligomeric_state == 1 && n == 0 ||
                error("Mechanism: a mechanism with no allosteric states has oligomeric " *
                      "state 1 and no regulatory site; got oligomeric state " *
                      "$(allosteric.oligomeric_state) and $n regulatory " *
                      (n == 1 ? "site" : "sites"))
        else
            _check_tags(states, length(catalytic.groups), allosteric.oligomeric_state)
        end
        sites = sort(allosteric.regulatory_sites; by = _site_key)
        _check_ligands(sites, catalytic)
        _labels_in_order(catalytic.steps) &&
            return new(catalytic,
                       Allostery(copy(states), allosteric.oligomeric_state, sites))
        _, group_of = _number_groups(Int[e.group for e in catalytic.steps],
                                     length(catalytic.groups))
        new(_relabelled(catalytic, group_of),
            Allostery(isempty(states) ? Symbol[] : _tags_by_group(group_of, states),
                      allosteric.oligomeric_state, sites))
    end
end

"""Whether the kinetic group labels of `steps` are numbered by their first step: the first
step of each group comes after the first step of each group with a smaller label."""
function _labels_in_order(steps::Vector{Step})
    seen = 0
    for e in steps
        e.group > seen + 1 && return false
        seen = max(seen, e.group)
    end
    true
end

"""A site's place in the order of a mechanism's sites: its ligands, its multiplicity,
then its ligands' tags."""
_site_key(s::RegulatorySite) = (s.ligands, s.multiplicity, s.allosteric_states)

# ─── Content order ──────────────────────────────────────────────────────

"""Nodes by label, then bound metabolites, residual additions and residual subtractions,
each list in lexicographic order, a list coming before every longer list it begins."""
function Base.isless(f::EnzymeForm, g::EnzymeForm)
    f.conformation === g.conformation || return isless(f.conformation, g.conformation)
    c = cmp(f.bound, g.bound)
    c == 0 && (c = cmp(f.added, g.added))
    c == 0 && (c = cmp(f.subtracted, g.subtracted))
    c < 0
end

"""Steps by node indices, then taken-up and given-off metabolites as lists (as
`isless(::EnzymeForm)` orders them), then the flag, steady state first; the group takes no
part. Tuples and `argmin` over steps compare by it."""
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
        -> (; scheme, edge_of, flipped, group_of, node_of)

The canonical `CatalyticScheme` that `nodes` and `steps` write, its kinetic groups
numbered by their first step, with no representation check:
the nodes in any order, a node any number of times, each step between two of them in
either direction, its `group` a written kinetic group label (1, 2, …, each used), and
every list of a node or a step sorted. Each step is stored in the direction
`_orientation` gives. `edge_of[j]` is the step written step `j` became, `flipped[j]`
whether that step runs against step `j` as written, `group_of[l]` the kinetic group
written label `l` became, and `node_of[i]` the node written node `i` became. The checked
constructors build on it. The nodes and steps share the lists of `nodes` and `steps`.
"""
function _canonicalize(substrates::Vector{Symbol}, products::Vector{Symbol},
                       nodes::Vector{EnzymeForm}, steps::Vector{Step})
    distinct, node_of = _sorted_nodes(nodes)
    written = [Step(node_of[e.from], node_of[e.to], e.taken_up, e.given_off,
                    e.rapid_equilibrium, e.group) for e in steps]
    substrates, products = sort(substrates), sort(products)
    forms = [_form_name(n, substrates, products) for n in distinct]
    fwd = _orientation(written, forms, substrates, products)
    stored = [fwd[j] ? written[j] : _reversed(written[j]) for j in eachindex(written)]
    p = sortperm(stored)
    group, group_of = _number_groups([stored[i].group for i in p],
                                     maximum(e -> e.group, steps; init = 0))
    canonical = [Step(e.from, e.to, e.taken_up, e.given_off, e.rapid_equilibrium, group[k])
                 for (k, e) in enumerate(view(stored, p))]
    scheme = _unchecked_scheme(substrates, products, distinct, canonical, forms)
    (; scheme, edge_of = invperm(p), flipped = .!fwd, group_of, node_of)
end

"""The distinct nodes of `written` in canonical order, free E (label `E`, nothing bound,
no residual) first, and the index each written node became. A search keeps every
mechanism it emits, so the nodes, built in a vector grown past its length, are stored as
a copy at its length."""
function _sorted_nodes(written::Vector{EnzymeForm})
    nodes, node_of = EnzymeForm[], Vector{Int}(undef, length(written))
    for t in sortperm(written; by = n -> (!_is_free_E(n), n))
        (isempty(nodes) || nodes[end] != written[t]) && push!(nodes, written[t])
        node_of[t] = length(nodes)
    end
    copy(nodes), node_of
end

"""The side order `_runs_toward_products` reads: substrates, bound plus free, then fewer
products. `free` lists the metabolites a step takes up or gives off at node `n`; a
modifier counts on neither side, and every other name that is no substrate counts as a
product."""
function _side_score(n::EnzymeForm, free::Vector{Metabolite}, substrates::Vector{Symbol})
    s = p = 0
    for list in (n.bound, free), x in list
        x.modifier && continue
        x.name in substrates ? (s += 1) : (p += 1)
    end
    (s, -p)
end

"""
    _orientation(written, forms, substrates, products) -> BitVector

Whether to store each step of `written` as written, from the steps alone: their ends,
which index `forms`, the forms' names with free E first, and the names they take up and
give off. The first rule that decides gives a step's direction:

1. A step that only takes up, or only gives off, metabolites runs as the uptake.
2. A step that takes up some and gives off others runs the way that takes up more
   substrates and gives off more products, net.
3. Any other step, an isomerization or one that trades one modifier for another, runs
   the way the shortest cycle through it that turns the reaction over runs
   (`_turnover_directions!`).
4. With no such cycle, or one each way of equal length, it runs away from free E, from
   the form fewer steps from it.
5. With both ends equally far, it runs from the form whose name comes first, free E
   before every other.

A metabolite takes its role from its name, so a reactant's modifier copy counts as the
reactant. The steps of a kinetic group share their constants, so they turn together,
after the member nearest free E among those rules 1 to 3 decide, or among all of them
when they decide none (`_turn_groups_together!`).
"""
function _orientation(written::Vector{Step}, forms::Vector{String},
                      substrates::Vector{Symbol}, products::Vector{Symbol})
    nf = length(forms)
    net(list::Vector{Metabolite}) =
        count(x -> x.name in substrates, list) - count(x -> x.name in products, list)
    # Rules 1 and 2: positive to store a step as written, negative to turn it, zero when
    # neither decides, as for an isomerization.
    verdict = [isempty(e.taken_up) || isempty(e.given_off) ?
               length(e.taken_up) - length(e.given_off) :
               net(e.taken_up) - net(e.given_off) for e in written]
    at = _steps_at_forms(written, nf)
    _turnover_directions!(verdict, written, at, substrates[1])
    distance = _distance_from_free_E(written, at)
    # Each form's place: free E first, then the others by name. Rule 5 and the group rule
    # break their ties by it, so a direction follows the forms' names and not the order
    # the forms are stored in.
    place = nf == 0 ? Int[] : invperm([1; sortperm(view(forms, 2:nf)) .+ 1])
    fwd, decided = trues(length(written)), trues(length(written))
    for (k, e) in enumerate(written)
        decided[k] = verdict[k] != 0
        a, b = distance[e.from], distance[e.to]
        fwd[k] = decided[k] ? verdict[k] > 0 : a == b ? place[e.from] < place[e.to] : a < b
    end
    _turn_groups_together!(fwd, decided, written, distance, place)
    fwd
end

"""The steps of `written` at each of its `nf` forms: those that touch form `u` are
`at.steps[at.start[u]:(at.start[u + 1] - 1)]`, a step joining a form to itself twice."""
function _steps_at_forms(written::Vector{Step}, nf::Int)
    start = zeros(Int, nf + 1)
    for e in written
        start[e.from + 1] += 1
        start[e.to + 1] += 1
    end
    start[1] = 1
    for u in 1:nf
        start[u + 1] += start[u]
    end
    steps, filled = Vector{Int}(undef, 2 * length(written)), start[1:nf]
    for (k, e) in enumerate(written), u in (e.from, e.to)
        steps[filled[u]] = k
        filled[u] += 1
    end
    (; start, steps)
end

"""The form at the other end of step `e` from form `u`."""
_other_end(e::Step, u::Int) = e.from == u ? e.to : e.from

"""How many steps of `written` separate each of its forms from free E, form 1, along
steps taken either way (`at`, from `_steps_at_forms`); `typemax(Int)` for a form no
steps join to it."""
function _distance_from_free_E(written::Vector{Step}, at)
    distance = fill(typemax(Int), length(at.start) - 1)
    isempty(distance) && return distance
    distance[1] = 0
    queue, head = [1], 0
    while head < length(queue)
        u = queue[head += 1]
        for i in at.start[u]:(at.start[u + 1] - 1)
            v = _other_end(written[at.steps[i]], u)
            distance[v] == typemax(Int) || continue
            distance[v] = distance[u] + 1
            push!(queue, v)
        end
    end
    distance
end

"""
    _turnover_directions!(verdict, written, at, s1) -> verdict

Decide each step `k` of `written` (`at` lists its steps by form, `_steps_at_forms`) that
`verdict` leaves at `0`: `verdict[k]` becomes `1` when the shortest closed walk that takes
step `k` once, as written, and takes up one of the first substrate `s1` net (a modifier
copy counting as it) is shorter than any that takes the step the other way, `-1` when it
is longer, and stays `0` when there is none either way or both are equally long. The walk
may pass a form twice, so a step on a cycle that turns nothing over finds a walk each way
of one length, round a turnover cycle in either direction, and is decided by neither.
"""
function _turnover_directions!(verdict::Vector{Int}, written::Vector{Step}, at, s1::Symbol)
    any(iszero, verdict) || return verdict
    nf = length(at.start) - 1
    uptake = [count(x -> x.name === s1, e.taken_up) - count(x -> x.name === s1, e.given_off)
              for e in written]
    reached, forms = falses(nf), Int[]
    # The search runs breadth first over (form, net uptake of s1) and keeps the net uptake
    # within ±span, as many as there are forms, so that it ends: a walk round a turnover
    # cycle takes up more on every turn. Each queued state carries its walk's length.
    span = nf
    seen = falses(nf, 2span + 1)
    queue = NTuple{3, Int}[]
    for k in eachindex(written)
        verdict[k] == 0 || continue
        from, to = written[k].from, written[k].to
        # A step whose ends no other steps join lies on no cycle. Looking first keeps the
        # search below from visiting every (form, net uptake) pair for such a step.
        fill!(reached, false)
        empty!(forms)
        reached[to] = true
        push!(forms, to)
        head = 0
        while head < length(forms) && !reached[from]
            u = forms[head += 1]
            for i in at.start[u]:(at.start[u + 1] - 1)
                at.steps[i] == k && continue
                v = _other_end(written[at.steps[i]], u)
                reached[v] || (reached[v] = true; push!(forms, v))
            end
        end
        reached[from] || continue
        # The walk back from `to` to `from`, without step `k`, takes up `1 - uptake[k]`
        # net when step `k` as written turns the reaction over, and `-1 - uptake[k]` when
        # step `k` taken the other way does.
        ahead, behind = 1 - uptake[k], -1 - uptake[k]
        abs(ahead) <= span && abs(behind) <= span || continue
        fill!(seen, false)
        empty!(queue)
        seen[to, span + 1] = true
        push!(queue, (to, 0, 0))
        a = b = typemax(Int)
        head = 0
        # Once the search reaches a state as far out as the shortest walk found, every
        # walk of that length has been found.
        while head < length(queue) && queue[head + 1][3] < min(a, b)
            u, n, d = queue[head += 1]
            for i in at.start[u]:(at.start[u + 1] - 1)
                j = at.steps[i]
                j == k && continue
                v = _other_end(written[j], u)
                x = written[j].from == u ? uptake[j] : -uptake[j]
                (abs(n + x) > span || seen[v, n + x + span + 1]) && continue
                seen[v, n + x + span + 1] = true
                v == from && n + x == ahead && (a = d + 1)
                v == from && n + x == behind && (b = d + 1)
                push!(queue, (v, n + x, d + 1))
            end
        end
        verdict[k] = a == b ? 0 : a < b ? 1 : -1
    end
    verdict
end

"""
Turn the steps of each kinetic group (each written step's `group` label) with the
group's lead, its member nearest free E among those `decided` marks, or among all of
them when none is marked: nearest by the distances of its two ends from free E
(`distance`), the nearer end first, then by the places of its two forms (`place`), the
earlier first. A member that takes up or gives off metabolites takes the direction in
which it takes up what the lead takes up; an isomerization keeps the direction it was
written in relative to the lead. `fwd[k]` says whether written step `k` is stored as
written.
"""
function _turn_groups_together!(fwd::BitVector, decided::BitVector, written::Vector{Step},
                                distance::Vector{Int}, place::Vector{Int})
    near(k::Int) = (e = written[k];
                    (minmax(distance[e.from], distance[e.to])...,
                     minmax(place[e.from], place[e.to])...))
    for g in 1:maximum(e -> e.group, written; init = 0)
        in_group(k::Int) = written[k].group == g
        count(in_group, eachindex(written)) < 2 && continue
        any_decided = any(k -> in_group(k) && decided[k], eachindex(written))
        lead = 0
        for k in eachindex(written)
            in_group(k) && (decided[k] || !any_decided) || continue
            (lead == 0 || near(k) < near(lead)) && (lead = k)
        end
        taken_up = fwd[lead] ? written[lead].taken_up : written[lead].given_off
        for k in eachindex(written)
            in_group(k) || continue
            e = written[k]
            if _is_iso(e)
                fwd[k] = fwd[lead]
            elseif e.taken_up == taken_up || e.given_off == taken_up
                fwd[k] = e.taken_up == taken_up
            end
        end
    end
    fwd
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

function CatalyticScheme(substrates::Vector{Symbol}, products::Vector{Symbol},
                         nodes::Vector{EnzymeForm}, steps::Vector{Step})
    for (k, e) in enumerate(steps)
        e.from in eachindex(nodes) && e.to in eachindex(nodes) && continue
        error("Step $k: its ends are nodes $(e.from) and $(e.to), but there are " *
              "$(length(nodes)) nodes; a step's ends index the nodes")
    end
    nodes = [EnzymeForm(n.conformation, sort(n.bound), sort(n.added), sort(n.subtracted))
             for n in nodes]
    steps = [Step(e.from, e.to, sort(e.taken_up), sort(e.given_off), e.rapid_equilibrium,
                  e.group) for e in steps]
    # Every name a node or a step binds as a modifier is a modifier.
    modifiers = _modifiers(nodes, steps)
    _check_declarations(substrates, products, modifiers, Int[e.group for e in steps])
    c = _canonicalize(substrates, products, nodes, steps)
    # The checks name each kinetic group by the number a mechanism gives it; the scheme
    # keeps the written labels, which `Mechanism` reads to give each group its state.
    _check_representation(c.scheme, modifiers)
    c.group_of == eachindex(c.group_of) ? c.scheme :
        _relabelled(c.scheme, invperm(c.group_of))
end

"""The tags by canonical group: `allosteric_states[l]`, the tag of written label `l`,
moved to group `group_of[l]`."""
function _tags_by_group(group_of::Vector{Int}, allosteric_states::Vector{Symbol})
    tags = similar(allosteric_states)
    for (l, g) in enumerate(group_of)
        tags[g] = allosteric_states[l]
    end
    tags
end

# ─── Checks ─────────────────────────────────────────────────────────────

"""Error unless `allosteric_states` holds one valid tag (`_ALLOSTERIC_STATES`) for each
of `G` kinetic groups and the oligomeric state is at least 1."""
function _check_tags(allosteric_states::Vector{Symbol}, G::Int, oligomeric_state::Int)
    length(allosteric_states) == G ||
        error("Mechanism: allosteric_states length " *
              "$(length(allosteric_states)) must match the number of kinetic groups $G")
    oligomeric_state ≥ 1 ||
        error("Mechanism: oligomeric_state must be ≥ 1, " *
              "got $oligomeric_state")
    bad = findfirst(∉(_ALLOSTERIC_STATES), allosteric_states)
    bad === nothing ||
        error("Mechanism: catalytic group $bad has invalid allosteric state " *
              "$(allosteric_states[bad]) (must be one of $_ALLOSTERIC_STATES)")
    nothing
end

"""Error when a site ligand cannot name a metabolite (`_check_metabolite_names`) or is a
conformation label of `catalytic`, or when a ligand sits on two sites: each ligand's
constant is named after the ligand alone."""
function _check_ligands(sites::Vector{RegulatorySite}, catalytic::CatalyticScheme)
    ligs = Symbol[l for s in sites for l in s.ligands]
    _check_metabolite_names("Mechanism", ligs)
    for l in ligs
        any(n -> n.conformation === l, catalytic.nodes) &&
            _label_collision_error(l, "Mechanism", "ligand")
    end
    dup = findfirst(i -> ligs[i] in view(ligs, 1:(i - 1)), eachindex(ligs))
    dup === nothing ||
        error("Mechanism: ligand $(ligs[dup]) appears in two " *
              "distinct regulatory sites; rendered regulator parameter names would " *
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
Error unless the canonical scheme `m` represents a mechanism: conformation labels shaped as
labels (`check_conformation_label`), none a metabolite's name; every name declared (in
`modifiers`, the names bound as modifiers, or among its substrates and products), a plain
entry only for a reactant, residuals that add substrates and subtract products; no step
joining a node to itself, taking up and giving off one metabolite, or binding a
metabolite into a node that holds it in another role; no two nodes that render one name;
each kinetic group one kind of step with one flag; each reaction once; free E; one
connected graph; and a step that exchanges each declared name and each name of
`modifiers`, a modifier as a modifier.
"""
function _check_representation(m::CatalyticScheme, modifiers::Vector{Symbol})
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
    _check_distinct_form_names(m)
    _check_uniform_groups(m, m.groups)
    _check_each_reaction_once(m, m.groups)
    _check_free_E(m)
    _check_connected(m)
    _check_declared_used(m, modifiers)
end

"""Error unless node 1 of `m` is free E: label `E`, nothing bound, no residual. Canonical
form puts free E first wherever a mechanism holds it."""
function _check_free_E(m::CatalyticScheme)
    !isempty(m.nodes) && _is_free_E(m.nodes[1]) ||
        error("Mechanism: no free enzyme form E; every mechanism holds the free form E, " *
              "with nothing bound and no residual")
end

"""Error unless the steps of `m` join every node to free E, node 1: nodes that never
meet free E split the enzyme between parts with no single steady state."""
function _check_connected(m::CatalyticScheme)
    f = findfirst(!, _joined_to_free_E(m, falses(length(m.steps))))
    f === nothing ||
        error("Mechanism: the steps do not form one connected graph; no steps lead from " *
              "free E to $(_form_name(m, f))")
end

"""Error unless some step of `m` takes up or gives off each declared substrate and
product, and some step takes up or gives off each modifier `modifiers` declares as a
modifier: a declared name that no step exchanges takes no part in the rate law."""
function _check_declared_used(m::CatalyticScheme, modifiers::Vector{Symbol})
    used(x::Metabolite) = any(e -> x in e.taken_up || x in e.given_off, m.steps)
    for (names, kind) in ((m.substrates, "substrate"), (m.products, "product"))
        k = findfirst(n -> !used(Metabolite(n, false)), names)
        k === nothing ||
            error("Mechanism: $kind $(names[k]) is declared, but no step uses it")
    end
    k = findfirst(n -> !used(Metabolite(n, true)), modifiers)
    k === nothing ||
        error("Mechanism: modifier $(modifiers[k]) is bound at a node, " *
              "but no step binds it")
end

"""Error when two nodes of `m` render one name: every constant is named after the nodes
its step joins (`_form_name`)."""
function _check_distinct_form_names(m::CatalyticScheme)
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
function _check_uniform_groups(m::CatalyticScheme, groups::Vector{Vector{Int}})
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
function _check_each_reaction_once(m::CatalyticScheme, groups::Vector{Vector{Int}})
    seen = Dict{NTuple{2, Tuple{Int, Vector{Metabolite}}}, Tuple{Int, Int}}()
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
# seeded with the type's name. A catalytic scheme's stored facts take no part: they are
# functions of its first four fields. A struct's default hash mixes in its type's
# identity, which changes between precompile sessions; a value's hash reads only
# symbols, integers and flags, so every process agrees on it. `Metabolite` keeps Julia's
# `==`, which compares its two fields, and has its own `hash`.
for T in (RegulatorySite, Reactant, Regulator, EnzymeReaction, EnzymeForm, Step,
          CatalyticScheme, Allostery, Mechanism)
    fs = T === CatalyticScheme ? (:substrates, :products, :nodes, :steps) : fieldnames(T)
    @eval Base.:(==)(a::$T, b::$T) =
        $(reduce((l, r) -> :($l && $r), [:(a.$f == b.$f) for f in fs]))
    @eval Base.hash(x::$T, h::UInt) =
        $(foldl((acc, f) -> :(hash(x.$f, $acc)), fs;
                init = :(hash($(QuoteNode(nameof(T))), h))))
end

# ─── Reads ──────────────────────────────────────────────────────────────

"""
Every metabolite entry `m` holds, sorted by (name, modifier) and once each: the
metabolites its nodes bind and its steps take up and give off, and its residual names as
non-modifier entries. A reactant and its modifier copy are two entries.
"""
function _table(m::CatalyticScheme)
    table = Metabolite[]
    for n in m.nodes
        append!(table, n.bound)
        foreach(x -> push!(table, Metabolite(x, false)), n.added)
        foreach(x -> push!(table, Metabolite(x, false)), n.subtracted)
    end
    for e in m.steps
        append!(table, e.taken_up)
        append!(table, e.given_off)
    end
    unique!(sort!(table))
end
"""Every name `nodes` and `steps` bind as a modifier, sorted, each once."""
function _modifiers(nodes::Vector{EnzymeForm}, steps::Vector{Step})
    names = Symbol[]
    add(list::Vector{Metabolite}) =
        foreach(x -> x.modifier && !(x.name in names) && push!(names, x.name), list)
    foreach(n -> add(n.bound), nodes)
    foreach(e -> (add(e.taken_up); add(e.given_off)), steps)
    sort!(names)
end
"""
    metabolites(m) -> Vector{Symbol}

The names of the metabolites of mechanism `m`, each once: its substrates, then its
products, then the names its steps bind as modifiers only, each sorted by name, then the
ligands of an allosteric mechanism's regulatory sites, site by site. These are the
concentrations its rate law reads, in the order `metabolites(rate)` of its compiled law
([`compile_rate_equation`](@ref EnzymeRates.Derivation.compile_rate_equation)) gives
them as a tuple.

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
2-element Vector{Symbol}:
 :S
 :P
```
"""
metabolites(m::Mechanism) =
    unique!([metabolites(m.catalytic);
             [l for s in m.allosteric.regulatory_sites for l in s.ligands]])
"""The names of the metabolites of catalytic scheme `c`, each once: its substrates, then
its products, then the names its steps bind as modifiers only, each sorted by name."""
metabolites(c::CatalyticScheme) = unique!([c.substrates; c.products; c.modifiers])
"""Whether `x` names a substrate or a product of `m`."""
_is_reactant(m::CatalyticScheme, x::Symbol) =
    insorted(x, m.substrates) || insorted(x, m.products)
_has_residual(n::EnzymeForm) = !isempty(n.added) || !isempty(n.subtracted)
"""Whether step `e` takes up one metabolite and gives off none."""
_is_binding(e::Step) = isempty(e.given_off) && length(e.taken_up) == 1
"""Whether step `e` takes up and gives off nothing: an isomerization."""
_is_iso(e::Step) = isempty(e.taken_up) && isempty(e.given_off)
"""Whether binding `e` joins two nodes of `m` that carry one residual."""
function _keeps_residual(m::CatalyticScheme, e::Step)
    f, t = m.nodes[e.from], m.nodes[e.to]
    f.added == t.added && f.subtracted == t.subtracted
end
"""Whether the sorted list `t` is the sorted list `f` with one more `x`."""
function _holds_one_more(t::Vector{Metabolite}, f::Vector{Metabolite}, x::Metabolite)
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
_is_plain_binding(m::CatalyticScheme, e::Step) =
    _is_binding(e) && _keeps_residual(m, e) &&
    _holds_one_more(m.nodes[e.to].bound, m.nodes[e.from].bound, only(e.taken_up))
"""Whether binding `e` takes up its metabolite into a node that holds its source's
metabolites plus that one by name, with the same residual: a plain binding, or one that
binds a reactant where its modifier copy sits, or the reverse."""
function _binds_by_name(m::CatalyticScheme, e::Step)
    _is_binding(e) && _keeps_residual(m, e) || return false
    names(n::Int) = Symbol[x.name for x in m.nodes[n].bound]
    names(e.to) == sort!(push!(names(e.from), only(e.taken_up).name))
end
"""
Whether step `e` of `m` runs, as stored, toward products: from the side that holds more
substrates, bound or free, then fewer products (`_side_score`). A step that takes up or
gives off one metabolite is stored as its uptake, so a half-reaction that releases a
product, `E(A) <--> E(; residual = A - P) + P`, runs toward products against its stored
direction. A step whose two sides tie, such as an isomerization of free forms, runs as
stored.
"""
_runs_toward_products(m::CatalyticScheme, e::Step) =
    _side_score(m.nodes[e.from], e.taken_up, m.substrates) >=
    _side_score(m.nodes[e.to], e.given_off, m.substrates)
"""Whether `n` is a free node: nothing bound and no residual, whatever its label."""
_is_free(n::EnzymeForm) = isempty(n.bound) && !_has_residual(n)
"""Whether `n` is free E: label `E`, nothing bound and no residual."""
_is_free_E(n::EnzymeForm) = n.conformation === :E && _is_free(n)
