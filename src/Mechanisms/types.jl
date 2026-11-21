# ABOUTME: Value types: EnzymeReaction, and the mechanism as a reaction graph (its catalytic
# ABOUTME: scheme in canonical form) and its allostery, with their checks and identity.

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
    Step(from, to, taken_up, given_off, rapid_equilibrium, kinetic_group)

A step between the enzyme forms `from` and `to` of its catalytic scheme, indices into its
`enzyme_forms`: going `from → to` it takes up the metabolites `taken_up` and gives off
`given_off`, indices into its `metabolites`, each list sorted and holding a metabolite as
often as the step moves it. `rapid_equilibrium` marks a rapid-equilibrium step, and
`kinetic_group` is the step's kinetic group, whose steps share their constants.
"""
struct Step
    from::Int
    to::Int
    taken_up::Vector{Int}
    given_off::Vector{Int}
    rapid_equilibrium::Bool
    kinetic_group::Int
end

"""
    CatalyticScheme(metabolites, enzyme_forms, substrates, products, steps)

The catalytic scheme of a mechanism: a reaction graph. `metabolites` names the metabolites
its steps take up or give off, and `substrates` and `products` index them; every other
metabolite is a modifier. `enzyme_forms` names its enzyme forms, the free enzyme `E` among
them; a name is opaque, and only `E` means something. Each step of `steps` joins two forms
and takes up and gives off metabolites by index (`Step`); its `kinetic_group` labels its
kinetic group, the labels running 1, 2, … with each used. The lists may come in any order
and a step may be written in either direction.

The constructor checks the representation and stores the scheme in canonical form: the
metabolites as substrates, then products, then modifiers, each by name; `E` first and the
other forms by name; each step in the direction the graph gives it (`_orientation`), the
steps sorted. Each step keeps the kinetic group label it is written with;
[`Mechanism`](@ref) reads the labels to give each group its allosteric state, then numbers
the groups by their first step. It throws an `ArgumentError` when no form is named `E`,
and an `ErrorException` when an index is out of range, a name is listed twice or cannot be
written, a step joins a form to itself or takes up and gives off one metabolite, a form is
not joined to `E`, a kinetic group mixes kinds of step, a reaction appears twice, or a
metabolite takes part in no step.
"""
struct CatalyticScheme
    metabolites::Vector{Symbol}
    enzyme_forms::Vector{Symbol}
    substrates::Vector{Int}
    products::Vector{Int}
    steps::Vector{Step}
    function CatalyticScheme(metabolites::Vector{Symbol}, enzyme_forms::Vector{Symbol},
                             substrates::Vector{Int}, products::Vector{Int},
                             steps::Vector{Step})
        _check_written(metabolites, enzyme_forms, substrates, products, steps)
        c = _canonical(metabolites, enzyme_forms, substrates, products, steps)
        _check_graph(c)
        c
    end
    # The scheme with these fields, unchecked and as given: the canonical scheme is built
    # through it, and a mechanism renumbers its kinetic groups through it.
    global _unchecked_scheme(metabolites::Vector{Symbol}, enzyme_forms::Vector{Symbol},
                             substrates::Vector{Int}, products::Vector{Int},
                             steps::Vector{Step}) =
        new(metabolites, enzyme_forms, substrates, products, steps)
end

# ─── Canonical form ─────────────────────────────────────────────────────

"""
    _canonical(metabolites, enzyme_forms, substrates, products, steps) -> CatalyticScheme

The canonical scheme that the arguments write, unchecked but for what it needs: indices
in range and a form named `E` (`_check_written`). The metabolites are stored as
substrates, then products, then modifiers, each by name; `E` first and the other forms by
name; each step's lists sorted, the step in the direction `_orientation` gives, and the
steps sorted by their ends, then their lists, then their flag. Each step keeps its kinetic
group label.
"""
function _canonical(metabolites::Vector{Symbol}, enzyme_forms::Vector{Symbol},
                    substrates::Vector{Int}, products::Vector{Int}, steps::Vector{Step})
    role(i) = i in substrates ? 1 : i in products ? 2 : 3
    met_order = sortperm(eachindex(metabolites); by = i -> (role(i), metabolites[i]))
    met_of = invperm(met_order)
    form_order = sortperm(eachindex(enzyme_forms);
                          by = f -> (enzyme_forms[f] !== :E, enzyme_forms[f]))
    form_of = invperm(form_order)
    ns, np = length(substrates), length(products)
    written = [Step(form_of[e.from], form_of[e.to], sort!(met_of[e.taken_up]),
                    sort!(met_of[e.given_off]), e.rapid_equilibrium, e.kinetic_group)
               for e in steps]
    fwd = _orientation(written, length(enzyme_forms), ns, np)
    stored = [fwd[k] ? e : _reversed(e) for (k, e) in enumerate(written)]
    sort!(stored; by = e -> (e.from, e.to, e.taken_up, e.given_off, e.rapid_equilibrium))
    _unchecked_scheme(metabolites[met_order], enzyme_forms[form_order], collect(1:ns),
                      collect((ns + 1):(ns + np)), stored)
end

"""`e` stored in the opposite direction: forms swapped, taken-up and given-off swapped."""
_reversed(e::Step) =
    Step(e.to, e.from, e.given_off, e.taken_up, e.rapid_equilibrium, e.kinetic_group)

"""Whether step `e` takes up and gives off nothing: an isomerization."""
_is_iso(e::Step) = isempty(e.taken_up) && isempty(e.given_off)

"""
    _orientation(written, nf, ns, np) -> BitVector

Whether to store each step of `written` as written, from the graph alone: its `nf` forms
stored in canonical order, `E` first and the others by name, and its metabolites the
`ns` substrates, then the `np` products, then the modifiers. The first rule that decides
gives a step's direction:

1. A step that only takes up, or only gives off, metabolites runs as the uptake.
2. A step that takes up some and gives off others runs the way that takes up more
   substrates and gives off more products, net.
3. Any other step, an isomerization or one that trades one modifier for another, runs
   the way the shortest cycle through it that turns the reaction over runs
   (`_turnover_directions!`).
4. With no such cycle, or one each way of equal length, it runs away from `E`, from the
   form fewer steps from it.
5. With both ends equally far, it runs from the form stored first.

The steps of a kinetic group share their constants, so they turn together, after the
member nearest `E` among those rules 1 to 3 decide, or among all of them when they decide
none (`_turn_groups_together!`).
"""
function _orientation(written::Vector{Step}, nf::Int, ns::Int, np::Int)
    net(list::Vector{Int}) = count(<=(ns), list) - count(i -> ns < i <= ns + np, list)
    # Rules 1 and 2: positive to store a step as written, negative to turn it, zero when
    # neither decides, as for an isomerization.
    verdict = [isempty(e.taken_up) || isempty(e.given_off) ?
               length(e.taken_up) - length(e.given_off) :
               net(e.taken_up) - net(e.given_off) for e in written]
    at = _steps_at_forms(written, nf)
    # The first substrate, metabolite 1, counts turnover.
    _turnover_directions!(verdict, written, at, 1)
    distance = _distance_from_free_E(written, at)
    fwd, decided = trues(length(written)), trues(length(written))
    for (k, e) in enumerate(written)
        decided[k] = verdict[k] != 0
        a, b = distance[e.from], distance[e.to]
        fwd[k] = decided[k] ? verdict[k] > 0 : a == b ? e.from < e.to : a < b
    end
    _turn_groups_together!(fwd, decided, written, distance)
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

"""How many steps of `written` separate each of its forms from `E`, form 1, along
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
step `k` once, as written, and takes up one of the first substrate, metabolite `s1`, net
is shorter than any that takes the step the other way, `-1` when it
is longer, and stays `0` when there is none either way or both are equally long. The walk
may pass a form twice, so a step on a cycle that turns nothing over finds a walk each way
of one length, round a turnover cycle in either direction, and is decided by neither.
"""
function _turnover_directions!(verdict::Vector{Int}, written::Vector{Step}, at, s1::Int)
    any(iszero, verdict) || return verdict
    nf = length(at.start) - 1
    uptake = [count(==(s1), e.taken_up) - count(==(s1), e.given_off) for e in written]
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
        # A step whose two ends the other steps do not join lies on no cycle. Looking
        # first keeps the search below from visiting every (form, net uptake) pair for
        # such a step.
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
Turn the steps of each kinetic group (each written step's `kinetic_group` label) with the
group's lead, its member nearest `E` among those `decided` marks, or among all of them
when none is marked: nearest by the distances of its two ends from `E` (`distance`), the
nearer end first, then by its two forms, the one stored first first. A member that takes
up or gives off metabolites takes the direction in which it takes up what the lead takes
up; an isomerization keeps the direction it was written in relative to the lead. `fwd[k]`
says whether written step `k` is stored as written.
"""
function _turn_groups_together!(fwd::BitVector, decided::BitVector, written::Vector{Step},
                                distance::Vector{Int})
    near(k::Int) = (e = written[k];
                    (minmax(distance[e.from], distance[e.to])..., minmax(e.from, e.to)...))
    for g in 1:maximum(e -> e.kinetic_group, written; init = 0)
        in_group(k::Int) = written[k].kinetic_group == g
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

# ─── Checks ─────────────────────────────────────────────────────────────

"""
Error unless the scheme as written can be stored: every index in range; `substrates` and
`products` not empty; each metabolite name once, a substrate's not a product's and a
modifier's not a reactant's, and a name a metabolite can take (`_check_metabolite_names`);
the form names distinct, names a text can write (`_check_identifiers`) and no metabolite's;
a form named `E` (an `ArgumentError`); and the kinetic group labels running from 1 with
each used.
"""
function _check_written(metabolites::Vector{Symbol}, enzyme_forms::Vector{Symbol},
                        substrates::Vector{Int}, products::Vector{Int},
                        steps::Vector{Step})
    nm, nf = length(metabolites), length(enzyme_forms)
    for (k, e) in enumerate(steps)
        e.from in 1:nf && e.to in 1:nf ||
            error("Mechanism: step $k joins forms $(e.from) and $(e.to), but there are " *
                  "$nf forms; a step's ends index the enzyme forms")
        all(in(1:nm), e.taken_up) && all(in(1:nm), e.given_off) ||
            error("Mechanism: step $k takes up $(e.taken_up) and gives off " *
                  "$(e.given_off), but there are $nm metabolites; a step's lists index " *
                  "the metabolites")
    end
    for (list, kind) in ((substrates, "substrate"), (products, "product"))
        isempty(list) && error("Mechanism: $(kind)s must not be empty")
        all(in(1:nm), list) ||
            error("Mechanism: the $(kind)s $list index $nm metabolites")
        allunique(metabolites[list]) || error("Mechanism: duplicate $kind names")
    end
    both = findfirst(in(metabolites[products]), metabolites[substrates])
    both === nothing ||
        error("Mechanism: $(metabolites[substrates[both]]) named as both a substrate " *
              "and a product; concentrations and constants are keyed by name")
    modifiers = metabolites[[i for i in 1:nm if !(i in substrates || i in products)]]
    allunique(modifiers) || error("Mechanism: duplicate modifier names")
    reactant = findfirst(in(metabolites[[substrates; products]]), modifiers)
    reactant === nothing ||
        error("Mechanism: $(modifiers[reactant]) is declared both as a reactant and as a " *
              "modifier; a step that binds a reactant as a modifier simply takes it up, " *
              "as `ES + S ⇌ ESS` does")
    _check_metabolite_names("Mechanism", metabolites)
    _check_identifiers("Mechanism", "enzyme form", enzyme_forms)
    dup = findfirst(i -> enzyme_forms[i] in view(enzyme_forms, 1:(i - 1)),
                    eachindex(enzyme_forms))
    dup === nothing || error("Mechanism: two enzyme forms are named $(enzyme_forms[dup])")
    shared = findfirst(in(metabolites), enzyme_forms)
    shared === nothing ||
        error("Mechanism: $(enzyme_forms[shared]) names both an enzyme form and a " *
              "metabolite; a text could not tell them apart")
    :E in enzyme_forms ||
        throw(ArgumentError("Mechanism: no enzyme form is named `E`; the free enzyme is " *
                            "the form named `E`"))
    labels = Int[e.kinetic_group for e in steps]
    all(≥(1), labels) ||
        error("Mechanism: each step takes one kinetic group label from 1 to the number " *
              "of groups")
    empty = findfirst(∉(labels), 1:maximum(labels; init = 0))
    empty === nothing ||
        error("Mechanism: kinetic group $empty is empty; a kinetic group holds at least " *
              "one step")
    nothing
end

"""
Error unless the canonical scheme `c` is a graph a mechanism can have: no step joins a
form to itself or takes up and gives off one metabolite, the steps join every form to
`E`, the steps of each kinetic group take up and give off the same metabolites with one
flag, each reaction appears once, and every metabolite takes part in some step.
"""
function _check_graph(c::CatalyticScheme)
    for e in c.steps
        e.from == e.to &&
            error("Step: both ends are the form $(c.enzyme_forms[e.from]); a step joins " *
                  "two forms")
        shared = unique!([i for i in e.taken_up if i in e.given_off])
        isempty(shared) ||
            error("Step $(c.enzyme_forms[e.from]) → $(c.enzyme_forms[e.to]): " *
                  "$(join(c.metabolites[shared], ", ")) both consumed and released")
    end
    cut_off = findall(!, _joined_to_free_E(c, falses(length(c.steps))))
    isempty(cut_off) ||
        error("Mechanism: the steps do not form one connected graph; no steps lead from " *
              "free E to " * join(c.enzyme_forms[cut_off], ", "))
    # The checks name each kinetic group by the number a mechanism gives it, by its first
    # step; the scheme keeps the written labels, which `Mechanism` reads.
    written = Int[e.kinetic_group for e in c.steps]
    labels, _ = _number_groups(written, maximum(written; init = 0))
    groups = [findall(==(g), labels) for g in 1:maximum(labels; init = 0)]
    _check_uniform_groups(c, groups)
    _check_each_reaction_once(c, groups)
    used = falses(length(c.metabolites))
    for e in c.steps, i in Iterators.flatten((e.taken_up, e.given_off))
        used[i] = true
    end
    k = findfirst(!, used)
    k === nothing && return nothing
    k in c.substrates || k in c.products ||
        error("Mechanism: modifier $(c.metabolites[k]) is declared, but no step binds it")
    error("Mechanism: $(k in c.substrates ? "substrate" : "product") " *
          "$(c.metabolites[k]) is declared, but no step uses it")
end

"""Error unless the steps of each kinetic group (`groups`) of `c` take up and give off the
same metabolites and carry one flag: they share their constants."""
function _check_uniform_groups(c::CatalyticScheme, groups::Vector{Vector{Int}})
    label(e::Step) = join(_sides(c, e), " → ") * (e.rapid_equilibrium ? " (RE)" : " (SS)")
    for members in groups
        e1 = c.steps[first(members)]
        for k in members
            e = c.steps[k]
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
function _check_each_reaction_once(c::CatalyticScheme, groups::Vector{Vector{Int}})
    seen = Dict{NTuple{2, Tuple{Int, Vector{Int}}}, Tuple{Int, Int}}()
    flag(e::Step) = e.rapid_equilibrium ? "RE" : "SS"
    for (g, members) in enumerate(groups), k in members
        e = c.steps[k]
        key = minmax((e.from, e.taken_up), (e.to, e.given_off))
        if haskey(seen, key)
            g0, k0 = seen[key]
            e0 = c.steps[k0]
            text = join(_sides(c, e0), " ⇌ ")
            g0 == g && error("Mechanism: kinetic group $g holds the reaction $text " *
                             "twice; a mechanism holds each reaction once")
            flags = flag(e0) == flag(e) ? "" :
                " ($(flag(e0)) in group $g0, $(flag(e)) in group $g)"
            error("Mechanism: kinetic groups $g0 and $g both hold the reaction " *
                  "$text$flags; a reaction belongs to one kinetic group, whose " *
                  "constants are named after it")
        end
        seen[key] = (g, k)
    end
end

"""The two sides of step `e` of `c` in its stored direction, as names: each side its form's
name and the metabolites the step takes up or gives off there, joined by `_` (`EA_B`)."""
_sides(c::CatalyticScheme, e::Step) =
    (join([c.enzyme_forms[e.from]; c.metabolites[e.taken_up]], "_"),
     join([c.enzyme_forms[e.to]; c.metabolites[e.given_off]], "_"))

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

An enzyme mechanism: its catalytic scheme (`catalytic`), a reaction graph, and its
allostery (`allosteric`). It is plain when `allosteric.allosteric_states` is empty, and
allosteric (MWC) otherwise. `@enzyme_mechanism` builds one from its text, and `show`
prints that text.

`allosteric.allosteric_states[l]` is the allosteric state of the steps of `catalytic`
whose kinetic group label is `l`. The constructor stores the mechanism in canonical form,
so two writings of one mechanism are `==` and hash alike: the kinetic groups numbered by
their first step, the states moved with them, and the regulatory sites sorted. A plain
mechanism and an allosteric one with the same catalytic scheme differ. It throws an
`ErrorException` when an allosteric mechanism does not give each kinetic group one known
state or has an oligomeric state below 1, when a plain one has an oligomeric state other
than 1 or a regulatory site, or when a ligand cannot name a metabolite, names an enzyme
form or sits at two sites. Whether its rate law is defined is
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
            _check_tags(states, maximum(e -> e.kinetic_group, catalytic.steps; init = 0),
                        allosteric.oligomeric_state)
        end
        sites = sort(allosteric.regulatory_sites; by = _site_key)
        _check_ligands(sites, catalytic)
        _labels_in_order(catalytic.steps) &&
            return new(catalytic,
                       Allostery(copy(states), allosteric.oligomeric_state, sites))
        group, group_of = _number_groups(Int[e.kinetic_group for e in catalytic.steps],
                                         maximum(e -> e.kinetic_group, catalytic.steps))
        steps = [Step(e.from, e.to, e.taken_up, e.given_off, e.rapid_equilibrium, group[k])
                 for (k, e) in enumerate(catalytic.steps)]
        new(_unchecked_scheme(catalytic.metabolites, catalytic.enzyme_forms,
                              catalytic.substrates, catalytic.products, steps),
            Allostery(isempty(states) ? Symbol[] : _tags_by_group(group_of, states),
                      allosteric.oligomeric_state, sites))
    end
end

"""Whether the kinetic group labels of `steps` are numbered by their first step: the first
step of each group comes after the first step of each group with a smaller label."""
function _labels_in_order(steps::Vector{Step})
    seen = 0
    for e in steps
        e.kinetic_group > seen + 1 && return false
        seen = max(seen, e.kinetic_group)
    end
    true
end

"""A site's place in the order of a mechanism's sites: its ligands, its multiplicity, then
its ligands' states."""
_site_key(s::RegulatorySite) = (s.ligands, s.multiplicity, s.allosteric_states)

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

"""The states by canonical group: `allosteric_states[l]`, the state of written label `l`,
moved to group `group_of[l]`."""
function _tags_by_group(group_of::Vector{Int}, allosteric_states::Vector{Symbol})
    tags = similar(allosteric_states)
    for (l, g) in enumerate(group_of)
        tags[g] = allosteric_states[l]
    end
    tags
end

"""Error unless `allosteric_states` holds one valid state (`_ALLOSTERIC_STATES`) for each
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

"""Error when a site ligand cannot name a metabolite (`_check_metabolite_names`) or names
an enzyme form of `catalytic`, which a text would read as the form, or when a ligand sits
on two sites: each ligand's constant is named after the ligand alone."""
function _check_ligands(sites::Vector{RegulatorySite}, catalytic::CatalyticScheme)
    ligs = Symbol[l for s in sites for l in s.ligands]
    _check_metabolite_names("Mechanism", ligs)
    for l in ligs
        l in catalytic.enzyme_forms &&
            error("Mechanism: ligand $l is also the name of an enzyme form; a text could " *
                  "not tell them apart")
    end
    dup = findfirst(i -> ligs[i] in view(ligs, 1:(i - 1)), eachindex(ligs))
    dup === nothing ||
        error("Mechanism: ligand $(ligs[dup]) appears in two " *
              "distinct regulatory sites; rendered regulator parameter names would " *
              "collide. Same-ligand-two-sites is not enumerated.")
    nothing
end

# ─── Identity ───────────────────────────────────────────────────────────

# `==` and `hash` by content for every value type whose fields hold vectors or structs,
# whose default `==` would compare identity: each field in declaration order, `hash`
# seeded with the type's name. A struct's default hash mixes in its type's identity, which
# changes between precompile sessions; a value's hash reads only symbols, integers and
# flags, so every process agrees on it.
for T in (RegulatorySite, Reactant, Regulator, EnzymeReaction, Step, CatalyticScheme,
          Allostery, Mechanism)
    fs = fieldnames(T)
    @eval Base.:(==)(a::$T, b::$T) =
        $(reduce((l, r) -> :($l && $r), [:(a.$f == b.$f) for f in fs]))
    @eval Base.hash(x::$T, h::UInt) =
        $(foldl((acc, f) -> :(hash(x.$f, $acc)), fs;
                init = :(hash($(QuoteNode(nameof(T))), h))))
end

# ─── Reads ──────────────────────────────────────────────────────────────

"""
    metabolites(m) -> Vector{Symbol}

The names of the metabolites of mechanism `m`, each once: those of its catalytic scheme,
its substrates, then its products, then its modifiers, each sorted by name, then the
ligands of an allosteric mechanism's regulatory sites that are not among them, site by
site. These are the concentrations its rate law reads, in the order `metabolites(rate)` of
its compiled law
([`compile_rate_equation`](@ref EnzymeRates.Derivation.compile_rate_equation)) gives
them as a tuple.

```jldoctest
julia> using EnzymeRates

julia> m = @enzyme_mechanism begin
           substrates: S
           products: P
           steps: begin
               E + S <--> ES
               ES <--> E + P
           end
       end;

julia> metabolites(m)
2-element Vector{Symbol}:
 :S
 :P
```
"""
metabolites(m::Mechanism) =
    unique!([m.catalytic.metabolites;
             [l for s in m.allosteric.regulatory_sites for l in s.ligands]])
