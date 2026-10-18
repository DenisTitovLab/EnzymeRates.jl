# ABOUTME: The mechanism value types as a reaction graph of enzyme forms and steps in
# ABOUTME: canonical form: raw input to checked graph, content order, identity, Sig.
module ReactionGraph

using ..EnzymeRates: Met, RawForm, RawStep, RawMechanism, RawSite, RawAllostericMechanism,
                     RegulatorySite, _VALID_CAT_ALLO_STATES, _check_metabolite_names

# ─── Types ──────────────────────────────────────────────────────────────

"""
Copies of each metabolite-table entry, four bits per entry: entry `i` holds bits
`4(i - 1)` to `4i - 1`. A mechanism therefore holds at most 32 table entries and at most
15 copies of one entry on a form or a step side.
"""
const Counts = UInt128

"""How many copies of table entry `i` `c` counts."""
copies(c::Counts, i::Int) = Int((c >> (4 * (i - 1))) & 0xf)

"""The number of copies `c` counts over every entry."""
function total(c::Counts)
    n = 0
    while c != 0
        n += Int(c & 0xf)
        c >>= 4
    end
    n
end

"""The table entries `c` counts, ascending, each as often as counted."""
function entries(c::Counts)
    out = Int[]
    i = 1
    while c != 0
        for _ in 1:Int(c & 0xf)
            push!(out, i)
        end
        c >>= 4
        i += 1
    end
    out
end

"""
An enzyme form: its conformation label (an index into the mechanism's sorted
`conformations`), the table entries it binds, and its residual: the substrates it adds
and the products it subtracts.
"""
struct Form
    conformation::Int
    bound::Counts
    added::Counts
    subtracted::Counts
end

"""
A step between forms `from` and `to`, stored in its physical-forward direction: going
`from → to` it takes up `consumed` and gives off `released`. `rapid_equilibrium` marks
a rapid-equilibrium step.
"""
struct Edge
    from::Int
    to::Int
    consumed::Counts
    released::Counts
    rapid_equilibrium::Bool
end

"""
    Mechanism

A mechanism as its reaction graph, in canonical form: two writings of one mechanism give
equal fields. `substrates`, `products` and `conformations` are sorted; `table` holds
each metabolite entry once, sorted by (name, modifier), a reactant and its modifier copy
being two entries; `forms` hold free E first, then the other forms in content order;
`edges` run physical-forward, in content order; `group[k]` is edge `k`'s kinetic group,
the groups numbered by their first edge. `Mechanism(raw::RawMechanism)` builds one and
checks it; the field constructor checks nothing and serves canonical form and the Sig
decoder.
"""
struct Mechanism
    substrates::Vector{Symbol}
    products::Vector{Symbol}
    table::Vector{Met}
    conformations::Vector{Symbol}
    forms::Vector{Form}
    edges::Vector{Edge}
    group::Vector{Int}
end

"""
    AllostericMechanism

An MWC enzyme: its catalytic `Mechanism`, one allosteric tag per kinetic group of it
(`cat_allo_states[g]`), the catalytic multiplicity and the regulatory sites, sorted. A
future model with its own derivation wraps the same `Mechanism` in a type of its own.
"""
struct AllostericMechanism
    catalytic::Mechanism
    cat_allo_states::Vector{Symbol}
    catalytic_multiplicity::Int
    regulatory_sites::Vector{RegulatorySite}
    function AllostericMechanism(catalytic::Mechanism, cat_allo_states::Vector{Symbol},
                                 catalytic_multiplicity::Int,
                                 regulatory_sites::Vector{RegulatorySite})
        _check_tags(cat_allo_states, ngroups(catalytic), catalytic_multiplicity)
        sites = sort(regulatory_sites; by = _site_key)
        _check_ligands(sites)
        new(catalytic, copy(cat_allo_states), catalytic_multiplicity, sites)
    end
end

"""A site's place in the order of a mechanism's sites: its ligands, its multiplicity,
then its ligands' tags."""
_site_key(s::RegulatorySite) = (s.ligands, s.multiplicity, s.allo_states)

abstract type AbstractEnzymeMechanism end

"""A plain mechanism as a type: `Sig`, its only parameter, is the mechanism encoded as
integers (`_sig_of`), so equal mechanisms share one type and the derivation compiles once
per mechanism."""
struct EnzymeMechanism{Sig} <: AbstractEnzymeMechanism end

"""An allosteric mechanism as a type: `Sig`, its only parameter, is the mechanism
encoded as integers (`_sig_of`)."""
struct AllostericEnzymeMechanism{Sig} <: AbstractEnzymeMechanism end

# ─── Content order ──────────────────────────────────────────────────────

"""
Whether the entries `a` counts, listed ascending with repeats, come before those of `b`
in lexicographic order, a list coming before every longer list it begins. With the table
sorted by (name, modifier), this is the order of the two metabolite lists by name, then
by modifier.
"""
function _seq_less(a::Counts, b::Counts)
    a == b && return false
    i = trailing_zeros(a ⊻ b) >> 2           # first entry, from 0, whose copies differ
    ca, cb = (a >> (4i)) & 0xf, (b >> (4i)) & 0xf
    # Both lists go on with min(ca, cb) copies of entry i. The one with more copies goes
    # on with entry i; the other goes on with a later entry, or ends.
    ca > cb ? (b >> (4i + 4)) != 0 : (a >> (4i + 4)) == 0
end

"""Forms by label, then bound entries, residual additions and residual subtractions,
each entry list in lexicographic order (`_seq_less`)."""
function Base.isless(f::Form, g::Form)
    f.conformation != g.conformation && return f.conformation < g.conformation
    f.bound != g.bound && return _seq_less(f.bound, g.bound)
    f.added != g.added && return _seq_less(f.added, g.added)
    _seq_less(f.subtracted, g.subtracted)
end

"""Edges by form indices, then consumed and released entries as lists (`_seq_less`), then
the flag, steady state first. Tuples and `argmin` over edges compare by it, never by the
packed numbers."""
function Base.isless(e::Edge, f::Edge)
    e.from != f.from && return e.from < f.from
    e.to != f.to && return e.to < f.to
    e.consumed != f.consumed && return _seq_less(e.consumed, f.consumed)
    e.released != f.released && return _seq_less(e.released, f.released)
    e.rapid_equilibrium < f.rapid_equilibrium
end

"""`e` stored in the opposite direction: forms swapped, consumed and released swapped."""
_reversed(e::Edge) = Edge(e.to, e.from, e.released, e.consumed, e.rapid_equilibrium)

# ─── Canonical form ─────────────────────────────────────────────────────

"""
    _canonicalize(raw::RawMechanism)
        -> (; mechanism, edge_of, flipped, group_of, form_of, misaligned)

The canonical `Mechanism` that `raw` writes, with no representation check. `edge_of[j]`
is the edge written step `j` became, `flipped[j]` whether that edge runs against step
`j` as written, `group_of[l]` the kinetic group written label `l` became,
`form_of[2j - 1]` and `form_of[2j]` the forms step `j`'s two written ends became, and
`misaligned` the lead and member of a tied kinetic group that no orientation aligns, or
no edge. Alignment stops at that group, so a result with a misaligned pair is not
canonical. Errors only where no encoding exists: more than 32 table entries or more than
15 copies of one entry on a form or a step side. The checked constructors build on it.
"""
function _canonicalize(raw::RawMechanism)
    substrates, products = sort(raw.substrates), sort(raw.products)
    table = _metabolite_table(raw)
    conformations = sort!(unique!([f.conformation for s in raw.steps
                                   for f in (s.from, s.to)]))
    n = length(raw.steps)
    ends = Vector{Form}(undef, 2n)
    for (j, s) in enumerate(raw.steps)
        ends[2j - 1] = _form(table, conformations, s.from)
        ends[2j] = _form(table, conformations, s.to)
    end
    forms, form_of = _sorted_forms(ends, conformations)
    written = [Edge(form_of[2j - 1], form_of[2j], _counts(table, s.consumed),
                    _counts(table, s.released), s.rapid_equilibrium)
               for (j, s) in enumerate(raw.steps)]
    substrate_mask, product_mask = _role_masks(table, substrates)
    fwd, misaligned = _orientation(forms, written, raw.group, substrate_mask, product_mask)
    stored = [fwd[j] ? written[j] : _reversed(written[j]) for j in 1:n]
    p = sortperm(stored)
    group, group_of = _number_groups(raw.group[p], maximum(raw.group; init = 0))
    m = Mechanism(substrates, products, table, conformations, forms, stored[p], group)
    (; mechanism = m, edge_of = invperm(p), flipped = .!fwd, group_of, form_of,
     misaligned)
end

"""Every metabolite entry `raw` holds, sorted and once each: bound, consumed and released
metabolites, and the residual names as non-modifier entries."""
function _metabolite_table(raw::RawMechanism)
    seen = Met[]
    for s in raw.steps
        append!(seen, s.consumed, s.released)
        for f in (s.from, s.to)
            append!(seen, f.bound)
            foreach(x -> push!(seen, Met(x, false)), f.added)
            foreach(x -> push!(seen, Met(x, false)), f.subtracted)
        end
    end
    table = unique!(sort!(seen))
    length(table) <= 32 || error("Mechanism: a mechanism holds at most 32 metabolite " *
        "entries (a name and its modifier copy count apart); this one has " *
        "$(length(table))")
    table
end

"""The counts of `mets`, entries of the sorted `table`."""
function _counts(table::Vector{Met}, mets)
    c = Counts(0)
    for x in mets
        i = searchsortedfirst(table, x)
        copies(c, i) == 15 && error("Mechanism: a form or a step side holds $(x.name) " *
            "more than 15 times; a mechanism holds at most 15 copies of one metabolite " *
            "there")
        c += Counts(1) << (4 * (i - 1))
    end
    c
end

"""Form `f` as entries of the sorted `table` and `conformations`."""
_form(table::Vector{Met}, conformations::Vector{Symbol}, f::RawForm) =
    Form(searchsortedfirst(conformations, f.conformation), _counts(table, f.bound),
         _counts(table, (Met(x, false) for x in f.added)),
         _counts(table, (Met(x, false) for x in f.subtracted)))

"""The distinct forms of `ends` in canonical order, free E (label `E`, nothing bound, no
residual) first, and the index each end became."""
function _sorted_forms(ends::Vector{Form}, conformations::Vector{Symbol})
    free_E = Form(something(findfirst(==(:E), conformations), 0), 0, 0, 0)
    forms, form_of = Form[], Vector{Int}(undef, length(ends))
    for t in sortperm(ends; by = f -> (f != free_E, f))
        (isempty(forms) || forms[end] != ends[t]) && push!(forms, ends[t])
        form_of[t] = length(forms)
    end
    forms, form_of
end

"""Masks of the table's substrate entries and of its product entries; a modifier entry
is in neither, so it counts on neither side of a reaction."""
function _role_masks(table::Vector{Met}, substrates::Vector{Symbol})
    substrate_mask = product_mask = Counts(0)
    for (i, x) in enumerate(table)
        x.modifier && continue
        nibble = Counts(0xf) << (4 * (i - 1))
        insorted(x.name, substrates) ? (substrate_mask |= nibble) :
                                       (product_mask |= nibble)
    end
    substrate_mask, product_mask
end

"""
    _orientation(forms, written, labels, substrate_mask, product_mask)
        -> (fwd::BitVector, misaligned::Vector{Edge})

Whether each written edge runs in its physical-forward direction. A step with one
metabolite in all takes it up (tier 0). Otherwise the side holding more substrates,
bound plus free, then fewer products, comes first (tier 1); otherwise the step runs from
a form where only products enter or leave to one where only substrates do (tier 2);
otherwise from the earlier form (tier 3). A modifier entry counts on neither side. Entry
sides are read over every edge as written; reversing an edge swaps its ends and lists
together, so they do not depend on any written direction. Tied groups then turn as one
(`_align_tied_groups!`).
"""
function _orientation(forms::Vector{Form}, written::Vector{Edge}, labels::Vector{Int},
                      substrate_mask::Counts, product_mask::Counts)
    side(c::Counts) = UInt8((c & substrate_mask != 0) + 2 * (c & product_mask != 0))
    enters = zeros(UInt8, length(forms))
    for e in written
        enters[e.from] |= side(e.consumed)
        enters[e.to] |= side(e.released)
    end
    score(f::Int, free::Counts) =
        (total(forms[f].bound & substrate_mask) + total(free & substrate_mask),
         -total(forms[f].bound & product_mask) - total(free & product_mask))
    n = length(written)
    fwd, tier = trues(n), zeros(Int, n)
    for (k, e) in enumerate(written)
        if total(e.consumed) + total(e.released) == 1
            fwd[k] = e.consumed != 0
        elseif (a = score(e.from, e.consumed)) != (b = score(e.to, e.released))
            fwd[k], tier[k] = a > b, 1
        elseif (enters[e.from], enters[e.to]) in ((0x02, 0x01), (0x01, 0x02))
            fwd[k], tier[k] = enters[e.from] == 0x02, 2
        else
            fwd[k], tier[k] = e.from < e.to, 3
        end
    end
    misaligned = _align_tied_groups!(fwd, tier, forms, written, labels)
    fwd, misaligned
end

"""
Turn the steps tiers 2 and 3 decide in each kinetic group (`labels`) with the group's
lead: the tied step that changes conformation first, then by tier, then by its oriented
edge. A group's steps share their constants, so a conformational change and its mirror
at other bound forms must run the same way. Every other tied step that takes up and
gives off the lead's metabolites takes the orientation that does, and, if it changes
conformation, makes the lead's change; a step that fits both ways keeps its own. Returns
the lead and the first step that fits the lead's metabolites but not its change, or no
edge.
"""
function _align_tied_groups!(fwd::BitVector, tier::Vector{Int}, forms::Vector{Form},
                             written::Vector{Edge}, labels::Vector{Int})
    label_pair(e::Edge) = (forms[e.from].conformation, forms[e.to].conformation)
    turns(e::Edge) = label_pair(e)[1] != label_pair(e)[2]
    oriented(k::Int) = fwd[k] ? written[k] : _reversed(written[k])
    for l in unique(labels)
        tied = [k for k in eachindex(written) if labels[k] == l && tier[k] > 1]
        length(tied) < 2 && continue
        lead = oriented(argmin(k -> (!turns(oriented(k)), tier[k], oriented(k)), tied))
        for k in tied
            o = oriented(k)
            same = [x for x in (o, _reversed(o))
                    if (x.consumed, x.released) == (lead.consumed, lead.released)]
            isempty(same) && continue
            i = findfirst(x -> !turns(x) || label_pair(x) == label_pair(lead), same)
            i === nothing && return [lead, same[1]]
            fwd[k] = same[i] == written[k]
        end
    end
    Edge[]
end

"""Group labels in edge order renumbered by first occurrence, and the number each
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
    Mechanism(raw::RawMechanism)

The mechanism `raw` writes, in canonical form, checked: its declarations
(`_check_declarations`), then its representation (`_check_representation`). Whether its
rate law is defined is the boundary check's question, not the constructor's.
"""
function Mechanism(raw::RawMechanism)
    _check_declarations(raw)
    c = _canonicalize(raw)
    _check_representation(c.mechanism, c.misaligned, raw.modifiers)
    c.mechanism
end

"""
    AllostericMechanism(raw::RawAllostericMechanism)

The allosteric mechanism `raw` writes: its regulatory sites, its tags (each written
group label's), then its catalytic part as `Mechanism(raw.catalytic)` checks it, each
tag moved to the canonical group its label became.
"""
function AllostericMechanism(raw::RawAllostericMechanism)
    sites = [RegulatorySite(s.ligands, s.multiplicity, s.allo_states) for s in raw.sites]
    _check_tags(raw.cat_allo_states, maximum(raw.catalytic.group; init = 0),
                raw.catalytic_multiplicity)
    _check_declarations(raw.catalytic)
    c = _canonicalize(raw.catalytic)
    _check_representation(c.mechanism, c.misaligned, raw.catalytic.modifiers)
    tags = similar(raw.cat_allo_states)
    for (l, g) in enumerate(c.group_of)
        tags[g] = raw.cat_allo_states[l]
    end
    AllostericMechanism(c.mechanism, tags, raw.catalytic_multiplicity, sites)
end

# ─── Checks ─────────────────────────────────────────────────────────────

"""Error unless `cat_allo_states` holds one valid tag for each of `G` kinetic groups and
the catalytic multiplicity is at least 1. No catalytic group is `:OnlyI`: the active
conformation is the one that catalyzes."""
function _check_tags(cat_allo_states::Vector{Symbol}, G::Int, catalytic_multiplicity::Int)
    length(cat_allo_states) == G ||
        error("AllostericMechanism: cat_allo_states length $(length(cat_allo_states)) " *
              "must match the number of kinetic groups $G")
    catalytic_multiplicity ≥ 1 ||
        error("AllostericMechanism: catalytic_multiplicity must be ≥ 1, " *
              "got $catalytic_multiplicity")
    bad = findfirst(∉(_VALID_CAT_ALLO_STATES), cat_allo_states)
    bad === nothing ||
        error("AllostericMechanism: catalytic group $bad has invalid " *
              "allo state $(cat_allo_states[bad]) (must be one of " *
              "$_VALID_CAT_ALLO_STATES); :OnlyI is rejected for " *
              "catalytic groups (active-state-active convention)")
    nothing
end

"""Error when a site ligand's name is reserved or holds `_`, or when a ligand sits on two
sites: each ligand's constant is named after the ligand alone."""
function _check_ligands(sites::Vector{RegulatorySite})
    ligs = [l for s in sites for l in s.ligands]
    _check_metabolite_names("AllostericMechanism", ligs)
    dup = findfirst(i -> ligs[i] in view(ligs, 1:(i - 1)), eachindex(ligs))
    dup === nothing ||
        error("AllostericMechanism: ligand $(ligs[dup]) appears in two " *
              "distinct regulatory sites; rendered Kreg names would " *
              "collide. Same-ligand-two-sites is not enumerated.")
    nothing
end

"""Error unless `raw` declares substrates and products, each name once, no name both a
substrate and a product, and gives each step one kinetic group label from 1 to the
number of groups, every label used."""
function _check_declarations(raw::RawMechanism)
    isempty(raw.substrates) && error("Mechanism: substrates must not be empty")
    isempty(raw.products) && error("Mechanism: products must not be empty")
    for (names, kind) in ((raw.substrates, "substrate"), (raw.products, "product"),
                          (raw.modifiers, "modifier"))
        allunique(names) || error("Mechanism: duplicate $kind names")
    end
    both = findfirst(in(raw.products), raw.substrates)
    both === nothing ||
        error("Mechanism: $(raw.substrates[both]) named as both a substrate and a " *
              "product; concentrations and constants are keyed by name")
    length(raw.group) == length(raw.steps) && all(≥(1), raw.group) ||
        error("Mechanism: each step takes one kinetic group label from 1 to the number " *
              "of groups")
    empty = findfirst(∉(raw.group), 1:maximum(raw.group; init = 0))
    empty === nothing ||
        error("Mechanism: kinetic group $empty is empty; a kinetic group holds at least " *
              "one step")
    nothing
end

"""
Error unless the canonical mechanism `m` represents one: every name declared (in
`modifiers`, the declared modifiers, or among its substrates and products), a modifier
entry only for a declared modifier and a plain entry only for a reactant, residuals that
add substrates and subtract products; no step joining a form to itself, taking up and
giving off one metabolite, or binding a metabolite into a form that holds it in another
role; tied steps that make one conformational change (`misaligned`, from
`_canonicalize`); no two forms that render one name; each kinetic group one kind of step
with one flag; each reaction once.
"""
function _check_representation(m::Mechanism, misaligned::Vector{Edge},
                               modifiers::Vector{Symbol})
    for x in m.table
        _is_reactant(m, x.name) || x.name in modifiers ||
            error("Mechanism: $(x.name) is not declared; declare it under substrates:, " *
                  "products: or modifiers:")
        x.modifier && !(x.name in modifiers) &&
            error("Mechanism: $(x.name) binds as a modifier, but modifiers: does not " *
                  "declare it")
        x.modifier || _is_reactant(m, x.name) ||
            error("Mechanism: $(x.name) is a modifier, so it binds only as a modifier")
    end
    substrate_mask, product_mask = _role_masks(m)
    for f in m.forms, (bad, verb) in ((f.added & product_mask, "adds"),
                                      (f.subtracted & substrate_mask, "subtracts"))
        bad == 0 ||
            error("Mechanism: a residual adds substrates and subtracts products; one " *
                  "$verb $(m.table[first(entries(bad))].name)")
    end
    for e in m.edges
        e.from == e.to &&
            error("Step: both ends are the form $(_form_name(m, e.from)); a step joins " *
                  "two forms")
        shared = _occupied(e.consumed) & _occupied(e.released)
        shared == 0 ||
            error("Step $(_form_name(m, e.from)) → $(_form_name(m, e.to)): " *
                  "$(join([m.table[i].name for i in entries(shared)], ", ")) both " *
                  "consumed and released")
        _binds_by_name(m, e) && !_is_plain_binding(m, e) &&
            error("Step $(_form_name(m, e.from)) → $(_form_name(m, e.to)) changes a " *
                  "metabolite's role: $(_form_name(m, e.to)) holds " *
                  "$(_form_name(m, e.from))'s metabolites plus " *
                  "$(m.table[_bound_entry(e)].name) by name, not by role; write a " *
                  "modifier copy as X::Mod on both sides")
    end
    isempty(misaligned) ||
        error("Mechanism: a kinetic group holds $(join(_sides(m, misaligned[1]), " → ")) " *
              "and $(join(_sides(m, misaligned[2]), " → ")), which make different " *
              "conformational changes; the steps of a kinetic group share their " *
              "constants, so they must make the same conformational change")
    _check_distinct_form_names(m)
    groups = _group_members(m)
    _check_uniform_groups(m, groups)
    _check_each_reaction_once(m, groups)
end

"""A count of 1 for every entry `c` counts at least once."""
function _occupied(c::Counts)
    x = c | (c >> 1)
    x |= x >> 2
    x & Counts(0x11111111111111111111111111111111)
end

"""Error when two forms of `m` render one name: every constant is named after the forms
its step joins (`_form_name`)."""
function _check_distinct_form_names(m::Mechanism)
    seen = Dict{String, Int}()
    for f in eachindex(m.forms)
        other = get!(seen, _form_name(m, f), f)
        other == f ||
            error("Mechanism: the enzyme forms $(_form_text(m, other)) and " *
                  "$(_form_text(m, f)) both render the name $(_form_name(m, f)); every " *
                  "constant is named after the forms its step joins, so rename a " *
                  "metabolite to tell the two forms apart")
    end
end

"""Error unless the edges of each kinetic group (`groups`) take up and give off the same
entries and carry one flag: they share their constants."""
function _check_uniform_groups(m::Mechanism, groups::Vector{Vector{Int}})
    label(e::Edge) = join(_sides(m, e), " → ") * (e.rapid_equilibrium ? " (RE)" : " (SS)")
    for members in groups
        e1 = m.edges[first(members)]
        for k in members
            e = m.edges[k]
            (e.consumed, e.released, e.rapid_equilibrium) ==
                (e1.consumed, e1.released, e1.rapid_equilibrium) && continue
            error("Mechanism: a kinetic group holds $(label(e1)) and $(label(e)), " *
                  "which differ in kind or in RE/SS; the steps of a kinetic group " *
                  "share their constants, so they must be the same kind of step " *
                  "with the same flag")
        end
    end
end

"""Error when one reaction, the unordered pair of an edge's two sides, appears in two
edges: a reaction belongs to one kinetic group, whose constants are named after it, and
one group holding it twice would count its edge twice. Reads the groups (`groups`) in
order, each in edge order, and names the reaction as the edge seen first stores it."""
function _check_each_reaction_once(m::Mechanism, groups::Vector{Vector{Int}})
    seen = Dict{NTuple{2, Tuple{Int, Counts}}, Tuple{Int, Int}}()
    flag(e::Edge) = e.rapid_equilibrium ? "RE" : "SS"
    for (g, members) in enumerate(groups), k in members
        e = m.edges[k]
        key = minmax((e.from, e.consumed), (e.to, e.released))
        if haskey(seen, key)
            g0, k0 = seen[key]
            text = join(_sides(m, k0), " ⇌ ")
            g0 == g && error("Mechanism: kinetic group $g holds the reaction $text " *
                             "twice; a mechanism holds each reaction once")
            e0 = m.edges[k0]
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

# `hash` by content for every graph type, and `==` by content for those whose fields hold
# vectors: each field in declaration order, seeded with the type's name. A struct's
# default hash mixes in its type's identity, which changes between precompile sessions;
# these hashes read only symbols, integers and flags, so every process agrees on them.
# `Form` and `Edge` keep Julia's `==`, which compares their fields.
for T in (Form, Edge, Mechanism, AllostericMechanism)
    @eval Base.hash(x::$T, h::UInt) =
        $(foldl((acc, f) -> :(hash(x.$f, $acc)), fieldnames(T);
                init = :(hash($(QuoteNode(nameof(T))), h))))
end
for T in (Mechanism, AllostericMechanism)
    @eval Base.:(==)(a::$T, b::$T) =
        $(reduce((l, r) -> :($l && $r), [:(a.$f == b.$f) for f in fieldnames(T)]))
end

# ─── Sig ────────────────────────────────────────────────────────────────
#
# A mechanism's Sig, every entry an `Int`:
#   plain       the numbers of substrates, products, table entries, conformation labels,
#               forms and edges; the substrate and product names; each table entry's
#               name and modifier flag; the conformation labels; each form (label,
#               bound, added, subtracted); each edge (from, to, consumed, released,
#               flag); each edge's group. A packed count takes two entries, high half
#               first.
#   name        its number of characters, then each character's code point.
#   allosteric  the plain Sig of its catalytic mechanism; one tag code per group; the
#               catalytic multiplicity; the number of sites, then per site its
#               multiplicity, its number of ligands and per ligand its name and tag code.

"""The allosteric tags, coded by position: catalytic groups and site ligands share it."""
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
_sig_counts!(out::Vector{Int}, c::Counts) = push!(out, (c >> 64) % Int, c % Int)

function _sig_mechanism!(out::Vector{Int}, m::Mechanism)
    push!(out, length(m.substrates), length(m.products), length(m.table),
          length(m.conformations), length(m.forms), length(m.edges))
    foreach(s -> _sig_name!(out, s), m.substrates)
    foreach(s -> _sig_name!(out, s), m.products)
    for x in m.table
        _sig_name!(out, x.name)
        push!(out, x.modifier)
    end
    foreach(s -> _sig_name!(out, s), m.conformations)
    for f in m.forms
        push!(out, f.conformation)
        _sig_counts!(out, f.bound)
        _sig_counts!(out, f.added)
        _sig_counts!(out, f.subtracted)
    end
    for e in m.edges
        push!(out, e.from, e.to)
        _sig_counts!(out, e.consumed)
        _sig_counts!(out, e.released)
        push!(out, e.rapid_equilibrium)
    end
    append!(out, m.group)
end

"""`m` as one homogeneous integer tuple, the only type parameter of its compiled type."""
_sig_of(m::Mechanism) = Tuple(_sig_mechanism!(Int[], m))
function _sig_of(am::AllostericMechanism)
    out = _sig_mechanism!(Int[], am.catalytic)
    foreach(t -> push!(out, _tag_code(t)), am.cat_allo_states)
    push!(out, am.catalytic_multiplicity, length(am.regulatory_sites))
    for s in am.regulatory_sites
        push!(out, s.multiplicity, length(s.ligands))
        for (l, t) in zip(s.ligands, s.allo_states)
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
function _next_counts!(r::_SigReader)
    hi = _next!(r) % UInt64
    lo = _next!(r) % UInt64
    (Counts(hi) << 64) | Counts(lo)
end

function _read_mechanism!(r::_SigReader)
    ns = _next!(r); np = _next!(r); nt = _next!(r)
    nc = _next!(r); nf = _next!(r); ne = _next!(r)
    substrates = [_next_name!(r) for _ in 1:ns]
    products = [_next_name!(r) for _ in 1:np]
    table = [Met(_next_name!(r), _next!(r) == 1) for _ in 1:nt]
    conformations = [_next_name!(r) for _ in 1:nc]
    forms = [Form(_next!(r), _next_counts!(r), _next_counts!(r), _next_counts!(r))
             for _ in 1:nf]
    edges = [Edge(_next!(r), _next!(r), _next_counts!(r), _next_counts!(r),
                  _next!(r) == 1) for _ in 1:ne]
    group = [_next!(r) for _ in 1:ne]
    Mechanism(substrates, products, table, conformations, forms, edges, group)
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
once for every Sig length, as long as callers splat a Sig whose length their code leaves
open."""
_mechanism_from_sig(sig::Int...) = _read_mechanism!(_SigReader(collect(sig), 0))
"""The allosteric mechanism whose Sig is `sig`, splatted; only the composition
constructor's tag and site checks run."""
function _allosteric_mechanism_from_sig(sig::Int...)
    r = _SigReader(collect(sig), 0)
    m = _read_mechanism!(r)
    tags = [_TAG_CODES[_next!(r)] for _ in 1:ngroups(m)]
    multiplicity = _next!(r)
    sites = [_read_site!(r) for _ in 1:_next!(r)]
    AllostericMechanism(m, tags, multiplicity, sites)
end

EnzymeMechanism(m::Mechanism) = EnzymeMechanism{_sig_of(m)}()
AllostericEnzymeMechanism(am::AllostericMechanism) =
    AllostericEnzymeMechanism{_sig_of(am)}()
compile_mechanism(m::Mechanism) = EnzymeMechanism(m)
compile_mechanism(am::AllostericMechanism) = AllostericEnzymeMechanism(am)

# The lifts read the Sig at run time, so they compile once and no method specializes on
# a Sig but the derivation's own.
Mechanism(@nospecialize(em::EnzymeMechanism)) =
    _mechanism_from_sig(typeof(em).parameters[1]...)
AllostericMechanism(@nospecialize(aem::AllostericEnzymeMechanism)) =
    _allosteric_mechanism_from_sig(typeof(aem).parameters[1]...)
_concrete(@nospecialize(em::EnzymeMechanism)) = Mechanism(em)
_concrete(@nospecialize(aem::AllostericEnzymeMechanism)) = AllostericMechanism(aem)
catalytic_mechanism(@nospecialize(aem::AllostericEnzymeMechanism)) =
    EnzymeMechanism(AllostericMechanism(aem).catalytic)
catalytic_multiplicity(@nospecialize(aem::AllostericEnzymeMechanism)) =
    AllostericMechanism(aem).catalytic_multiplicity

# ─── Reads ──────────────────────────────────────────────────────────────

substrates(m::Mechanism) = m.substrates
products(m::Mechanism) = m.products
"""Every name bound as a modifier, reactant copies included, sorted."""
modifiers(m::Mechanism) = unique!([x.name for x in m.table if x.modifier])
# The concentration names in the order the rate law reads them: substrates, products,
# then the modifier-only names, each by name; an allosteric mechanism appends each site's
# ligands in site order, once each.
metabolites(m::Mechanism) = unique!([m.substrates; m.products; modifiers(m)])
metabolites(am::AllostericMechanism) =
    unique!([metabolites(am.catalytic);
             [l for s in am.regulatory_sites for l in s.ligands]])
"""Whether `x` names a substrate or a product of `m`."""
_is_reactant(m::Mechanism, x::Symbol) = insorted(x, m.substrates) || insorted(x, m.products)
_role_masks(m::Mechanism) = _role_masks(m.table, m.substrates)

catalytic(m::Mechanism) = m
catalytic(am::AllostericMechanism) = am.catalytic
cat_allo_states(am::AllostericMechanism) = am.cat_allo_states
cat_allo_state(am::AllostericMechanism, g::Int) = am.cat_allo_states[g]
catalytic_multiplicity(am::AllostericMechanism) = am.catalytic_multiplicity
regulatory_sites(am::AllostericMechanism) = am.regulatory_sites
"""The number of kinetic groups of `m`."""
ngroups(m::Mechanism) = maximum(m.group; init = 0)
"""Each kinetic group's edges, in edge order."""
function _group_members(m::Mechanism)
    members = [Int[] for _ in 1:ngroups(m)]
    for (k, g) in enumerate(m.group)
        push!(members[g], k)
    end
    members
end
_has_residual(f::Form) = f.added != 0 || f.subtracted != 0
_is_binding(e::Edge) = e.released == 0 && total(e.consumed) == 1
"""The entry a binding takes up; 0 for any other edge."""
_bound_entry(e::Edge) = _is_binding(e) ? (trailing_zeros(e.consumed) >> 2) + 1 : 0
"""Whether binding `e` takes up its entry into a form that holds its source's entries plus
that one, with the same residual; the label may change."""
function _is_plain_binding(m::Mechanism, e::Edge)
    i = _bound_entry(e)
    f, t = m.forms[e.from], m.forms[e.to]
    i != 0 && (f.added, f.subtracted) == (t.added, t.subtracted) &&
        copies(f.bound, i) < 15 && t.bound == f.bound + (Counts(1) << (4 * (i - 1)))
end
"""Whether binding `e` takes up its entry into a form that holds its source's metabolites
plus that one by name, with the same residual: a plain binding, or one that binds a
reactant where its modifier copy sits, or the reverse."""
function _binds_by_name(m::Mechanism, e::Edge)
    i = _bound_entry(e)
    f, t = m.forms[e.from], m.forms[e.to]
    i != 0 && (f.added, f.subtracted) == (t.added, t.subtracted) || return false
    names(c::Counts) = [m.table[j].name for j in entries(c)]
    names(t.bound) == sort!(push!(names(f.bound), m.table[i].name))
end

"""`m` as raw input: every edge a step in its stored direction, labelled with its group.
Every vector is fresh, so a move that edits the raw value never edits `m`."""
function _raw(m::Mechanism)
    form(f::Int) = (x = m.forms[f];
        RawForm(m.conformations[x.conformation], [m.table[i] for i in entries(x.bound)],
                [m.table[i].name for i in entries(x.added)],
                [m.table[i].name for i in entries(x.subtracted)]))
    steps = [RawStep(form(e.from), form(e.to), [m.table[i] for i in entries(e.consumed)],
                     [m.table[i] for i in entries(e.released)], e.rapid_equilibrium)
             for e in m.edges]
    RawMechanism(copy(m.substrates), copy(m.products), modifiers(m), steps,
                 copy(m.group))
end
_raw(am::AllostericMechanism) = RawAllostericMechanism(_raw(am.catalytic),
    copy(am.cat_allo_states), am.catalytic_multiplicity,
    [RawSite(copy(s.ligands), s.multiplicity, copy(s.allo_states))
     for s in am.regulatory_sites])

# ─── Names ──────────────────────────────────────────────────────────────

"""Entry `i` as names write it: its name, with `mod` on a reactant's modifier copy; a
modifier-only name stays bare."""
_entry_label(m::Mechanism, i::Int) =
    m.table[i].modifier && _is_reactant(m, m.table[i].name) ?
        string(m.table[i].name, "mod") : String(m.table[i].name)

"""Form `f`'s name: its label, its bound entries' labels in table order, then any
residual as `_res`, `_plus<name>` per substrate it adds and `_minus<name>` per product it
subtracts (`E`, `EAB`, `EPmodQ`, `EB_res_plusA_minusP`)."""
function _form_name(m::Mechanism, f::Int)
    x = m.forms[f]
    io = IOBuffer()
    print(io, m.conformations[x.conformation])
    foreach(i -> print(io, _entry_label(m, i)), entries(x.bound))
    if _has_residual(x)
        print(io, "_res")
        foreach(i -> print(io, "_plus", m.table[i].name), entries(x.added))
        foreach(i -> print(io, "_minus", m.table[i].name), entries(x.subtracted))
    end
    String(take!(io))
end

"""Form `f` as a mechanism macro writes it: `E(NAD, P)`, `E(A::Mod)`, `E(I)`,
`E(B; residual = A - P)`."""
function _form_text(m::Mechanism, f::Int)
    x = m.forms[f]
    mets = [_is_reactant(m, m.table[i].name) && m.table[i].modifier ?
            "$(m.table[i].name)::Mod" : "$(m.table[i].name)" for i in entries(x.bound)]
    terms = join([m.table[i].name for i in entries(x.added)], " + ") *
            join(" - $(m.table[i].name)" for i in entries(x.subtracted))
    res = _has_residual(x) ? "; residual = " * lstrip(terms) : ""
    "$(m.conformations[x.conformation])(" * join(mets, ", ") * res * ")"
end

"""A step side: form `f`'s name and the labels of the free entries `c`, joined by `_`
(`EA_B` for `E(A) + B`)."""
_side_name(m::Mechanism, f::Int, c::Counts) =
    join([_form_name(m, f); [_entry_label(m, i) for i in entries(c)]], "_")
"""An edge's two sides in its stored direction; every step constant is named after
them."""
_sides(m::Mechanism, e::Edge) =
    (_side_name(m, e.from, e.consumed), _side_name(m, e.to, e.released))
_sides(m::Mechanism, k::Int) = _sides(m, m.edges[k])

end
