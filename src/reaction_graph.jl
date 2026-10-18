# ABOUTME: The mechanism value types as a reaction graph of enzyme forms and steps in
# ABOUTME: canonical form: raw input to canonical graph, its content order and identity.
module ReactionGraph

using ..EnzymeRates: Met, RawForm, RawMechanism, RegulatorySite, _VALID_CAT_ALLO_STATES,
                     _check_metabolite_names

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
the groups numbered by their first edge. The field constructor checks nothing.
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
15 copies of one entry on a form or a step side.
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

# ─── Reads ──────────────────────────────────────────────────────────────

"""The number of kinetic groups of `m`."""
ngroups(m::Mechanism) = maximum(m.group; init = 0)

end
