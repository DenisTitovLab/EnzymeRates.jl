# ABOUTME: The graph accessors of a mechanism: its species, steps, kinetic groups, the
# ABOUTME: metabolites a form binds or a step exchanges, and an allosteric mechanism's tags.

# ─── Metabolites ────────────────────────────────────────────────────────

"""
    Metabolites

The metabolites an enzyme form binds or a step side exchanges, as a multiset over its
mechanism's metabolite table. A reactant and its modifier copy are two entries of the
table, read under one name.
"""
struct Metabolites
    counts::Counts
    table::Vector{Met}
end

"""The number of copies of metabolite `name` in `s`, a reactant's modifier copies
included."""
function copies(s::Metabolites, name::Symbol)
    n = 0
    for (i, x) in enumerate(s.table)
        x.name === name && (n += copies(s.counts, i))
    end
    n
end

"""The names of the metabolites `s` holds, sorted, each as often as held."""
metabolite_names(s::Metabolites) = Symbol[s.table[i].name for i in entries(s.counts)]

Base.isempty(s::Metabolites) = s.counts == 0

"""The number of metabolites `s` holds, each copy counted."""
Base.length(s::Metabolites) = total(s.counts)

"""The copies `s` holds as modifiers: a modifier's, and a reactant's modifier copy, but not
the reactant itself."""
function modifier_copies(s::Metabolites)
    mask = Counts(0)
    for (i, x) in enumerate(s.table)
        x.modifier && (mask |= Counts(0xf) << (4 * (i - 1)))
    end
    Metabolites(s.counts & mask, s.table)
end

"""The multiset holding each metabolite as often as the fuller of `s` and `t` does."""
Base.union(s::Metabolites, t::Metabolites) = _combine(max, s, t)

"""The multiset holding each metabolite as often as the emptier of `s` and `t` does."""
Base.intersect(s::Metabolites, t::Metabolites) = _combine(min, s, t)

"""The multiset holding each table entry `pick(copies in s, copies in t)` times."""
function _combine(pick, s::Metabolites, t::Metabolites)
    s.table == t.table ||
        throw(ArgumentError("Metabolites: ∪ and ∩ combine the metabolites of one " *
                            "mechanism, over one metabolite table"))
    c = Counts(0)
    for i in eachindex(s.table)
        c |= Counts(pick(copies(s.counts, i), copies(t.counts, i))) << (4 * (i - 1))
    end
    Metabolites(c, s.table)
end

Base.:(==)(s::Metabolites, t::Metabolites) = s.counts == t.counts && s.table == t.table
Base.hash(s::Metabolites, h::UInt) = hash(s.table, hash(s.counts, hash(:Metabolites, h)))

# ─── Species ────────────────────────────────────────────────────────────

"""The enzyme forms of `m`, as form numbers."""
species(m::Mechanism) = 1:length(m.forms)

"""The metabolites form `x` of `m` binds."""
bound(m::Mechanism, x::Int) = Metabolites(m.forms[x].bound, m.table)

"""The residual of form `x` of `m`, as `(added, subtracted)`: the substrates it adds and
the products it subtracts."""
residual(m::Mechanism, x::Int) =
    (Metabolites(m.forms[x].added, m.table), Metabolites(m.forms[x].subtracted, m.table))

"""The conformation label of form `x` of `m`."""
conformation(m::Mechanism, x::Int) = m.conformations[m.forms[x].conformation]

# ─── Steps ──────────────────────────────────────────────────────────────

"""The steps of `m`, as step numbers."""
steps(m::Mechanism) = 1:length(m.edges)

"""The forms step `k` of `m` joins, in its stored direction: physical-forward."""
ends(m::Mechanism, k::Int) = (m.edges[k].from, m.edges[k].to)

"""The metabolites step `k` of `m` takes up, in its stored direction."""
taken_up(m::Mechanism, k::Int) = Metabolites(m.edges[k].consumed, m.table)

"""The metabolites step `k` of `m` gives off, in its stored direction."""
given_off(m::Mechanism, k::Int) = Metabolites(m.edges[k].released, m.table)

"""Whether step `k` of `m` is at rapid equilibrium."""
is_rapid_equilibrium(m::Mechanism, k::Int) = m.edges[k].rapid_equilibrium

"""The steps of `m` at rapid equilibrium, one flag per step (`is_rapid_equilibrium`)."""
rapid_equilibrium_steps(m::Mechanism) =
    BitVector([is_rapid_equilibrium(m, k) for k in steps(m)])

"""
What step `k` of `m` does: `:isomerization` when it takes up and gives off nothing
(`_is_iso`), `:binding` when it is a plain binding (`_is_chemistry` fails), else
`:chemistry`.
"""
function step_kind(m::Mechanism, k::Int)
    e = m.edges[k]
    _is_iso(e) ? :isomerization : _is_chemistry(m, e) ? :chemistry : :binding
end

"""Whether step `k` of `m` is a binding: it takes up one metabolite and gives off none,
whether or not its two forms differ by that metabolite alone (a fused binding is of kind
`:chemistry`, `step_kind`)."""
is_binding(m::Mechanism, k::Int) = _is_binding(m.edges[k])

# ─── Kinetic groups ─────────────────────────────────────────────────────

"""The kinetic group of step `k` of `m`."""
group(m::Mechanism, k::Int) = m.group[k]

"""The steps of kinetic group `g` of `m`, in step order."""
members(m::Mechanism, g::Int) = findall(==(g), m.group)

# ─── Allosteric ─────────────────────────────────────────────────────────

"""The allosteric tag of kinetic group `g` of `am`'s catalytic graph."""
tag(am::AllostericMechanism, g::Int) = am.mwc_tags[g]
