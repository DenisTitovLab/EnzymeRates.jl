# ABOUTME: The graph accessors of a mechanism: its species, steps, kinetic groups, the
# ABOUTME: metabolites a form binds or a step exchanges, and an allosteric mechanism's tags.

# ─── Metabolites ────────────────────────────────────────────────────────

"""
    Metabolites

The metabolites an enzyme form binds or a step side exchanges, as a multiset: one sorted
list, a metabolite as often as held. A reactant and its modifier copy are two entries,
read under one name.
"""
struct Metabolites
    mets::Vector{Met}
end

"""The number of copies of metabolite `name` in `s`, a reactant's modifier copies
included."""
copies(s::Metabolites, name::Symbol) = count(x -> x.name === name, s.mets)

"""The names of the metabolites `s` holds, sorted, each as often as held."""
metabolite_names(s::Metabolites) = Symbol[x.name for x in s.mets]

Base.isempty(s::Metabolites) = isempty(s.mets)

"""The number of metabolites `s` holds, each copy counted."""
Base.length(s::Metabolites) = length(s.mets)

"""The copies `s` holds as modifiers: a modifier's, and a reactant's modifier copy, but not
the reactant itself."""
modifier_copies(s::Metabolites) = Metabolites(filter(x -> x.modifier, s.mets))

"""The multiset holding each metabolite as often as the fuller of `s` and `t` does."""
Base.union(s::Metabolites, t::Metabolites) = _combine(max, s, t)

"""The multiset holding each metabolite as often as the emptier of `s` and `t` does."""
Base.intersect(s::Metabolites, t::Metabolites) = _combine(min, s, t)

"""The multiset holding each metabolite `pick(copies in s, copies in t)` times, merged
from the two sorted lists."""
function _combine(pick, s::Metabolites, t::Metabolites)
    a, b, out = s.mets, t.mets, Met[]
    i = j = 1
    while i <= length(a) || j <= length(b)
        x = j > length(b) || (i <= length(a) && isless(a[i], b[j])) ? a[i] : b[j]
        na = nb = 0
        while i <= length(a) && a[i] == x
            i += 1
            na += 1
        end
        while j <= length(b) && b[j] == x
            j += 1
            nb += 1
        end
        for _ in 1:pick(na, nb)
            push!(out, x)
        end
    end
    Metabolites(out)
end

Base.:(==)(s::Metabolites, t::Metabolites) = s.mets == t.mets
Base.hash(s::Metabolites, h::UInt) = hash(s.mets, hash(:Metabolites, h))

# ─── Species ────────────────────────────────────────────────────────────

"""The enzyme forms of `m`, as form numbers."""
species(m::Mechanism) = 1:length(m.nodes)

"""The metabolites form `x` of `m` binds."""
bound(m::Mechanism, x::Int) = Metabolites(m.nodes[x].bound)

"""The residual of form `x` of `m`, as `(added, subtracted)`: the substrates it adds and
the products it subtracts."""
residual(m::Mechanism, x::Int) =
    (Metabolites(Met[Met(n, false) for n in m.nodes[x].added]),
     Metabolites(Met[Met(n, false) for n in m.nodes[x].subtracted]))

"""The conformation label of form `x` of `m`."""
conformation(m::Mechanism, x::Int) = m.nodes[x].conformation

# ─── Steps ──────────────────────────────────────────────────────────────

"""The steps of `m`, as step numbers."""
steps(m::Mechanism) = 1:length(m.steps)

"""The forms step `k` of `m` joins, in its stored direction: physical-forward."""
ends(m::Mechanism, k::Int) = (m.steps[k].from, m.steps[k].to)

"""The metabolites step `k` of `m` takes up, in its stored direction."""
taken_up(m::Mechanism, k::Int) = Metabolites(m.steps[k].taken_up)

"""The metabolites step `k` of `m` gives off, in its stored direction."""
given_off(m::Mechanism, k::Int) = Metabolites(m.steps[k].given_off)

"""Whether step `k` of `m` is at rapid equilibrium."""
is_rapid_equilibrium(m::Mechanism, k::Int) = m.steps[k].rapid_equilibrium

"""The steps of `m` at rapid equilibrium, one flag per step (`is_rapid_equilibrium`)."""
rapid_equilibrium_steps(m::Mechanism) =
    BitVector([is_rapid_equilibrium(m, k) for k in steps(m)])

"""
What step `k` of `m` does: `:isomerization` when it takes up and gives off nothing,
`:binding` when it is a plain binding (`_is_plain_binding`), `:fused_binding` when it
takes up one metabolite and gives off none but its two forms differ by more than that
metabolite, else `:chemistry`.
"""
function step_kind(m::Mechanism, k::Int)
    e = m.steps[k]
    _is_iso(e) ? :isomerization : _is_plain_binding(m, e) ? :binding :
        _is_binding(e) ? :fused_binding : :chemistry
end

# ─── Kinetic groups ─────────────────────────────────────────────────────

"""The kinetic group of step `k` of `m`."""
group(m::Mechanism, k::Int) = m.steps[k].group

"""Each kinetic group's steps, in step order."""
function groups(m::Mechanism)
    members = [Int[] for _ in 1:ngroups(m)]
    for (k, e) in enumerate(m.steps)
        push!(members[e.group], k)
    end
    members
end

"""The steps of kinetic group `g` of `m`, in step order."""
members(m::Mechanism, g::Int) = findall(e -> e.group == g, m.steps)

# ─── Allosteric ─────────────────────────────────────────────────────────

"""The allosteric tag of kinetic group `g` of `am`'s catalytic graph."""
tag(am::AllostericMechanism, g::Int) = am.mwc_tags[g]
