# ABOUTME: The guard against a wrong rate law, check_valid_rate_equation: its verdicts,
# ABOUTME: networks and faces, classes 1 to 3 with the limit test, and its error messages.

# ─── check_valid_rate_equation ──────────────────────────────────────────

"""
    check_valid_rate_equation(m) -> true
    check_valid_rate_equation(Bool, m) -> Bool

Whether every cycle of `m`'s steps takes up and gives off each metabolite as the reaction
does, `m`'s rate law is the unique limit of mass-action kinetics, as its rapid-equilibrium
steps speed up and its tagged constants vanish, on every face where data sit, and each of
its constants has a name of its own: the guard (`_guard`) finds no fault. The one-argument
form throws an error whose message (`_guard_message`) names each step at fault by its DSL
line, or the names its constants share, and returns `true` otherwise; the `Bool` form
returns the verdict (`false` on a fault) instead of throwing. It reads the reaction graph
through masks and builds nothing. The mechanism macro and the parser call the one-argument
form, and every derivation entry runs its check; both forms read the verdict through the
memo (`_Memo`), and a compiled law's generated methods take it afresh (`_derivation_input`).
The constructors and the Sig decoder do not call it, so code that builds a mechanism
internally never trips a check meant for the mechanisms it returns.
"""
function check_valid_rate_equation(m::Mechanisms.Mechanism)
    _check_verdict(m, _memo_verdict(m))
    true
end
check_valid_rate_equation(::Type{Bool}, m::Mechanisms.Mechanism) =
    _memo_verdict(m).fault === _NO_FAULT

"""Throw the guard's message (`_guard_message`) unless `m`'s verdict `v` holds no fault."""
function _check_verdict(m::Mechanisms.Mechanism, v)
    v.fault === _NO_FAULT || error(_guard_message(m))
    nothing
end

# ─── The guard: verdicts, networks, faces, classes 1 and 2 ─────────────

"""
What the guard finds wrong with a mechanism: nothing (`_NO_FAULT`); two constants of its
rate law share a name (`_SHARED_NAME`, `_shared_constant_names`); a cycle of steps does
not balance a metabolite (`_UNBALANCED_CYCLE`, `_unbalanced_cycle`); no cycle turns the
reaction over (`_NO_CYCLE`, class 1); a cycle of rapid-equilibrium steps turns it over
(`_RE_TURNOVER`, class 2); on a face (class 3), a second closed class (`_NONUNIQUE`), a
rapid-equilibrium segment with no lowest form (`_BOTTOMLESS`), or a law the limit test
rejects where the face's structure (`_face_stage`) stopped at a trap fed on the face
(`_TRAP`), a segment whose lowest forms split (`_SPLIT`), a law segment chain with several
closed classes (`_LAW_MULTI`) or an aggregated chain with several (`_MULTILEVEL`); in
admission's interior, each excluded step switched off as a reading says (`_Reading`), a
rapid-equilibrium step that a conformation excludes and the reading cannot switch off or
stops (`_FAST_STEP_PINNED`; the guard reports it only for a step inside one island, which
no reading switches off), steps that lose a direction and, each pointing uphill, close a
directed cycle or point into free E's component (`_ARC_CYCLE`), forms free E's component
does not reach along those steps (`_ISLAND`), or an active conformation that no longer
turns the reaction over (`_NO_ACTIVE_TURNOVER`); on a face, in admission's face condition,
a second closed class of the ε-network (`_EPS_NONUNIQUE`), a form that a fast step's
vanishing direction fills one way faster than it empties (`_FILLED`), a stranded basin
entered by, or holding, an alive step whose reverse is dead, or entered from outside it
and the law's forms (`_PUMPED`), or tagged steps that shorten the fed forms' cheapest
spanning in-tree where no potential absorbs their turnover (`_BRIDGE`); or, no reading the
guard tried being consistent and the active conformation turning the reaction over, more
readings than it searches (`_TOO_MANY_READINGS`). A reading is consistent when
admission's interior finds none of its first three faults, so the guard reports neither
`_ARC_CYCLE` nor `_ISLAND` (`_guard_with_reading`).
"""
@enum _Fault::UInt8 begin
    _NO_FAULT
    _NO_CYCLE
    _RE_TURNOVER
    _NONUNIQUE
    _BOTTOMLESS
    _TRAP
    _SPLIT
    _LAW_MULTI
    _MULTILEVEL
    _FAST_STEP_PINNED
    _ARC_CYCLE
    _ISLAND
    _EPS_NONUNIQUE
    _FILLED
    _PUMPED
    _BRIDGE
    _NO_ACTIVE_TURNOVER
    _TOO_MANY_READINGS
    _SHARED_NAME
    _UNBALANCED_CYCLE
end

"""The faults of admission the guard reports; `_admission_message` writes their text."""
const _ADMISSION_FAULTS = (_FAST_STEP_PINNED, _NO_ACTIVE_TURNOVER, _EPS_NONUNIQUE, _FILLED,
                           _PUMPED, _BRIDGE, _TOO_MANY_READINGS)

"""The faults under one consistent reading that another may avoid: the face condition's,
which read which transitions vanish. Every other fault holds under every consistent
reading: classes 1 and 3 and the active conformation's turnover read no reading, the
ε-network's closed classes read no rate, and class 2 reads the fast transitions, which are
every rapid-equilibrium step's under every consistent reading."""
const _FAULTS_A_READING_CHANGES = (_FILLED, _PUMPED, _BRIDGE)

"""The most ways of switching off the steps one conformation excludes that the guard
searches (`_consistent_columns`); when a conformation has more, it tries its built
readings alone."""
const _MAX_CHOICES_PER_CONFORMATION = 20_000

"""The most pairs of consistent ways, one per conformation, that the guard searches
(`_guard_with_reading`); when there are more, it tries its built readings alone."""
const _MAX_READING_PAIRS = 5_000

"""
The limit test's verdict on a face (class 3): not run (`_LIMIT_UNTESTED`); the fast-RE
limit of the steady state equals the law (`_LIMIT_EQUAL`), differs from it
(`_LIMIT_DIFFER`) or moves with the relative speeds of the rapid-equilibrium steps
(`_LIMIT_SPEED`); the law has no value there (`_LIMIT_LAW_NAN`).
"""
@enum _LimitVerdict::UInt8 begin
    _LIMIT_UNTESTED
    _LIMIT_EQUAL
    _LIMIT_DIFFER
    _LIMIT_SPEED
    _LIMIT_LAW_NAN
end

"""
The guard's verdict on a mechanism, an isbits value, so a loop over verdicts keeps its
types concrete: its `fault`, the limit test's verdict where it decided (`limit`), the face
(an index into `_faces`; 0 when the fault belongs to no face), the conformation it lies in
(1 active, 2 inactive; 0 for none, or when the forms it names lie in both) and the edge its
message leads with (0 none).
"""
struct _Verdict
    fault::_Fault
    limit::_LimitVerdict
    face::Int
    conformation::Int
    edge::Int
end

"""The verdict on a mechanism that passes every check."""
const _WELL_POSED = _Verdict(_NO_FAULT, _LIMIT_UNTESTED, 0, 0, 0)

"""
A network the guard reads: `ncopies` copies of a mechanism's `nforms` forms, node
`(c - 1) * nforms + f` for form `f` in copy `c` (`_node`), `present` marking the nodes it
holds. Transitions come in pairs: `2j - 1` runs an edge (`edge`; 0 for the flip at free
E) in its stored direction in copy `conformation` (0 for the flip), and `2j` runs it in
reverse; `forward` marks the first of each pair. `fast` marks a rapid-equilibrium
transition and `vanishing` one whose rate goes to zero in the ε-network. `change[n, t]`
counts the copies of metabolite `n` (an index into `m.metabolites`) that transition `t`
takes up less those it gives off; bit `n - 1` of `takes[t]` is set when it takes up
metabolite `n`.
"""
struct _Network
    nforms::Int
    ncopies::Int
    present::BitVector
    from::Vector{Int}
    to::Vector{Int}
    edge::Vector{Int}
    conformation::Vector{Int}
    forward::BitVector
    fast::BitVector
    vanishing::BitVector
    change::Matrix{Int}
    takes::Vector{UInt64}
end

_node(net::_Network, f::Int, c::Int) = (c - 1) * net.nforms + f

"""Each node's transitions of `net` that `keep` flags and that leave it, in transition
order."""
function _out_transitions(net::_Network, keep::BitVector)
    out = [Int[] for _ in 1:(net.nforms * net.ncopies)]
    for t in eachindex(net.from)
        keep[t] && push!(out[net.from[t]], t)
    end
    out
end

"""Each node's successors in `net` along the transitions `keep` flags, in transition
order; its predecessors when `backward`."""
function _successors(net::_Network, keep::BitVector; backward::Bool = false)
    a, b = backward ? (net.to, net.from) : (net.from, net.to)
    next = [Int[] for _ in 1:(net.nforms * net.ncopies)]
    for t in eachindex(a)
        keep[t] && push!(next[a[t]], b[t])
    end
    next
end

"""The nodes that the nodes `start` reach along the successor lists `next`, `start`
included."""
function _reached(next::Vector{Vector{Int}}, start::Vector{Int})
    seen = falses(length(next))
    seen[start] .= true
    todo = copy(start)
    while !isempty(todo)
        for v in next[pop!(todo)]
            seen[v] || (seen[v] = true; push!(todo, v))
        end
    end
    seen
end

"""
The network over the `present` nodes holding both directions of edge `k` of `m` in copy
`c` wherever `use[k, c]` holds, fast where `fast[k, c]` does, the stored direction
vanishing where `vanish_forward[k, c]` does and the reverse where `vanish_reverse[k, c]`
does; with two copies, the fast flip between the copies of free E, which takes up and
gives off nothing.
"""
function _network(m::Mechanisms.CatalyticScheme, present::BitVector, use::BitMatrix,
                  fast::BitMatrix, vanish_forward::BitMatrix, vanish_reverse::BitMatrix)
    nf, ncopies = length(m.enzyme_forms), size(use, 2)
    nt = 2 * count(use) + (ncopies == 2 ? 2 : 0)
    from, to, edge, conformation = zeros(Int, nt), zeros(Int, nt), zeros(Int, nt),
                                   zeros(Int, nt)
    forward, isfast, vanishing = falses(nt), falses(nt), falses(nt)
    change, takes = zeros(Int, length(m.metabolites), nt), zeros(UInt64, nt)
    t = 0
    for c in 1:ncopies, k in eachindex(m.steps)
        use[k, c] || continue
        a, b = (c - 1) * nf .+ (m.steps[k].from, m.steps[k].to)
        up, off = m.steps[k].taken_up, m.steps[k].given_off
        for (u, v, gain, loss, gone) in ((a, b, up, off, vanish_forward[k, c]),
                                         (b, a, off, up, vanish_reverse[k, c]))
            t += 1
            from[t], to[t], edge[t], conformation[t] = u, v, k, c
            forward[t], isfast[t], vanishing[t] = u == a, fast[k, c], gone
            for n in gain
                change[n, t] += 1
                takes[t] |= UInt64(1) << (n - 1)
            end
            for n in loss
                change[n, t] -= 1
            end
        end
    end
    if ncopies == 2
        from[t + 1], to[t + 1], forward[t + 1], isfast[t + 1] = 1, nf + 1, true, true
        from[t + 2], to[t + 2], isfast[t + 2] = nf + 1, 1, true
    end
    _Network(nf, ncopies, present, from, to, edge, conformation, forward, isfast,
             vanishing, change, takes)
end

"""
The law network L of `m`: a plain mechanism's forms and edges as they are; an allosteric
mechanism's two copies, copy `c` holding the steps of its conformation
(`steps_running_in_allosteric_conformations`) and the forms they touch, free E always,
joined by the flip at free E.
"""
function _law_network(am::Mechanisms.Mechanism)
    m = am.catalytic
    re = _re_flags(m)
    n, nf = length(m.steps), length(m.enzyme_forms)
    _is_allosteric(am) ||
        return _network(m, trues(nf), trues(n, 1), reshape(re, n, 1), falses(n, 1),
                        falses(n, 1))
    use, fast, present = falses(n, 2), falses(n, 2), falses(2 * nf)
    runs = Mechanisms.steps_running_in_allosteric_conformations(am)
    for (c, mask) in enumerate((runs.active, runs.inactive))
        use[:, c] .= mask
        fast[:, c] .= mask .& re
        shift = (c - 1) * nf
        present[shift + 1] = true
        for (k, e) in enumerate(m.steps)
            mask[k] && (present[shift + e.from] = present[shift + e.to] = true)
        end
    end
    _network(m, present, use, fast, falses(n, 2), falses(n, 2))
end

"""
The ε-network N of `mech`, its rapid-equilibrium flags overridden by `re`: every form and
edge in one copy per conformation of `reading`, each edge switched off as `reading` says
there. An edge that runs keeps its kind and both rates; one that loses a direction keeps
its kind while that direction vanishes and the other keeps its rate; any other is slow and
vanishes both ways (admission rejects a rapid-equilibrium one, `_FAST_STEP_PINNED`). A
plain mechanism's one conformation runs every step (`_built_reading`), so its ε-network is
its law network.
"""
function _eps_network(mech::Mechanisms.Mechanism, re::BitVector, reading::_Reading)
    m = mech.catalytic
    dims = size(reading)
    fast, vanish_forward, vanish_reverse = falses(dims), falses(dims), falses(dims)
    for c in axes(reading, 2), k in axes(reading, 1)
        o = reading[k, c]
        fast[k, c] = re[k] && o in (_RUNS, _LOSES_FORWARD, _LOSES_REVERSE)
        vanish_forward[k, c] = o !== _RUNS && o !== _LOSES_REVERSE
        vanish_reverse[k, c] = o !== _RUNS && o !== _LOSES_FORWARD
    end
    _network(m, trues(dims[2] * length(m.enzyme_forms)), trues(dims), fast,
             vanish_forward, vanish_reverse)
end

"""
The faces the guard checks of `m`, each a mask of the metabolite names it zeroes (bit
`n - 1` for metabolite `n` of `m.metabolites`), in order: the interior (0), every product,
every substrate, each substrate and then each product alone on a side of two or more; then
each of these with every modifier too. A regulatory-site ligand is zeroed only when it
also binds the catalytic site. No face zeroes a substrate together with a product.
"""
function _faces(m::Mechanisms.CatalyticScheme)
    substrates, products = m.substrates, m.products
    bit(n::Int) = UInt64(1) << (n - 1)
    side(ns) = reduce(|, (bit(n) for n in ns); init = UInt64(0))
    base = UInt64[0, side(products), side(substrates)]
    length(substrates) >= 2 && append!(base, bit.(substrates))
    length(products) >= 2 && append!(base, bit.(products))
    mods = side(n for n in eachindex(m.metabolites)
                if !(n in substrates || n in products))
    mods == 0 ? base : [base; base .| mods]
end

_zeroes(face::UInt64, n::Int) = (face >> (n - 1)) & 1 == 1

"""The metabolite names that face `face` of `m` zeroes (`_faces`), sorted."""
_face_names(m::Mechanisms.CatalyticScheme, face::UInt64) =
    sort!([x for (n, x) in enumerate(m.metabolites) if _zeroes(face, n)])

"""The tree transitions from node `x` up to its root in the forest `parent` (each node's
tree transition, 0 at a root)."""
function _tree_path(net::_Network, parent::Vector{Int}, x::Int)
    path = Int[]
    while parent[x] != 0
        t = parent[x]
        push!(path, t)
        x = net.from[t] == x ? net.to[t] : net.from[t]
    end
    path
end

"""The transitions of the cycle that the non-tree transition `t` between nodes `u` and
`v` closes in the forest `parent`, in cycle order from `u`."""
function _tree_cycle(net::_Network, parent::Vector{Int}, u::Int, v::Int, t::Int)
    up, down = _tree_path(net, parent, u), _tree_path(net, parent, v)
    while !isempty(up) && !isempty(down) && up[end] == down[end]
        pop!(up)
        pop!(down)
    end
    [up; reverse(down); t]
end

"""
    _fast_segments(net) -> (; segment, segments, offset, cycle)

The rapid-equilibrium segments of `net`, walked breadth-first along its fast transitions
from each present node no earlier walk reached, in node order: `segment[u]`, node `u`'s
(0 for a node `net` does not hold); `segments`, each one's nodes in walk order;
`offset[n, u]`, the copies of metabolite `n` (an index into `m.metabolites`) the walk
takes up from the
segment's first node to `u`; and `cycle`, the transitions of the first fast cycle that
takes up some name net, which the first transition whose offsets disagree closes with its
tree path, in cycle order.
When `cycle` is empty, each offset is a potential: it does not depend on the walk.
"""
function _fast_segments(net::_Network)
    nn = net.nforms * net.ncopies
    out = _out_transitions(net, net.fast)
    segment, segments, parent = zeros(Int, nn), Vector{Int}[], zeros(Int, nn)
    offset, disagreeing = zeros(Int, size(net.change, 1), nn), 0
    for root in 1:nn
        net.present[root] && segment[root] == 0 || continue
        members = [root]
        push!(segments, members)
        segment[root] = length(segments)
        head = 0
        while head < length(members)
            head += 1
            u = members[head]
            for t in out[u]
                v = net.to[t]
                if segment[v] == 0
                    segment[v], parent[v] = length(segments), t
                    push!(members, v)
                    for n in axes(offset, 1)
                        offset[n, v] = offset[n, u] + net.change[n, t]
                    end
                elseif disagreeing == 0 &&
                       any(n -> offset[n, v] != offset[n, u] + net.change[n, t],
                           axes(offset, 1))
                    disagreeing = t
                end
            end
        end
    end
    d = disagreeing
    cycle = d == 0 ? Int[] : _tree_cycle(net, parent, net.from[d], net.to[d], d)
    (; segment, segments, offset, cycle)
end

"""The transitions of the first cycle of fast transitions of `net` that takes up some
metabolite name net (class 2: the rate is infinite), in cycle order; empty when none does
(`_fast_segments`). Class 2 runs once, at the interior, since a fast cycle that runs on a
face runs at the interior too."""
_re_turnover_transitions(net::_Network) = _fast_segments(net).cycle

"""
    _unbalanced_cycle(m) -> (; steps, change)

A cycle of the steps of `m` that does not balance some metabolite, and its net change. A
cycle turns the reaction over as many times as it takes up the first substrate net, and
balances when it takes up each substrate and gives off each product that many times net,
and gives off each modifier as often as it takes it up. `steps` lists the cycle's steps in
cycle order. `change[x]` counts the copies of metabolite `x` (an index into
`m.metabolites`) that the cycle takes up less those it gives off, read in the direction
that takes up the first metabolite it changes, in stored order: the first substrate,
whenever the cycle changes it. Both are empty when every cycle balances. A step that
leaves out a metabolite, writes one on the wrong side or names the wrong form makes such a
cycle, and the joint solve has no law for it (`_fundamental_cycles`). The search borrows
class 2's walk (`_fast_segments`) on a network of every step, each marked fast.
"""
function _unbalanced_cycle(m::Mechanisms.CatalyticScheme)
    n, nf = length(m.steps), length(m.enzyme_forms)
    net = _network(m, trues(nf), trues(n, 1), trues(n, 1), falses(n, 1), falses(n, 1))
    # A transition turns the reaction over as many times as it takes up the first
    # substrate. With the reaction's uptake over those turns taken off each reactant's
    # row, a balanced cycle takes up no metabolite net, so the walk of `_fast_segments`,
    # all of whose transitions are fast here, finds a cycle that does not balance. `turns`
    # is a copy, as the loop zeroes the first substrate's own row.
    turns = net.change[first(m.substrates), :]
    for x in m.substrates
        net.change[x, :] .-= turns
    end
    for x in m.products
        net.change[x, :] .+= turns
    end
    cycle = _fast_segments(net).cycle
    isempty(cycle) && return (; steps = Int[], change = Int[])
    # The cycle is a closed walk from the start of its last transition. Count each step's
    # metabolites as written, in the direction the walk runs the step.
    change, at = zeros(Int, length(m.metabolites)), net.from[last(cycle)]
    for t in cycle
        along = net.from[t] == at
        at = along ? net.to[t] : net.from[t]
        dir = along == net.forward[t] ? 1 : -1
        for x in m.steps[net.edge[t]].taken_up
            change[x] += dir
        end
        for x in m.steps[net.edge[t]].given_off
            change[x] -= dir
        end
    end
    change[findfirst(!=(0), change)::Int] < 0 && (change .= .-change)
    (; steps = net.edge[cycle], change)
end

"""The first fault `check(i)` returns over the indices `i` of `faces`, in order;
`_WELL_POSED` when every face passes."""
function _first_face_fault(check, faces::Vector{UInt64})
    for i in eachindex(faces)
        v = check(i)
        v.fault === _NO_FAULT || return v
    end
    _WELL_POSED
end

"""
    _guard(m) -> _Verdict

The well-posedness guard on `m`, plain or allosteric: its verdict
(`_guard_with_reading`), `_WELL_POSED` when every check passes.
"""
_guard(mech::Mechanisms.Mechanism) = _guard_with_reading(mech).verdict

"""
    _guard_with_reading(m) -> (; verdict, reading)

The well-posedness guard on `m`, plain or allosteric, and the reading of how its
conformations switch off its excluded steps (`_Reading`) that its `verdict` was judged
under. First that each constant has a name of its own (`_shared_constant_names`), then
that every cycle balances each metabolite (`_unbalanced_cycle`), then class 1 (a cycle
turns the reaction over, `Mechanisms.has_turnover_cycle`); then, under a reading
(`_verdict_under`), class 2, admission's interior, and on every face class 3 and
admission's face condition. A plain mechanism's one reading runs every step.

An allosteric mechanism's verdict follows the energy-level rule: each conformation cuts off
every island but free E's, and every island it cuts off lies above free E's in energy and
holds none of that conformation's enzyme. A reading is consistent with the rule when
admission's interior finds none of `_interior_verdict`'s faults (`_FAST_STEP_PINNED`,
`_ARC_CYCLE`, `_ISLAND`) under it; the interior's last check, `_NO_ACTIVE_TURNOVER`, reads
no reading. The guard admits the mechanism under the first consistent reading that passes
every check: its three built readings first (`_built_reading`), then every other pair of
consistent columns (`_consistent_columns`), conformation 1's outer, each in its order.
When a conformation has more than `_MAX_CHOICES_PER_CONFORMATION` ways, or the pairs
number more than `_MAX_READING_PAIRS`, it searches none of them, so its verdict depends on
no order. A rejection reports, in order: a fast cycle that turns the reaction over under
the first built reading (`_RE_TURNOVER`), fast under every consistent reading too; an
excluded rapid-equilibrium step inside one island (`_FAST_STEP_PINNED`), which no reading
switches off; the verdict of the first consistent reading tried; or, when it tried none,
an active conformation that turns nothing over (`_NO_ACTIVE_TURNOVER`), else that there
are too many readings to search (`_TOO_MANY_READINGS`).
"""
function _guard_with_reading(am::Mechanisms.Mechanism)
    m = am.catalytic
    first_reading = _built_reading(am)
    isempty(_shared_constant_names(am)) ||
        return (; verdict = _Verdict(_SHARED_NAME, _LIMIT_UNTESTED, 0, 0, 0),
                reading = first_reading)
    isempty(_unbalanced_cycle(m).steps) ||
        return (; verdict = _Verdict(_UNBALANCED_CYCLE, _LIMIT_UNTESTED, 0, 0, 0),
                reading = first_reading)
    Mechanisms.has_turnover_cycle(m, trues(length(m.steps))) ||
        return (; verdict = _Verdict(_NO_CYCLE, _LIMIT_UNTESTED, 0, 0, 0),
                reading = first_reading)
    # A plain mechanism excludes no step, so its one reading runs every step.
    _is_allosteric(am) ||
        return (; verdict = _verdict_under(am, first_reading), reading = first_reading)
    # The built readings, in order. The first consistent reading tried, and its verdict,
    # are the ones a rejection reports.
    tried = _Reading[]
    reported, verdict = _Reading(undef, 0, 0), _WELL_POSED
    for r in instances(_SlowUptakeReading)
        reading = _built_reading(am, r)
        reading in tried && continue
        v = _verdict_under(am, reading)
        v === _WELL_POSED && return (; verdict = v, reading)
        # The built readings keep fast every rapid-equilibrium step but those they cannot
        # switch off, so a cycle fast under them is fast under every consistent reading.
        v.fault === _RE_TURNOVER && return (; verdict = v, reading)
        push!(tried, reading)
        # Class 2 passed, so `v` is admission's interior verdict whenever the interior
        # fails: the reading is consistent exactly when `v` holds none of its first three
        # faults.
        if isempty(reported) && v.fault ∉ (_FAST_STEP_PINNED, _ARC_CYCLE, _ISLAND)
            reported, verdict = reading, v
            verdict.fault in _FAULTS_A_READING_CHANGES ||
                return (; verdict, reading = reported)
        end
    end
    active = _consistent_columns(am, 1)
    active.pinned == 0 ||
        return (; verdict = _Verdict(_FAST_STEP_PINNED, _LIMIT_UNTESTED, 0, 1,
                                     active.pinned), reading = first_reading)
    inactive = _consistent_columns(am, 2)
    inactive.pinned == 0 ||
        return (; verdict = _Verdict(_FAST_STEP_PINNED, _LIMIT_UNTESTED, 0, 2,
                                     inactive.pinned), reading = first_reading)
    if active.too_many || inactive.too_many ||
       length(active.columns) * length(inactive.columns) > _MAX_READING_PAIRS
        isempty(reported) || return (; verdict, reading = reported)
        # No reading lets an active conformation without turnover catalyze, so that fault
        # is the one to report: the user can fix it, while too many readings to search
        # would only say the guard did not look.
        v = _active_turnover_verdict(am)
        v === _WELL_POSED || return (; verdict = v, reading = first_reading)
        c = active.too_many ? 1 : inactive.too_many ? 2 : 0
        return (; verdict = _Verdict(_TOO_MANY_READINGS, _LIMIT_UNTESTED, 0, c, 0),
                reading = first_reading)
    end
    # With no built reading consistent, the search's first reading is the one to report.
    # Each conformation has a consistent column: rank its islands by their distance from
    # free E's, ties by label, and let each excluded step between two islands lose toward
    # the later one.
    if isempty(reported)
        reported = hcat(active.columns[1], inactive.columns[1])
        verdict = _verdict_under(am, reported)
        verdict === _WELL_POSED && return (; verdict, reading = reported)
        verdict.fault in _FAULTS_A_READING_CHANGES || return (; verdict, reading = reported)
        push!(tried, reported)
    end
    # Every consistent reading passes class 2 and the active conformation's turnover, as
    # the reported one does. Class 3 reads no reading: a face it rejects rejects them all.
    law, faces = _law_network(am), _faces(m)
    _first_face_fault(i -> _check_face(am, law, faces, i), faces) === _WELL_POSED ||
        return (; verdict, reading = reported)
    re = _re_flags(m)
    pair = similar(first_reading)
    for column1 in active.columns, column2 in inactive.columns
        pair[:, 1] .= column1
        pair[:, 2] .= column2
        pair in tried && continue
        eps = _eps_network(am, re, pair)
        _first_face_fault(i -> _admission_face(am, eps, law, faces, i), faces) ===
            _WELL_POSED && return (; verdict = _WELL_POSED, reading = pair)
    end
    (; verdict, reading = reported)
end

"""
    _verdict_under(m, reading) -> _Verdict

The guard's verdict on `m` but class 1, its excluded steps switched off as `reading` says:
class 2 on the ε-network of that reading (`_re_turnover_verdict`), admission's interior
(`_admission_interior`), then class 3 and admission's face condition on every face
(`_faces_verdict`), the first fault winning; `_WELL_POSED` when every check passes.
"""
function _verdict_under(mech::Mechanisms.Mechanism, reading::_Reading)
    eps = _eps_network(mech, _re_flags(mech.catalytic), reading)
    v = _re_turnover_verdict(eps)
    v === _WELL_POSED || return v
    v = _admission_interior(mech, reading)
    v === _WELL_POSED || return v
    _faces_verdict(mech, eps)
end

"""Class 2's verdict on the ε-network `net`: `_RE_TURNOVER` in the conformation of the
first fast cycle that turns the reaction over, led by its first edge; else
`_WELL_POSED`."""
function _re_turnover_verdict(net::_Network)
    cycle = _re_turnover_transitions(net)
    isempty(cycle) && return _WELL_POSED
    _Verdict(_RE_TURNOVER, _LIMIT_UNTESTED, 0, net.conformation[cycle[1]],
             net.edge[cycle[1]])
end

# ─── Class 3: the structure of every face ──────────────────────────────

"""Each node's strongly connected component in the digraph with successor lists `next`,
searched from the nodes `nodes` marks, numbered in the order Tarjan's algorithm closes
them; 0 for a node no search reaches."""
function _strong_components(next::Vector{Vector{Int}}, nodes::BitVector)
    n = length(next)
    index, low, comp = zeros(Int, n), zeros(Int, n), zeros(Int, n)
    onstack, stack, work = falses(n), Int[], Tuple{Int, Int}[]
    counter = ncomp = 0
    for s in 1:n
        nodes[s] && index[s] == 0 || continue
        counter += 1
        index[s] = low[s] = counter
        push!(stack, s)
        onstack[s] = true
        push!(work, (s, 1))
        while !isempty(work)
            u, i = work[end]
            if i <= length(next[u])
                work[end] = (u, i + 1)
                v = next[u][i]
                if index[v] == 0
                    counter += 1
                    index[v] = low[v] = counter
                    push!(stack, v)
                    onstack[v] = true
                    push!(work, (v, 1))
                elseif onstack[v]
                    low[u] = min(low[u], index[v])
                end
            else
                pop!(work)
                isempty(work) || (low[work[end][1]] = min(low[work[end][1]], low[u]))
                if low[u] == index[u]
                    ncomp += 1
                    while true
                        v = pop!(stack)
                        onstack[v] = false
                        comp[v] = ncomp
                        v == u && break
                    end
                end
            end
        end
    end
    comp
end

"""The closed classes (bottom strongly connected components, `_strong_components`) of the
digraph with successor lists `next`, over the nodes `nodes` marks, each as its sorted
nodes; a node with no successor is a class of its own."""
function _bottom_sccs(next::Vector{Vector{Int}}, nodes::BitVector)
    comp = _strong_components(next, nodes)
    ncomp = maximum(comp; init = 0)
    leaves = falses(ncomp)
    for u in eachindex(next), v in next[u]
        nodes[u] && comp[u] != comp[v] && (leaves[comp[u]] = true)
    end
    classes = [Int[] for _ in 1:ncomp]
    for u in eachindex(next)
        nodes[u] && push!(classes[comp[u]], u)
    end
    [classes[k] for k in 1:ncomp if !leaves[k]]
end

"""The closed classes of the transitions of `net` that `keep` marks, over its present
nodes."""
_closed_classes(net::_Network, keep::BitVector) =
    _bottom_sccs(_successors(net, keep), net.present)

"""The fast classes (`class_of`) that the fast transitions (`fast_next`) reach from node
`w`, its own included."""
function _classes_reached(fast_next::Vector{Vector{Int}}, class_of::Vector{Int}, w::Int)
    seen = falses(length(fast_next))
    seen[w] = true
    todo, found = [w], Int[]
    while !isempty(todo)
        x = pop!(todo)
        class_of[x] > 0 && !(class_of[x] in found) && push!(found, class_of[x])
        for y in fast_next[x]
            seen[y] || (seen[y] = true; push!(todo, y))
        end
    end
    found
end

"""
    _aggregated_chain(net, keep) -> (; fast_classes, class_of, next)

The chain over the fast classes of the transitions of `net` that `keep` flags:
`fast_classes`, the closed classes of the fast ones among them (`_closed_classes`; a node
with no fast way out is a class of its own); `class_of`, each node's class (0 for a node
that the fast transitions carry into another class); and `next`, for each class the other
classes that a slow transition out of it reaches, followed along the fast ones
(`_classes_reached`). A slow transition out of a node in no fast class does not count: its
share vanishes as the fast steps speed up.
"""
function _aggregated_chain(net::_Network, keep::BitVector)
    fast_keep = keep .& net.fast
    fast_classes = _closed_classes(net, fast_keep)
    class_of = zeros(Int, net.nforms * net.ncopies)
    for (k, C) in enumerate(fast_classes), u in C
        class_of[u] = k
    end
    fast_next = _successors(net, fast_keep)
    next = [Int[] for _ in fast_classes]
    for t in eachindex(net.from)
        q = class_of[net.from[t]]
        keep[t] && !net.fast[t] && q > 0 || continue
        for k in _classes_reached(fast_next, class_of, net.to[t])
            k == q || push!(next[q], k)
        end
    end
    (; fast_classes, class_of, next)
end

"""
    _face_structure(net, face) -> NamedTuple

The structure that `_face_stage` reads on face `face` of the law network
`net`: `alive`, the transitions that take up no zeroed name; `closed`, the closed classes
of the alive ones; the rapid-equilibrium segments (`segment` per node, `segments`
breadth-first) and each node's uptake of every name from its segment's first node
(`offset`, `_fast_segments`); `bottomless`, the first segment with no lowest node (0 when
none), a lowest node being one at the segment's componentwise least uptake of the zeroed
names, whose weight stays positive on the face; the fast classes of the alive transitions
(`_aggregated_chain`), `class_of` each node's (0 for a node that drains elsewhere) and
`trap`, the classes of no lowest node; `aggregated`, the closed classes of their
aggregated chain; and `law_classes`, the closed classes of the law's segment chain, its
alive slow transitions out of lowest nodes.
"""
function _face_structure(net::_Network, face::UInt64)
    nn, nt = net.nforms * net.ncopies, length(net.from)
    alive = BitVector([net.takes[t] & face == 0 for t in 1:nt])
    closed = _closed_classes(net, alive)
    zeroed = [n for n in axes(net.change, 1) if _zeroes(face, n)]
    (; segment, segments, offset) = _fast_segments(net)
    lowest, bottomless = falses(nn), 0
    for (g, members) in enumerate(segments)
        least = [minimum(offset[n, u] for u in members) for n in zeroed]
        for u in members
            lowest[u] = all(j -> offset[zeroed[j], u] == least[j], eachindex(zeroed))
        end
        bottomless == 0 && !any(u -> lowest[u], members) && (bottomless = g)
    end
    (; fast_classes, class_of, next) = _aggregated_chain(net, alive)
    all(C -> all(u -> lowest[u] == lowest[C[1]], C), fast_classes) ||
        error("_face_structure: a fast class holds lowest and other forms")
    trap = BitVector([!lowest[C[1]] for C in fast_classes])
    law_next = [Int[] for _ in segments]
    for t in 1:nt
        u, v = net.from[t], net.to[t]
        alive[t] && !net.fast[t] && lowest[u] && segment[u] != segment[v] &&
            push!(law_next[segment[u]], segment[v])
    end
    aggregated = _bottom_sccs(next, trues(length(fast_classes)))
    law_classes = _bottom_sccs(law_next, trues(length(segments)))
    (; alive, closed, segment, segments, offset, bottomless, fast_classes, class_of, trap,
     aggregated, law_classes)
end

"""
What the structure `s` of a face (`_face_structure`) settles: a second closed class
(`_NONUNIQUE`); a segment with no lowest form (`_BOTTOMLESS`); `_NO_FAULT` when every
segment holds one fast class and the law's segment chain one closed class; else the stage
at which the limit test decides: the law's chain has several closed classes
(`_LAW_MULTI`), the aggregated chain has (`_MULTILEVEL`), its closed class holds a trap
(`_TRAP`), or none of these (`_SPLIT`, where a segment it meets may split).
"""
function _face_stage(s)
    length(s.closed) > 1 && return _NONUNIQUE
    s.bottomless > 0 && return _BOTTOMLESS
    classes = zeros(Int, length(s.segments))
    for C in s.fast_classes
        classes[s.segment[C[1]]] += 1
    end
    one_law = length(s.law_classes) == 1
    one_law && all(==(1), classes) && return _NO_FAULT
    !one_law ? _LAW_MULTI : length(s.aggregated) > 1 ? _MULTILEVEL :
        any(k -> s.trap[k], s.aggregated[1]) ? _TRAP : _SPLIT
end

"""The copy (1 active, 2 inactive) of node `u` of `net`."""
_copy_of(net::_Network, u::Int) = (u - 1) ÷ net.nforms + 1
_form_of(net::_Network, u::Int) = (u - 1) % net.nforms + 1
"""The copy of `net` that every node of `nodes` lies in; 0 when they lie in both."""
_common_copy(net::_Network, nodes::Vector{Int}) =
    (c = _copy_of(net, nodes[1]); all(u -> _copy_of(net, u) == c, nodes) ? c : 0)

"""The closed class of face structure `s` that holds no active free E (node 1), or its
first when each holds none."""
_away_class(s) = s.closed[something(findfirst(C -> !insorted(1, C), s.closed), 1)]

"""The edges whose transitions leave the nodes `C` and take up a zeroed name on the face
of `s`, in transition order, each once."""
_exits(net::_Network, s, C::Vector{Int}) =
    unique!([net.edge[t] for t in eachindex(net.from)
             if insorted(net.from[t], C) && !insorted(net.to[t], C) && !s.alive[t] &&
                net.edge[t] != 0])

"""
The trap class 3 found on a face (stage `_TRAP`) and the rapid-equilibrium edge that feeds
it or joins it to its segment, as `(trap, edge, joined)`, the trap an index into
`s.fast_classes`. For an alive slow transition out of a class of the aggregated chain's
closed class, the edge is the last fast transition of an alive fast path from its target
into a trap of that class (`joined` false); when the slow transition lands in the trap, or
no transition feeds one (the closed class is a trap alone), it is the rapid-equilibrium
edge that joins the trap to the rest of its segment, whose way out takes up a zeroed name
(`joined` true).
"""
function _trap_feeder(net::_Network, s)
    nn = net.nforms * net.ncopies
    in_R = falses(length(s.fast_classes))
    in_R[s.aggregated[1]] .= true
    trap_of(u) = (q = s.class_of[u]; q > 0 && in_R[q] && s.trap[q] ? q : 0)
    function joining(q)
        for t in eachindex(net.from)
            net.fast[t] && s.class_of[net.from[t]] == q && s.class_of[net.to[t]] != q &&
                return (q, net.edge[t], true)
        end
        (q, 0, true)
    end
    fast_out = _out_transitions(net, s.alive .& net.fast)
    for t in eachindex(net.from)
        k = s.class_of[net.from[t]]
        s.alive[t] && !net.fast[t] && k > 0 && in_R[k] || continue
        v = net.to[t]
        landed = trap_of(v)
        landed > 0 && return joining(landed)
        seen = falses(nn)
        seen[v] = true
        todo = [v]
        while !isempty(todo)
            x = pop!(todo)
            for t2 in fast_out[x]
                y = net.to[t2]
                fed = trap_of(y)
                fed > 0 && s.class_of[x] != fed && return (fed, net.edge[t2], false)
                seen[y] || (seen[y] = true; push!(todo, y))
            end
        end
    end
    joining(s.aggregated[1][findfirst(k -> s.trap[k], s.aggregated[1])])
end

"""
The split class 3 found on a face (stage `_SPLIT`): of the first rapid-equilibrium segment
that the aggregated chain's closed class meets and that holds more than one lowest fast
class, those classes (`parts`), and the rapid-equilibrium edges whose transitions leave
one of them and take up a zeroed name (`joins`), the only edges that could join them.
"""
function _split_parts(net::_Network, s)
    lowest_classes = [Int[] for _ in s.segments]
    for (k, C) in enumerate(s.fast_classes)
        s.trap[k] || push!(lowest_classes[s.segment[C[1]]], k)
    end
    segment_of(k) = s.segment[s.fast_classes[k][1]]
    R = s.aggregated[1]
    first_split = findfirst(k -> length(lowest_classes[segment_of(k)]) != 1, R)
    segment = segment_of(R[something(first_split, 1)])
    parts = [s.fast_classes[k] for k in lowest_classes[segment]]
    joins = unique!([net.edge[t] for t in eachindex(net.from)
                     if net.fast[t] && !s.alive[t] && net.edge[t] != 0 &&
                        any(C -> insorted(net.from[t], C), parts)])
    (; parts, joins)
end

"""
    _check_face(mech, net, faces, i) -> _Verdict

Class 3 on face `faces[i]` of `mech`, whose law network is `net`. The face's structure
(`_face_stage`) rejects a second closed class or a bottomless segment and accepts a face
where it fixes the limit; otherwise the limit test (`_face_limit_check`, seeded by the
content hash and the face) decides, and a rejection carries the stage it stopped at, the
edge that feeds the trap or joins it to its segment (`_trap_feeder`) and the split's first
joining edge. A bottomless segment, a split or a trap carries the copy its forms share
(`_common_copy`). A face the limit test accepts reports `_LIMIT_EQUAL`.
"""
function _check_face(mech::Mechanisms.Mechanism, net::_Network, faces::Vector{UInt64},
                     i::Int)
    s = _face_structure(net, faces[i])
    stage = _face_stage(s)
    stage === _NO_FAULT && return _Verdict(_NO_FAULT, _LIMIT_UNTESTED, i, 0, 0)
    if stage === _NONUNIQUE
        C = _away_class(s)
        exits = _exits(net, s, C)
        return _Verdict(_NONUNIQUE, _LIMIT_UNTESTED, i, _copy_of(net, C[1]),
                        isempty(exits) ? 0 : exits[1])
    end
    stage === _BOTTOMLESS &&
        return _Verdict(_BOTTOMLESS, _LIMIT_UNTESTED, i,
                        _common_copy(net, s.segments[s.bottomless]), 0)
    limit = _face_limit_check(mech, net, s, faces[i], hash(i, hash(mech)))
    limit === _LIMIT_EQUAL && return _Verdict(_NO_FAULT, _LIMIT_EQUAL, i, 0, 0)
    if stage === _TRAP
        trap, edge = _trap_feeder(net, s)
        return _Verdict(_TRAP, limit, i, _common_copy(net, s.fast_classes[trap]), edge)
    elseif stage === _SPLIT
        split = _split_parts(net, s)
        return _Verdict(_SPLIT, limit, i, _common_copy(net, reduce(vcat, split.parts)),
                        isempty(split.joins) ? 0 : split.joins[1])
    end
    _Verdict(stage, limit, i, 0, 0)
end

"""
    _faces_verdict(mech, eps) -> _Verdict

The first fault of `mech` face by face, in the order of `_faces`, from its ε-network `eps`
(`_eps_network`): class 3 on its law network (`_check_face`), which for a plain mechanism
is `eps`; then, for an allosteric mechanism, admission's face condition on `eps`
(`_admission_face`). `_WELL_POSED` when every face passes.
"""
function _faces_verdict(am::Mechanisms.Mechanism, eps::_Network)
    faces = _faces(am.catalytic)
    _is_allosteric(am) ||
        return _first_face_fault(i -> _check_face(am, eps, faces, i), faces)
    law = _law_network(am)
    _first_face_fault(faces) do i
        v = _check_face(am, law, faces, i)
        v.fault === _NO_FAULT ? _admission_face(am, eps, law, faces, i) : v
    end
end

# ─── Class 3: the limit test mod 2^61 − 1 ───────────────────────────────

"""The prime 2^61 − 1: the limit test computes every rate, weight and limit modulo it."""
const _P61 = (Int64(1) << 61) - 1
"""A limit that grows without bound, outside the residues 0…p − 1."""
const _INFINITE = Int64(-1)
"""A law that has no value on a face, outside the residues 0…p − 1."""
const _UNDEFINED = Int64(-2)

_mulp(a::Int64, b::Int64) = Int64(mod(widemul(a, b), _P61))
_addp(a::Int64, b::Int64) = (s = a + b; s >= _P61 ? s - _P61 : s)
_subp(a::Int64, b::Int64) = (s = a - b; s < 0 ? s + _P61 : s)
_modp(x::Integer) = Int64(mod(x, _P61))

"""`a^k` mod p; a negative `k` raises the inverse."""
function _powp(a::Int64, k::Integer)
    k < 0 && return _powp(_invp(a), -k)
    r = Int64(1)
    while k > 0
        isodd(k) && (r = _mulp(r, a))
        a = _mulp(a, a)
        k >>= 1
    end
    r
end

_invp(a::Int64) = (a == 0 && throw(DivideError()); _powp(a, _P61 - 2))

"""A deterministic stream of residues in 2…p − 2 (splitmix64), seeded by a mechanism's
content hash and a face, so every process draws the same points."""
mutable struct _ResidueStream
    state::UInt64
end

function _residue!(r::_ResidueStream)
    r.state += 0x9e3779b97f4a7c15
    z = r.state
    z = (z ⊻ (z >> 30)) * 0xbf58476d1ce4e5b9
    z = (z ⊻ (z >> 27)) * 0x94d049bb133111eb
    Int64(2 + (z ⊻ (z >> 31)) % UInt64(_P61 - 3))
end

"""`(A \\ B, det(A))` mod p for square `A`, by Gauss–Jordan elimination: the determinant is
the product of the pivots, its sign flipped by each row swap; `A \\ B` is an empty matrix,
and the determinant 0, when `A` is singular."""
function _solvep(A::Matrix{Int64}, B::Matrix{Int64})
    n = size(A, 1)
    M = hcat(A, B)
    d = Int64(1)
    for j in 1:n
        k = findfirst(i -> M[i, j] != 0, j:n)
        k === nothing && return Matrix{Int64}(undef, 0, 0), Int64(0)
        k += j - 1
        if k != j
            M[[j, k], :] = M[[k, j], :]
            d = _subp(Int64(0), d)
        end
        d = _mulp(d, M[j, j])
        pivot = _invp(M[j, j])
        for c in j:size(M, 2)
            M[j, c] = _mulp(M[j, c], pivot)
        end
        for i in 1:n
            (i == j || M[i, j] == 0) && continue
            f = M[i, j]
            for c in j:size(M, 2)
                M[i, c] = _subp(M[i, c], _mulp(f, M[j, c]))
            end
        end
    end
    M[:, (n + 1):end], d
end

"""The stationary distribution mod p of the chain with off-diagonal rates `R[i, j]`
(`i → j`) and one closed class covering every state, and the determinant of the system that
gives it (the generator's transpose, its last row ones); empty, and 0, when singular."""
function _stationaryp(R::Matrix{Int64})
    n = size(R, 1)
    n == 1 && return Int64[1], Int64(1)
    A = zeros(Int64, n, n)
    for i in 1:n, j in 1:n
        i == j && continue
        A[j, i] = _addp(A[j, i], R[i, j])
        A[i, i] = _subp(A[i, i], R[i, j])
    end
    A[n, :] .= 1
    b = zeros(Int64, n, 1)
    b[n] = 1
    x, d = _solvep(A, b)
    vec(x), d
end

# Polynomials mod p: coefficient vectors, constant term first.
_ptrimp!(a::Vector{Int64}) = (while length(a) > 1 && a[end] == 0; pop!(a); end; a)
_degp(a::Vector{Int64}) =
    (b = _ptrimp!(copy(a)); length(b) == 1 && b[1] == 0 ? -1 : length(b) - 1)
_paddp(a::Vector{Int64}, b::Vector{Int64}) =
    _ptrimp!([_addp(get(a, i, Int64(0)), get(b, i, Int64(0)))
              for i in 1:max(length(a), length(b))])
_psubp(a::Vector{Int64}, b::Vector{Int64}) =
    _ptrimp!([_subp(get(a, i, Int64(0)), get(b, i, Int64(0)))
              for i in 1:max(length(a), length(b))])
function _pmulp(a::Vector{Int64}, b::Vector{Int64})
    c = zeros(Int64, length(a) + length(b) - 1)
    for i in eachindex(a), j in eachindex(b)
        c[i + j - 1] = _addp(c[i + j - 1], _mulp(a[i], b[j]))
    end
    _ptrimp!(c)
end
function _pdivremp(a::Vector{Int64}, b::Vector{Int64})
    a, b = _ptrimp!(copy(a)), _ptrimp!(copy(b))
    db = _degp(b)
    db < 0 && throw(DivideError())
    q = zeros(Int64, max(_degp(a) - db + 1, 1))
    lead = _invp(b[end])
    while _degp(a) >= db
        k = _degp(a) - db
        c = _mulp(a[end], lead)
        q[k + 1] = c
        for i in 0:db
            a[i + k + 1] = _subp(a[i + k + 1], _mulp(c, b[i + 1]))
        end
        _ptrimp!(a)
    end
    _ptrimp!(q), a
end
_pevalp(a::Vector{Int64}, x::Int64) =
    foldr((c, acc) -> _addp(c, _mulp(acc, x)), a; init = Int64(0))
"""The lowest-order term of `a` as `(order, coefficient)`; `(typemax(Int), 0)` for 0."""
_plowp(a::Vector{Int64}) =
    (k = findfirst(!=(0), a); k === nothing ? (typemax(Int), Int64(0)) : (k - 1, a[k]))

"""Newton interpolation mod p through the points `(x, y)`."""
function _interpp(x::Vector{Int64}, y::Vector{Int64})
    n = length(x)
    c = copy(y)
    for j in 2:n, i in n:-1:j
        c[i] = _mulp(_subp(c[i], c[i - 1]), _invp(_subp(x[i], x[i - j + 1])))
    end
    p = Int64[c[n]]
    for i in (n - 1):-1:1
        p = _paddp(_pmulp(p, Int64[_subp(Int64(0), x[i]), 1]), Int64[c[i]])
    end
    p
end

"""`(N, D)` with degrees at most `dmax` and `N/D` through the points: the extended
Euclidean algorithm on the interpolant (it needs at least `2 dmax + 1` points)."""
function _ratrecp(x::Vector{Int64}, y::Vector{Int64}, dmax::Int)
    M = Int64[1]
    for xi in x
        M = _pmulp(M, Int64[_subp(Int64(0), xi), 1])
    end
    r0, r1 = M, _interpp(x, y)
    t0, t1 = Int64[0], Int64[1]
    while _degp(r1) > dmax
        q, r = _pdivremp(r0, r1)
        r0, r1 = r1, r
        t0, t1 = t1, _psubp(t0, _pmulp(q, t1))
    end
    r1, t1
end

"""The rational function of one variable that `f` evaluates mod p, reconstructed from
random points with a degree bound doubled until two more points agree; empty vectors past
degree 1024."""
function _reconstructp(f, rng::_ResidueStream)
    d = 8
    while d <= 1024
        xs = [_residue!(rng) for _ in 1:(2d + 3)]
        ys = [f(x) for x in xs]
        N, D = _ratrecp(xs[1:(end - 2)], ys[1:(end - 2)], d)
        all(k -> (dv = _pevalp(D, xs[k]); dv != 0 &&
                  _mulp(_pevalp(N, xs[k]), _invp(dv)) == ys[k]),
            (length(xs) - 1):length(xs)) && return N, D
        d *= 2
    end
    Int64[], Int64[]
end

"""
A random point mod p of the model of `mech`'s catalytic graph `m` for `ncopies`
conformations on face `face`: log-free energies of every form in each conformation and of
every metabolite name, an integer null-space basis of the kinetic-group ties (members
equal in each conformation that keeps them, an `:EqualAI` member equal across) raised to
random residues; from them each edge's equilibrium constant `K[k, c]` and `L`; a
steady-state speed per kinetic group and conformation (one for both of an `:EqualAI`
group), a rapid-equilibrium speed per edge and conformation, the flip's speed;
concentrations `conc` per name, zero on the face.
"""
function _model_point(mech::Mechanisms.Mechanism, ncopies::Int, face::UInt64,
                      rng::_ResidueStream)
    m = mech.catalytic
    nf, nnames = length(m.enzyme_forms), length(m.metabolites)
    nv = ncopies * nf + nnames
    excluded = [_excluded_steps(mech, c) for c in 1:ncopies]
    equal_ai(g) = _is_allosteric(mech) && mech.allosteric.allosteric_states[g] === :EqualAI
    function drop!(row, k, c, dir)
        a, b = m.steps[k].from, m.steps[k].to
        row[(c - 1) * nf + a] += dir
        row[(c - 1) * nf + b] -= dir
        for n in m.steps[k].taken_up
            row[ncopies * nf + n] += dir
        end
        for n in m.steps[k].given_off
            row[ncopies * nf + n] -= dir
        end
        row
    end
    ties = Vector{Int}[]
    groups = _kinetic_groups(m)
    for (g, members) in enumerate(groups), c in 1:ncopies
        lead = first(members)
        excluded[c][lead] && continue
        for k in members[2:end]
            push!(ties, drop!(drop!(zeros(Int, nv), k, c, 1), lead, c, -1))
        end
        if c == 2 && equal_ai(g)
            for k in members
                push!(ties, drop!(drop!(zeros(Int, nv), k, 2, 1), k, 1, -1))
            end
        end
    end
    B = _integer_nullspace(isempty(ties) ? zeros(Int, 0, nv) :
                           permutedims(reduce(hcat, ties)))
    t = [_residue!(rng) for _ in axes(B, 2)]
    x = [reduce(_mulp, (_powp(t[j], B[i, j]) for j in axes(B, 2) if B[i, j] != 0);
                init = Int64(1)) for i in 1:nv]
    K = ones(Int64, length(m.steps), ncopies)
    xnames = view(x, (ncopies * nf + 1):nv)
    for k in eachindex(m.steps), c in 1:ncopies
        a, b = m.steps[k].from, m.steps[k].to
        num = _mul_held(x[(c - 1) * nf + a], m.steps[k].taken_up, xnames)
        den = _mul_held(x[(c - 1) * nf + b], m.steps[k].given_off, xnames)
        K[k, c] = _mulp(num, _invp(den))
    end
    ngroups = length(groups)
    ss = [_residue!(rng) for _ in 1:ngroups, _ in 1:ncopies]
    for g in 1:ngroups
        ncopies == 2 && equal_ai(g) && (ss[g, 2] = ss[g, 1])
    end
    re = [_residue!(rng) for _ in eachindex(m.steps), _ in 1:ncopies]
    flip = _residue!(rng)
    L = ncopies == 2 ? _mulp(x[1], _invp(x[nf + 1])) : Int64(1)
    conc = [_zeroes(face, n) ? Int64(0) : _residue!(rng) for n in 1:nnames]
    (; K, ss, re, flip, L, conc)
end

"""`r` times `value[n]` once for each copy of metabolite `n` (an index into a scheme's
`metabolites`) that `held` holds, mod p."""
function _mul_held(r::Int64, held::Vector{Int}, value::AbstractVector{Int64})
    for n in held
        r = _mulp(r, value[n])
    end
    r
end

"""`pt` with new rapid-equilibrium and flip speeds: another speed pattern."""
_respeed_point(pt, rng::_ResidueStream) =
    merge(pt, (re = [_residue!(rng) for _ in axes(pt.re, 1), _ in axes(pt.re, 2)],
               flip = _residue!(rng)))

"""The rate mod p of every transition of `net` at the concentrations `conc` (per name):
its speed times `K` in the stored direction, its speed alone in reverse, times each
metabolite it takes up, a fast one times `Λ` unless `Λ` is 0."""
function _rates_p(m::Mechanisms.CatalyticScheme, net::_Network, pt, conc::Vector{Int64},
                  Λ::Int64)
    r = zeros(Int64, length(net.from))
    for t in eachindex(net.from)
        k, c = net.edge[t], net.conformation[t]
        if k == 0
            r[t] = net.forward[t] ? _mulp(pt.flip, pt.L) : pt.flip
        else
            speed = net.fast[t] ? pt.re[k, c] : pt.ss[m.steps[k].kinetic_group, c]
            r[t] = net.forward[t] ? _mulp(speed, pt.K[k, c]) : speed
            held = net.forward[t] ? m.steps[k].taken_up : m.steps[k].given_off
            r[t] = _mul_held(r[t], held, conc)
        end
        net.fast[t] && Λ != 0 && (r[t] = _mulp(r[t], Λ))
    end
    r
end

"""ρ(u) mod p: over the slow transitions out of node `u`, their rate `r` times their
uptake of `x1` plus the fall of its potential `h` along them."""
function _slow_rho(net::_Network, r::Vector{Int64}, h::Vector{Int}, x1::Int)
    rho = zeros(Int64, net.nforms * net.ncopies)
    for t in eachindex(net.from)
        (net.fast[t] || r[t] == 0) && continue
        u = net.from[t]
        rho[u] = _addp(rho[u],
                       _mulp(r[t], _modp(net.change[x1, t] + h[u] - h[net.to[t]])))
    end
    rho
end

"""Equilibrium weights mod p within each part of `parts` over the fast transitions inside
it, from its first node, at rates `r`."""
function _eq_weights_p(net::_Network, parts::Vector{Vector{Int}}, r::Vector{Int64})
    nn = net.nforms * net.ncopies
    part = zeros(Int, nn)
    for (q, P) in enumerate(parts), u in P
        part[u] = q
    end
    next = [Int[] for _ in 1:nn]
    for t in eachindex(net.from)
        u = net.from[t]
        net.fast[t] && part[u] > 0 && part[u] == part[net.to[t]] && push!(next[u], t)
    end
    w = zeros(Int64, nn)
    for P in parts
        w[P[1]] = 1
        done = falses(nn)
        done[P[1]] = true
        todo = [P[1]]
        while !isempty(todo)
            u = pop!(todo)
            for t in next[u]
                v = net.to[t]
                done[v] && continue
                back = isodd(t) ? t + 1 : t - 1
                r[t] != 0 && r[back] != 0 ||
                    error("_face_limit_check: a fast transition inside a part runs one way")
                w[v] = _mulp(w[u], _mulp(r[t], _invp(r[back])))
                done[v] = true
                push!(todo, v)
            end
        end
        all(u -> done[u], P) || error("_face_limit_check: a part is not fast-connected")
    end
    w
end

"""The face chain's limit mod p as every rapid-equilibrium speed grows, by interpolating
the rate's numerator and its tree-sum denominator in Λ (degree below the closed class's
size) on the face chain's closed class; `_INFINITE` when the rate grows with Λ."""
function _limit_by_interpolation(m::Mechanisms.CatalyticScheme, net::_Network, s, pt,
                                 h::Vector{Int}, x1::Int, rng::_ResidueStream)
    C = s.closed[1]
    n = length(C)
    pos = zeros(Int, net.nforms * net.ncopies)
    pos[C] .= 1:n
    xs, Ds, Ns = Int64[], Int64[], Int64[]
    for _ in 1:(n + 1)
        Λ = _residue!(rng)
        r = _rates_p(m, net, pt, pt.conc, Λ)
        Q = zeros(Int64, n, n)
        for t in eachindex(net.from)
            a, b = pos[net.from[t]], pos[net.to[t]]
            s.alive[t] && r[t] != 0 && a > 0 && b > 0 && a != b || continue
            Q[a, b] = _addp(Q[a, b], r[t])
        end
        x, D = _stationaryp(Q)
        isempty(x) && error("_face_limit_check: singular at a random Λ")
        rho = _slow_rho(net, r, h, x1)
        v = reduce(_addp, (_mulp(x[pos[u]], rho[u]) for u in C); init = Int64(0))
        push!(xs, Λ)
        push!(Ds, D)
        push!(Ns, _mulp(D, v))
    end
    cD, cN = _interpp(xs[1:n], Ds[1:n]), _interpp(xs[1:n], Ns[1:n])
    _pevalp(cD, xs[end]) == Ds[end] && _pevalp(cN, xs[end]) == Ns[end] ||
        error("_face_limit_check: the degree bound failed")
    dD, dN = _degp(cD), _degp(cN)
    dN > dD && return _INFINITE
    dN < dD && return Int64(0)
    _mulp(cN[dD + 1], _invp(cD[dD + 1]))
end

"""The law mod p at concentrations `conc` with every name present: segment weights and
the slow segment chain over every segment."""
function _law_inside(m::Mechanisms.CatalyticScheme, net::_Network, s, pt,
                     conc::Vector{Int64}, h::Vector{Int}, x1::Int)
    r = _rates_p(m, net, pt, conc, Int64(0))
    w = _eq_weights_p(net, s.segments, r)
    W = [_invp(reduce(_addp, w[S]; init = Int64(0))) for S in s.segments]
    Q = zeros(Int64, length(s.segments), length(s.segments))
    for t in eachindex(net.from)
        net.fast[t] && continue
        a, b = s.segment[net.from[t]], s.segment[net.to[t]]
        a == b && continue
        Q[a, b] = _addp(Q[a, b], _mulp(_mulp(w[net.from[t]], W[a]), r[t]))
    end
    mu = first(_stationaryp(Q))
    isempty(mu) && error("_face_limit_check: singular segment chain")
    rho = _slow_rho(net, r, h, x1)
    v = Int64(0)
    for (g, S) in enumerate(s.segments), u in S
        v = _addp(v, _mulp(_mulp(mu[g], _mulp(w[u], W[g])), rho[u]))
    end
    v
end

"""The law's value mod p on the face: its limit as the zeroed concentrations go to zero
along two random directions, by rational reconstruction in their scale; `_UNDEFINED` when
a direction meets a pole or the two differ."""
function _law_by_reconstruction(m::Mechanisms.CatalyticScheme, net::_Network, s, pt,
                                h::Vector{Int}, x1::Int, face::UInt64,
                                rng::_ResidueStream)
    found = Int64[]
    for _ in 1:2
        direction = [_residue!(rng) for _ in pt.conc]
        law_at(δ) = _law_inside(m, net, s, pt,
                                [_zeroes(face, n) ? _mulp(direction[n], δ) : pt.conc[n]
                                 for n in eachindex(pt.conc)], h, x1)
        num, den = _reconstructp(law_at, rng)
        isempty(den) && error("_face_limit_check: no reconstruction under the degree cap")
        (a, lead_num), (b, lead_den) = _plowp(num), _plowp(den)
        push!(found, a > b ? Int64(0) : a < b ? _INFINITE :
                     _mulp(lead_num, _invp(lead_den)))
    end
    found[1] == _INFINITE || found[1] != found[2] ? _UNDEFINED : found[1]
end

"""
    _face_limit_check(mech, net, s, face, seed) -> _LimitVerdict

The limit test of class 3 on face `face` of `mech`, whose law
network is `net` and face structure `s`: at two random points of the model mod 2^61 − 1,
drawn from `seed`, the fast-RE limit of the face chain by interpolation in Λ against the
law's face value by reconstruction as the zeroed concentrations vanish. On a difference a
third limit at new speeds tells `_LIMIT_DIFFER` from `_LIMIT_SPEED`. A false
`_LIMIT_EQUAL` has probability at most twice the degree over p.
"""
function _face_limit_check(mech::Mechanisms.Mechanism, net::_Network, s, face::UInt64,
                           seed::UInt64)
    m = mech.catalytic
    rng = _ResidueStream(seed)
    x1 = first(m.substrates)
    # Class 2, which runs first, leaves no fast cycle that takes `x1` up: a potential.
    h = s.offset[x1, :]
    limit(pt) = _limit_by_interpolation(m, net, s, pt, h, x1, rng)
    for _ in 1:2
        pt = _model_point(mech, net.ncopies, face, rng)
        law = _law_by_reconstruction(m, net, s, pt, h, x1, face, rng)
        law == _UNDEFINED && return _LIMIT_LAW_NAN
        v = limit(pt)
        # The face chain's fast transitions are the law network's, a subset of the
        # ε-network's, where class 2, which runs first, finds no turnover cycle.
        v == _INFINITE && error("_face_limit_check: the limit grows without bound")
        v == law && continue
        return limit(_respeed_point(pt, rng)) == v ? _LIMIT_DIFFER : _LIMIT_SPEED
    end
    _LIMIT_EQUAL
end

# ─── Cycle-space helpers of the limit test ─────────────────────────────────

"""
Reduced row echelon form over `Rational{BigInt}`. Returns the pivot and
free column indices (pivot_cols in row-pivot order, so pivot_cols[i] is the
pivot at reduced-matrix row i) plus the reduced matrix R. Used by
`_rational_nullspace` (nullspace basis).
"""
function _rref_partition(A::AbstractMatrix)
    m, n = size(A)
    R = Matrix{Rational{BigInt}}(A)
    pivot_cols = Int[]
    row = 1
    for col in 1:n
        piv = findfirst(r -> R[r, col] != 0, row:m)
        piv === nothing && continue
        piv += row - 1
        R[row, :], R[piv, :] = R[piv, :], R[row, :]
        R[row, :] ./= R[row, col]
        for r in 1:m
            r != row && R[r, col] != 0 && (R[r, :] .-= R[r, col] .* R[row, :])
        end
        push!(pivot_cols, col); row += 1
    end
    free_cols = setdiff(1:n, pivot_cols)
    return pivot_cols, free_cols, R
end

"""
Basis of `{x : M·x = 0}` as the COLUMNS of the returned matrix, exact over the
rationals: each free column of `M`'s reduced row echelon form yields one basis
vector, 1 at its free coordinate.
"""
function _rational_nullspace(M::AbstractMatrix)
    pivots, free, R = _rref_partition(M)
    N = zeros(Rational{BigInt}, size(M, 2), length(free))
    for (k, f) in enumerate(free)
        N[f, k] = 1
        for (r, c) in enumerate(pivots); N[c, k] = -R[r, f]; end
    end
    N
end

"""
`_rational_nullspace(A)` with each basis column scaled to coprime integers whose
first nonzero entry is positive.
"""
function _integer_nullspace(A::Matrix{Int})
    N = _rational_nullspace(A)
    result = zeros(Int, size(N))
    for k in axes(N, 2)
        col = Int.(N[:, k] .* lcm(denominator.(N[:, k])))
        col .÷= gcd(col)
        result[:, k] .= col[findfirst(!iszero, col)] < 0 ? -col : col
    end
    result
end

# ─── Guard messages ────────────────────────────────────────────────────

"""
    _guard_message(mech) -> String

The error for the guard's verdict on `mech`:
`Mechanism: <cause, naming each step by its line>. <fix>.`, the
cause read under the reading the verdict was judged under (`_guard_with_reading`). Class
2's fix is printed when each step's variant leaves no fast turnover cycle
(`_every_step_cures`); every other fix only when its variant passes `_guard`. When two
constants share a name (`_SHARED_NAME`), the error names them and asks for a renaming;
for a cycle that does not balance a metabolite (`_UNBALANCED_CYCLE`), it names the
cycle's steps and its net change (`_unbalanced_cycle`) and the slips that make one.
An error when the guard admits `mech`.
"""
function _guard_message(mech::Mechanisms.Mechanism)
    (; verdict, reading) = _guard_with_reading(mech)
    verdict.fault === _NO_FAULT && error("_guard_message: verdict $verdict holds no fault")
    verdict.fault === _SHARED_NAME &&
        return "Mechanism: constants of this mechanism's rate equation share a name (" *
               join(_shared_constant_names(mech), ", ") * "). Rename a metabolite or an " *
               "enzyme form so that each constant has a name of its own."
    if verdict.fault === _UNBALANCED_CYCLE
        m = mech.catalytic
        (; steps, change) = _unbalanced_cycle(m)
        lines = join([Mechanisms.step_text(m, k) for k in steps], ", ", " and ")
        # The metabolites on one side of the net change, each with its count when above 1.
        function side(sgn)
            listed = [(abs(c) == 1 ? "" : "$(abs(c)) ") * String(m.metabolites[x])
                      for (x, c) in enumerate(change) if sgn * c > 0]
            isempty(listed) ? "nothing" : join(listed, ", ", " and ")
        end
        return "Mechanism: the steps $lines close a cycle that, net, takes up $(side(1)) " *
               "and gives off $(side(-1)). Round every cycle, the substrates go in and " *
               "the products come out as often as the reaction runs, and each modifier " *
               "comes off as often as it binds. Check whether one of these steps leaves " *
               "out a metabolite, writes one on the wrong side or names the wrong form."
    end
    verdict.fault === _NO_CYCLE &&
        return "Mechanism: no cycle of steps turns substrates into products, so the rate " *
               "is zero everywhere."
    verdict.fault === _RE_TURNOVER && return _re_turnover_message(mech, reading)
    # Only an allosteric mechanism excludes steps, so only its verdict is admission's.
    verdict.fault in _ADMISSION_FAULTS && return _admission_message(mech, verdict, reading)
    _face_message(mech, verdict)
end

"""Whether making each step of `cycle`, transitions of `mech`'s ε-network `net` under
`reading`, steady-state with its kinetic group leaves no fast cycle that turns the
reaction over: then class 2's fix holds whichever step of the cycle is picked."""
function _every_step_cures(mech::Mechanisms.Mechanism, net::_Network, cycle::Vector{Int},
                           reading::_Reading)
    m = mech.catalytic
    groups = _kinetic_groups(m)
    for t in cycle
        re = _re_flags(m)
        re[groups[m.steps[net.edge[t]].kinetic_group]] .= false
        isempty(_re_turnover_transitions(_eps_network(mech, re, reading))) || return false
    end
    true
end

"""Class 2's message, the excluded steps switched off as `reading` says: the steps of the
first fast cycle that turns the reaction over, each excluded step among them marked with
the conformation where it loses a direction (its affinity, when it takes up one
metabolite and loses that uptake), and the fix where each of those steps made
steady-state leaves no fast cycle that turns the reaction over (`_every_step_cures`)."""
function _re_turnover_message(mech::Mechanisms.Mechanism, reading::_Reading)
    m = mech.catalytic
    net = _eps_network(mech, _re_flags(m), reading)
    cycle = _re_turnover_transitions(net)
    c = net.conformation[cycle[1]]
    state = c == 2 ? "inactive" : "active"
    lines = String[]
    for t in cycle
        k = net.edge[t]
        line = Mechanisms.step_text(m, k)
        if _excluded_steps(mech, c)[k]
            e = m.steps[k]
            affinity = reading[k, c] === _LOSES_FORWARD && length(e.taken_up) == 1 &&
                       isempty(e.given_off)
            line *= affinity ? " (whose affinity vanishes in the $state conformation)" :
                    " (whose rate one way vanishes in the $state conformation)"
        end
        push!(lines, line)
    end
    # Only an allosteric mechanism's network has an inactive copy.
    place = c == 2 ? " in the inactive conformation" : ""
    msg = "Mechanism: the rapid-equilibrium steps $(join(lines, ", ", " and ")) close a " *
          "cycle$place that turns the reaction over, so the mechanism has no finite rate."
    fix = " Make one step of the cycle steady-state."
    _every_step_cures(mech, net, cycle, reading) ? msg * fix : msg
end

"""Where a class-3 fault of `mech` lies: nothing for a plain mechanism; the conformation
`c` of an allosteric one, nothing when `c` is 0."""
_place(mech::Mechanisms.Mechanism, c::Int) =
    !_is_allosteric(mech) ? "" : c == 2 ? " in the inactive conformation" :
    c == 1 ? " in the active conformation" : ""

"""
The message of a class-3 verdict `v` on `mech`: the face, the forms and steps at fault,
and a fix only where its variant passes the guard (`_steady_state_passes`): for a trap,
the kinetic group of the step that feeds it or joins it to its segment made steady-state;
for a bottomless segment, the kinetic group of the first of its rapid-equilibrium bindings
whose variant passes. A verdict in no one conformation labels each form with its own
(`_conformation_suffix`).
"""
function _face_message(mech::Mechanisms.Mechanism, v::_Verdict)
    m = mech.catalytic
    net = _law_network(mech)
    face = _faces(m)[v.face]
    s = _face_structure(net, face)
    at = _face_phrase(m, face, _place(mech, v.conformation))
    # Forms that lie in different conformations name each form's conformation.
    label(u) = v.conformation == 0 ? _conformation_suffix(net, u) : ""
    forms(nodes) =
        "{" * join([String(m.enzyme_forms[_form_of(net, u)]) * label(u) for u in nodes],
                   ", ") * "}"
    lines(edges) = join([Mechanisms.step_text(m, k) for k in edges], ", ", " and ")
    if v.fault === _NONUNIQUE
        C = _away_class(s)
        return "Mechanism: $(at)the forms $(forms(C)) keep enzyme away from free E: " *
               "every step that leaves them takes up a zeroed metabolite " *
               "($(lines(_exits(net, s, C)))), so the steady state is not unique."
    elseif v.fault === _BOTTOMLESS
        binding(k) = length(m.steps[k].taken_up) + length(m.steps[k].given_off) == 1
        bindings = unique!([net.edge[t] for t in eachindex(net.from)
                            if net.fast[t] && net.edge[t] != 0 &&
                               s.segment[net.from[t]] == s.bottomless &&
                               binding(net.edge[t])])
        k = findfirst(k -> _steady_state_passes(mech, k), bindings)
        fix = k === nothing ? "" :
              " Make the kinetic group of $(Mechanisms.step_text(m, bindings[k])) " *
              "steady-state."
        return "Mechanism: $(at)the rapid-equilibrium segment " *
               "$(forms(s.segments[s.bottomless])) has no form whose weight keeps the " *
               "segment's lowest power of each zeroed metabolite, so the forms' shares " *
               "depend on how those concentrations vanish and the rate is 0/0 there." * fix
    elseif v.fault === _TRAP
        trap, edge, joined = _trap_feeder(net, s)
        fix = _steady_state_passes(mech, edge) ?
              " Make that step's kinetic group steady-state." : ""
        joined &&
            return "Mechanism: $(at)the forms $(forms(s.fast_classes[trap])) hold " *
                   "enzyme, but the rapid-equilibrium step " *
                   "$(Mechanisms.step_text(m, edge)) that joins them to the rest of " *
                   "their segment takes up a zeroed metabolite, so the rate law " *
                   "misweights them." * fix
        return "Mechanism: $(at)the rapid-equilibrium step " *
               "$(Mechanisms.step_text(m, edge)) fills the forms " *
               "$(forms(s.fast_classes[trap])), which no rapid-equilibrium step running " *
               "there empties, so the rate law misweights them." * fix
    elseif v.fault === _SPLIT
        split = _split_parts(net, s)
        return "Mechanism: $(at)the forms $(join(forms.(split.parts), " and ")) of one " *
               "rapid-equilibrium segment are joined only by steps that take up a " *
               "zeroed metabolite ($(lines(split.joins))), so the rate law cannot weight " *
               "one against the other."
    elseif v.fault === _LAW_MULTI || v.fault === _MULTILEVEL
        why = v.limit === _LIMIT_LAW_NAN ?
                  "the rate law has no value (0/0), while the steady state does" :
              v.limit === _LIMIT_SPEED ?
                  "the steady state depends on how fast the rapid-equilibrium steps " *
                  "run, which the rate law does not carry" :
                  "the steady state's fast-equilibrium limit differs from the rate law"
        return "Mechanism: $(at)$why."
    end
    error("_face_message: no message for $(v.fault)")
end

# ─── Message helpers that do work for every mechanism ───────────────────

"""
Whether `mech` with the kinetic group of edge `k` made steady-state passes the guard
(`_guard`). Every fix that makes a group steady-state in the guard's messages, but class
2's, is checked this way.
"""
function _steady_state_passes(mech::Mechanisms.Mechanism, k::Int)
    m = mech.catalytic
    g = m.steps[k].kinetic_group
    steps = [Mechanisms.Step(s.from, s.to, s.taken_up, s.given_off,
                             s.rapid_equilibrium && s.kinetic_group != g, s.kinetic_group)
             for s in m.steps]
    _guard(_with_steps(mech, steps)) === _WELL_POSED
end

"""`m` with the catalytic steps `steps`, each in the kinetic group of `m` its
`kinetic_group` labels, built by the constructors; it keeps the allostery of `m`."""
function _with_steps(m::Mechanisms.Mechanism, steps::Vector{Mechanisms.Step})
    c = m.catalytic
    Mechanisms.Mechanism(Mechanisms.CatalyticScheme(c.metabolites, c.enzyme_forms,
                                                    c.substrates, c.products, steps),
                         m.allosteric)
end

"""How a message names face `face` of `m`: nothing for the interior, else `at zero `, its
zeroed names (`_face_names`) and `place`, then a comma."""
_face_phrase(m::Mechanisms.CatalyticScheme, face::UInt64, place::String = "") =
    face == 0 ? "" :
    "at zero " * join(String.(_face_names(m, face)), ", ", " and ") * place * ", "
