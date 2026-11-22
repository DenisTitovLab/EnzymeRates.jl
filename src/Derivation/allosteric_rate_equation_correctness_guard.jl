# ABOUTME: Admission, the guard's check that an allosteric law is the limit as the excluded
# ABOUTME: steps are switched off: the readings, its interior and face conditions, messages.

# ─── How each conformation switches off the steps it excludes ───────────

"""
    exclusion_components(m::CatalyticScheme, excluded::BitVector) -> Vector{Int}

Each form's component under the edges of the catalytic scheme `m` that a conformation
keeps, all but those `excluded` flags (`_excluded_steps`), every form included, labelled
by its smallest form index, so free E's label is 1. The components are the conformation's
islands: a form outside free E's is cut off in that conformation.
"""
function exclusion_components(m::Mechanisms.CatalyticScheme, excluded::BitVector)
    parent = collect(eachindex(m.enzyme_forms))
    for (k, e) in enumerate(m.steps)
        excluded[k] && continue
        a, b = _union_find_root(parent, e.from), _union_find_root(parent, e.to)
        a == b || (parent[max(a, b)] = min(a, b))
    end
    [_union_find_root(parent, f) for f in eachindex(m.enzyme_forms)]
end

"""
Which excluded steady-state uptakes into a cut-off form (`_built_reading` defines them)
lose their affinity in a built reading of the guard; the others stop:
- `_ALL_LOSE_AFFINITY`: every one;
- `_OUTWARD_LOSE_AFFINITY`: those that lead away from free E;
- `_ALL_STOPPED`: none.
The guard tries the built readings in this order, before it searches the other readings
(`_guard_with_reading`).

The steps a conformation keeps join the forms into islands (`exclusion_components`). An
island's distance is the fewest excluded uptakes into a cut-off form, rapid-equilibrium or
steady-state, each followed from its unbound side to its bound side, that lead to it from
free E's island; an island they do not reach is infinitely far. Such an uptake leads away
from free E when its bound side's island is farther than its unbound side's.
"""
@enum _SlowUptakeReading::UInt8 _ALL_LOSE_AFFINITY _OUTWARD_LOSE_AFFINITY _ALL_STOPPED

"""
    _built_reading(am, r = _ALL_LOSE_AFFINITY) -> _Reading

How the conformations of `am` switch off its catalytic steps in the guard's built reading
`r`, read from the reaction graph alone. A conformation keeps a step it does not exclude.
An excluded uptake into a cut-off form is an excluded step that takes up one metabolite
and gives off none, whose two sides lie in different islands (`exclusion_components`) and
whose bound side lies outside free E's island. At rapid equilibrium such a step loses its
affinity: its stored direction, the uptake, vanishes (`_LOSES_FORWARD`). At steady state it
loses its affinity or stops, as `r` says (`_SlowUptakeReading`). Any other excluded step
keeps its kind, rapid equilibrium or steady state, as every step does in both
conformations. At steady state it stops (`_STOPPED`). At rapid equilibrium the built
reading has no way to switch it off (`_CANNOT_SWITCH_OFF`): the ε-network treats it as
slow, vanishing both ways, which makes the reading inconsistent. A searched reading
(`_consistent_columns`) may still switch off such a step between two islands by letting it
lose a direction. The developer page, "Well-posedness guard", says why no reading switches
off one inside an island. A plain mechanism has one conformation, which runs every step.
"""
function _built_reading(am::Mechanisms.Mechanism,
                        r::_SlowUptakeReading = _ALL_LOSE_AFFINITY)
    steps = am.catalytic.steps
    _is_allosteric(am) || return fill(_RUNS, length(steps), 1)
    reading = fill(_RUNS, length(steps), 2)
    for c in 1:2
        excluded = _excluded_steps(am, c)
        any(excluded) || continue
        comp = exclusion_components(am.catalytic, excluded)
        uptake = [excluded[k] && length(e.taken_up) == 1 && isempty(e.given_off) &&
                  comp[e.from] != comp[e.to] && comp[e.to] != comp[1]
                  for (k, e) in enumerate(steps)]
        # Each island's distance from free E's, breadth-first along the uptakes; an
        # island's label is its smallest form index.
        distance = fill(typemax(Int), length(comp))
        distance[comp[1]] = 0
        queue, head = [comp[1]], 0
        while head < length(queue)
            head += 1
            u = queue[head]
            for (k, e) in enumerate(steps)
                uptake[k] && comp[e.from] == u && distance[comp[e.to]] == typemax(Int) ||
                    continue
                distance[comp[e.to]] = distance[u] + 1
                push!(queue, comp[e.to])
            end
        end
        for (k, e) in enumerate(steps)
            excluded[k] || continue
            loses = uptake[k] &&
                    (e.rapid_equilibrium || r === _ALL_LOSE_AFFINITY ||
                     r === _OUTWARD_LOSE_AFFINITY &&
                     distance[comp[e.to]] > distance[comp[e.from]])
            reading[k, c] = loses ? _LOSES_FORWARD :
                            e.rapid_equilibrium ? _CANNOT_SWITCH_OFF : _STOPPED
        end
    end
    reading
end

# ─── Admission: the interior ────────────────────────────────────────────

"""Admission's arcs in conformation `c` of `reading`: the edges of `am`'s catalytic graph
that lose a direction there, in edge order, each as `(lower side's form, higher side's
form, edge)`, so that it points uphill, along the direction that vanishes."""
function _arcs(am::Mechanisms.Mechanism, c::Int, reading::_Reading)
    arcs = NTuple{3, Int}[]
    for (k, e) in enumerate(am.catalytic.steps)
        reading[k, c] === _LOSES_FORWARD && push!(arcs, (e.from, e.to, k))
        reading[k, c] === _LOSES_REVERSE && push!(arcs, (e.to, e.from, k))
    end
    arcs
end

"""
The arcs (`_arcs`) of one directed cycle that `arcs` close over the components `comp`,
each arc read from its lower side's component to its higher side's, in cycle order; empty
when they close none. An arc inside one component closes a cycle alone.
"""
function _arc_cycle(comp::Vector{Int}, arcs::Vector{NTuple{3, Int}})
    # Kahn's algorithm removes every component no cycle feeds; a walk back along the
    # remaining arcs then lands on a cycle.
    tail(i) = comp[arcs[i][1]]
    head(i) = comp[arcs[i][2]]
    indeg = zeros(Int, length(comp))
    for i in eachindex(arcs)
        indeg[head(i)] += 1
    end
    todo = [r for r in unique(comp) if indeg[r] == 0]
    while !isempty(todo)
        u = pop!(todo)
        for i in eachindex(arcs)
            tail(i) == u && (indeg[head(i)] -= 1) == 0 && push!(todo, head(i))
        end
    end
    start = 0
    for i in eachindex(arcs)
        indeg[tail(i)] > 0 && indeg[head(i)] > 0 && (start = i; break)
    end
    start == 0 && return NTuple{3, Int}[]
    function back(x)
        for i in eachindex(arcs)
            head(i) == x && indeg[tail(i)] > 0 && return i
        end
        0
    end
    x = head(start)
    for _ in eachindex(comp)
        x = tail(back(x))
    end
    cycle = [back(x)]
    while tail(cycle[end]) != x
        push!(cycle, back(tail(cycle[end])))
    end
    arcs[reverse!(cycle)]
end

"""The components (`comp` labels) that free E's reaches along `arcs` (`_arcs`), each read
from its lower side's component to its higher side's."""
function _reached_components(comp::Vector{Int}, arcs::Vector{NTuple{3, Int}})
    next = [Int[] for _ in comp]
    for (low, high, _) in arcs
        push!(next[comp[low]], comp[high])
    end
    _reached(next, [comp[1]])
end

"""
    _interior_verdict(am, c, excluded, comp, reading) -> _Verdict

Admission's interior checks on conformation `c` of `am`, its excluded edges (`excluded`,
`_excluded_steps`) switched off as `reading` says. `comp` holds the components of the forms
its other edges join (`exclusion_components`), and the edges that lose a direction are
arcs (`_arcs`), each pointing uphill. The checks, in order:
- kind (`_FAST_STEP_PINNED`, the first such edge): `reading` stops no excluded
  rapid-equilibrium edge and leaves none unswitched (`_CANNOT_SWITCH_OFF`), since a step
  keeps its kind in both conformations;
- existence (`_ARC_CYCLE`, led by the cycle's first edge, else by the first arc into free
  E's component): the arcs close no directed cycle, and none points into free E's
  component. With every other constant fixed, the exponents of a direction's equilibrium
  ratio are differences of component potentials. The arcs must rise, and every cut-off
  component must lie above free E's: it holds no enzyme, while forms below free E's would
  hold all of it. Such potentials exist exactly when this check passes;
- uniqueness (`_ISLAND`, the first excluded edge between a component free E's reaches and
  one it does not): free E's component reaches every component along the arcs, so every
  island lies behind an arc on its higher side and the limit is deletion plus pruning.
The verdict carries the conformation and the edge (face 0); `_WELL_POSED` when all hold.
"""
function _interior_verdict(am::Mechanisms.Mechanism, c::Int, excluded::BitVector,
                           comp::Vector{Int}, reading::_Reading)
    for (k, e) in enumerate(am.catalytic.steps)
        o = reading[k, c]
        (o === _CANNOT_SWITCH_OFF || o === _STOPPED && e.rapid_equilibrium) &&
            return _Verdict(_FAST_STEP_PINNED, _LIMIT_UNTESTED, 0, c, k)
    end
    arcs = _arcs(am, c, reading)
    cycle = _arc_cycle(comp, arcs)
    isempty(cycle) || return _Verdict(_ARC_CYCLE, _LIMIT_UNTESTED, 0, c, cycle[1][3])
    for (_, high, k) in arcs
        comp[high] == comp[1] && return _Verdict(_ARC_CYCLE, _LIMIT_UNTESTED, 0, c, k)
    end
    reached = _reached_components(comp, arcs)
    for (k, e) in enumerate(am.catalytic.steps)
        excluded[k] && reached[comp[e.from]] != reached[comp[e.to]] &&
            return _Verdict(_ISLAND, _LIMIT_UNTESTED, 0, c, k)
    end
    _WELL_POSED
end

"""
    _admission_interior(am, reading) -> _Verdict

Admission's interior checks, each excluded edge switched off as `reading` says: for each
conformation whose tags exclude an edge (`_excluded_steps`), the checks of
`_interior_verdict`; then the active conformation's turnover (`_active_turnover_verdict`).
A reading is consistent when none of `_interior_verdict`'s faults arises. A plain
mechanism excludes no edge, so its interior passes.
"""
function _admission_interior(am::Mechanisms.Mechanism, reading::_Reading)
    for c in 1:2
        excluded = _excluded_steps(am, c)
        any(excluded) || continue
        v = _interior_verdict(am, c, excluded, exclusion_components(am.catalytic, excluded),
                              reading)
        v === _WELL_POSED || return v
    end
    _active_turnover_verdict(am)
end

"""Admission's turnover check, which reads no reading: when tags exclude edges of the
active conformation of `am`, its law network still holds a cycle that turns the reaction
over (`_NO_ACTIVE_TURNOVER`, led by the first edge it excludes); `_WELL_POSED` otherwise."""
function _active_turnover_verdict(am::Mechanisms.Mechanism)
    k = something(findfirst(_excluded_steps(am, 1)), 0)
    k == 0 && return _WELL_POSED
    runs = Mechanisms.steps_running_in_allosteric_conformations(am)
    Mechanisms.has_turnover_cycle(am.catalytic, runs.active) ||
        return _Verdict(_NO_ACTIVE_TURNOVER, _LIMIT_UNTESTED, 0, 1, k)
    _WELL_POSED
end

"""
    _consistent_columns(am, c) -> (; pinned, too_many, columns)

The ways conformation `c` of `am` may switch off the steps it excludes that admission's
interior accepts (`_interior_verdict`), each a column of a `_Reading`, in the search's
order. The steps `c` keeps join the forms into islands (`exclusion_components`), and every
island but free E's lies higher. An excluded step inside one island stops at steady state;
at rapid equilibrium no way switches it off, so `pinned` is the first such step and
`columns` is empty. An excluded step between two islands has up to three ways, in this
order: it loses its stored direction unless that points into free E's island, loses its
reverse unless that does, or, at steady state, stops. The columns run through these
choices with the last step changing fastest. With more than
`_MAX_CHOICES_PER_CONFORMATION` choices, `too_many` holds and `columns` is empty.
"""
function _consistent_columns(am::Mechanisms.Mechanism, c::Int)
    steps = am.catalytic.steps
    excluded = _excluded_steps(am, c)
    comp = exclusion_components(am.catalytic, excluded)
    reading = fill(_RUNS, length(steps), 2)
    between, options = Int[], Vector{_SwitchOff}[]
    columns = Vector{_SwitchOff}[]
    choices = 1
    for (k, e) in enumerate(steps)
        excluded[k] || continue
        if comp[e.from] == comp[e.to]
            e.rapid_equilibrium && return (; pinned = k, too_many = false, columns)
            reading[k, c] = _STOPPED
            continue
        end
        ways = _SwitchOff[]
        comp[e.to] == comp[1] || push!(ways, _LOSES_FORWARD)
        comp[e.from] == comp[1] || push!(ways, _LOSES_REVERSE)
        e.rapid_equilibrium || push!(ways, _STOPPED)
        push!(between, k)
        push!(options, ways)
        choices = min(choices * length(ways), _MAX_CHOICES_PER_CONFORMATION + 1)
    end
    choices > _MAX_CHOICES_PER_CONFORMATION &&
        return (; pinned = 0, too_many = true, columns)
    pick = ones(Int, length(between))
    while true
        for (j, k) in enumerate(between)
            reading[k, c] = options[j][pick[j]]
        end
        _interior_verdict(am, c, excluded, comp, reading) === _WELL_POSED &&
            push!(columns, reading[:, c])
        # The next choice: the last step that has a later way takes it, and every step
        # after it starts again at its first.
        d = length(between)
        while d > 0 && pick[d] == length(options[d])
            pick[d] = 1
            d -= 1
        end
        d == 0 && return (; pinned = 0, too_many = false, columns)
        pick[d] += 1
    end
end

# ─── Admission: the face condition ──────────────────────────────────────

"""
    _face_sets(net, law_nodes, face) -> (; alive, fed, n0, stranded)

What admission's face condition reads on face `face` of the ε-network `net`, whose law
network holds the nodes `law_nodes`:
- `alive`, the transitions alive there, which take up no zeroed name;
- `fed`, the nodes free E's active copy reaches along them;
- `n0`, the alive transitions that do not vanish;
- `stranded`, the stranded basin: fed nodes outside the law network that reach it along
  no chain of alive transitions that do not vanish.
"""
function _face_sets(net::_Network, law_nodes::BitVector, face::UInt64)
    alive = BitVector([(net.takes[t] & face) == 0 for t in eachindex(net.takes)])
    n0 = alive .& .!net.vanishing
    fed = _reached(_successors(net, alive), [_node(net, 1, 1)])
    drains = _reached(_successors(net, n0; backward = true),
                      findall(law_nodes .& net.present))
    (; alive, fed, n0, stranded = fed .& .!law_nodes .& .!drains)
end

"""Copy `c`'s alive transitions with both ends in `nodes`: `van`, those that vanish with
their tags, and `held`, those `_potential_on_sccs` reads (not vanishing, or fast)."""
function _copy_transitions(net::_Network, alive::BitVector, nodes::BitVector, c::Int)
    van, held = falses(length(net.from)), falses(length(net.from))
    for t in eachindex(net.from)
        alive[t] && net.conformation[t] == c && nodes[net.from[t]] && nodes[net.to[t]] ||
            continue
        van[t] = net.vanishing[t]
        held[t] = !net.vanishing[t] || net.fast[t]
    end
    (; van, held)
end

"""
The cost of a cheapest spanning out-arborescence of the nodes `1:n` rooted at `root` over
`arcs`, each `(from, to, cost)`; -1 when some node has no entering arc.
"""
function _min_arborescence(n::Int, arcs::Vector{NTuple{3, Int}}, root::Int)
    # Chu–Liu/Edmonds: take each node's cheapest entering arc, contract the cycles they
    # close and repeat.
    INF = typemax(Int) ÷ 4
    cost = 0
    while true
        inw, pre = fill(INF, n), zeros(Int, n)
        for (u, v, w) in arcs
            u != v && w < inw[v] && (inw[v] = w; pre[v] = u)
        end
        for v in 1:n
            v != root && inw[v] == INF && return -1
        end
        inw[root] = 0
        id, vis = zeros(Int, n), zeros(Int, n)
        cnt = 0
        for v in 1:n
            cost += inw[v]
            x = v
            while vis[x] != v && id[x] == 0 && x != root
                vis[x] = v
                x = pre[x]
            end
            if x != root && id[x] == 0
                cnt += 1
                y = pre[x]
                while y != x
                    id[y] = cnt
                    y = pre[y]
                end
                id[x] = cnt
            end
        end
        cnt == 0 && return cost
        for v in 1:n
            id[v] == 0 && (cnt += 1; id[v] = cnt)
        end
        next = NTuple{3, Int}[]
        for (u, v, w) in arcs
            id[u] != id[v] && push!(next, (id[u], id[v], w - inw[v]))
        end
        arcs, n, root = next, cnt, id[root]
    end
end

"""
    _fewest_slow(net, nodes, keep) -> Int

The fewest slow transitions in a spanning in-tree of the nodes `nodes` over the
transitions `keep` flags (cost 1 slow, 0 fast), any root; -1 when no in-tree spans them.
"""
function _fewest_slow(net::_Network, nodes::BitVector, keep::BitVector)
    # A virtual root joined to every node at a cost above any tree's turns the in-tree into
    # Chu–Liu/Edmonds's out-arborescence (`_min_arborescence`) on reversed transitions; a
    # cost of two such joins means no single root serves.
    at = zeros(Int, length(nodes))
    n = 0
    for u in eachindex(nodes)
        nodes[u] && (n += 1; at[u] = n)
    end
    arcs = NTuple{3, Int}[]
    for t in eachindex(net.from)
        a, b = at[net.from[t]], at[net.to[t]]
        keep[t] && a != 0 && b != 0 && push!(arcs, (b, a, net.fast[t] ? 0 : 1))
    end
    for u in 1:n
        push!(arcs, (n + 1, u, n + 1))
    end
    c = _min_arborescence(n + 1, arcs, n + 1)
    c < 0 || c >= 2 * (n + 1) ? -1 : c - (n + 1)
end

"""
    _potential_on_sccs(net, keep, x) -> Bool

Whether the uptake of metabolite `x` (a `change` row) is a potential on every strongly
connected part (`_strong_components`) of the transitions `keep` flags, so no cycle of them
turns the reaction over.
"""
function _potential_on_sccs(net::_Network, keep::BitVector, x::Int)
    comp = _strong_components(_successors(net, keep), net.present)
    inside = BitVector([keep[t] && comp[net.from[t]] == comp[net.to[t]]
                        for t in eachindex(net.from)])
    n = length(net.present)
    !Mechanisms.has_nonzero_weight_cycle(net.from, net.to, view(net.change, x, :), inside,
                                         zeros(Int, n), zeros(Int, n))
end

"""Whether the transitions `extra` shorten the cheapest spanning in-tree of `nodes` over
`base` (`_fewest_slow`)."""
function _shortens_in_tree(net::_Network, nodes::BitVector, base::BitVector,
                           extra::BitVector)
    without = _fewest_slow(net, nodes, base)
    with = _fewest_slow(net, nodes, base .| extra)
    with >= 0 && (without < 0 || with < without)
end

"""
The edges of the vanishing transitions on a cheapest spanning in-tree of `nodes` over
`base` and `van` (`_BRIDGE`), in transition order: an inclusion-minimal set of the
transitions `van` with which `base` still gives that tree's cost (`_fewest_slow`), so
every such tree uses each of them.
"""
function _bridging_edges(net::_Network, nodes::BitVector, base::BitVector, van::BitVector)
    best = _fewest_slow(net, nodes, base .| van)
    kept = copy(van)
    for t in eachindex(van)
        kept[t] || continue
        kept[t] = false
        _fewest_slow(net, nodes, base .| kept) == best || (kept[t] = true)
    end
    unique!([net.edge[t] for t in findall(kept)])
end

"""
Of two or more closed classes of the ε-network, the first that holds no node `law_nodes`
marks: the class an `_EPS_NONUNIQUE` verdict names. Admission reads a face only once class
3 has accepted it, and then such a class exists, because a closed class of the ε-network
that holds a law node holds a closed class of the law network, and class 3 leaves that
network only one.
"""
_class_outside_law(closed::Vector{Vector{Int}}, law_nodes::BitVector) =
    closed[findfirst(cls -> !any(view(law_nodes, cls)), closed)::Int]

"""
    _fast_law_reach(net, law_nodes, s) -> (; from_law, to_law)

On the face whose sets `s` (`_face_sets`) admission's face condition reads, of the
ε-network `net` whose law network holds the nodes `law_nodes`:
- `from_law`, the nodes the law network reaches along alive rapid-equilibrium
  transitions, vanishing or not, so the flow into them can grow with the speed of the
  fast steps;
- `to_law`, the nodes that reach it along alive rapid-equilibrium transitions that do not
  vanish, so they empty as fast as that.
"""
function _fast_law_reach(net::_Network, law_nodes::BitVector, s)
    law = findall(law_nodes .& net.present)
    (; from_law = _reached(_successors(net, s.alive .& net.fast), law),
     to_law = _reached(_successors(net, s.n0 .& net.fast; backward = true), law))
end

"""
    _face_condition(net, law_nodes, face, x) -> _Verdict

Admission's face condition, on face `face` of the ε-network `net` of an allosteric
mechanism, whose law network holds the nodes `law_nodes`; `x` is the `change` row of its
first substrate (`first(m.substrates)`). It refuses, in order:
- `_EPS_NONUNIQUE`: the alive transitions leave two closed classes, so at every ε the
  steady state depends on where the enzyme starts. The verdict names the closed class
  outside the law network (`_class_outside_law`).
- `_FILLED`, naming the step and conformation of a rapid-equilibrium transition that
  vanishes and meets all three conditions below, so that its target fills faster than it
  empties:
  - the transition is alive while its reverse, which does not vanish, is dead;
  - its source is fed and reached from the law network along alive rapid-equilibrium
    transitions, vanishing or not;
  - its target reaches the law network along no alive rapid-equilibrium transitions that
    do not vanish.
- `_PUMPED`: a stranded basin (`_face_sets`) is entered by, or holds, an alive transition
  whose reverse is dead, or is entered from outside it and the law network, so admission
  cannot show that its weight does not depend on how fast its exits vanish. The verdict
  names a feed whose reverse is dead when there is one. A basin entered only reversibly
  from the law network is sealed: its weight vanishes, and it leaves the fed nodes.
- `_BRIDGE`, in either of two cases:
  - in some conformation both of these hold, and the verdict leads with the first of that
    conformation's vanishing steps on a cheapest tree (`_bridging_edges`):
    - its vanishing transitions shorten the cheapest spanning in-tree of the fed nodes
      (`_shortens_in_tree`);
    - the first substrate's uptake is no potential on the strongly connected parts of its
      alive transitions that do not vanish or are fast (`_potential_on_sccs`);
  - no conformation's vanishing transitions shorten the in-tree alone, but both
    conformations' together do; the verdict names no conformation or step.
The tests label these conditions. R1: a conformation's vanishing transitions do not shorten
the in-tree. Joint R1: both conformations' vanishing transitions together do not. T0+: the
first substrate's uptake is a potential on those strongly connected parts. A face whose fed
nodes no in-tree spans passes: class 3 decides it. The verdict leaves the face index 0
(`_admission_face` sets it). The developer page, "Well-posedness guard", argues for the
filling check and for the set-aside rule (`_admission_face`), and gives the numerical
evidence for both.
"""
function _face_condition(net::_Network, law_nodes::BitVector, face::UInt64, x::Int)
    s = _face_sets(net, law_nodes, face)
    closed = _closed_classes(net, s.alive)
    length(closed) > 1 &&
        return _Verdict(_EPS_NONUNIQUE, _LIMIT_UNTESTED, 0,
                        _copy_of(net, _class_outside_law(closed, law_nodes)[1]), 0)
    # Why a one-way fast transition must not fill its form, and the numerical evidence for
    # this check: the developer page, "Well-posedness guard", on the filling check. No
    # built reading has such a transition, since the direction it keeps is a release, alive
    # on every face, so under them the searches of `_fast_law_reach` never run here.
    one_way(t) = (r = isodd(t) ? t + 1 : t - 1;
                  s.alive[t] && !s.alive[r] && net.vanishing[t] && !net.vanishing[r] &&
                  net.fast[t] && s.fed[net.from[t]])
    if any(one_way, eachindex(net.from))
        reach = _fast_law_reach(net, law_nodes, s)
        for t in eachindex(net.from)
            one_way(t) && reach.from_law[net.from[t]] && !reach.to_law[net.to[t]] &&
                return _Verdict(_FILLED, _LIMIT_UNTESTED, 0, net.conformation[t],
                                net.edge[t])
        end
    end
    # The verdict names a feed whose reverse is dead when there is one, and else the first
    # feed whose reverse runs back but that enters from outside the basin and the law
    # network, so that a message saying nothing runs back is true where it can be.
    runs_back = 0
    for t in eachindex(net.from)
        u, w = net.from[t], net.to[t]
        s.alive[t] && s.fed[u] && s.stranded[w] || continue
        s.alive[isodd(t) ? t + 1 : t - 1] ||
            return _Verdict(_PUMPED, _LIMIT_UNTESTED, 0, net.conformation[t], net.edge[t])
        runs_back == 0 && !(s.stranded[u] || law_nodes[u]) && (runs_back = t)
    end
    runs_back == 0 ||
        return _Verdict(_PUMPED, _LIMIT_UNTESTED, 0, net.conformation[runs_back],
                        net.edge[runs_back])
    nodes = s.fed .& .!s.stranded
    _fewest_slow(net, nodes, s.n0) < 0 && return _WELL_POSED
    passed = 0
    for c in 1:net.ncopies
        (; van, held) = _copy_transitions(net, s.alive, nodes, c)
        any(van) || continue
        if !_shortens_in_tree(net, nodes, s.n0, van)
            passed += 1
        elseif !_potential_on_sccs(net, held, x)
            return _Verdict(_BRIDGE, _LIMIT_UNTESTED, 0, c,
                            first(_bridging_edges(net, nodes, s.n0, van)))
        end
    end
    passed == 2 && _shortens_in_tree(net, nodes, s.n0, s.alive .& net.vanishing) &&
        return _Verdict(_BRIDGE, _LIMIT_UNTESTED, 0, 0, 0)
    _WELL_POSED
end

"""
    _without_set_aside_groups(net, law_nodes, face) -> _Network

The ε-network `net`, whose law network holds the nodes `law_nodes`, without the groups
that cannot matter on face `face` and every transition at them; `net` itself when no group
qualifies. A group is the fed nodes outside the law network (`_face_sets`) that alive
transitions join to one another, with every pair of transitions at them that has an alive
direction. It is set aside when
- the law nodes its pairs touch lie in one set that the law network's rapid-equilibrium
  pairs, alive both ways, join (one touching node counts);
- its fast classes drain: in the aggregated chain of the alive transitions that do not
  vanish (`_aggregated_chain`), each of its fast classes leads to a class that holds a law
  node. The chain is the whole network's; it is the group's own, joined to the law's, when
  no transition that does not vanish runs from a law node into a group, as under a
  consistent reading;
- each of its pairs is safe:
  - both directions are alive; or
  - one direction alone is alive, and it does not vanish, so it only empties its source;
    or
  - the pair joins two nodes outside the law network, its vanishing direction alone is
    alive, and that direction does one of these (`_fast_law_reach`):
    - it is slow;
    - it leaves a node the law network does not reach along alive rapid-equilibrium
      transitions;
    - it enters a node that reaches the law network along alive rapid-equilibrium
      transitions that do not vanish.
A node outside the law network that free E's active copy does not reach on the face and
whose alive transitions all lead to nodes set aside goes with them: left alone it would
read as a second closed class. The developer page, "Well-posedness guard", says why such a
group cannot change the limit.
"""
function _without_set_aside_groups(net::_Network, law_nodes::BitVector, face::UInt64)
    s = _face_sets(net, law_nodes, face)
    reach = _fast_law_reach(net, law_nodes, s)
    nn, npairs = length(net.present), length(net.from) ÷ 2
    cut = s.fed .& .!law_nodes
    function join!(parent, a, b)
        ra, rb = _union_find_root(parent, a), _union_find_root(parent, b)
        ra == rb || (parent[max(ra, rb)] = min(ra, rb))
        nothing
    end
    # Each node's group, and each law node's set of nodes that the law network's
    # rapid-equilibrium pairs hold in equilibrium, as union–find parents.
    group, held = collect(1:nn), collect(1:nn)
    for j in 1:npairs
        t, r = 2j - 1, 2j
        a, b = net.from[t], net.to[t]
        cut[a] && cut[b] && (s.alive[t] || s.alive[r]) && join!(group, a, b)
        law_nodes[a] && law_nodes[b] && s.alive[t] && s.alive[r] && net.fast[t] &&
            !net.vanishing[t] && !net.vanishing[r] && join!(held, a, b)
    end
    # The fast classes that drain: those that lead, in the aggregated chain, to a class
    # holding a law node. A slow way out of a node in no fast class, one that the fast
    # steps carry into another class, does not count, or a fast trap inside a group would
    # read as draining.
    chain = _aggregated_chain(net, s.n0)
    back = [Int[] for _ in chain.fast_classes]
    for (q, ks) in enumerate(chain.next), k in ks
        push!(back[k], q)
    end
    drained = _reached(back, [q for (q, C) in enumerate(chain.fast_classes)
                              if any(u -> law_nodes[u], C)])
    # Per group root: whether it may be set aside, and the equilibrium set of the law nodes
    # its pairs touch (0 none yet, -1 more than one).
    qualifies, touched = trues(nn), zeros(Int, nn)
    for u in 1:nn
        q = chain.class_of[u]
        cut[u] && q > 0 && !drained[q] && (qualifies[_union_find_root(group, u)] = false)
    end
    for j in 1:npairs
        t, r = 2j - 1, 2j
        a, b = net.from[t], net.to[t]
        (s.alive[t] || s.alive[r]) && (cut[a] || cut[b]) || continue
        g = _union_find_root(group, cut[a] ? a : b)
        for u in (a, b)
            law_nodes[u] || continue
            h = _union_find_root(held, u)
            touched[g] = touched[g] in (0, h) ? h : -1
        end
        k = s.alive[t] ? t : r
        safe = s.alive[t] && s.alive[r] || !net.vanishing[k] ||
               !law_nodes[a] && !law_nodes[b] &&
               (!net.fast[k] || !reach.from_law[net.from[k]] || reach.to_law[net.to[k]])
        safe || (qualifies[g] = false)
    end
    gone = falses(nn)
    for u in 1:nn
        g = _union_find_root(group, u)
        gone[u] = cut[u] && qualifies[g] && touched[g] != -1
    end
    any(gone) || return net
    out = _out_transitions(net, s.alive)
    grew = true
    while grew
        grew = false
        for u in 1:nn
            !gone[u] && !s.fed[u] && !law_nodes[u] && !isempty(out[u]) &&
                all(t -> gone[net.to[t]], out[u]) && (gone[u] = grew = true)
        end
    end
    keep = BitVector([!gone[net.from[t]] && !gone[net.to[t]] for t in eachindex(net.from)])
    _Network(net.nforms, net.ncopies, net.present .& .!gone, net.from[keep], net.to[keep],
             net.edge[keep], net.conformation[keep], net.forward[keep], net.fast[keep],
             net.vanishing[keep], net.change[:, keep], net.takes[keep])
end

"""
    _admission_face(am, eps, law, faces, i) -> _Verdict

Admission's face condition (`_face_condition`) on face `faces[i]` of `am`, from its
ε-network `eps` and law network `law`; a rejection carries the face index `i`. A face the
condition refuses with `_BRIDGE` passes when the condition finds no fault once more on
`eps` without the groups that cannot matter there (`_without_set_aside_groups`). A face
refused with `_PUMPED` is not run again: its stranded form's group cannot drain, so it is
never set aside, and free E still feeds it through the transition that pumps it.
"""
function _admission_face(am::Mechanisms.Mechanism, eps::_Network, law::_Network,
                         faces::Vector{UInt64}, i::Int)
    x = first(am.catalytic.substrates)
    v = _face_condition(eps, law.present, faces[i], x)
    if v.fault === _BRIDGE
        # Why a group set aside cannot change the limit, and the numerical evidence for
        # it: the developer page, "Well-posedness guard", on the set-aside rule.
        rest = _without_set_aside_groups(eps, law.present, faces[i])
        rest !== eps && _face_condition(rest, law.present, faces[i], x) === _WELL_POSED &&
            return _WELL_POSED
    end
    v.fault === _NO_FAULT ? v : _Verdict(v.fault, v.limit, i, v.conformation, v.edge)
end

# ─── Admission's messages and the checks behind their fixes ─────────────

"""The edges of a shortest path from form `a` to form `b` over the edges of `am`'s
catalytic graph that conformation `c` keeps, in path order; empty when `a == b`. Callers
pass two forms those edges join."""
function _kept_path(am::Mechanisms.Mechanism, c::Int, a::Int, b::Int)
    m = am.catalytic
    nforms = length(m.enzyme_forms)
    via = zeros(Int, nforms)
    seen = falses(nforms)
    seen[a] = true
    queue = [a]
    head = 0
    excluded = _excluded_steps(am, c)
    while !seen[b]
        head += 1
        u = queue[head]
        for k in eachindex(m.steps)
            excluded[k] && continue
            from, to = m.steps[k].from, m.steps[k].to
            v = from == u ? to : to == u ? from : 0
            (v == 0 || seen[v]) && continue
            seen[v] = true
            via[v] = k
            push!(queue, v)
        end
    end
    path = Int[]
    f = b
    while f != a
        k = via[f]
        pushfirst!(path, k)
        from, to = m.steps[k].from, m.steps[k].to
        f = from == f ? to : from
    end
    path
end

"""
Whether `am` with the kinetic group of each edge in `ks` tagged `tag` passes the guard
(`_guard`). Every fix that retags steps in admission's messages is checked this way.
"""
function _retag_passes(am::Mechanisms.Mechanism, ks, tag::Symbol)
    states = copy(am.allosteric.allosteric_states)
    states[[am.catalytic.steps[k].kinetic_group for k in ks]] .= tag
    retagged = Mechanisms.Mechanism(am.catalytic,
                                    Mechanisms.Allostery(states,
                                                         am.allosteric.oligomeric_state,
                                                         am.allosteric.regulatory_sites))
    _guard(retagged) === _WELL_POSED
end

"""
The message of admission's verdict `v` on `am`, judged with its excluded steps switched
off as `reading` says: the cause, naming each step by its DSL line (`step_text`). It
offers each fix only when its variant passes the guard: a retag through `_retag_passes`, a
steady-state fix through `_steady_state_passes`. A fault of the active conformation is
offered the fixes of the same fault in the inactive conformation, with `:OnlyI` in place
of `:OnlyA`.
"""
function _admission_message(am::Mechanisms.Mechanism, v::_Verdict,
                            reading::_Reading)
    m = am.catalytic
    c, k = v.conformation, v.edge
    conf = c == 1 ? "active" : "inactive"
    # The state that excludes steps from conformation `c`, and its text.
    excluding = c == 1 ? :OnlyI : :OnlyA
    tag = ":$excluding"
    line(j) = Mechanisms.step_text(m, j)
    lines(js) = join((Mechanisms.step_text(m, j) for j in js), "; ")
    # The steps `js` as "step …" or "steps …", and the clause that says they vanish.
    steps(js) = (length(js) == 1 ? "step " : "steps ") * lines(js)
    vanish(js) = length(js) == 1 ? "that step vanishes" : "those steps vanish"
    # The fixes `found` as one sentence of alternatives; nothing when there are none.
    offer(found) = isempty(found) ? "" : " " * uppercasefirst(join(found, ", or ")) * "."
    if v.fault === _FAST_STEP_PINNED
        # The guard reports the fault only for a step inside one island, whose sides the
        # steps conformation `c` keeps join.
        joins = _kept_path(am, c, m.steps[k].from, m.steps[k].to)
        msg = "Mechanism: in the $conf conformation the rapid-equilibrium step " *
              "$(line(k)) is tagged $tag, but the steps that conformation keeps still " *
              "join its two sides. Thermodynamics then fixes its equilibrium constant, " *
              "and a step keeps its kind in both conformations, so it cannot be switched " *
              "off."
        offered = String[]
        _steady_state_passes(am, k) && push!(offered, "make its kinetic group steady-state")
        _retag_passes(am, joins, excluding) &&
            push!(offered, "tag $tag the steps that join its sides: $(lines(joins))")
        return msg * offer(offered)
    elseif v.fault === _TOO_MANY_READINGS
        who = c == 0 ? "the two conformations together" : "the $conf conformation"
        most = c == 0 ? _MAX_READING_PAIRS : _MAX_CHOICES_PER_CONFORMATION
        return "Mechanism: $who can switch off the tagged steps in more " *
               "than $most ways, too many for the guard to search, and none of the three " *
               "ways it tries first fits the rule for switching them off, so it cannot " *
               "show that the rate law is the limit."
    elseif v.fault in (_EPS_NONUNIQUE, _FILLED, _PUMPED, _BRIDGE)
        eps = _eps_network(am, _re_flags(m), reading)
        law = _law_network(am)
        faces = _faces(m)
        face = faces[v.face]
        phrase = _face_phrase(m, face)
        s = _face_sets(eps, law.present, face)
        if v.fault === _EPS_NONUNIQUE
            outside = _class_outside_law(_closed_classes(eps, s.alive), law.present)
            forms = [String(m.enzyme_forms[_form_of(eps, u)]) *
                     _conformation_suffix(eps, u) for u in outside]
            return "Mechanism: $(phrase)the forms $(join(forms, ", ")) " *
                   "hold the enzyme they receive and are reached only through tagged " *
                   "steps, so the steady state depends on where the enzyme starts."
        elseif v.fault === _FILLED
            # The step's vanishing direction fills its target; the forms the target reaches
            # along fast transitions that do not vanish leave only through slower ones.
            filling = only(t for t in eachindex(eps.from)
                           if eps.edge[t] == k && eps.conformation[t] == c &&
                              eps.vanishing[t])
            filled = eps.to[filling]
            pool = _reached(_successors(eps, s.n0 .& eps.fast), [filled])
            ways = unique!([eps.edge[t] for t in eachindex(eps.from)
                            if s.alive[t] && pool[eps.from[t]] && !pool[eps.to[t]]])
            name = String(m.enzyme_forms[_form_of(eps, filled)])
            pooled = [String(m.enzyme_forms[_form_of(eps, u)]) for u in findall(pool)]
            leave = length(pooled) == 1 ? "$name leaves" :
                    "the forms $(join(pooled, ", ")) leave"
            msg = "Mechanism: $(phrase)in the $conf conformation the $tag " *
                  "rapid-equilibrium step $(line(k)) runs one way into $name, at a rate " *
                  "that vanishes yet grows with the rapid-equilibrium speed, and $leave " *
                  "only through the slower $(steps(ways)), so $name fills faster than it " *
                  "empties and its weight depends on how fast the tagged step vanishes."
            return _steady_state_passes(am, k) ?
                   msg * " Make the kinetic group of $(line(k)) steady-state." : msg
        elseif v.fault === _PUMPED
            basin = [u for u in findall(s.stranded) if _copy_of(eps, u) == c]
            exits = unique!([eps.edge[t] for t in eachindex(eps.from)
                             if s.alive[t] && eps.vanishing[t] && eps.from[t] in basin &&
                                !s.stranded[eps.to[t]]])
            out = "which leave only through the vanishing $(steps(exits)) while " *
                  "nothing runs back the way it came, so admission cannot show that " *
                  "their weight does not depend on how fast $(vanish(exits))"
            from, to = m.steps[k].from, m.steps[k].to
            inside = s.stranded[_node(eps, from, c)] && s.stranded[_node(eps, to, c)]
            return "Mechanism: $(phrase)$(line(k)) " *
                   (inside ? "runs one way among" : "feeds") * " the forms " *
                   "$(join([String(m.enzyme_forms[_form_of(eps, u)])
                            for u in basin], ", ")) of " *
                   "the $conf conformation, $out."
        elseif c == 0
            ks = _bridging_edges(eps, s.fed .& .!s.stranded, s.n0, s.alive .& eps.vanishing)
            return "Mechanism: $(phrase)the tagged steps $(lines(ks)) of both " *
                   "conformations together bridge forms that each conformation's other " *
                   "steps keep apart, so admission cannot show that the limit does not " *
                   "depend on how fast those steps vanish."
        end
        nodes = s.fed .& .!s.stranded
        (; van) = _copy_transitions(eps, s.alive, nodes, c)
        ks = _bridging_edges(eps, nodes, s.n0, van)
        bridge = steps(ks) * (length(ks) == 1 ? " bridges" : " bridge")
        msg = "Mechanism: $(phrase)the $tag $bridge forms of the $conf " *
              "conformation that its other steps keep apart, so admission cannot show " *
              "that the limit does not depend on how fast $(vanish(ks))."
        return _steady_state_passes(am, k) ?
               msg * " Make the kinetic group of $(line(k)) steady-state." : msg
    elseif v.fault === _NO_ACTIVE_TURNOVER
        ks = findall(_excluded_steps(am, 1))
        msg = "Mechanism: without its :OnlyI steps $(lines(ks)) the active " *
              "conformation turns no substrate into product, but it must catalyze the " *
              "reaction."
        # A closure over `ks`, which other branches assign too, would box it.
        untag = _retag_passes(am, ks, :EqualAI) && _retag_passes(am, ks, :NonequalAI)
        return untag ? msg * " Tag their kinetic groups :EqualAI or :NonequalAI." : msg
    else
        error("_admission_message: $(v.fault) is no admission fault")
    end
end

# ─── A form's conformation in a message ─────────────────────────────────

"""The suffix a message writes after the form of node `u` of the two-copy network `net`:
` (active)` or ` (inactive)`."""
_conformation_suffix(net::_Network, u::Int) =
    _copy_of(net, u) == 1 ? " (active)" : " (inactive)"
