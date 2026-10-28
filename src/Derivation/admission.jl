# ABOUTME: Admission of allosteric mechanisms: the interior (excluded bindings, arc cycles,
# ABOUTME: islands, active turnover) and the face condition on each face of the ε-network.

# ─── Admission: the interior (spec section 7) ──────────────────────────

"""The plain bindings of `am`'s catalytic graph that conformation `c` excludes, in edge
order: admission's arcs, each stored as its uptake, unbound side first."""
function _arcs(am::Mechanisms.AllostericMechanism, c::Int)
    m = Mechanisms.catalytic(am)
    [k for k in Mechanisms.steps(m)
     if Mechanisms._excludes(am, k, c) && Mechanisms.step_kind(m, k) === :binding]
end

"""The first edge of `am`'s catalytic graph that conformation `c` excludes; 0 when it
excludes none."""
_first_excluded(am::Mechanisms.AllostericMechanism, c::Int) =
    something(findfirst(k -> Mechanisms._excludes(am, k, c),
                        Mechanisms.steps(Mechanisms.catalytic(am))), 0)

"""
The arcs of one directed cycle that `arcs` close over the components `comp` (each arc an
edge of `m` read from its unbound side's component to its bound side's), in cycle order;
empty when they close none. Kahn's algorithm removes every component no cycle feeds; a
walk back along the remaining arcs then lands on a cycle.
"""
function _arc_cycle(m::Mechanisms.Mechanism, comp::Vector{Int}, arcs::Vector{Int})
    tail(k) = comp[Mechanisms.ends(m, k)[1]]
    head(k) = comp[Mechanisms.ends(m, k)[2]]
    indeg = zeros(Int, length(comp))
    for k in arcs
        indeg[head(k)] += 1
    end
    todo = [r for r in unique(comp) if indeg[r] == 0]
    while !isempty(todo)
        u = pop!(todo)
        for k in arcs
            tail(k) == u && (indeg[head(k)] -= 1) == 0 && push!(todo, head(k))
        end
    end
    start = 0
    for k in arcs
        indeg[tail(k)] > 0 && indeg[head(k)] > 0 && (start = k; break)
    end
    start == 0 && return Int[]
    function back(x)
        for k in arcs
            head(k) == x && indeg[tail(k)] > 0 && return k
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
    reverse!(cycle)
end

"""The components (`comp` labels) that free E's reaches along `arcs`, edges of `m` read
from their unbound side's component to their bound side's."""
function _reached_components(m::Mechanisms.Mechanism, comp::Vector{Int}, arcs::Vector{Int})
    next = [Int[] for _ in comp]
    for k in arcs
        a, b = Mechanisms.ends(m, k)
        push!(next[comp[a]], comp[b])
    end
    _reached(next, [comp[1]])
end

"""
The edge an island's message names, 0 when free E's component reaches every component
(`reached`): of the edges conformation `c` excludes that join a reached component to an
unreached one, the first chemistry edge (a `step_kind` other than `:binding`) that, run
toward products (`_runs_toward_products`), enters the unreached one; else the first
chemistry edge; else the first edge.
"""
function _island_edge(am::Mechanisms.AllostericMechanism, c::Int, comp::Vector{Int},
                      reached::BitVector)
    m = Mechanisms.catalytic(am)
    best, rank = 0, 4
    for k in Mechanisms.steps(m)
        Mechanisms._excludes(am, k, c) || continue
        from, to = Mechanisms.ends(m, k)
        a, b = reached[comp[from]], reached[comp[to]]
        a == b && continue
        r = Mechanisms.step_kind(m, k) === :binding ? 3 :
            a == Mechanisms._runs_toward_products(m, k) ? 1 : 2
        if r < rank
            best, rank = k, r
        end
    end
    best
end

"""
    _admission_interior(am) -> _Verdict

Admission's interior checks (spec section 7), for each conformation `c` whose tags exclude
an edge (`_excludes`), on the components of the forms its other edges join
(`_exclusion_components`) with its tagged plain bindings as arcs (`_arcs`):
- existence: no arc stays inside one component (`_NO_LIMIT`) and the arcs close no
  directed cycle (`_ARC_CYCLE`); equivalent to the exact cone test, since with every other
  constant fixed a direction's equilibrium-ratio exponents are differences of component
  potentials and only the arcs need a sign;
- uniqueness: free E's component reaches every component along the arcs, so every island
  lies behind a tagged plain binding on its bound side and the limit is deletion plus
  pruning (`_ISLAND`);
- turnover: when tags exclude edges of the active conformation, its law network still
  holds a cycle that turns the reaction over (`_NO_ACTIVE_TURNOVER`).
The verdict carries the conformation and the edge its message leads with (face 0).
"""
function _admission_interior(am::Mechanisms.AllostericMechanism)
    m = Mechanisms.catalytic(am)
    for c in 1:2
        _first_excluded(am, c) == 0 && continue
        comp = Mechanisms._exclusion_components(am, c)
        arcs = _arcs(am, c)
        for k in arcs
            a, b = Mechanisms.ends(m, k)
            comp[a] == comp[b] && return _Verdict(_NO_LIMIT, _LIMIT_UNTESTED, 0, c, k)
        end
        cycle = _arc_cycle(m, comp, arcs)
        isempty(cycle) || return _Verdict(_ARC_CYCLE, _LIMIT_UNTESTED, 0, c, cycle[1])
        k = _island_edge(am, c, comp, _reached_components(m, comp, arcs))
        k == 0 || return _Verdict(_ISLAND, _LIMIT_UNTESTED, 0, c, k)
    end
    k = _first_excluded(am, 1)
    k != 0 &&
        !Mechanisms.has_turnover_cycle(m, Mechanisms.conformation_steps(am, :active)) &&
        return _Verdict(_NO_ACTIVE_TURNOVER, _LIMIT_UNTESTED, 0, 1, k)
    _WELL_POSED
end

"""The edges of a shortest path from form `a` to form `b` over the edges of `am`'s
catalytic graph that conformation `c` keeps, in path order; empty when `a == b`. Callers
pass two forms those edges join."""
function _kept_path(am::Mechanisms.AllostericMechanism, c::Int, a::Int, b::Int)
    m = Mechanisms.catalytic(am)
    nforms = length(Mechanisms.species(m))
    via = zeros(Int, nforms)
    seen = falses(nforms)
    seen[a] = true
    queue = [a]
    head = 0
    while !seen[b]
        head += 1
        u = queue[head]
        for k in Mechanisms.steps(m)
            Mechanisms._excludes(am, k, c) && continue
            from, to = Mechanisms.ends(m, k)
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
        from, to = Mechanisms.ends(m, k)
        f = from == f ? to : from
    end
    path
end

"""
Whether `am` with the kinetic group of each edge in `ks` tagged `tag` passes the guard
(`_guard`). Every fix that retags steps in admission's messages is checked this way.
"""
function _retag_passes(am::Mechanisms.AllostericMechanism, ks, tag::Symbol)
    states = copy(Mechanisms.cat_allo_states(am))
    states[[Mechanisms.group(Mechanisms.catalytic(am), k) for k in ks]] .= tag
    _guard(Mechanisms._with(am; states)) === _WELL_POSED
end

"""
Whether `mech` with the kinetic group of edge `k` made steady-state passes the guard
(`_guard`). Every fix that makes a group steady-state in the guard's messages, but class
2's, is checked this way.
"""
function _steady_state_passes(mech::Mechanisms._AnyMech, k::Int)
    m = Mechanisms.catalytic(mech)
    re = Mechanisms._flags(m)
    re[Mechanisms.members(m, Mechanisms.group(m, k))] .= false
    _guard(Mechanisms._build(Mechanisms._raw(mech, re))) === _WELL_POSED
end

"""
The message of admission's verdict `v` on `am`: the cause, naming each step by its DSL
line (`_step_text`), and each fix only when its variant passes the guard: a retag through
`_retag_passes`, a steady-state fix through `_steady_state_passes`. In the active
conformation the only retag offered untags the `:OnlyI` groups that leave it no turnover,
since a step tagged `:OnlyI` must bind modifiers alone.
"""
function _admission_message(am::Mechanisms.AllostericMechanism, v::_Verdict)
    m = Mechanisms.catalytic(am)
    c, k = v.conformation, v.edge
    conf = c == 1 ? "active" : "inactive"
    tag = c == 1 ? ":OnlyI" : ":OnlyA"
    line(j) = Mechanisms._step_text(m, j)
    lines(js) = join((Mechanisms._step_text(m, j) for j in js), "; ")
    # The steps `js` as "step …" or "steps …", and the clause that says they vanish.
    steps(js) = (length(js) == 1 ? "step " : "steps ") * lines(js)
    vanish(js) = length(js) == 1 ? "that step vanishes" : "those steps vanish"
    untag_passes(ks) = all(t -> _retag_passes(am, ks, t), (:EqualAI, :NonequalAI))
    # The fixes whose variants pass: the groups of `untag` (`what`) tagged :EqualAI or
    # :NonequalAI, and the steps `joins` that join `sides` tagged :OnlyA.
    function fixes(untag::Vector{Int}, what::String, joins::Vector{Int}, sides::String)
        found = String[]
        untag_passes(untag) && push!(found, "tag $what :EqualAI or :NonequalAI")
        _retag_passes(am, joins, :OnlyA) &&
            push!(found, "tag :OnlyA the steps that join $sides: $(lines(joins))")
        isempty(found) ? "" : " " * uppercasefirst(join(found, ", or ")) * "."
    end
    if v.fault === _NO_LIMIT
        msg = "AllostericMechanism: in the $conf conformation, steps outside $tag " *
              "groups still join the two sides of the $tag binding $(line(k)), so " *
              "thermodynamics keeps its affinity finite and the binding has no limit " *
              "there."
        c == 2 || return msg
        return msg * fixes([k], "its group", _kept_path(am, c, Mechanisms.ends(m, k)...),
                           "its sides")
    elseif v.fault === _ARC_CYCLE
        cycle = _arc_cycle(m, Mechanisms._exclusion_components(am, c), _arcs(am, c))
        msg = "AllostericMechanism: the $tag bindings $(lines(cycle)) close a cycle in " *
              "the $conf conformation: steps outside $tag groups join each one's bound " *
              "side to the next one's unbound side, so thermodynamics forbids their " *
              "affinities to vanish together."
        c == 2 || return msg
        joins = Int[]
        for (i, j) in enumerate(cycle)
            next = cycle[mod1(i + 1, length(cycle))]
            append!(joins, _kept_path(am, c, Mechanisms.ends(m, j)[2],
                                      Mechanisms.ends(m, next)[1]))
        end
        return msg * fixes(cycle, "them", unique!(joins), "them")
    elseif v.fault === _ISLAND
        comp = Mechanisms._exclusion_components(am, c)
        reached = _reached_components(m, comp, _arcs(am, c))
        from, to = Mechanisms.ends(m, k)
        lost = reached[comp[from]] ? comp[to] : comp[from]
        fixable = c == 2 && untag_passes([k])
        msg = "AllostericMechanism: in the $conf conformation the forms " *
              join([Mechanisms._form_text(m, f)
                    for f in Mechanisms.species(m) if comp[f] == lost], ", ") *
              " join free E only through $tag steps that are not bindings toward them, " *
              "such as $(line(k)), so their weight depends on how fast those steps vanish."
        msg = c == 2 && fixable ?
              msg * " Tag the group of $(line(k)) :EqualAI or :NonequalAI." : msg
        # A ping-pong whose inactive conformation runs neither half-reaction strands the
        # covalent forms between the halves; the roadmap lists what would hold them.
        neither = c == 2 && all(j -> Mechanisms.step_kind(m, j) === :binding ||
                                     Mechanisms._excludes(am, j, 2),
                                Mechanisms.steps(m))
        covalent = any(f -> comp[f] == lost && !all(isempty, Mechanisms.residual(m, f)),
                       Mechanisms.species(m))
        return neither && covalent ?
               msg * " An inactive conformation that runs neither half-reaction of a " *
                     "ping-pong is not supported yet (see the Roadmap page)." : msg
    elseif v.fault in (_EPS_NONUNIQUE, _PUMPED, _BRIDGE)
        eps, law = _eps_network(am), _law_network(am)
        faces = _faces(m)
        face = faces[v.face]
        phrase = _face_phrase(m, face)
        s = _face_sets(eps, law.present, face)
        if v.fault === _EPS_NONUNIQUE
            outside = _class_outside_law(_closed_classes(eps, s.alive), law.present)
            forms = [Mechanisms._form_text(m, _form_of(eps, u)) *
                     _conformation_suffix(eps, u) for u in outside]
            return "AllostericMechanism: $(phrase)the forms $(join(forms, ", ")) " *
                   "hold the enzyme they receive and are reached only through tagged " *
                   "steps, so the steady state depends on where the enzyme starts."
        elseif v.fault === _PUMPED
            basin = [u for u in findall(s.stranded) if _copy_of(eps, u) == c]
            exits = unique!([eps.edge[t] for t in eachindex(eps.from)
                             if s.alive[t] && eps.vanishing[t] && eps.from[t] in basin &&
                                !s.stranded[eps.to[t]]])
            out = "which leave only through the vanishing $(steps(exits)) while " *
                  "nothing runs back the way it came, so admission cannot show that " *
                  "their weight does not depend on how fast $(vanish(exits))"
            from, to = Mechanisms.ends(m, k)
            inside = s.stranded[_node(eps, from, c)] && s.stranded[_node(eps, to, c)]
            return "AllostericMechanism: $(phrase)$(line(k)) " *
                   (inside ? "runs one way among" : "feeds") * " the forms " *
                   "$(join([Mechanisms._form_text(m, _form_of(eps, u))
                            for u in basin], ", ")) of " *
                   "the $conf conformation, $out."
        elseif c == 0
            ks = _bridging_edges(eps, s.fed .& .!s.stranded, s.n0, s.alive .& eps.vanishing)
            return "AllostericMechanism: $(phrase)the tagged steps $(lines(ks)) of both " *
                   "conformations together bridge forms that each conformation's other " *
                   "steps keep apart, so admission cannot show that the limit does not " *
                   "depend on how fast those steps vanish."
        end
        nodes = s.fed .& .!s.stranded
        (; van) = _copy_transitions(eps, s.alive, nodes, c)
        ks = _bridging_edges(eps, nodes, s.n0, van)
        bridge = steps(ks) * (length(ks) == 1 ? " bridges" : " bridge")
        msg = "AllostericMechanism: $(phrase)the $tag $bridge forms of the $conf " *
              "conformation that its other steps keep apart, so admission cannot show " *
              "that the limit does not depend on how fast $(vanish(ks))."
        # Offer the steady-state fix only where its variant passes the guard: R2funnel
        # keeps its bridge with its tagged S binding steady-state (R2funnel_SS).
        return _steady_state_passes(am, k) ?
               msg * " Make the group of $(line(k)) steady-state." : msg
    elseif v.fault === _NO_ACTIVE_TURNOVER
        ks = [j for j in Mechanisms.steps(m) if Mechanisms._excludes(am, j, 1)]
        msg = "AllostericMechanism: without its :OnlyI steps $(lines(ks)) the active " *
              "conformation turns no substrate into product, and only the active " *
              "conformation may catalyze."
        return untag_passes(ks) ? msg * " Tag their groups :EqualAI or :NonequalAI." : msg
    else
        error("_admission_message: $(v.fault) is no admission fault")
    end
end

# ─── Admission: the face condition ─────────────────────────────────────

"""
    _face_sets(net, law_nodes, face) -> (; alive, fed, n0, stranded)

What admission's face condition reads on face `face` of the ε-network `net`, whose law
network holds the nodes `law_nodes`: the transitions alive there (taking up no zeroed
name), the nodes free E's active copy reaches along them, the alive transitions that do
not vanish, and the stranded basin: fed nodes outside the law network that reach it along
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
their tags, and `held`, those T0+ reads (not vanishing, or fast)."""
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
`arcs`, each `(from, to, cost)`, by Chu–Liu/Edmonds: take each node's cheapest entering
arc, contract the cycles they close and repeat; -1 when some node has no entering arc.
"""
function _min_arborescence(n::Int, arcs::Vector{NTuple{3, Int}}, root::Int)
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

MINSLOW: the fewest slow transitions in a spanning in-tree of the nodes `nodes` over the
transitions `keep` flags (cost 1 slow, 0 fast), any root; -1 when no in-tree spans them. A
virtual root joined to every node at a cost above any tree's turns the in-tree into
Chu–Liu/Edmonds's out-arborescence (`_min_arborescence`) on reversed transitions; a cost of
two such joins means no single root serves.
"""
function _fewest_slow(net::_Network, nodes::BitVector, keep::BitVector)
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

T0+: whether the uptake of name `x` (a `change` row, `_turnover_name`) is a potential on
every strongly connected part (`_strong_components`) of the transitions `keep` flags, so
no cycle of them turns the reaction over.
"""
function _potential_on_sccs(net::_Network, keep::BitVector, x::Int)
    comp = _strong_components(_successors(net, keep), net.present)
    inside = BitVector([keep[t] && comp[net.from[t]] == comp[net.to[t]]
                        for t in eachindex(net.from)])
    n = length(net.present)
    !Mechanisms._has_unbalanced_cycle(net.from, net.to, view(net.change, x, :), inside,
                                      zeros(Int, n), zeros(Int, n))
end

"""Whether the transitions `extra` shorten the cheapest spanning in-tree of `nodes` over
`base` (`_fewest_slow`): R1's failure for one conformation's vanishing transitions, joint
R1's for both conformations'."""
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
marks, else the second: the class an `_EPS_NONUNIQUE` verdict names. Whenever class 3 has
accepted the face such a class exists, because a closed class of the ε-network that holds
a law node holds a closed class of the law network, and class 3 leaves that network only
one. The fallback arises only from `_admission` on a face that class 3 rejects.
"""
function _class_outside_law(closed::Vector{Vector{Int}}, law_nodes::BitVector)
    for cls in closed
        any(view(law_nodes, cls)) || return cls
    end
    closed[2]
end

"""How a message names face `face` of `m`: nothing for the interior, else `at zero `, its
zeroed names (`_face_names`) and `place`, then a comma."""
_face_phrase(m::Mechanisms.Mechanism, face::UInt64, place::String = "") =
    face == 0 ? "" :
    "at zero " * join(String.(_face_names(m, face)), ", ", " and ") * place * ", "

"""
    _face_condition(net, law_nodes, face, x) -> _Verdict

Admission's face condition, on face `face` of the ε-network `net` of an allosteric
mechanism, whose law network holds the nodes `law_nodes`; `x` is the `change` row of its
first substrate (`_turnover_name`). In order:
- `_EPS_NONUNIQUE`: the alive transitions leave two closed classes, so at every ε the
  steady state depends on where the enzyme starts; the verdict names the closed class
  outside the law network (`_class_outside_law`);
- `_PUMPED`: a stranded basin (`_face_sets`) is entered by, or holds, an alive transition
  whose reverse is dead, or is entered from outside it and the law network, so admission
  cannot show that its weight does not depend on how fast its exits vanish; a basin
  entered only reversibly from the law network is sealed: its weight vanishes, and it
  leaves the fed nodes;
- `_BRIDGE`: in some conformation the vanishing transitions shorten the cheapest spanning
  in-tree of the fed nodes (R1 fails, `_shortens_in_tree`) and the first substrate's uptake
  is no potential on the strongly connected parts of its alive transitions that are not
  vanishing or are fast (T0+ fails, `_potential_on_sccs`), led by the first of its
  vanishing steps on a cheapest tree (`_bridging_edges`); or both conformations pass R1 and
  their vanishing transitions together shorten it (joint R1, no conformation or step).
A face whose fed nodes no in-tree spans passes: class 3 decides it. The verdict leaves the
face index 0 (`_admission_face` sets it).
"""
function _face_condition(net::_Network, law_nodes::BitVector, face::UInt64, x::Int)
    s = _face_sets(net, law_nodes, face)
    closed = _closed_classes(net, s.alive)
    length(closed) > 1 &&
        return _Verdict(_EPS_NONUNIQUE, _LIMIT_UNTESTED, 0,
                        _copy_of(net, _class_outside_law(closed, law_nodes)[1]), 0)
    for t in eachindex(net.from)
        u, w = net.from[t], net.to[t]
        s.alive[t] && s.fed[u] && s.stranded[w] || continue
        s.alive[isodd(t) ? t + 1 : t - 1] && (s.stranded[u] || law_nodes[u]) && continue
        return _Verdict(_PUMPED, _LIMIT_UNTESTED, 0, net.conformation[t], net.edge[t])
    end
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
    _admission_face(am, eps, law, faces, i) -> _Verdict

Admission's face condition (`_face_condition`) on face `faces[i]` of `am`, from its
ε-network `eps` and law network `law`; a rejection carries the face index `i`.
"""
function _admission_face(am::Mechanisms.AllostericMechanism, eps::_Network, law::_Network,
                         faces::Vector{UInt64}, i::Int)
    v = _face_condition(eps, law.present, faces[i],
                        _turnover_name(Mechanisms.catalytic(am)))
    v.fault === _NO_FAULT ? v : _Verdict(v.fault, v.limit, i, v.conformation, v.edge)
end

"""
    _admission(am) -> _Verdict

Admission alone (spec section 7): the interior checks (`_admission_interior`), then the
face condition on every face (`_admission_face`), class 3 not run.
"""
function _admission(am::Mechanisms.AllostericMechanism)
    v = _admission_interior(am)
    v.fault === _NO_FAULT || return v
    eps, law = _eps_network(am), _law_network(am)
    faces = _faces(Mechanisms.catalytic(am))
    _first_face_fault(i -> _admission_face(am, eps, law, faces, i), faces)
end
