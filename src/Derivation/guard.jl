# ABOUTME: The well-posedness guard the macros, parser, exits and derivation entries run:
# ABOUTME: its verdicts, reaction networks and faces, classes 1 and 2, and class 3 per face.

"""
    check_valid_rate_equation(m; throw = true) -> Bool

Whether `m`'s rate law is the unique limit of mass-action kinetics, as its rapid-equilibrium
steps speed up and its tagged constants vanish, on every face where data sit (spec
section 8): the guard (`_guard`) finds no fault. With `throw = true` a fault is an error
whose message (`_guard_message`) names each step at fault by its DSL line, and the result
is `true`; with `throw = false` the result is the verdict. It reads the reaction graph
through masks and builds nothing. The mechanism macros and the parser call it, and every
derivation entry runs its check (`_derivation_input`, `_well_posed`). It reads the verdict
through the memo (`_Memo`). Identify's base tier, the enumerator's exits and its flip,
split and tagging searches keep what it accepts with `throw = false`. The constructors and
the Sig decoder do not call it, so code that builds a mechanism internally never trips a
check meant for the mechanisms it returns. The keyword shadows `Base.throw` in the body.
"""
function check_valid_rate_equation(m::_AnyMech; throw::Bool = true)
    v = _memo_verdict(m)
    throw && _check_verdict(m, v)
    v.fault === _NO_FAULT
end

"""Error with the guard's message (`_guard_message`) unless `m`'s verdict `v` holds no
fault."""
function _check_verdict(m::_AnyMech, v)
    v.fault === _NO_FAULT || error(_guard_message(m, v))
    nothing
end

"""`m`, once its guard verdict, taken afresh rather than through the memo (`_Memo`),
holds no fault (`_check_verdict`)."""
_well_posed(m) = (_check_verdict(m, _guard(m)); m)

"""The concrete mechanism the compiled type `M` encodes (`Mechanism(em)`,
`AllostericMechanism(aem)`), checked as `check_valid_rate_equation` checks it, without the
memo (`_well_posed`): the input of every derivation entry."""
_derivation_input(@nospecialize(M::Type{<:Mechanisms.EnzymeMechanism})) =
    _well_posed(Mechanisms.Mechanism(M()))
_derivation_input(@nospecialize(M::Type{<:Mechanisms.AllostericEnzymeMechanism})) =
    _well_posed(Mechanisms.AllostericMechanism(M()))

# ─── The guard: verdicts, networks, faces, classes 1 and 2 ─────────────

"""
What the guard finds wrong with a mechanism (spec section 8): nothing (`_NO_FAULT`); no
cycle turns the reaction over (`_NO_CYCLE`, class 1); a cycle of rapid-equilibrium steps
turns it over (`_RE_TURNOVER`, class 2); on a face (class 3), a second closed class
(`_NONUNIQUE`), a rapid-equilibrium segment with no lowest form (`_BOTTOMLESS`), or a law
the limit test rejects where procedure A stopped at a trap fed on the face (`_TRAP`), a
segment whose lowest forms split (`_SPLIT`), a law segment chain with several closed
classes (`_LAW_MULTI`) or an aggregated chain with several (`_MULTILEVEL`); in admission's
interior (class 4, spec section 7), a tagged binding whose two sides the conformation's
other steps join (`_NO_LIMIT`), tagged bindings that close a directed cycle
(`_ARC_CYCLE`), forms free E's component does not reach along the tagged bindings
(`_ISLAND`), or an active conformation that no longer turns the reaction over
(`_NO_ACTIVE_TURNOVER`); on a face, in admission's face condition, a second closed class
of the ε-network (`_EPS_NONUNIQUE`), a stranded basin entered by, or holding, an alive step
whose reverse is dead, or entered from outside it and the law's forms (`_PUMPED`), or
tagged steps that shorten the fed forms' cheapest spanning in-tree where no potential
absorbs their turnover (`_BRIDGE`).
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
    _NO_LIMIT
    _ARC_CYCLE
    _ISLAND
    _EPS_NONUNIQUE
    _PUMPED
    _BRIDGE
    _NO_ACTIVE_TURNOVER
end

"""The faults admission reports (spec section 7); `_admission_message` writes their text."""
const _ADMISSION_FAULTS = (_NO_LIMIT, _ARC_CYCLE, _ISLAND, _NO_ACTIVE_TURNOVER,
                           _EPS_NONUNIQUE, _PUMPED, _BRIDGE)

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
counts the copies of metabolite name `n` (`_guard_names`) that transition `t` takes up
less those it gives off; bit `n - 1` of `takes[t]` is set when it takes up name `n`.
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

"""The metabolite names of `m`, sorted: the guard's name `n` (a row of a network's
`change`, a bit of its `takes` and of a face) is the `n`th. A reactant's modifier copy
reads as its reactant, as they share a concentration."""
_guard_names(m::Mechanisms.Mechanism) = sort!(Mechanisms.metabolites(m))

"""The guard's name (`_guard_names`) of `m`'s first substrate, whose uptake counts
turnover."""
_turnover_name(m::Mechanisms.Mechanism) =
    searchsortedfirst(_guard_names(m), first(Mechanisms.substrates(m)))

"""The node of form `f` in copy `c` of `net`."""
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
function _network(m::Mechanisms.Mechanism, present::BitVector, use::BitMatrix,
                  fast::BitMatrix, vanish_forward::BitMatrix, vanish_reverse::BitMatrix)
    nf, ncopies = length(Mechanisms.species(m)), size(use, 2)
    names = _guard_names(m)
    nt = 2 * count(use) + (ncopies == 2 ? 2 : 0)
    from, to, edge, conformation = zeros(Int, nt), zeros(Int, nt), zeros(Int, nt),
                                   zeros(Int, nt)
    forward, isfast, vanishing = falses(nt), falses(nt), falses(nt)
    change, takes = zeros(Int, length(names), nt), zeros(UInt64, nt)
    t = 0
    for c in 1:ncopies, k in Mechanisms.steps(m)
        use[k, c] || continue
        a, b = (c - 1) * nf .+ Mechanisms.ends(m, k)
        up, off = Mechanisms.taken_up(m, k), Mechanisms.given_off(m, k)
        for (u, v, gain, loss, gone) in ((a, b, up, off, vanish_forward[k, c]),
                                         (b, a, off, up, vanish_reverse[k, c]))
            t += 1
            from[t], to[t], edge[t], conformation[t] = u, v, k, c
            forward[t], isfast[t], vanishing[t] = u == a, fast[k, c], gone
            for (n, x) in enumerate(names)
                taken = Mechanisms.copies(gain, x)
                change[n, t] = taken - Mechanisms.copies(loss, x)
                taken > 0 && (takes[t] |= UInt64(1) << (n - 1))
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
The law network L of `m`: a plain mechanism's forms and edges as they are, its
rapid-equilibrium flags overridden by `re`; an allosteric mechanism's two copies, copy `c`
holding the steps of its conformation (`conformation_steps`) and the forms they touch
(`present_species`), free E always, joined by the flip at free E.
"""
function _law_network(m::Mechanisms.Mechanism,
                      re::BitVector = Mechanisms.rapid_equilibrium_steps(m))
    n = length(Mechanisms.steps(m))
    fast = falses(n, 1)
    fast[:, 1] .= re
    _network(m, trues(length(Mechanisms.species(m))), trues(n, 1), fast, falses(n, 1),
             falses(n, 1))
end
function _law_network(am::Mechanisms.AllostericMechanism)
    m = Mechanisms.catalytic(am)
    re = Mechanisms.rapid_equilibrium_steps(m)
    n, nf = length(Mechanisms.steps(m)), length(Mechanisms.species(m))
    use, fast, present = falses(n, 2), falses(n, 2), falses(2 * nf)
    for (c, conformation) in enumerate(_CONFORMATIONS)
        mask = Mechanisms.conformation_steps(am, conformation)
        use[:, c] .= mask
        fast[:, c] .= mask .& re
        nodes = (c - 1) * nf .+ (1:nf)
        present[nodes] .= Mechanisms.present_species(m, mask)
        present[nodes[1]] = true
    end
    _network(m, present, use, fast, falses(n, 2), falses(n, 2))
end

"""
The ε-network N of `m`, its rapid-equilibrium flags overridden by `re`: every form and
edge in each copy. In the copy that excludes an edge (`excludes`), a rapid-equilibrium
plain binding stays fast while its uptake vanishes and its release keeps its rate; any
other excluded edge is slow and vanishes both ways. A plain mechanism's ε-network is its
law network.
"""
_eps_network(m::Mechanisms.Mechanism,
             re::BitVector = Mechanisms.rapid_equilibrium_steps(m)) = _law_network(m, re)
function _eps_network(am::Mechanisms.AllostericMechanism,
                      re::BitVector =
                          Mechanisms.rapid_equilibrium_steps(Mechanisms.catalytic(am)))
    m = Mechanisms.catalytic(am)
    n = length(Mechanisms.steps(m))
    fast, vanish_forward, vanish_reverse = falses(n, 2), falses(n, 2), falses(n, 2)
    for c in 1:2, k in Mechanisms.steps(m)
        excluded = Mechanisms.excludes(am, k, _CONFORMATIONS[c])
        plain = Mechanisms.step_kind(m, k) === :binding
        fast[k, c] = re[k] && (!excluded || plain)
        vanish_forward[k, c] = excluded
        vanish_reverse[k, c] = excluded && !plain
    end
    _network(m, trues(2 * length(Mechanisms.species(m))), trues(n, 2), fast, vanish_forward,
             vanish_reverse)
end

"""
The faces the guard checks (spec section 8) of `m`, each a mask of the metabolite names it
zeroes (bit `n - 1` for name `n` of `_guard_names`), in order: the interior (0), every
product, every substrate, each substrate and then each product alone on a side of two or
more; then each of these with every modifier-only name too. A reactant's `::Mod` copy
shares its reactant's name and so is zeroed with it; a regulatory-site ligand is zeroed
only when it also binds the catalytic site. No face zeroes a substrate together with a
product.
"""
function _faces(m::Mechanisms.Mechanism)
    names = _guard_names(m)
    substrates, products = Mechanisms.substrates(m), Mechanisms.products(m)
    bit(x::Symbol) = UInt64(1) << (searchsortedfirst(names, x) - 1)
    side(xs::Vector{Symbol}) = reduce(|, (bit(x) for x in xs); init = UInt64(0))
    base = UInt64[0, side(products), side(substrates)]
    length(substrates) >= 2 && append!(base, bit.(substrates))
    length(products) >= 2 && append!(base, bit.(products))
    mods = side(filter(x -> !(insorted(x, substrates) || insorted(x, products)), names))
    mods == 0 ? base : [base; base .| mods]
end

"""Whether face `face` (`_faces`) zeroes metabolite name `n` of `_guard_names`."""
_zeroes(face::UInt64, n::Int) = (face >> (n - 1)) & 1 == 1

"""The metabolite names that face `face` of `m` zeroes (`_faces`), sorted."""
_face_names(m::Mechanisms.Mechanism, face::UInt64) =
    [x for (n, x) in enumerate(_guard_names(m)) if _zeroes(face, n)]

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
(`segment_walk`) from each present node no earlier walk reached: `segment[u]`, node
`u`'s (0 for a node `net` does not hold); `segments`, each one's nodes in walk order;
`offset[n, u]`, the copies of name `n` (`_guard_names`) the walk takes up from the
segment's first node to `u`; and `cycle`, the transitions of the first fast cycle that
takes up some name net, which the first transition whose offsets disagree closes with its
tree path, in cycle order.
When `cycle` is empty, each offset is a potential: it does not depend on the walk.
"""
function _fast_segments(net::_Network)
    w = Mechanisms.segment_walk(net.present, _out_transitions(net, net.fast), net.to,
                                net.change)
    t = w.disagreeing
    cycle = t == 0 ? Int[] : _tree_cycle(net, w.parent, net.from[t], net.to[t], t)
    (; segment = w.segment_of, w.segments, offset = w.offsets, cycle)
end

"""The transitions of the first cycle of fast transitions of `net` that takes up some
metabolite name net (class 2: the rate is infinite), in cycle order; empty when none does
(`_fast_segments`). Class 2 runs once, at the interior, since a fast cycle that runs on a
face runs at the interior too."""
_re_turnover_transitions(net::_Network) = _fast_segments(net).cycle

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

The well-posedness guard (spec section 8) on `m`, plain or allosteric: class 1 (a cycle
turns the reaction over, `has_turnover_cycle`); class 2 (no cycle of fast transitions
does, `_re_turnover_transitions`, on the ε-network, where an excluded rapid-equilibrium
binding stays fast while its affinity vanishes); for an allosteric mechanism, admission's
interior (`_admission_interior`); then on every face class 3 (`_check_face`, on the law
network) and, for an allosteric mechanism, admission's face condition (`_admission_face`,
on the ε-network), the first fault winning. `_WELL_POSED` when every check passes. It
reads the reaction graph through masks and builds nothing.
"""
function _guard(mech::_AnyMech)
    m = Mechanisms.catalytic(mech)
    Mechanisms.has_turnover_cycle(m, trues(length(Mechanisms.steps(m)))) ||
        return _Verdict(_NO_CYCLE, _LIMIT_UNTESTED, 0, 0, 0)
    eps = _eps_network(mech)
    v = _re_turnover_verdict(eps)
    v === _WELL_POSED || return v
    # A plain mechanism's ε-network is its law network, and admission reads no tag of it.
    mech isa Mechanisms.AllostericMechanism || return _faces_verdict(mech, eps)
    v = _admission_interior(mech)
    v === _WELL_POSED || return v
    law, faces = _law_network(mech), _faces(m)
    _first_face_fault(faces) do i
        face_verdict = _check_face(mech, law, faces, i)
        face_verdict.fault === _NO_FAULT ? _admission_face(mech, eps, law, faces, i) :
            face_verdict
    end
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

# ─── Class 3: procedure A on every face ────────────────────────────────

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

"""The fast classes (`class_of`) that the alive fast transitions (`fast_next`) reach from
node `w`, its own included."""
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
    _face_structure(net, face) -> NamedTuple

The structure procedure A reads on face `face` of the law network
`net`: `alive`, the transitions that take up no zeroed name; `closed`, the closed classes
of the alive ones; the rapid-equilibrium segments (`segment` per node, `segments`
breadth-first) and each node's uptake of every name from its segment's first node
(`offset`, `_fast_segments`); `bottomless`, the first segment with no lowest node (0 when
none), a lowest node being one at the segment's componentwise least uptake of the zeroed
names, whose weight stays positive on the face; the fast classes (closed classes of the
alive fast transitions, a node with no alive fast exit its own), `class_of` each node's
(0 for a node that drains elsewhere) and `trap`, the classes of no lowest node;
`aggregated`, the closed classes of the chain over fast classes that an alive slow
transition and then alive fast paths make; and `law_classes`, the closed classes of the
law's segment chain, its alive slow transitions out of lowest nodes.
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
    fast_alive = alive .& net.fast
    fast_classes = _closed_classes(net, fast_alive)
    class_of = zeros(Int, nn)
    for (k, C) in enumerate(fast_classes), u in C
        class_of[u] = k
    end
    all(C -> all(u -> lowest[u] == lowest[C[1]], C), fast_classes) ||
        error("_face_structure: a fast class holds lowest and other forms")
    trap = BitVector([!lowest[C[1]] for C in fast_classes])
    fast_next = _successors(net, fast_alive)
    aggregated_next = [Int[] for _ in fast_classes]
    law_next = [Int[] for _ in segments]
    for t in 1:nt
        alive[t] && !net.fast[t] || continue
        u, v = net.from[t], net.to[t]
        if class_of[u] > 0
            for k in _classes_reached(fast_next, class_of, v)
                k == class_of[u] || push!(aggregated_next[class_of[u]], k)
            end
        end
        lowest[u] && segment[u] != segment[v] && push!(law_next[segment[u]], segment[v])
    end
    aggregated = _bottom_sccs(aggregated_next, trues(length(fast_classes)))
    law_classes = _bottom_sccs(law_next, trues(length(segments)))
    (; alive, closed, segment, segments, offset, bottomless, fast_classes, class_of, trap,
     aggregated, law_classes)
end

"""
Where procedure A settles a face from its structure `s`
(`_face_structure`): a second closed class (`_NONUNIQUE`); a segment with no lowest form
(`_BOTTOMLESS`); `_NO_FAULT` by Case I (every segment holds one fast class, the law's
segment chain one closed class); else the stage at which the limit test decides: the
law's chain has several closed classes (`_LAW_MULTI`), the aggregated chain has
(`_MULTILEVEL`), its closed class holds a trap (`_TRAP`), or none of these (`_SPLIT`,
where a segment it meets may split).
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
"""The form of node `u` of `net`."""
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

Class 3 (spec section 8) on face `faces[i]` of `mech`, whose law
network is `net`: procedure A rejects a second closed class or a bottomless segment and
accepts a face whose structure fixes the limit (`_face_stage`); otherwise the limit test
(`_face_limit_check`, seeded by the content hash and the face) decides, and a rejection
carries the stage it stopped at, the edge that feeds the trap or joins it to its segment
(`_trap_feeder`) and the split's first joining edge. A bottomless segment, a split or a
trap carries the copy its forms share (`_common_copy`). A face the limit test accepts
reports `_LIMIT_EQUAL`.
"""
function _check_face(mech::_AnyMech, net::_Network, faces::Vector{UInt64},
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

"""The first class-3 fault of `mech` on its law network `net`, face by face in the order
of `_faces`; `_WELL_POSED` when every face passes."""
function _faces_verdict(mech::_AnyMech, net::_Network)
    faces = _faces(Mechanisms.catalytic(mech))
    _first_face_fault(i -> _check_face(mech, net, faces, i), faces)
end

