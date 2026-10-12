# ABOUTME: Two-conformation mass-action oracle for allosteric mechanisms in which only free E
# ABOUTME: flips; :OnlyA steps deleted in the inactive state, or kept there scaled by ε.
using EnzymeRates, LinearAlgebra, Random, Printf
const ER = EnzymeRates
setprecision(BigFloat, 512)

"""The mechanism the DSL text `src` builds, or the error message it raises."""
function build(src::String)
    try
        Core.eval(Main, Meta.parse(src))
    catch err
        err isa LoadError && (err = err.error)
        sprint(showerror, err)
    end
end

"""Evaluate a dependent-parameter right-hand side (Symbol, number or call Expr)."""
function evalx(x, vals)
    x isa Symbol && return vals[x]
    x isa Number && return BigFloat(x)
    x isa Expr && x.head === :call || error("unexpected expression $x")
    f, a = x.args[1], [evalx(y, vals) for y in x.args[2:end]]
    f === :* && return prod(a)
    f === :+ && return sum(a)
    f === :- && return length(a) == 1 ? -a[1] : a[1] - a[2]
    f === :/ && return a[1] / a[2]
    f === :// && return a[1] / a[2]
    f === :^ && return a[1]^a[2]
    f === :sqrt && return sqrt(a[1])
    f === :inv && return inv(a[1])
    error("unknown operation $f in $x")
end

"""Every parameter value: the fitted ones, Keq, and each dependent of the package's
combined solve (`_dependent_param_exprs`), resolved until none is left."""
function all_values(am, fitted::Dict{Symbol, BigFloat}, Keq::BigFloat)
    dep, _ = ER._dependent_param_exprs(am)
    vals = copy(fitted)
    vals[:Keq] = Keq
    todo = collect(dep)
    while !isempty(todo)
        left = similar(todo, 0)
        for (k, x) in todo
            try
                vals[k] = evalx(x, vals)
            catch err
                err isa KeyError || rethrow()
                push!(left, k => x)
            end
        end
        length(left) == length(todo) && error("unresolved dependents: $(first.(left))")
        todo = left
    end
    vals
end

"""Identity of a step regardless of its stored direction: its two sides and RE flag."""
ukey(s) = (Set([(ER.name(ER.from_species(s)), sort(String.(ER.name.(ER.consumed(s))))),
                (ER.name(ER.to_species(s)), sort(String.(ER.name.(ER.released(s)))))]),
           ER.is_equilibrium(s))

plain_binding(s) = ER.is_binding(s) && !ER._is_chemistry(s)

"""A step of the oracle network: `a + cons → b + rel`. `k` is `[K]` (package convention: a
binding's dissociation constant, another RE step's forward constant) or `[kf, kr]`."""
struct OStep
    a::Symbol
    b::Symbol
    cons::Vector{Symbol}
    rel::Vector{Symbol}
    re::Bool
    binding::Bool
    k::Vector{BigFloat}
end
OStep(s::ER.Step, k) = OStep(ER.name(ER.from_species(s)), ER.name(ER.to_species(s)),
                             ER.name.(ER.consumed(s)), ER.name.(ER.released(s)),
                             ER.is_equilibrium(s), ER.is_binding(s), k)

"""Constant part of the forward equilibrium ratio [b]·Πrel / ([a]·Πcons)."""
eqratio(e::OStep) = e.re ? (e.binding ? 1 / e.k[1] : e.k[1]) : e.k[1] / e.k[2]

"""Per step of conformation `state`'s package mechanism (`_state_parts`): ukey => (step as
stored there, constant names as the derivation reads them)."""
function state_names(am, state)
    cm, sp = ER._state_parts(am, state)
    rename = ER._build_wegscheider_rename_map(cm; step_params = sp)
    out = Dict{Any, Tuple{ER.Step, Vector{Symbol}}}()
    for (j, (s, _)) in enumerate(ER._flat_steps(cm))
        nm = Symbol[ER.name(p, cm) for p in sp[j]]
        ER.is_equilibrium(s) && (nm = [get(rename, nm[1], nm[1])])
        out[ukey(s)] = (s, nm)
    end
    out
end

tag_of(am) = Dict(ukey(s) => ER.cat_allo_states(am)[g]
                  for (g, group) in enumerate(ER.steps(am)) for s in group)

"""Components of the forms joined by the steps outside `:OnlyA` groups (all forms), the
arcs (each `:OnlyA` plain binding, unbound side to bound side), the potentials (longest
path along the arcs) and the island-rule verdict: G1 (an arc inside a component or a
directed cycle of arcs) and the islands (components free E's does not reach along arcs)."""
function island_analysis(am)
    tags = ER.cat_allo_states(am)
    forms = unique(ER.name(f) for g in ER.steps(am) for s in g
                   for f in (ER.from_species(s), ER.to_species(s)))
    parent = Dict(f => f for f in forms)
    top(f) = (while parent[f] != f; f = parent[f]; end; f)
    for (g, group) in enumerate(ER.steps(am)), s in group
        tags[g] === :OnlyA && continue
        a, b = top(ER.name(ER.from_species(s))), top(ER.name(ER.to_species(s)))
        a == b || (parent[a] = b)
    end
    comps = Dict{Symbol, Vector{Symbol}}()
    for f in forms; push!(get!(comps, top(f), Symbol[]), f); end
    arcs = [(top(ER.name(ER.from_species(s))), top(ER.name(ER.to_species(s))))
            for (g, group) in enumerate(ER.steps(am)) for s in group
            if tags[g] === :OnlyA && plain_binding(s)]
    g1 = any(((u, v),) -> u == v, arcs)
    phi = Dict(c => 0 for c in keys(comps))
    for _ in 1:(length(comps) + 1), (u, v) in arcs
        u == v && continue
        phi[v] = max(phi[v], phi[u] + 1)
    end
    g1 |= any(v -> v > length(comps), values(phi))    # a directed cycle of arcs
    root = top(:E)
    reached = Set([root]); changed = true
    while changed
        changed = false
        for (u, v) in arcs
            u in reached && !(v in reached) && (push!(reached, v); changed = true)
        end
    end
    islands = [sort(comps[c]) for c in keys(comps) if !(c in reached)]
    admitted = [sort(comps[c]) for c in reached if c != root]
    # Breadth-first distance from free E's component over the components that any
    # :OnlyA step joins (the potential of mode :away).
    links = [(top(ER.name(ER.from_species(s))), top(ER.name(ER.to_species(s))))
             for (g, group) in enumerate(ER.steps(am)) for s in group if tags[g] === :OnlyA]
    dist = Dict(root => 0); frontier = [root]
    while !isempty(frontier)
        nxt = Symbol[]
        for c in frontier, (u, v) in links, (x, y) in ((u, v), (v, u))
            x == c && !haskey(dist, y) && (dist[y] = dist[c] + 1; push!(nxt, y))
        end
        frontier = nxt
    end
    (; comp = Dict(f => top(f) for f in forms), comps, arcs, phi, g1, islands, admitted, dist)
end

"""The oracle's steps of both conformations. Active: the package's active constants for
every step. Inactive: the active constants for an `:EqualAI` step, the package's
inactive constants for a `:NonequalAI` step, and for an `:OnlyA` step nothing
(`mode = :delete`) or an ε-scaled copy of its active constants (times `base` on the
equilibrium ratio): `:literal` scales a plain binding's uptake by ε (an RE one keeps RE
with K/ε) and both rates of any other step by ε^a (an RE one becomes steady state with
rates ε^a·(ratio, 1)); `:phi` keeps every ε thermodynamically consistent, scaling each
equilibrium ratio by ε to the potential difference (`island_analysis`) and the slower
rate of a non-binding step by ε^(a + |Δφ|); `:away` treats every :OnlyA step between
components like a binding: its rate away from free E vanishes, its rate back is kept."""
function oracle_steps(am, vals; mode = :delete, ε = big(0), a = 1, base = big(1))
    nA, nI, tags = state_names(am, :A), state_names(am, :I), tag_of(am)
    val(names) = BigFloat[vals[x] for x in names]
    stepsA = [OStep(s, val(nm)) for (s, nm) in values(nA)]
    stepsI, missingI = OStep[], String[]
    ia = mode in (:phi, :away) ? island_analysis(am) : nothing
    for (key, (sA, nmA)) in nA
        t = tags[key]
        if t === :EqualAI
            push!(stepsI, OStep(sA, val(nmA)))
        elseif t === :NonequalAI
            if haskey(nI, key)
                sI, nmI = nI[key]
                push!(stepsI, OStep(sI, val(nmI)))
            else
                push!(missingI, "$(ER.name(ER.from_species(sA)))-$(ER.name(ER.to_species(sA)))")
            end
        elseif mode !== :delete
            e = OStep(sA, val(nmA))
            r = eqratio(e) * base
            kr = e.re ? big(1) : e.k[2]
            if mode === :literal
                if e.binding && plain_binding(sA)
                    push!(stepsI, e.re ? OStep(e.a, e.b, e.cons, e.rel, true, true, [1 / (r * ε)]) :
                                         OStep(e.a, e.b, e.cons, e.rel, false, true, [kr * r * ε, kr]))
                else
                    push!(stepsI, OStep(e.a, e.b, e.cons, e.rel, false, e.binding,
                                        [kr * r * ε^a, kr * ε^a]))
                end
            elseif mode === :phi
                dphi = ia.phi[ia.comp[e.b]] - ia.phi[ia.comp[e.a]]
                rn = r * ε^dphi
                if e.binding && plain_binding(sA)
                    dphi >= 1 || error("arc without a positive potential difference")
                    push!(stepsI, e.re ? OStep(e.a, e.b, e.cons, e.rel, true, true, [1 / rn]) :
                                         OStep(e.a, e.b, e.cons, e.rel, false, true, [kr * rn, kr]))
                else
                    s = a + max(0, -dphi)
                    push!(stepsI, OStep(e.a, e.b, e.cons, e.rel, false, e.binding,
                                        [kr * rn * ε^s, kr * ε^s]))
                end
            elseif mode === :away
                # Every :OnlyA step joining two components loses its rate away from free E
                # (its equilibrium ratio scales by ε to the distance difference) and keeps
                # its rate toward E; one inside a component is kinetically blocked.
                d = ia.dist[ia.comp[e.b]] - ia.dist[ia.comp[e.a]]
                rn = r * ε^d
                if d == 0
                    push!(stepsI, OStep(e.a, e.b, e.cons, e.rel, false, e.binding,
                                        [kr * r * ε, kr * ε]))
                elseif e.re
                    push!(stepsI, OStep(e.a, e.b, e.cons, e.rel, true, e.binding,
                                        [e.binding ? 1 / rn : rn]))
                else
                    push!(stepsI, OStep(e.a, e.b, e.cons, e.rel, false, e.binding,
                                        d > 0 ? [kr * rn, kr] : [kr * r, kr * ε^(-d)]))
                end
            else
                error("unknown mode $mode")
            end
        end
    end
    stepsA, stepsI, missingI
end

"""Forms that `root` (free E by default) reaches along `steps`."""
function reach_from_E(steps, root = :E)
    reach = Set([root]); changed = true
    while changed
        changed = false
        for e in steps
            (e.a in reach) != (e.b in reach) && (push!(reach, e.a, e.b); changed = true)
        end
    end
    reach
end

"""Steady-state velocity per active site (net uptake of the first substrate `x1`, total
enzyme 1) of the concerted network: `n` protomers move independently within a
conformation, and only the oligomer with every protomer at free E flips, fast, with
[E…E_I]/[E…E_A] = L (`flip` names another flip form). RE steps run at rates
Λ·(ratio·Πcons, Πrel). Inactive forms the flip form cannot reach are left out (with
`strict`, an error)."""
function oracle_v(stepsA, stepsI, concs, L, x1; n = 1, Λ = big(10)^50, strict = false,
                  flip = :E)
    c(ms) = prod((concs[m] for m in ms); init = big(1))
    trans(steps) = [(e.a, e.b,
                     e.re ? Λ * eqratio(e) * c(e.cons) : e.k[1] * c(e.cons),
                     e.re ? Λ * c(e.rel) : e.k[2] * c(e.rel),
                     count(==(x1), e.cons) - count(==(x1), e.rel)) for e in steps]
    reachI = reach_from_E(stepsI, flip)
    formsI_all = Set(f for e in stepsI for f in (e.a, e.b))
    left_out = sort!(collect(setdiff(formsI_all, reachI)))
    strict && !isempty(left_out) && error("inactive forms unreachable from E: $left_out")
    tA = trans(stepsA)
    tI = trans(filter(e -> e.a in reachI && e.b in reachI, stepsI))
    formsA = sort!(unique([f for t in tA for f in (t[1], t[2])]))
    formsI = sort!(collect(reachI))
    joint(forms) = vec([collect(t) for t in Iterators.product(ntuple(_ -> forms, n)...)])
    nodes = [[(st, :A) for st in joint(formsA)]; [(st, :I) for st in joint(formsI)]]
    idx = Dict(nd => i for (i, nd) in enumerate(nodes))
    N = length(nodes)
    M = zeros(BigFloat, N, N)
    flux = Tuple{Int, Int, BigFloat, BigFloat, Int}[]
    for (conf, tr) in ((:A, tA), (:I, tI)), (st, cf) in nodes
        cf === conf || continue
        for i in 1:n, (fa, fb, rf, rr, w) in tr
            st[i] === fa || continue
            st2 = copy(st); st2[i] = fb
            p, q = idx[(st, conf)], idx[(st2, conf)]
            M[q, p] += rf; M[p, p] -= rf; M[p, q] += rr; M[q, q] -= rr
            w == 0 || push!(flux, (p, q, rf, rr, w))
        end
    end
    allE = fill(flip, n)
    if haskey(idx, (allE, :I))
        p, q = idx[(allE, :A)], idx[(allE, :I)]
        M[q, p] += Λ * L; M[p, p] -= Λ * L; M[p, q] += Λ; M[q, q] -= Λ
    end
    M[N, :] .= 1
    rhs = zeros(BigFloat, N); rhs[N] = 1
    x = M \ rhs
    v = sum((w * (rf * x[p] - rr * x[q]) for (p, q, rf, rr, w) in flux); init = big(0))
    mass_I = sum(x[i] for (i, nd) in enumerate(nodes) if nd[2] === :I; init = big(0))
    v / n, left_out, mass_I
end

"""Random fitted values: log-uniform on [0.1, 10], `L` on [0.2, 5]."""
draw_fitted(rng, fp) = Dict{Symbol, BigFloat}(
    p => (p === :L ? exp(BigFloat(log(0.2) + log(25) * rand(rng))) :
                     exp(BigFloat(log(0.1) + log(100) * rand(rng)))) for p in fp)
draw_conc(rng) = exp(BigFloat(log(0.1) + log(100) * rand(rng)))

"""The package's rate (BigFloat inputs, E_total = 1)."""
function pkg_rate(m, vals, concs)
    names = (ER.fitted_params(m)..., :Keq, :E_total)
    prm = NamedTuple{names}(Tuple(x === :E_total ? big(1) : vals[x] for x in names))
    mets = Tuple(ER.metabolites(m))
    cn = NamedTuple{mets}(Tuple(concs[x] for x in mets))
    BigFloat(ER.rate_equation(m, cn, prm))
end

relerr(a, b) = abs(a - b) / max(abs(b), big(10)^-300)

"""The stored steps of `am` with their tags, one line each."""
function step_lines(am)
    tags = ER.cat_allo_states(am)
    side(f, ms) = join([string(ER.name(f)); [string(ER.name(m)) for m in ms]], " + ")
    [side(ER.from_species(s), ER.consumed(s)) * (ER.is_equilibrium(s) ? " ⇌ " : " <--> ") *
     side(ER.to_species(s), ER.released(s)) * " :: " * string(tags[g])
     for (g, group) in enumerate(ER.steps(am)) for s in group]
end

"""Forms the package keeps in the inactive conformation (`_inactive_groups`)."""
inactive_forms(am) = sort!(unique([ER.name(f) for grp in ER._inactive_groups(am) for s in grp
                                   for f in (ER.from_species(s), ER.to_species(s))]))

"""Concentrations at equilibrium: random, with the last product set so ΠP = Keq·ΠS."""
function eq_concs(rng, am, Keq, base)
    rxn = ER.reaction(am)
    S = ER.name.(ER.substrates(rxn)); P = ER.name.(ER.products(rxn))
    c = copy(base)
    for x in S; c[x] = draw_conc(rng); end
    for x in P[1:end-1]; c[x] = draw_conc(rng); end
    c[P[end]] = Keq * prod(c[x] for x in S) / prod((c[x] for x in P[1:end-1]); init = big(1))
    c
end

"""DSL text of an allosteric mechanism from step lines `(text, tag)`."""
function dsl(subs, prods, steps; n = 1, inhibitors = String[])
    inh = isempty(inhibitors) ? "" : "    catalytic_inhibitors: " * join(inhibitors, ", ") * "\n"
    "@allosteric_mechanism begin\n    substrates: $(join(subs, ", "))\n" *
    "    products: $(join(prods, ", "))\n$inh    catalytic_multiplicity: $n\n" *
    "    catalytic_steps: begin\n" *
    join(["        $s :: $t" for (s, t) in steps], "\n") * "\n    end\nend"
end
