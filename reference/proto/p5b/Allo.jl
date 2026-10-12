# ABOUTME: Single-graph allosteric (MWC) rate law of a Pop.PMech: both conformations from one step
# ABOUTME: list, one joint log-linear solve, whole-mechanism names (probe P5b); stdlib only.
isdefined(Main, :Engine) || include(joinpath(@__DIR__, "..", "common", "Engine.jl"))
isdefined(Main, :Thermo) || include(joinpath(@__DIR__, "..", "common", "Thermo.jl"))

"""
Both conformations of an allosteric mechanism come from its one step list.

- **Inactive graph** (today's `:OnlyA` semantics, kept for this probe): the steps of the
  non-`:OnlyA` groups whose two forms are reachable from a free form (unbound,
  residual-free) over those steps.
- **Columns**: one per bare stored constant of each conformation, `(state, step, slot)`:
  slot 1 is an RE step's `K` or an SS step's `kf`, slot 2 an SS step's `kr`, both in the
  step's stored direction.
- **Rows**: the fundamental cycles of each conformation's graph (Haldane: `c·log Keq`,
  Wegscheider: 0), kinetic-group equalities within each conformation, and the `:EqualAI`
  ties `I = A`.
- **Priority**: every inactive column before every active one, then P2's content key
  `(kind, occupied, ss, against)` and the content name; one exact RREF picks the dependents.
- **Names**: a column renders its group's whole-mechanism name for its state
  (`step_names_A` / `step_names_I`), so group members, and an `:EqualAI` group in both
  states, share one name. A dependent whose right side is one name with exponent 1 is a
  fold (a rename applied to the law); every other dependent is an assignment.
- **Law**: each conformation's bare law (Engine) with its step variables renamed, so
  `v = E_t·(N_A·Q_A^(n-1)·W_A + L·N_I·Q_I^(n-1)·W_I) / (Q_A^n·W_A + L·Q_I^n·W_I)` with each
  state's `N`, `Q` divided by its free-enzyme weight `d`.
"""
module Allo

import ..Pop, ..Engine, ..Thermo
using ..Poly: Poly, Vars, var!, Mono, mono, mmul, P, subst, lowest_shift, nterms

const Q = Rational{Int64}
const PI = P{Int}
const Col = Tuple{Symbol, Int, Int}

# ─── The two graphs ────────────────────────────────────────────────────

"""Steps of the inactive conformation: non-`:OnlyA` steps whose forms are reachable from a
free form over the non-`:OnlyA` steps."""
function inactive_kept(m::Pop.PMech)
    g = Pop.group_of_step(m)
    live = [m.group_tags[g[j]] !== :OnlyA for j in eachindex(m.steps)]
    reach = Set(f for (j, s) in enumerate(m.steps) if live[j] for f in (s.from, s.to)
                if Pop.is_free(m.species[f]))
    grown = true
    while grown
        grown = false
        for (j, s) in enumerate(m.steps)
            live[j] && (s.from in reach) != (s.to in reach) || continue
            push!(reach, s.from, s.to); grown = true
        end
    end
    [live[j] && s.from in reach && s.to in reach for (j, s) in enumerate(m.steps)]
end

"""The plain mechanism of the steps `keep` and the forms they touch, with each of its steps'
index in `m`."""
function submech(m::Pop.PMech, keep)
    js = findall(keep)
    forms = sort!(unique([f for j in js for f in (m.steps[j].from, m.steps[j].to)]))
    loc = Dict(f => i for (i, f) in enumerate(forms))
    steps = [Pop.PStep(loc[s.from], loc[s.to], s.consumed, s.released, s.re) for s in m.steps[js]]
    pm = Pop.PMech(; id = m.id, source = m.source, substrates = m.substrates,
                   products = m.products, inhibitors = m.inhibitors,
                   metabolites = m.metabolites, species = m.species[forms], steps,
                   groups = [[k] for k in eachindex(steps)])
    pm, js
end

# ─── The joint solve ───────────────────────────────────────────────────

nslots(s::Pop.PStep) = s.re ? 1 : 2

"""Fundamental cycles of the multigraph `edges[k] = (a, b)` as ±1 vectors over the edges
(+1 traversed a → b), from a breadth-first spanning forest."""
function fundamental_cycles(edges)
    adj = Dict{Int, Vector{NTuple{3, Int}}}()
    for (k, (a, b)) in enumerate(edges)
        push!(get!(adj, a, NTuple{3, Int}[]), (b, k, 1))
        push!(get!(adj, b, NTuple{3, Int}[]), (a, k, -1))
    end
    parent, depth, tree = Dict{Int, NTuple{3, Int}}(), Dict{Int, Int}(), falses(length(edges))
    for r in sort!(collect(keys(adj)))
        haskey(depth, r) && continue
        depth[r] = 0; queue = [r]
        while !isempty(queue)
            u = popfirst!(queue)
            for (v, k, sg) in adj[u]
                haskey(depth, v) && continue
                depth[v] = depth[u] + 1; parent[v] = (u, k, sg); tree[k] = true
                push!(queue, v)
            end
        end
    end
    map(findall(!, tree)) do k
        z = zeros(Int, length(edges)); z[k] = 1
        x, y = edges[k][2], edges[k][1]                 # close b → a through the tree
        while x != y
            if depth[x] >= depth[y]
                p, kk, sg = parent[x]; z[kk] -= sg; x = p
            else
                p, kk, sg = parent[y]; z[kk] += sg; y = p
            end
        end
        z
    end
end

"""The multiple of the net reaction that cycle `z` over the steps `js` of `m` runs."""
function cycle_multiple(m::Pop.PMech, js, z)
    tot = Dict{Symbol, Int}()
    for (k, zk) in enumerate(z), (ms, sg) in ((m.steps[js[k]].consumed, -1),
                                              (m.steps[js[k]].released, 1)), x in ms
        tot[x.name] = get(tot, x.name, 0) + sg * zk
    end
    filter!(p -> p.second != 0, tot)
    isempty(tot) && return 0
    net = Dict{Symbol, Int}()
    for x in m.substrates; net[x] = get(net, x, 0) - 1; end
    for x in m.products; net[x] = get(net, x, 0) + 1; end
    c = get(tot, first(m.substrates), 0) // net[first(m.substrates)]
    isinteger(c) && keys(tot) ⊆ keys(net) && all(get(tot, x, 0) == c * v for (x, v) in net) ||
        error("$(m.id): a cycle's net change $tot is no multiple of the net reaction")
    Int(c)
end

"""Elimination key of a column: the inactive state first, then P2's content key and name."""
function priority(cs, (st, j, r)::Col)
    c = cs[j]
    cr = c.re ? 1 : (c.flipped ? 3 - r : r)
    against = c.re ? 0 : Int((cr == 1) == (c.fwd != 1))
    ((st === :I ? 1 : 0, c.kind, c.occ, c.re ? 0 : 1, against), c.names[cr])
end

"""Every column as `Keq^q · Π independent^e` (an independent column is itself)."""
const Expr_ = Tuple{Q, Vector{Pair{Int, Q}}}

"""
    joint_solve(m, keepI) -> (cols, indep, expr, (; rows, maxint))

One RREF of the cycle rows of both conformations, the kinetic-group rows of each and the
`:EqualAI` ties, over the columns in elimination order.
"""
function joint_solve(m::Pop.PMech, keepI)
    cs = Thermo.canonical_steps(first(submech(m, trues(length(m.steps)))))
    cols = [(st, j, r) for st in (:A, :I) for (j, s) in enumerate(m.steps)
            if st === :A || keepI[j] for r in 1:nslots(s)]
    sort!(cols; by = c -> priority(cs, c), rev = true)
    at = Dict(c => i for (i, c) in enumerate(cols))
    rows = Vector{Pair{Vector{Int}, Int}}()
    newrow() = zeros(Int, length(cols))
    for (st, keep) in ((:A, trues(length(m.steps))), (:I, keepI))
        js = findall(keep)
        for z in fundamental_cycles([(m.steps[j].from, m.steps[j].to) for j in js])
            row = newrow()
            for (k, zk) in enumerate(z)
                zk == 0 && continue
                j = js[k]; s = m.steps[j]
                if s.re
                    row[at[(st, j, 1)]] += Pop.is_binding(s) ? -zk : zk
                else
                    row[at[(st, j, 1)]] += zk; row[at[(st, j, 2)]] -= zk
                end
            end
            push!(rows, row => cycle_multiple(m, js, z))
        end
    end
    tie!(a, b) = (row = newrow(); row[at[a]] += 1; row[at[b]] -= 1; push!(rows, row => 0))
    for (g, grp) in enumerate(m.groups)
        for st in (:A, :I)
            mem = [j for j in grp if st === :A || keepI[j]]
            for j in mem[2:end], r in 1:nslots(m.steps[j])
                tie!((st, j, r), (st, mem[1], r))
            end
        end
        m.group_tags[g] === :EqualAI || continue
        for j in grp, r in 1:nslots(m.steps[j])
            keepI[j] && tie!((:I, j, r), (:A, j, r))
        end
    end
    n = length(cols)
    M = zeros(Q, length(rows), n + 1)
    for (i, (row, c)) in enumerate(rows)
        M[i, 1:n] .= row; M[i, n + 1] = c
    end
    piv, maxint = Thermo.rref!(M)
    (n + 1) in piv && error("$(m.id): the constraints reduce to 0 = log Keq")
    indep = setdiff(1:n, piv)
    expr = Expr_[(0 // 1, [i => 1 // 1]) for i in 1:n]
    for (r, p) in enumerate(piv)
        expr[p] = (M[r, n + 1], [f => -M[r, f] for f in indep if M[r, f] != 0])
    end
    cols, indep, expr, (rows = length(rows), maxint)
end

# ─── Names ─────────────────────────────────────────────────────────────

colname(m::Pop.PMech, (st, j, r)::Col) = (st === :A ? m.step_names_A : m.step_names_I)[j][r]

"""Regulator-site constants of conformation `st` (`:A` or `:I`): `nothing` for a ligand
absent from it; an `:EqualAI` ligand shares the active name."""
function regname(st, lig, tag)
    tag === (st === :A ? :OnlyI : :OnlyA) && return nothing
    Symbol("K_", st === :I && tag !== :EqualAI ? "I_" : "A_", lig, "reg")
end

"""
The solve at the level of names. `value[name]` is each name's right side
`(Keq exponent, [independent name => exponent])`; `fitted` the independent catalytic names
(sorted), then the regulator constants (active site by site, inactive site by site), then
`L`; `folds` the names whose right side is one other name; `deps` every other dependent with
its label `:haldane` (mentions Keq) or `:wegscheider`; `tie_kind` how each fold arose
(`:collapse` = the inactive copy of a group equals its active one, `:wegscheider` otherwise).
"""
function names_of(m::Pop.PMech, cols, indep, expr)
    nm = [colname(m, c) for c in cols]
    value = Dict{Symbol, Tuple{Q, Vector{Pair{Symbol, Q}}}}()
    for (i, (q, ts)) in enumerate(expr)
        v = (q, sort!([nm[f] => e for (f, e) in ts]; by = first))
        get!(value, nm[i], v) == v ||
            error("$(m.id): the columns named $(nm[i]) take different values")
    end
    cat = sort!([nm[i] for i in indep]; by = string)
    allunique(cat) || error("$(m.id): two independent columns share a name")
    regs = [Symbol[regname(st, l, t) for site in m.sites
                   for (l, t) in zip(site.ligands, site.states)
                   if (st === :A || t !== :EqualAI) && regname(st, l, t) !== nothing]
            for st in (:A, :I)]
    fitted = [cat; regs[1]; regs[2]; :L]
    g = Pop.group_of_step(m)
    folds, deps, tie_kind = Dict{Symbol, Symbol}(), Dict{Symbol, Symbol}(), Dict{Symbol, Symbol}()
    for (i, c) in enumerate(cols)
        q, ts = value[nm[i]]
        (i in indep || ts == [nm[i] => 1]) && continue        # the name itself is independent
        if q == 0 && length(ts) == 1 && ts[1].second == 1
            folds[nm[i]] = ts[1].first
            t = cols[only(f for (f, _) in expr[i][2])]
            tie_kind[nm[i]] = c[1] === :I && t[1] === :A && g[c[2]] == g[t[2]] ? :collapse :
                              :wegscheider
        else
            deps[nm[i]] = q != 0 ? :haldane : :wegscheider
        end
    end
    (; value, fitted, folds, deps, tie_kind)
end

# ─── The law ───────────────────────────────────────────────────────────

"""Per state `(num, den, dfree)` over the name table `T`: the Engine bare law of the state's
graph with each step constant renamed to its column's name (folds applied) and the
concentrations brought to lowest terms. An inactive state with no step is `(0, 1, 1)`."""
function named_law(m, L::Union{Nothing, Engine.Law}, js, st, at, nm, folds, value, T::Vars)
    L === nothing && return zero(PI), one(PI), one(PI)
    sub = Dict{Int, Mono}()
    # The segment weights follow one spanning tree of the RE steps; around an RE cycle a
    # dependent constant written through independents makes them independent of that tree.
    re = [k for k in eachindex(js) if m.steps[js[k]].re]
    oncycle = Set(re[i] for z in fundamental_cycles([(m.steps[js[k]].from, m.steps[js[k]].to)
                                                     for k in re]) for i in findall(!iszero, z))
    for (k, vs) in enumerate(L.stepvars), (r, v) in enumerate(vs)
        n = nm[at[(st, js[k], r)]]
        sub[v] = mono(var!(T, get(folds, n, n)) => 1)
        if k in oncycle && value[n] != (0, [n => 1]) && !haskey(folds, n)
            q, ts = value[n]
            sub[v] = mono((var!(T, x) => Int(e) for (x, e) in (q == 0 ? ts : [:Keq => q; ts]))...)
        end
    end
    for v in L.concvars
        sub[v] = mono(var!(T, L.vars.names[v]) => 1)
    end
    num, den, dfree = (subst(p, v -> sub[v]) for p in (L.num, L.den, L.dfree))
    shift = lowest_shift((num, den), [T.index[x] for x in m.metabolites if haskey(T.index, x)])
    num * shift, den * shift, dfree * shift
end

"""The derived allosteric model of `m`."""
struct Model
    m::Pop.PMech
    keepI::Vector{Bool}
    cols::Vector{Col}
    indep::Vector{Int}
    expr::Vector{Expr_}
    solve_stats::NamedTuple
    names::NamedTuple
    T::Vars
    A::NTuple{3, PI}
    I::NTuple{3, PI}
end

function derive(m::Pop.PMech)
    keepI = inactive_kept(m)
    cols, indep, expr, solve_stats = joint_solve(m, keepI)
    names = names_of(m, cols, indep, expr)
    at = Dict(c => i for (i, c) in enumerate(cols))
    nm = [colname(m, c) for c in cols]
    T = Vars()
    laws = map((:A, :I)) do st
        pm, js = submech(m, st === :A ? trues(length(m.steps)) : keepI)
        L = isempty(js) ? nothing : Engine.bare_rate_law(pm; check = false)
        named_law(m, L, js, st, at, nm, names.folds, names.value, T)
    end
    Model(m, keepI, cols, indep, expr, solve_stats, names, T, laws[1], laws[2])
end

"""Normalization case of the free-enzyme weights, as the package decides it: `:same`
(`d_A == d_I`), `:divide` (both metabolite-free monomials) or `:cross`."""
function norm_case(md::Model)
    dA, dI = md.A[3], md.I[3]
    dA == dI && return :same
    conc = Set(md.T.index[x] for x in md.m.metabolites if haskey(md.T.index, x))
    monofree(d) = nterms(d) == 1 && !any(v in conc for (v, _) in only(keys(d.t)))
    monofree(dA) && monofree(dI) ? :divide : :cross
end

# ─── Numbers ───────────────────────────────────────────────────────────

"""Values of every name of `md` from values of its fitted names and `Keq` (BigFloat)."""
function all_values(md::Model, fit::AbstractDict{Symbol})
    v = Dict{Symbol, BigFloat}(k => BigFloat(x) for (k, x) in fit)
    for (n, (q, ts)) in md.names.value
        haskey(v, n) && continue
        v[n] = BigFloat(v[:Keq])^BigFloat(q) *
               prod((BigFloat(v[s])^BigFloat(e) for (s, e) in ts); init = big(1))
    end
    v
end

"""Each regulatory site's factor `(1 + Σ [lig]/K)^multiplicity` in conformation `st`."""
reg_factor(m, st, vals, concs) =
    prod((big(1) + sum((BigFloat(concs[l]) / vals[regname(st, l, t)]
                        for (l, t) in zip(site.ligands, site.states)
                        if regname(st, l, t) !== nothing); init = big(0)))^site.multiplicity
         for site in m.sites; init = big(1))

"""`v / E_t` of `md` at the values `vals` of every name (`all_values`) and `concs`; `n`
overrides the catalytic multiplicity."""
function rate(md::Model, vals, concs; n = md.m.cat_multiplicity, regs = true)
    x = zeros(BigFloat, length(md.T.names))
    for (i, s) in enumerate(md.T.names)
        x[i] = haskey(vals, s) ? vals[s] : BigFloat(concs[s])
    end
    ev(p) = Poly.evaluate(p, x)
    nA, qA = ev(md.A[1]) / ev(md.A[3]), ev(md.A[2]) / ev(md.A[3])
    nI, qI = ev(md.I[1]) / ev(md.I[3]), ev(md.I[2]) / ev(md.I[3])
    WA = regs ? reg_factor(md.m, :A, vals, concs) : big(1)
    WI = regs ? reg_factor(md.m, :I, vals, concs) : big(1)
    L = vals[:L]
    (nA * qA^(n - 1) * WA + L * nI * qI^(n - 1) * WI) / (qA^n * WA + L * qI^n * WI)
end

# ─── Text ──────────────────────────────────────────────────────────────

"""Sort key of a name that ignores its conformation tag, so both conformations' terms come
out in the same order: (the name without `A_`/`I_`, the tag)."""
function rolekey(s::Symbol)
    t = string(s)
    mt = match(r"^([kK])_([AI])_(.*)$", t)
    mt === nothing ? (t, "") : (string(mt[1], "_", mt[3]), mt[2])
end

"""`p` as text: terms by (concentration degree, sign, role keys), parameters before
concentrations inside a term, negative exponents in a divisor."""
function poly_text(p::PI, T::Vars, isconc)
    isempty(p.t) && return "0"
    term_key((mo, c)) = (sum((e for (v, e) in mo if isconc(T.names[v])); init = 0), c < 0,
                         sort!([(rolekey(T.names[v]), e) for (v, e) in mo]))
    pos, neg = String[], String[]
    for (mo, c) in sort!(collect(p.t); by = term_key)
        sym = sort(mo; by = ((v, _),) -> (isconc(T.names[v]), rolekey(T.names[v])))
        pw(v, e) = abs(e) == 1 ? string(T.names[v]) : string(T.names[v], " ^ ", abs(e))
        up = [pw(v, e) for (v, e) in sym if e > 0]
        dn = [pw(v, e) for (v, e) in sym if e < 0]
        abs(c) != 1 && pushfirst!(up, string(abs(c)))
        t = isempty(up) ? "1" : join(up, " * ")
        isempty(dn) || (t *= length(dn) == 1 ? " / " * dn[1] : " / (" * join(dn, " * ") * ")")
        push!(c > 0 ? pos : neg, t)
    end
    s = join(pos, " + ")
    isempty(neg) && return s
    n = length(neg) == 1 ? neg[1] : "(" * join(neg, " + ") * ")"
    isempty(pos) ? "-" * n : s * " - " * n
end

paren(t) = occursin(r"[ +\-*/]", t) ? "(" * t * ")" : t

"""`base ^ k` as text, nothing for a factor of 1 or a power 0."""
function power_text(base::String, k::Int)
    (base == "1" || k == 0) && return nothing
    k == 1 ? base : string(paren(base), " ^ ", k)
end

"""A product of optional factors as text (`1` when every factor is absent); a sum is
parenthesized only next to another factor."""
function prod_text(fs)
    f = filter(!isnothing, fs)
    isempty(f) ? "1" : length(f) == 1 ? f[1] :
    join((occursin(r" [+-] ", x) && !startswith(x, "(") ? "(" * x * ")" : x for x in f), " * ")
end

"""
The MWC rate law as text: the `params` and `concs` destructures, the dependents under
`# Wegscheider constraints:` / `# Haldane constraints:`, the folds under `# Ties:`, and
`v = E_total * (num) / (den)`. A state's polynomials are divided by its free-enzyme weight
when both weights are metabolite-free monomials (`norm_case` `:divide`), and each state is
multiplied by the other's weight to the power `n` when they differ otherwise (`:cross`).
"""
function law_text(md::Model)
    m, T, n = md.m, md.T, md.m.cat_multiplicity
    isconc(s) = s in m.metabolites
    case = norm_case(md)
    pol(st) = (p = st === :A ? md.A : md.I;
               case === :divide ? (p[1] * inv_mono(p[3]), p[2] * inv_mono(p[3])) : (p[1], p[2]))
    sites(st) = [power_text(join(["1"; [string(l, " / ", regname(st, l, t))
                                        for (l, t) in zip(site.ligands, site.states)
                                        if regname(st, l, t) !== nothing]], " + "),
                            site.multiplicity) for site in m.sites]
    other(st) = case === :cross ? power_text(poly_text((st === :A ? md.I : md.A)[3], T, isconc), n) :
                nothing
    function branch(st, numerator)
        N, Q = pol(st)
        numerator && isempty(N.t) && return nothing
        fs = numerator ? [power_text(poly_text(N, T, isconc), 1); power_text(poly_text(Q, T, isconc), n - 1)] :
             [power_text(poly_text(Q, T, isconc), n)]
        prod_text([fs; sites(st); other(st)])
    end
    lterm(b) = b == "1" ? "L" : "L * " * (occursin(r" [+-] ", b) && !startswith(b, "(") ? "(" * b * ")" : b)
    num = join(filter(!isnothing, [branch(:A, true),
                                   (b = branch(:I, true); b === nothing ? nothing : lterm(b))]), " + ")
    den = branch(:A, false) * " + " * lterm(branch(:I, false))
    lines = ["(; " * join([md.names.fitted; :Keq; :E_total], ", ") * ") = params",
             "(; " * join(m.metabolites, ", ") * ") = concs"]
    ex(e) = isinteger(e) ? string(Int(abs(e))) : string(Float64(abs(e)))
    fac(s, e) = abs(e) == 1 ? string(s) : string(s, " ^ ", ex(e))
    function rhs((q, ts))
        fs = q == 0 ? ts : [:Keq => q; ts]
        up = [fac(s, e) for (s, e) in fs if e > 0]; dn = [fac(s, e) for (s, e) in fs if e < 0]
        t = isempty(up) ? "1" : join(up, " * ")
        isempty(dn) ? t : t * (length(dn) == 1 ? " / " * dn[1] : " / (" * join(dn, " * ") * ")")
    end
    for (label, head) in ((:wegscheider, "# Wegscheider constraints:"), (:haldane, "# Haldane constraints:"))
        ds = sort!([d for (d, l) in md.names.deps if l === label]; by = string)
        isempty(ds) || push!(lines, head, ("$d = $(rhs(md.names.value[d]))" for d in ds)...)
    end
    fs = sort!(collect(md.names.folds); by = string ∘ first)
    isempty(fs) || push!(lines, "# Ties:", ("$a = $b" for (a, b) in fs)...)
    push!(lines, "v = E_total * ($num) / ($den)")
    join(lines, "\n")
end

"""The inverse of a one-term polynomial with a unit coefficient."""
function inv_mono(p::PI)
    mo, c = only(p.t)
    abs(c) == 1 || error("non-unit monomial")
    Poly.term(PI, Mono([v => -e for (v, e) in mo]), c)
end

end # module
