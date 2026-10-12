# ABOUTME: Helpers shared by the integration-stage scripts: every population's records (the
# ABOUTME: converted sets and the ter-ter exports), law and oracle values, the timing protocol.
include(joinpath(@__DIR__, "..", "allo", "lib.jl"))   # records, edge values, both oracles
using SHA, InteractiveUtils

const IDATA = joinpath(@__DIR__, "data")
const CORE_SETS = [String(split(l, '\t')[1])
                   for l in readlines(joinpath(DATA, "MANIFEST.tsv"))[2:end]]
const TERTER_SETS = [["init_terter.part$k" for k in 1:11]; "expand_terter"]
const SHARED = ("specs", "handwritten", "init_uniuni", "init_unibi", "init_bibi",
                "init_bibi_pp", "expand_bibi", "allo_uniuni", "allo_bibi")

"""The records `(set, p, c)` of one set: a converted core set, or a ter-ter export converted
here (`pop_<set>.jls`)."""
records_of(set) = set in CORE_SETS ? load_records((set,)) :
    [(set = set, p = p, c = S.Converted(p, i))
     for (i, p) in enumerate(deserialize(joinpath(IDATA, "pop_$set.jls")))]

"""The population a set belongs to (the ter-ter parts are one population)."""
population(set) = startswith(set, "init_terter") ? "init_terter" : set

"""`v/E_total` of a `Law` or an `AlloLaw` at concentrations `conc` and parameter values `x`
(both read by name), in the number type `T`."""
function law_value(L::S.Law, conc, x, ::Type{T} = BigFloat) where {T}
    v = [L.vars.conc[i] ? T(conc[n]) : T(x[n]) for (i, n) in enumerate(L.vars.names)]
    S.evaluate(L.H, v) * S.evaluate(L.Q, v) / S.evaluate(L.D, v)
end
function law_value(A::S.AlloLaw, conc, x, ::Type{T} = BigFloat) where {T}
    v = [A.vars.conc[i] ? T(conc[n]) : T(x[n]) for (i, n) in enumerate(A.vars.names)]
    ev(p) = S.evaluate(p, v)
    (QA, DA), (QI, DI) = ev.((A.active.Q, A.active.D)), ev.(S.inactive_terms(A))
    fA, fI, n, L = ev(A.fA), ev(A.fI), A.n, v[A.vars.index[:L]]
    WA = prod((ev(s[1])^s[3] for s in A.sites); init = one(T))
    WI = prod((ev(s[2])^s[3] for s in A.sites); init = one(T))
    ev(A.H) * (fI^n * QA * DA^(n - 1) * WA + L * fA^n * QI * DI^(n - 1) * WI) /
        (fI^n * DA^n * WA + L * fA^n * DI^n * WI)
end

"""
    random_values(m, sol, rng) -> (x, cv)

Random log-uniform BigFloat values of `parameters(m)` by name (`x`, a Dict, so no new
NamedTuple type is compiled per mechanism) and the value of every column of the law
(`cv[i]`): an independent's, or a dependent's power product of independents and `Keq`.
"""
function random_values(m, sol, rng)
    names = S.column_names(m, sol)
    x = Dict{Symbol, BigFloat}(n => Oracle.logu(rng) for n in S.parameters(m))
    cv = Dict{Int, BigFloat}(i => x[names[i]] for i in sol.indep)
    for d in sol.deps
        cv[d.col] = prod((S.pow(cv[f], e) for (f, e) in d.terms); init = S.pow(x[:Keq], d.keq))
    end
    x, cv
end

"""Per Pop step constants of record `r` (oracle convention) in one conformation, from the
column values `cv`; an edge whose column is outside the law gets placeholder ones."""
function step_constants(r, sol, cv, inactive)
    m = r.c.mech
    S.pop_constants(r.c, [[get(cv, sol.index[S.Col(k, s, inactive)], big(1))
                           for s in 1:(e.re ? 1 : 2)] for (k, e) in enumerate(m.edges)])
end

"""The independent route to `v/E_total` of record `r` at parameter values `x` and column
values `cv` (`random_values`): the steady-state oracle for a plain mechanism, P5a's
two-conformation oracle (deletion semantics, flip at the slice's flip form) for an
allosteric one."""
function oracle_value(r, sol, x, cv, conc)
    m, p = r.c.mech, r.p
    S.is_allosteric(m) || return Oracle.steady_state_v(p, step_constants(r, sol, cv, false), conc)
    mask = S.inactive_mask(m)
    actI = [mask[r.c.edge_of[j]] for j in eachindex(p.steps)]
    flip = pop_forms(r.c, p)[1]
    keep = Set([flip; [f for (j, s) in enumerate(p.steps) if actI[j] for f in (s.from, s.to)]])
    rn = Dict(S.regulator_names(m))
    regA, regI = Dict{Tuple{Int, Symbol}, Any}(), Dict{Tuple{Int, Symbol}, Any}()
    for (si, site) in enumerate(p.sites), (lig, st) in zip(site.ligands, site.states)
        regA[(si, lig)] = st === :OnlyI ? nothing : x[rn[lig][1]]
        regI[(si, lig)] = st === :OnlyA ? nothing : x[rn[lig][end]]
    end
    v_flip(p, flip, step_constants(r, sol, cv, false), step_constants(r, sol, cv, true),
           [s.re for s in p.steps], actI, keep, x[:L], regA, regI, conc;
           n = p.cat_multiplicity)
end

"""A random point of record `r`: parameter and column values (`random_values`) and
concentrations whose mass-action ratio is `Keq` times 10^±(0.3 to 2), off equilibrium."""
function random_point(r, sol, rng)
    x, cv = random_values(r.c.mech, sol, rng)
    u = (0.3 + 1.7 * rand(rng)) * (rand(rng, Bool) ? 1 : -1)
    x, cv, Oracle.concs_at_ratio(r.p, rng, x[:Keq] * big(10)^u)
end

"""The spanning-tree count of the segment graph of `m`'s edges in `mask` that the engine
expands (simple edges), as `bare_law` builds it."""
function tree_count(m, mask = trues(length(m.edges)))
    node, G = S.segment_nodes(m, mask)
    pairs = unique([minmax(node[e.a], node[e.b]) for (k, e) in enumerate(m.edges)
                    if mask[k] && !e.re && node[e.a] != node[e.b]])
    length(S.spanning_trees(G, pairs))
end

"""The package's performance protocol through a `@noinline` barrier: the minimum over 5
batches of 10_000 calls, in seconds per call."""
@noinline barrier(em, c, p) = S.rate_equation(em, c, p)
function t1n_of(em, c, p)
    barrier(em, c, p)
    best, acc = Inf, 0.0
    for _ in 1:5
        acc = 0.0
        t = @elapsed for _ in 1:10_000
            acc += barrier(em, c, p)
        end
        best = min(best, t / 10_000)
    end
    isfinite(acc) || error("non-finite result")
    best
end

Base.cumulative_compile_timing(true)
"""Cumulative compile time (ns) and JIT code bytes so far."""
compile_snapshot() = (Int(Base.cumulative_compile_time_ns()[1]), Int(Base.jit_total_bytes()))

"""Number of calls in a generated body (its operation count)."""
countops(b) = (n = Ref(0); w(x) = (x isa Expr && (x.head === :call && (n[] += 1);
                                                   foreach(w, x.args))); w(b); n[])

q3(v) = isempty(v) ? "-" : @sprintf("median %.3g, p90 %.3g, max %.3g", median(v),
                                    quantile(v, 0.9), maximum(v))
