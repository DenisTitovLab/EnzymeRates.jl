# ABOUTME: Tests the gauge and the law on mechanisms the slice's checks reject or would admit if
# ABOUTME: they were dropped: RE/SS flips of enumerated parents, derived without the checks.
include(joinpath(@__DIR__, "vlib.jl"))

const PLAIN = ["init_bibi", "init_bibi_pp", "expand_bibi", "p1_enum_inh"]
const ALLO = ["allo_uniuni_inh", "allo_unibi", "allo_bibi", "allo_bibi_pp", "allo_terbi"]

"""`m` with the RE flag of each kinetic group in `gs` toggled (whole groups), or `nothing`."""
function toggle(m, gs)
    (subs, prods, steps, group), kw = S.raw(m)
    steps = copy(steps)
    for (j, s) in enumerate(steps)
        group[j] in gs || continue
        steps[j] = S.RawStep(s.from, s.to, s.consumed, s.released, !s.re)
    end
    try
        S.Mechanism(subs, prods, steps, collect(Int, group); tags = kw.tags,
                    multiplicity = kw.multiplicity, sites = kw.sites)
    catch
        nothing
    end
end

"""The current checks' verdict: :pass, :boundary (check_boundary error) or :admit."""
function verdict(m)
    try
        S.check_boundary(m)
    catch
        return :boundary
    end
    S.admit(m) === nothing ? :pass : :admit
end

"""Gauge at 3 parameter points: 50 vectors (all concentrations log-uniform 1e-3..1e3), the
first 10 with the products at zero. Returns per point (nan, inf, maxabs, maxabs over the
product-free vectors)."""
function gauge(L, nc, prods, rng)
    map(1:3) do _
        nl = NL(L, draw_params(L.vars, rng))
        vs = map(1:50) do i
            c = [10.0^(-3 + 6 * rand(rng)) for _ in 1:nc]
            i <= 10 && (c[prods] .= 0.0)
            rate(nl, c)
        end
        fin = filter(isfinite, vs)
        (nan = count(isnan, vs), inf = count(isinf, vs), maxabs = maximum(abs, fin; init = 0.0),
         maxabs0 = maximum(abs, filter(isfinite, vs[1:10]); init = 0.0))
    end
end

"""Law against the mass-action oracle (plain): at a positive point and on the forward face
(products zero), with the RE speed Λ = 1e20 and 1e30 and a second speed pattern."""
function oracle_check(m, sol, L, rng)
    names, subs, prods, mods = roles(m)
    pd = draw_params(L.vars, rng)
    kv = edge_constants(m, sol, pd)
    nl = NL(L, pd)
    w = [0.1 + 10 * rand(rng) for _ in m.edges]
    map((:positive, :face)) do where_
        c = big.([10.0^(-1 + 2 * rand(rng)) for _ in names])
        where_ === :face && (c[prods] .= 0)
        vl = rate(nl, c)
        o20, _ = oracle(m, kv, c; Λ = big(10)^20)
        o30, _ = oracle(m, kv, c; Λ = big(10)^30)
        o30w, _ = oracle(m, kv, c; Λ = big(10)^30, w)
        rel(a, b) = Float64(abs(a - b) / max(abs(b), big(1e-300)))
        (law = Float64(vl), o30 = Float64(o30), law_vs_oracle = rel(vl, o30),
         growth = Float64(abs(o30) / max(abs(o20), big(1e-300))), speed = rel(o30w, o30))
    end
end

function analyze(parent, m, how, rng)
    base = (parent = parent, how = how)
    v = verdict(m)
    ok, why, sol, L = derive_law(m; checks = false)
    ok || return merge(base, (verdict = v, law = false, why = why, allo = S.is_allosteric(m)))
    names, subs, prods, mods = roles(m)
    Qs = L isa S.Law ? [L.Q] : [L.active.Q, S.inactive_terms(L)[1]]
    zf = zero_facts(m, L)
    cls, _, _ = kcat_class(m, sol)
    g = gauge(L, length(names), prods, rng)
    oc = !S.is_allosteric(m) ? oracle_check(m, sol, L, rng) : nothing
    merge(base, (verdict = v, law = true, why = "", allo = S.is_allosteric(m),
                 Qzero = all(q -> isempty(q.terms), Qs), zf, cls, g, oc, lines = steps_text(m)))
end

t0 = time()
rng = Xoshiro(424242)
rows = Any[]
seen = Set{UInt}()
for (sets, nparents) in ((PLAIN, 300), (ALLO, 150))
    pool = [r for set in sets for r in core_records(set)]
    for i in randperm(rng, length(pool))[1:nparents]
        r = pool[i]
        G = S.ngroups(r.m)
        picks = [[rand(rng, 1:G)]]
        G >= 2 && push!(picks, randperm(rng, G)[1:2])
        for gs in picks
            c = toggle(r.m, gs)
            c === nothing && continue
            h = hash(c)
            h in seen && continue
            push!(seen, h)
            re_before = [r.m.edges[findfirst(==(g), r.m.group)].re for g in gs]
            dir = all(re_before) ? :re_to_ss : !any(re_before) ? :ss_to_re : :mixed
            push!(rows, merge(analyze("$(r.set)/$(r.index)", c, "toggle $gs", rng), (dir = dir,)))
        end
    end
    @printf("%s done: %d mutants at %.0f s\n", sets[1], length(rows), time() - t0)
    flush(stdout)
end
serialize(joinpath(VD, "data", "mutant_rows.jls"), rows)
@printf("done %d mutants in %.0f s, maxrss %.0f MB\n", length(rows), time() - t0, Sys.maxrss() / 2^20)
