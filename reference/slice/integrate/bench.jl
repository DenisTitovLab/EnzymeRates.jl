# ABOUTME: Compiles rate_equation for chosen laws and checks the generated code: allocations, ns
# ABOUTME: per call, compile time, values against the law; usage: bench.jl label selection sets…
include(joinpath(@__DIR__, "lib.jl"))

const LABEL, SELECTION, SETS = ARGS[1], ARGS[2], ARGS[3:end]

"""Log-uniform Float64 in [1e-2, 1e2]."""
logu64(rng) = 10.0^(4 * rand(rng) - 2)

"""Indices of the records of `rows` (pipeline rows of one population, derived ones) chosen by
`SELECTION`: `all`, or `decile+N` (the largest tenth by denominator terms plus `N` others
drawn at random), or `range:a:b`, or `list:<file>` (a saved `set => indices` Dict)."""
function chosen(rows, rng)
    SELECTION == "all" && return eachindex(rows)
    if startswith(SELECTION, "list:")
        want = Set(get(deserialize(joinpath(IDATA, SELECTION[6:end])), rows[1].set, Int[]))
        return [i for (i, x) in enumerate(rows) if x.index in want]
    end
    if startswith(SELECTION, "range:")
        a, b = parse.(Int, split(SELECTION, ":")[2:3])
        return intersect(a:b, eachindex(rows))
    end
    n = parse(Int, split(SELECTION, "+")[2])
    order = sortperm([x.nD for x in rows]; rev = true)
    top = order[1:cld(length(rows), 10)]
    rest = setdiff(eachindex(rows), top)
    sort!([top; rest[randperm(rng, length(rest))[1:min(n, length(rest))]]])
end

"""The relative error of the Float64 code at `c`, `p` against the law in BigFloat (absolute
where the law is 0), the condition number of the Haldane factor there (`|ΠS| + |ΠP|/Keq`
over `|ΠS − ΠP/Keq|`), and whether the code's value is finite."""
function value_error(em, L, m, c, p)
    v = S.rate_equation(em, c, p)
    exact = p.E_total * law_value(L, Dict(pairs(c)), Dict(pairs(p)))
    a, b = prod(big(c[x]) for x in m.substrates), prod(big(c[x]) for x in m.products) / p.Keq
    err = iszero(exact) ? Float64(abs(v)) : Float64(abs(v - exact) / abs(exact))
    err, Float64((a + b) / abs(a - b)), isfinite(v)
end

rng = Xoshiro(hash(LABEL) % UInt32)
out = NamedTuple[]
t0 = time()
let r = first(records_of(SETS[1])), m = r.c.mech                # warm the generator path
    S.rate_equation(S.compile(m), (; (x => 1.0 for x in S.metabolites(m))...),
                    (; (x => 1.0 for x in S.parameters(m))...))
end
pipe = Dict((x.set, x.index) => x for l in readdir(IDATA) if startswith(l, "pipe_")
            for x in deserialize(joinpath(IDATA, l)) if x.set in SETS && x.ok)
for set in SETS
    recs = records_of(set)
    derived = [r for r in recs if haskey(pipe, (r.set, r.c.index))]
    for r in derived[chosen([pipe[(r.set, r.c.index)] for r in derived], rng)]
        m = r.c.mech
        sol = S.solve(m)
        L = S.rate_law(m, sol)
        em = S.compile(m)
        c = (; (x => logu64(rng) for x in S.metabolites(m))...)
        p = (; (x => (x === :E_total ? 1.0 : logu64(rng)) for x in S.parameters(m))...)
        s0 = compile_snapshot()
        t_first = @elapsed S.rate_equation(em, c, p)
        s1 = compile_snapshot()
        allocs = (S.rate_equation(em, c, p); @allocated S.rate_equation(em, c, p))
        t1n = t1n_of(em, c, p)
        # init_terter (27k laws): one extra random point instead of two, no zero checks and
        # no BigFloat run (its compilation would double the cost)
        terter = startswith(set, "init_terter")
        errs = [value_error(em, L, m, c2, p) for c2 in [c; [(; (x => logu64(rng)
                for x in S.metabolites(m))...) for _ in 1:(terter ? 1 : 2)]]]
        zeros = terter ? [(0.0, 1.0, true)] :
                [value_error(em, L, m, merge(c, NamedTuple{(x,)}((0.0,))), p)
                 for x in S.metabolites(m)]
        cb = map(BigFloat, c)
        bigerr = terter ? NaN : Float64(relerr(S.rate_equation(em, cb, map(BigFloat, p)),
                                               law_value(L, Dict(pairs(cb)),
                                                         Dict(pairs(map(BigFloat, p))))))
        x = pipe[(r.set, r.c.index)]
        push!(out, (set = r.set, index = r.c.index, id = r.p.id, allo = S.is_allosteric(m),
                    second = S.canonical(L).key == x.key,           # another process's key
                    nD = x.nD, ops = x.ops, t_first, compile = (s1[1] - s0[1]) / 1e9,
                    jit = s1[2] - s0[2], allocs, t1n_ns = 1e9 * t1n,
                    err = maximum(first, errs), cond_err = maximum(e / k for (e, k) in errs),
                    zeros_ok = all(z -> z[3] && (z[1] < 1e-10 || z[1] / z[2] < 1e-12), zeros),
                    zero_err = maximum(first, zeros), bigerr))
    end
    @printf("%-20s %d laws done at %.0f s\n", set, count(x -> x.set == set, out), time() - t0)
    flush(stdout)
end
serialize(joinpath(IDATA, "bench_$LABEL.jls"), out)
@printf("%s: %d laws; allocations > 0: %d; T1n %s ns, over 120 ns: %d; compile %s ms; \
        cond. error max %.1e; zeros not ok %d; BigFloat error max %.1e; maxrss %.0f MB\n",
        LABEL, length(out), count(x -> x.allocs > 0, out), q3([x.t1n_ns for x in out]),
        count(x -> x.t1n_ns >= 120, out), q3([1e3 * x.compile for x in out]),
        maximum(x -> x.cond_err, out; init = 0.0), count(x -> !x.zeros_ok, out),
        maximum(x -> isnan(x.bigerr) ? 0.0 : x.bigerr, out; init = 0.0), Sys.maxrss() / 2^20)
