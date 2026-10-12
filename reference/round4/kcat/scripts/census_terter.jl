# ABOUTME: Every derivable init_terter law (the pipeline rows' ok flags): kcat class (core and full
# ABOUTME: law), per-substrate limits at zero products for kcat-zero and no-unique-limit laws, gauge.
include(joinpath(@__DIR__, "klib.jl"))

"""kcat's class on the whole law (every edge), plain mechanisms."""
function full_class(m, sol)
    try
        _, terms = S.kcat_terms(m, sol, trues(length(m.edges)))
        all(t -> iszero(t[1]), terms) ? :zero : :finite
    catch err
        msg = err isa ErrorException ? err.msg : sprint(showerror, err)
        occursin("never saturates", msg) ? :infinite : occursin("depends on how", msg) ? :nounique : :other
    end
end

out = NamedTuple[]
t0 = time()
for k in 1:11
    rows = deserialize(joinpath(IDATA, "pipe_terter$k.jls"))
    ps = deserialize(joinpath(IDATA, "pop_init_terter.part$k.jls"))
    for x in rows
        x.ok || continue
        m = S.Converted(ps[x.index], x.index).mech
        sol = S.solve(m)
        cls, msg, _, _ = kcat_class(m, sol)
        fcls = full_class(m, sol)
        L = S.law(m, sol)
        gap, gapfull = Int[], Int[]
        if cls in (:zero, :nounique)
            gap = zero_product_limits(m, S.law(m, sol, S.catalytic_core(m))).gap
            gapfull = zero_product_limits(m, L).gap
        end
        rng = Xoshiro(hash((k, x.index)))
        g = map(1:3) do _
            xx, _ = params_at(m, sol, rng)
            v = gauge_rates(NLaw(m, L, xx), rng)
            fin = filter(isfinite, v)
            (nnan = count(isnan, v), ninf = count(isinf, v), maxabs = maximum(abs, fin; init = 0.0),
             maxabs0 = maximum(abs, filter(isfinite, v[1:10]); init = 0.0))
        end
        push!(out, (part = k, index = x.index, nD = x.nD,
                    pipe = x.kcat isa String ? first(x.kcat, 60) : "", cls, fcls, gap, gapfull,
                    nnan = sum(p.nnan for p in g), ninf = sum(p.ninf for p in g),
                    maxabs = minimum(p.maxabs for p in g), maxabs0 = minimum(p.maxabs0 for p in g)))
    end
    ps = nothing
    GC.gc()
    @printf("part %d done: %d laws so far, %.0f s, maxrss %.0f MB\n", k, length(out), time() - t0,
            Sys.maxrss() / 2^20)
    flush(stdout)
    serialize(joinpath(KDATA, "census_terter.jls"), out)
end
cnt(f, xs = out) = (d = Dict{Any, Int}(); foreach(x -> (y = f(x); d[y] = get(d, y, 0) + 1), xs); sort!(collect(d); by = string))
println("init_terter derivable laws: ", length(out))
println("kcat class (core): ", cnt(x -> x.cls))
println("full-law class:    ", cnt(x -> x.fcls))
println("core => full:      ", cnt(x -> (x.cls, x.fcls)))
println("agreement with the pipeline's verdict (no unique limit or not): ",
        count(x -> (x.cls === :nounique) == occursin("depends on how", x.pipe), out), " of ", length(out))
for c in (:zero, :nounique)
    xs = filter(x -> x.cls === c, out)
    println(c, ": sign of deg_X(S*Q) - deg_X(D) per substrate (A, B, D), core law: ", cnt(x -> sign.(x.gap), xs))
    println(c, ": same, full law: ", cnt(x -> sign.(x.gapfull), xs))
end
println("gauge: laws with a NaN ", count(x -> x.nnan > 0, out), ", an Inf ", count(x -> x.ninf > 0, out),
        ", max|v| not finite > 0 at some parameter point ", count(x -> !(0 < x.maxabs < Inf), out),
        ", product-free max|v| not finite > 0 ", count(x -> !(0 < x.maxabs0 < Inf), out))
@printf("%.0f s\n", time() - t0)
