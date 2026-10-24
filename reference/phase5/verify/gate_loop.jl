# ABOUTME: Checks why the package's rate_equation timing loop can read under 0.5 ns per call: per-batch
# ABOUTME: times of its exact loop shape on a small inlinable law, against a @noinline barrier.
@inline law(c, p) = p.E_total * p.k1 * c.S / (p.K + c.S + c.P / p.Ki)
@noinline law_barrier(c, p) = law(c, p)
function helper_shape(f, c, p)                 # test_rate_eq_derivation.jl:214-231, batch times kept
    ts = Float64[]
    acc = 0.0
    for _ in 1:5
        acc = 0.0
        t = @elapsed for _ in 1:10_000
            acc += f(c, p)
        end
        push!(ts, t / 10_000)
    end
    isfinite(acc) || error("non-finite")
    ts
end
c = (S = 1.3, P = 0.4)
p = (E_total = 1.0, k1 = 2.0, K = 0.7, Ki = 3.0)
for f in (law, law_barrier)
    helper_shape(f, c, p)
    println(rpad(string(f), 12), " per-batch ns: ", join((round(1e9 * t; digits = 3) for t in helper_shape(f, c, p)), ", "))
end
