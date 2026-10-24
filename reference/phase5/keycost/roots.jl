# ABOUTME: Item 11: cost and accuracy of the root emissions codegen could use for x^(1/q), q > 2
# ABOUTME: (Float64 exponent, exponent in the input's type, cbrt, nested sqrt, exp/log), Float64 and BigFloat.
using Printf, Random

const EMIT = [
    ("x^(1/q) [Float64 exponent, slice]", (x, q) -> x^(1 / q)),
    ("x^(one(x)/q)", (x, q) -> x^(one(x) / q)),
    ("cbrt(x) [q = 3]", (x, q) -> cbrt(x)),
    ("sqrt(sqrt(x)) [q = 4]", (x, q) -> sqrt(sqrt(x))),
    ("exp(log(x)/q)", (x, q) -> exp(log(x) / q)),
]
applies(name, q) = !(occursin("q = 3", name) && q != 3) && !(occursin("q = 4", name) && q != 4)

"""Seconds per call of `f(x, q)` over `xs` through a `@noinline` barrier, minimum of 5."""
@noinline callf(f::F, x, q) where {F} = f(x, q)
function per_call(f::F, xs, q) where {F}
    acc, best = 0.0, Inf
    for _ in 1:5
        t = @elapsed for x in xs
            acc += callf(f, x, q)
        end
        best = min(best, t / length(xs))
    end
    best, acc
end

rng = Xoshiro(20261014)
xs = [10.0^(4 * rand(rng) - 2) for _ in 1:100_000]
setprecision(BigFloat, 256)
bx = [big(10.0)^(4 * big(rand(rng)) - 2) for _ in 1:2_000]
for q in 3:6
    for (name, f) in EMIT
        applies(name, q) || continue
        t, _ = per_call(f, xs, q)
        # Float64 accuracy against the 512-bit root, in units of eps(Float64)
        e64 = setprecision(BigFloat, 512) do
            maximum(abs(big(f(x, q)) / big(x)^(big(1) / q) - 1) for x in xs[1:2000]) / eps()
        end
        # 256-bit inputs: relative error against the 512-bit root of the same input
        ebig = maximum(bx) do x
            y = f(x, q)
            setprecision(BigFloat, 512) do
                abs(BigFloat(y; precision = 512) / BigFloat(x; precision = 512)^(big(1) / q) - 1)
            end
        end
        bt = @elapsed for x in bx
            f(x, q)
        end
        @printf("q=%d %-36s Float64 %6.2f ns, max err %.2f eps; BigFloat-256 max rel err %.1e, %.2f µs\n",
                q, name, 1e9 * t, e64, Float64(ebig), 1e6 * bt / length(bx))
    end
end
