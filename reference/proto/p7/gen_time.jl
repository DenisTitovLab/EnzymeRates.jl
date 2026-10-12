# ABOUTME: Julia-side generator time (decode + body emission, warm) per law, separate from compilation,
# ABOUTME: and the nesting depth of the emitted expressions; P7 probe.
include(joinpath(@__DIR__, "Encode.jl"))
include(joinpath(@__DIR__, "Select.jl"))
using Statistics, Printf
depth(x) = x isa Expr ? 1 + maximum(depth, x.args; init = 0) : 0
set = selected()
for (name, mk, opt) in (("flat", P7.flat_law, (false, false, 0)), ("flatr", P7.flat_law, (true, false, 0)),
                        ("facth", P7.factored_law, (true, true, 0)),
                        ("halh", m -> P7.factored_law(m; expand_den = true), (true, true, 0)))
    codes = [Tuple(code(mk(m), opt)) for (_, m) in set]
    Gen.body(Gen.decode(codes[1]))
    ts = [minimum(@elapsed(Gen.body(Gen.decode(c))) for _ in 1:3) for c in codes]
    ds = [depth(Gen.body(Gen.decode(c))) for c in codes]
    i = argmax(length.(codes))
    @printf("%-6s generator median %.2f ms, max %.2f ms (largest law %s: %.2f ms); expression depth median %d max %d\n",
            name, 1000median(ts), 1000maximum(ts), set[i][2].id, 1000ts[i], median(ds), maximum(ds))
end
