# ABOUTME: Bridge from P7.Law to Gen.RateLaw (the type-level positional bytecode) for the P7 probe.
# ABOUTME: `code(law, opt)` is the Int vector; `encode` wraps it as the RateLaw singleton.
isdefined(Main, :P7) || include(joinpath(@__DIR__, "P7.jl"))
isdefined(Main, :Gen) || include(joinpath(@__DIR__, "Gen.jl"))

"""The bytecode of `law` (see `Gen.RateLaw`) under codegen options `opt = (recip, horner, salt)`."""
function code(law::P7.Law, opt)
    c = Int[law.np, law.nc, length(law.deps), length(law.nodes), law.num, law.den, law.scale,
            opt[1], opt[2], opt[3]]
    pairs!(ps) = (push!(c, length(ps)); for (a, b) in ps; push!(c, a, b); end)
    foreach(pairs!, law.deps)
    for nd in law.nodes
        push!(c, length(nd))
        for t in nd
            push!(c, t.coef); pairs!(t.mono); pairs!(t.facs)
        end
    end
    c
end

encode(law::P7.Law, opt) = Gen.RateLaw{Tuple(code(law, opt))}()
