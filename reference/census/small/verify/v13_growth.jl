# ABOUTME: Follows the one growth-scan hit (bi_bi #2349) and the next largest slopes out to t = 1e30 in
# ABOUTME: BigFloat along every ray, to tell slow saturation from unbounded growth.
include(joinpath(@__DIR__, "vlib.jl"))
setprecision(BigFloat, 256)
rng = Xoshiro(131313)
for (name, g) in (("bi_bi_rxn", 2349), ("bi_bi_rxn.r2_100000", 13492), ("bi_bi_pp_rxn.r2_1200", 9915))
    d, mechs = loadrun(name)
    m = mechs[g]
    ms = mstruct(m)
    pl = plaw(m)
    println("\n== $name #$g ", d.recs[g][:text])
    s = ER.rate_equation_string(ER.compile_mechanism(m))
    println(first(s, 2500))
    for draw in 1:3
        θ = [big(lu10(rng, -1, 1)) for _ in pl.names]
        keq = big(lu10(rng, -1, 1))
        for (side, grow) in ((:fwd, ms.subs), (:rev, ms.prods))
            base = Dict(x => big(lu10(rng, -1, 1)) for x in grow)
            for ray in [Any[grow]; Any[[x] for x in grow]]
                vs = BigFloat[]
                for k in 0:3:30
                    c = Dict{Symbol, BigFloat}(x => (x in ray ? big(10)^k : big(1)) * base[x] for x in grow)
                    push!(vs, abs(real(prate(pl, θ, keq, c))))
                end
                slopes = [Float64(log10(vs[i + 1] / vs[i]) / 3) for i in 1:length(vs) - 1]
                println("  draw $draw $side ray $(join(ray, "+")): log-slopes per 3 decades from 1 to 1e30: ",
                        join((@sprintf("%.2g", x) for x in slopes), " "))
            end
        end
    end
end
