# ABOUTME: Re-examines candidates where the numeric probe saw unbounded growth but condition V
# ABOUTME: says a maximal rate exists: the slice law along a substrate (or product) ray, in BigFloat.
include(joinpath(@__DIR__, "lib.jl"))
for name in split(ARGS[1], ",")
    d = deserialize(joinpath(CEN, "data", "$name.jls"))
    mechs = deserialize(joinpath(CEN, "data", "$name.mechs.jls"))
    odd = [r for r in d.recs if r[:status] === :ok && r[:fitted] &&
           (get(r, :unb_fwd, false) || get(r, :unb_rev, false)) && r[:vmax_fwd] && r[:vmax_rev]]
    println("== $name: $(length(odd)) numeric-unbounded with Vmax both ways")
    for r in odd
        m = mechs[r[:gid]]
        sm = to_slice(m)
        L = S.rate_law(sm, S.solve(sm))
        br = law_branches(L)
        rng = Xoshiro(r[:gid])
        side = get(r, :unb_fwd, false) ? sm.substrates : sm.products
        println("  #$(r[:gid]) r$(r[:round]) $(r[:move]) unbounded $(get(r, :unb_fwd, false) ? "forward" : "reverse")\n    $(r[:text])")
        setprecision(BigFloat, 256) do
            for draw in 1:2
                vars = L.vars
                vals = Dict(x => logu(rng) for x in vars.names if !(x in S.metabolites(sm)))
                vs = BigFloat[]
                for t in (big(10)^3, big(10)^6, big(10)^9, big(10)^12, big(10)^18)
                    x = [vars.conc[i] ? (vars.names[i] in side ? t : big(0)) : vals[vars.names[i]]
                         for i in eachindex(vars.names)]
                    push!(vs, evalbranches(br.N, x) / evalbranches(br.D, x))
                end
                println("    slice law along the ray (1e3, 1e6, 1e9, 1e12, 1e18): ",
                        join((@sprintf("%.4g", Float64(v)) for v in vs), ", "))
            end
        end
        pl = pkg_law(ER.compile_mechanism(m))
        θ = [10.0^(2rand(rng) - 1) for _ in pl.indep]
        keq = 10.0^(2rand(rng) - 1)
        vv = [pkg_rate(pl, θ, keq, Dict(x => t for x in (ER.name.(get(r, :unb_fwd, false) ?
              ER.substrates(ER.reaction(m)) : ER.products(ER.reaction(m))))), identity)
              for t in (1e3, 1e6, 1e9, 1e12)]
        println("    package law in Float64 (1e3, 1e6, 1e9, 1e12): ", join((@sprintf("%.4g", v) for v in vv), ", "))
    end
end
