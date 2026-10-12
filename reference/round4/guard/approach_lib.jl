# ABOUTME: The law along random directions into a face (zs at zero) and the RE-limit oracle at
# ABOUTME: and just off the face: tells a removable 0/0 from a direction-dependent one.

"""Text: the law at ε = 1e-6, 1e-12 along three random directions into the face where `zs`
vanish, whether the limits agree, and the oracle at the face and at ε = 1e-10 (Λ = 1e30)."""
function approach(m, rb, zs, rng)
    setprecision(BigFloat, 512) do
        p, vals = random_params(m, rng)
        merge!(vals, Dict(pairs(p)))
        base = Dict{Symbol, BigFloat}(x => logu(rng, -1, 1) for x in S.metabolites(m))
        foreach(l -> base[l] = big(0.0), S.regulators(m))
        N1, D1 = n1_branches(m, rb)
        lims, rows = BigFloat[], String[]
        for _ in 1:3
            u = Dict(z => logu(rng, -1, 1) for z in zs)
            row = String[]
            for ε in (big(1e-6), big(1e-12))
                c = copy(base)
                foreach(z -> c[z] = ε * u[z], zs)
                x = varvals(rb.vars, c, vals)
                v = evalbranches(N1, x) / evalbranches(D1, x)
                push!(row, @sprintf("%.6g", Float64(v)))
                ε == big(1e-12) && push!(lims, v)
            end
            push!(rows, join(row, "→"))
        end
        sp = ones(BigFloat, length(m.edges) + 1)
        o(c) = try
            Float64(oracle_v(m, vals, c, big(10)^30, sp))
        catch e
            NaN
        end
        c0 = copy(base); foreach(z -> c0[z] = big(0.0), zs)
        c1 = copy(base); foreach(z -> c1[z] = big(1e-10), zs)
        spread = maximum(lims) - minimum(lims)
        kind = spread > 1e-6 * maximum(abs.(lims)) + 1e-14 ? "direction-dependent" : "one limit"
        @sprintf("%s; law %s; oracle at face %.6g, at 1e-10 %.6g", kind, join(rows, " | "),
                 o(c0), o(c1))
    end
end
