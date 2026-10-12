# ABOUTME: Is the law's face value the exact-face RE limit for some ordering of the RE speeds?
# ABOUTME: Enumerates every extreme speed pattern (each RE step 1e-6 or 1e6) on trap and split cases.
include(joinpath(@__DIR__, "voracle.jl"))
include(joinpath(@__DIR__, "v1_cases.jl"))

function orderings(cs, rng; npar = 3)
    println("\n=== ", cs.name)
    mets = unique!(vcat(cs.subs, cs.prods, [x for s in cs.steps for x in vcat(s[3], s[4])]))
    reidx = findall(s -> s[5], cs.steps)
    for ip in 1:npar
        steps, _ = consistent(cs.steps, rng; mets)
        base = Dict(x => lu10(rng, -1, 1) for x in mets)
        for zs in cs.faces[1:min(end, 2)]
            law = interior_limits(steps, base, zs, cs.ref, rng; ndir = 1)[1][end]
            c0 = copy(base); foreach(z -> c0[z] = big(0), zs)
            vals = BigFloat[]
            for bits in 0:(2^length(reidx) - 1)
                sp = ones(BigFloat, length(steps))
                for (j, k) in enumerate(reidx)
                    sp[k] = (bits >> (j - 1)) & 1 == 1 ? big(10)^6 : big(10)^-6
                end
                push!(vals, rate_at(steps, c0, big(10)^40, sp, cs.ref)[1])
            end
            lo, hi = minimum(vals), maximum(vals)
            near = minimum(abs.(vals .- law)) / abs(law)
            @printf("  [param %d] face %-8s law %-12s extreme orderings: min %-12s max %-12s; closest ordering to the law: rel. gap %s\n",
                    ip, join(zs, ","), fmt(law), fmt(lo), fmt(hi), fmt(near))
        end
    end
end

rng = Xoshiro(7)
for cs in CASES
    cs.name[1] in ('C', 'D', 'E', 'F', 'I') && orderings(cs, rng)
end
