# ABOUTME: Law-family comparison by fitting: a Float64 evaluator of a slice law and a
# ABOUTME: Levenberg–Marquardt fit of one law to another's rate values, both directions.

"""A Float64 evaluator of the slice law of `m`: `v(θ::Dict, keq, conc::Dict)`."""
function law_fn(m)
    L = S.rate_law(m)
    names = collect(S.fitted_params(m))
    function v(θ, keq, conc)
        x = [L.vars.conc[i] ? conc[n] : n === :Keq ? keq : θ[n]
             for (i, n) in enumerate(L.vars.names)]
        law_value(L, Float64.(x))
    end
    names, v
end

"""Fit the law `(names, v)` to the target values `ys` at the points `pts` by LM in log θ,
from `starts`; the best relative RMS residual."""
function lm_fit(names, v, keq, pts, ys, starts; iters = 150)
    np = length(names)
    best = Inf
    for θ0 in starts
        lx = [log(θ0[n]) for n in names]
        res(lx) = (θ = Dict(n => exp(l) for (n, l) in zip(names, lx));
                   [(v(θ, keq, p) - y) / abs(y) for (p, y) in zip(pts, ys)])
        r = res(lx)
        c = sum(abs2, r)
        isfinite(c) || continue
        λ = 1e-3
        for _ in 1:iters
            c < 1e-26 && break
            J = zeros(length(r), np)
            for j in 1:np
                e = zeros(np); e[j] = 1e-6
                J[:, j] = (res(lx .+ e) .- res(lx .- e)) ./ 2e-6
            end
            all(isfinite, J) || break
            A = J' * J
            g = J' * r
            improved = false
            for _ in 1:12
                step = (A + λ * Diagonal(diag(A) .+ 1e-12)) \ g
                r2 = res(lx .- step)
                c2 = sum(abs2, r2)
                if isfinite(c2) && c2 < c
                    lx, r, c = lx .- step, r2, c2
                    λ = max(λ / 5, 1e-15)
                    improved = true
                    break
                end
                λ *= 10
            end
            improved || break
        end
        best = min(best, sqrt(c / length(r)))
        best < 1e-11 && break
    end
    best
end

"""Both directions for child `c` and parent `p`: `p` fit to `c`'s values at a random `c`
point and `c` fit to `p`'s values at a random `p` point (30 concentration points with the
products present, a start from the shared names and 5 random starts)."""
function both_ways(mc, mp, rng)
    out = Float64[]
    for (src, dst) in ((mc, mp), (mp, mc))
        ns, vs = law_fn(src)
        nd, vd = law_fn(dst)
        θs = Dict(n => 10.0^(2rand(rng) - 1) for n in ns)
        keq = 10.0^(2rand(rng) - 1)
        mets = S.metabolites(src)
        pts = [Dict(x => 10.0^(2rand(rng) - 1) for x in mets) for _ in 1:30]
        ys = [vs(θs, keq, p) for p in pts]
        starts = [Dict(n => get(θs, n, 1.0) for n in nd);
                  [Dict(n => 10.0^(2rand(rng) - 1) for n in nd) for _ in 1:5]]
        push!(out, lm_fit(nd, vd, keq, pts, ys, starts))
    end
    out
end

