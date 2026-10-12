# ABOUTME: Cancels non-monomial factors common to every component of a unit of a Canon.Law (probe
# ABOUTME: P8): finds the factor's variables by modular univariate gcds, then divides exactly.
isdefined(Main, :Canon) || include(joinpath(@__DIR__, "Canon.jl"))

module Factor

import ..Canon
using Random

const IPoly = Dict{Canon.Mono, BigInt}

"""Variables that a common factor of the polynomials `ps` involves: those `x` for which the
gcd of the univariate images (every other name at a random value modulo a prime) has positive
degree. A false alarm has probability about degree/prime."""
function factor_vars(ps, rng)
    vars = sort!(unique(s for p in ps for m in keys(p) for (s, _) in m))
    filter(vars) do x
        vals = Dict(s => rand(rng, 2:Canon.PRIME - 1) for s in vars if s !== x)
        length(foldl(Canon.gcd_mod, (Canon.univariate(p, x, vals) for p in ps))) > 1
    end
end

"""Exponent vector of a monomial over the sorted variables `vs` (lex order key)."""
expvec(m, vs) = [(i = findfirst(p -> p.first === v, m); i === nothing ? 0 : m[i].second) for v in vs]

"""Exact quotient of `p` by `f` in Z[vars] (lex leading-term division), or `nothing` when `f`
does not divide `p`."""
function divide(p::IPoly, f::IPoly)
    vs = sort!(unique(s for q in (p, f) for m in keys(q) for (s, _) in m))
    lead(q) = argmax(m -> expvec(m, vs), collect(keys(q)))
    lf = lead(f); ef = expvec(lf, vs); cf = f[lf]
    r, out = copy(p), IPoly()
    while !isempty(r)
        lr = lead(r)
        e = expvec(lr, vs) .- ef
        (all(>=(0), e) && iszero(rem(r[lr], cf))) || return nothing
        t = [v => k for (v, k) in zip(vs, e) if k != 0]
        c = div(r[lr], cf)
        out[t] = c
        for (m, a) in f
            d = Dict(m)
            for (v, k) in t
                d[v] = get(d, v, 0) + k
            end
            k = sort!([s => x for (s, x) in d if x != 0]; by = first)
            nv = get(r, k, big(0)) - c * a
            iszero(nv) ? delete!(r, k) : (r[k] = nv)
        end
    end
    out
end

"""The coefficients of `p` over the monomials of the variables outside `vf`: each a polynomial in
the `vf` variables only."""
function coefficients(p::IPoly, vf)
    out = Dict{Canon.Mono, IPoly}()
    for (m, a) in p
        w = [x for x in m if !(x.first in vf)]
        get!(out, w, IPoly())[[x for x in m if x.first in vf]] = a
    end
    out
end

"""`p` divided by its monomial content and integer content."""
function primitive(p::IPoly)
    vs = unique(s for m in keys(p) for (s, _) in m)
    lo = Dict(v => minimum(get(Dict(m), v, 0) for m in keys(p)) for v in vs)
    g = foldl(gcd, values(p); init = big(0))
    IPoly(sort!([s => e - lo[s] for (s, e) in m if e != lo[s]]; by = first) => div(a, g)
          for (m, a) in p)
end

"""`law` with the non-monomial factors common to all nonzero components of each unit divided
out, and the factors removed (as `(unit, factor)`). A factor `f` in the variables `vf` divides
every coefficient of every component over the monomials of the other variables, so the
candidates are those coefficients' primitive parts, smallest first; a candidate that divides
every coefficient is divided out and the search repeats. `unresolved` lists units whose common
factor no candidate matched."""
function cancel(law::Canon.Law; rng = Xoshiro(7), tries = 10)
    comps = Canon.normalize(law)
    removed, unresolved = Tuple{Int, IPoly}[], Int[]
    for (ui, u) in enumerate(law.units)
        live = [i for i in u if !isempty(comps[i])]
        length(live) < 2 && continue
        while true
            vf = factor_vars([comps[i] for i in live], rng)
            isempty(vf) && break
            coef = Dict(i => coefficients(comps[i], vf) for i in live)
            cands = unique(primitive(c) for i in live for c in values(coef[i]))
            sort!(cands; by = length)
            hit = nothing
            for f in first(filter(f -> length(f) > 1, cands), tries)
                qs = Dict(i => Dict(w => divide(c, f) for (w, c) in coef[i]) for i in live)
                any(isnothing, (x for d in values(qs) for x in values(d))) && continue
                hit = f
                for i in live
                    comps[i] = IPoly(sort!([w; m]; by = first) => a
                                     for (w, c) in qs[i] for (m, a) in c)
                end
                break
            end
            hit === nothing && (push!(unresolved, ui); break)
            push!(removed, (ui, hit))
        end
    end
    Canon.Law([Canon.Poly(m => c // 1 for (m, c) in p) for p in comps], law.fixed, law.units,
              law.params, law.header), removed, unresolved
end

end # module
