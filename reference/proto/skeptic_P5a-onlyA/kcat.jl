# ABOUTME: kcat as a leading-order limit (P5a): substrates → ∞ after regulator/inhibitor corners,
# ABOUTME: on the package's own polynomials, compared with the package's _kcat_forward values.
include(joinpath(@__DIR__, "pkgvals.jl"))
using Printf

"""A Pop.PPoly as a Poly.P{Int} over table `T` (`rename` maps symbols first)."""
function topoly(p, T; rename = Dict{Symbol, Symbol}())
    d = Dict{Poly.Mono, Int}()
    for (mono, c) in p
        isinteger(c) || error("non-integer coefficient $c")
        mn = Poly.mono((Poly.var!(T, get(rename, s, s)) => e for (s, e) in mono)...)
        Poly.addto!(d, mn, Int(c))
    end
    Poly.P{Int}(d)
end

"""`p` without the monomials that hold a variable of `zero` (those concentrations set to 0)."""
drop(p::Poly.P{Int}, zero::Set{Int}) =
    Poly.P{Int}(Dict(mn => c for (mn, c) in p.t if !any(v -> v.first in zero, mn)))

"""Map an expression's leaves through `f`."""
mapexpr(e, f) = e[1] === :leaf ? (:leaf, f(e[2])) : (e[1], Any[mapexpr(x, f) for x in e[2]])

"""kcat by leading order: for every corner (each non-substrate, non-product metabolite or
regulator 0 or saturating), the limit as the saturating corner variables → ∞, then the
substrates → ∞, products 0; the maximum over corners. Returns the value and the per-corner
values."""
function kcat_lead(num, den, T, vals, subs, prods, others)
    vi(s) = T.index[s]
    corners = BigFloat[]
    for mask in 0:(2^length(others) - 1)
        sat = [others[k] for k in eachindex(others) if (mask >> (k - 1)) & 1 == 1]
        zero = Set(vi(s) for s in [prods; setdiff(others, sat)] if haskey(T.index, s))
        n2, d2 = mapexpr(num, p -> drop(p, zero)), mapexpr(den, p -> drop(p, zero))
        d1 = Dict(vi(s) => MWC.Q(-1) for s in sat if haskey(T.index, s))
        ds = Dict(vi(s) => MWC.Q(-1) for s in subs if haskey(T.index, s))
        r = MWC.lead_value(n2, d2, [d1, ds], vals)
        push!(corners, r.v)
    end
    maximum(filter(isfinite, corners); init = big(-Inf)), corners
end

"""The package's rule for a plain mechanism: group the forward (positive) numerator terms
and the denominator by their metabolite monomial, keep the product-free patterns present
in both, and take the largest ratio."""
function kcat_patterns(num::Poly.P{Int}, den::Poly.P{Int}, T, vals, metset, prods)
    key(mn) = [p for p in mn if T.names[p.first] in metset]
    gn, gd = Dict{Any, BigFloat}(), Dict{Any, BigFloat}()
    rest(mn) = [p for p in mn if !(T.names[p.first] in metset)]
    ev(mn, c) = setprecision(() -> c * prod((vals[v]^e for (v, e) in rest(mn)); init = big(1)),
                             BigFloat, MWC.PREC)
    for (mn, c) in num.t
        c > 0 || continue
        gn[key(mn)] = get(gn, key(mn), big(0)) + ev(mn, c)
    end
    for (mn, c) in den.t
        gd[key(mn)] = get(gd, key(mn), big(0)) + ev(mn, c)
    end
    ks = [k for k in keys(gn) if haskey(gd, k) && !any(T.names[p.first] in prods for p in k)]
    isempty(ks) && return big(NaN)
    maximum(gn[k] / gd[k] for k in ks)
end

rows = NamedTuple[]
for name in ("specs", "handwritten", "init_uniuni", "init_unibi", "init_bibi", "init_bibi_pp",
             "expand_bibi", "allo_uniuni", "allo_bibi")
    pop = Pop.load_population(name)
    for m in pop, (si, s) in enumerate(m.samples[1:1])
        (isnan(s.kcat) && !any(e -> e.first === :kcat, m.stage_errors)) && continue
        vals = pkg_values(m, s)
        T = Poly.Vars()
        if !m.allosteric
            num = MWC.leaf(topoly(m.pkg_num, T)); den = MWC.leaf(topoly(m.pkg_den, T))
            others = setdiff(m.metabolites, [m.substrates; m.products])
        else
            st(p) = topoly(p, T)
            lawA = st.(m.state_polys_A); lawI = st.(m.state_polys_I)
            regsym(lig) = Symbol(lig, "#reg")
            nA = Dict((si2, lig) => t === :OnlyI ? nothing : Symbol("K_A_", lig, "reg")
                      for (si2, site) in enumerate(m.sites) for (lig, t) in zip(site.ligands, site.states))
            nI = Dict((si2, lig) => t === :OnlyA ? nothing : t === :EqualAI ?
                      Symbol("K_A_", lig, "reg") : Symbol("K_I_", lig, "reg")
                      for (si2, site) in enumerate(m.sites) for (lig, t) in zip(site.ligands, site.states))
            # regulatory-site ligands are corner variables of their own, as in _kcat_forward
            msites = Pop.PMech(; id = m.id, source = m.source, substrates = m.substrates,
                products = m.products, species = m.species, steps = m.steps, groups = m.groups,
                sites = [Pop.PSite(regsym.(site.ligands), site.multiplicity, site.states)
                         for site in m.sites])
            nA2 = Dict((k[1], regsym(k[2])) => v for (k, v) in nA)
            nI2 = Dict((k[1], regsym(k[2])) => v for (k, v) in nI)
            tab = MWC.Table(T)
            WA = MWC.W_poly(msites, nA2, tab); WI = MWC.W_poly(msites, nI2, tab)
            Lp = Poly.term(Poly.P{Int}, Poly.mono(Poly.var!(T, :L) => 1))
            num, den = MWC.combine(Tuple(lawA), Tuple(lawI), WA, WI, Lp, m.cat_multiplicity)
            others = [setdiff(m.metabolites, [m.substrates; m.products; m.regulators]);
                      unique(regsym.(m.regulators))]
        end
        v = Vector{BigFloat}(undef, length(T.names))
        for (k, sym) in enumerate(T.names)
            v[k] = haskey(vals, sym) ? vals[sym] : MWC.big1024(1)   # concentrations: any value
        end
        kc, corners = kcat_lead(num, den, T, v, m.substrates, m.products, others)
        kp = m.allosteric ? big(NaN) :
             kcat_patterns(MWC.expand_expr(num), MWC.expand_expr(den), T, v,
                           Set(m.metabolites), Set(m.products))
        push!(rows, (id = m.id, allo = m.allosteric, pkg = s.kcat, lead = Float64(kc),
                     pat = Float64(kp), ncorner = length(corners),
                     err = any(e -> e.first === :kcat, m.stage_errors)))
    end
end

reld(a, b) = abs(a - b) / max(abs(a), abs(b), 1e-300)
println("mechanisms with a package kcat or a kcat error: ", length(rows))
for allo in (false, true)
    rs = filter(r -> r.allo == allo, rows)
    ok = filter(r -> !r.err, rs)
    agree = count(r -> reld(r.lead, r.pkg) < 1e-10, ok)
    @printf("%s: %d mechanisms; leading order == package kcat (1e-10): %d of %d\n",
            allo ? "allosteric" : "plain", length(rs), agree, length(ok))
    allo || @printf("   pattern rule reproduces the package: %d of %d\n",
                    count(r -> reld(r.pat, r.pkg) < 1e-10, ok), length(ok))
    for r in ok
        reld(r.lead, r.pkg) < 1e-10 && continue
        @printf("   differs: %-18s pkg %.6g lead %.6g (ratio %.4g)%s\n", r.id, r.pkg, r.lead,
                r.lead / r.pkg, allo ? "" : @sprintf(" pattern %.6g", r.pat))
    end
    for r in filter(r -> r.err, rs)
        @printf("   package error (no kcat components): %-18s lead %.6g\n", r.id, r.lead)
    end
end
