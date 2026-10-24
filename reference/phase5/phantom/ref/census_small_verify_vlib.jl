# ABOUTME: Independent verification library for the small census: package-law evaluation through
# ABOUTME: the compiled rate_equation, complex-step ranks, same-law fits, growth scans, admission, oracle.
using EnzymeRates, Serialization, Random, Printf, LinearAlgebra, Statistics
const ER = EnzymeRates
const SCR = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad"
const CEN = joinpath(SCR, "census", "small")
const VER = joinpath(CEN, "verify")

errmsg(e) = replace(first(sprint(showerror, e), 300), '\n' => ' ')

# ─── reactions (the census's test constants) ─────────────────────────────────────────────
const TESTFILE = joinpath(dirname(dirname(pathof(EnzymeRates))), "test", "test_mechanism_enumeration.jl")
const RXNS = (:uni_uni_rxn, :uni_uni_with_reg, :uni_uni_allo, :uni_uni_allo_reg,
              :uni_uni_allo_2reg, :uni_uni_reg_and_inhibitor, :uni_bi_rxn, :bi_bi_rxn,
              :bi_bi_pp_rxn)
let ast = Meta.parseall(read(TESTFILE, String))
    for ex in ast.args
        ex isa Expr && ex.head === :const && ex.args[1] isa Expr &&
            ex.args[1].args[1] in RXNS && Core.eval(Main, ex)
    end
end

loadrun(name) = (deserialize(joinpath(CEN, "data", "$name.jls")),
                 deserialize(joinpath(CEN, "data", "$name.mechs.jls")))

# ─── mechanism structure from the package types ─────────────────────────────────────────
mkey(x) = string(nameof(typeof(x)), ":", ER.name(x))
isfree(sp) = isempty(ER.bound(sp)) && isempty(ER.residual(sp))

"""Forms, flat steps (from, to, consumed, released, RE flag, group), tags and metabolites of a
package mechanism."""
function mstruct(m)
    gs = ER.steps(m)
    forms = ER.Species[]
    fidx = Dict{ER.Species, Int}()
    addf(s) = get!(() -> (push!(forms, s); length(forms)), fidx, s)
    steps = NamedTuple[]
    for (g, grp) in enumerate(gs), s in grp
        a, b = addf(ER.from_species(s)), addf(ER.to_species(s))
        push!(steps, (a = a, b = b, cons = Symbol[ER.name(x) for x in ER.consumed(s)],
                      rel = Symbol[ER.name(x) for x in ER.released(s)],
                      re = ER.is_equilibrium(s), g = g, step = s))
    end
    rxn = ER.reaction(m)
    allo = m isa ER.AllostericMechanism
    mets = unique!(Symbol[x for s in steps for x in [s.cons; s.rel]])
    frees = [i for i in eachindex(forms) if isfree(forms[i])]
    (forms = forms, steps = steps, ngroups = length(gs),
     tags = allo ? collect(Symbol, ER.cat_allo_states(m)) : Symbol[],
     n = allo ? ER.catalytic_multiplicity(m) : 1,
     subs = Symbol[ER.name(x) for x in ER.substrates(rxn)],
     prods = Symbol[ER.name(x) for x in ER.products(rxn)],
     mets = mets, frees = frees, allo = allo)
end

"""A plain binding: one metabolite consumed, none released, the bound form is the unbound one
plus that metabolite (same conformation, same residual)."""
function plain_binding(ms, s)
    length(s.cons) == 1 && isempty(s.rel) || return false
    fa, fb = ms.forms[s.a], ms.forms[s.b]
    st = s.step
    x = only(ER.consumed(st))
    ER.conformation(fa) == ER.conformation(fb) && ER.residual(fa) == ER.residual(fb) &&
        sort(mkey.(ER.bound(fb))) == sort([mkey.(ER.bound(fa)); mkey(x)])
end

# ─── admission, re-implemented from the definition (every violation, not the first) ─────
"""The :OnlyA admission violations of an allosteric mechanism: `nofree`, `secondfree`,
`G1_inside` (an :OnlyA plain binding whose sides the non-:OnlyA steps join), `G1_cycle` (the
bindings' arcs between components close a directed cycle), `island` (a form whose component
the free form's component does not reach along the arcs). Empty for a plain mechanism."""
function my_admit(ms)
    ms.allo || return String[]
    v = String[]
    nf = length(ms.forms)
    kept = [ms.tags[s.g] !== :OnlyA for s in ms.steps]
    isempty(ms.frees) && return ["nofree"]
    flip = ms.frees[1]
    uf(pairs) = begin
        p = collect(1:nf)
        f(i) = (while p[i] != i; i = p[i]; end; i)
        for (a, b) in pairs
            p[f(a)] = f(b)
        end
        [f(i) for i in 1:nf]
    end
    reseg = uf([(s.a, s.b) for (k, s) in enumerate(ms.steps) if kept[k] && s.re])
    any(i -> reseg[i] != reseg[flip], ms.frees) && push!(v, "secondfree")
    comp = uf([(s.a, s.b) for (k, s) in enumerate(ms.steps) if kept[k]])
    arcs = [(comp[s.a], comp[s.b]) for (k, s) in enumerate(ms.steps)
            if !kept[k] && plain_binding(ms, s)]
    any(a -> a[1] == a[2], arcs) && push!(v, "G1_inside")
    # directed cycle among components (self-loops excluded): Kahn's algorithm
    cs = unique(comp)
    ar = [(a, b) for (a, b) in arcs if a != b]
    indeg = Dict(c => 0 for c in cs)
    foreach(((a, b),) -> indeg[b] += 1, ar)
    todo = [c for c in cs if indeg[c] == 0]
    while !isempty(todo)
        u = pop!(todo)
        for (a, b) in ar
            a == u && (indeg[b] -= 1) == 0 && push!(todo, b)
        end
    end
    any(c -> indeg[c] > 0, cs) && push!(v, "G1_cycle")
    seen = Set([comp[flip]])
    todo = [comp[flip]]
    while !isempty(todo)
        u = pop!(todo)
        for (a, b) in ar
            a == u && !(b in seen) && (push!(seen, b); push!(todo, b))
        end
    end
    any(i -> !(comp[i] in seen), 1:nf) && push!(v, "island")
    v
end

# ─── the package law, through the compiled rate_equation ────────────────────────────────
struct PLaw
    em::Any
    names::Vector{Symbol}
    mets::Vector{Symbol}
end
function plaw(m)
    em = ER.compile_mechanism(m)
    PLaw(em, collect(Symbol, ER.fitted_params(em)), collect(Symbol, ER.metabolites(em)))
end

"""v at fitted values `θ` (in `pl.names` order), `keq`, concentrations `conc` (Dict by name,
missing names at zero), E_total = 1."""
function prate(pl::PLaw, θ::AbstractVector, keq, conc::AbstractDict)
    T = promote_type(eltype(θ), typeof(keq), valtype(conc))
    p = NamedTuple{(pl.names..., :Keq, :E_total)}((T.(θ)..., T(keq), one(T)))
    c = NamedTuple{(pl.mets...,)}(Tuple(T(get(conc, x, zero(T))) for x in pl.mets))
    ER.rate_equation(pl.em, c, p)
end

lu10(rng, lo, hi) = 10.0^(lo + (hi - lo) * rand(rng))
randconc(rng, mets; lo = -1.0, hi = 1.0) = Dict{Symbol, Float64}(x => lu10(rng, lo, hi) for x in mets)

"""Row of ∂ log v/∂ log θ at one point by complex step (exact to rounding)."""
function logjac_row(pl, θ, keq, conc; h = 1e-30)
    v = real(prate(pl, θ, keq, conc))
    row = zeros(length(θ))
    for j in eachindex(θ)
        θc = complex.(θ)
        θc[j] += im * θ[j] * h
        row[j] = imag(prate(pl, θc, keq, conc)) / h / v
    end
    row, v
end

"""The identifiable rank of the package law: singular values of the complex-step log-Jacobian
over `n + extra` random points (θ, Keq, concentrations 10^[lo, hi])."""
function cs_rank(pl, rng; extra = 10, lo = -2.0, hi = 2.0, tol = 1e-8)
    n = length(pl.names)
    n == 0 && return (rank = 0, n = 0, kept = 1.0, dropped = 0.0)
    θ = [lu10(rng, lo, hi) for _ in 1:n]
    keq = lu10(rng, lo, hi)
    J = zeros(n + extra, n)
    for r in 1:(n + extra)
        J[r, :] .= first(logjac_row(pl, θ, keq, randconc(rng, pl.mets; lo, hi)))
    end
    all(isfinite, J) || return (rank = -1, n = n, kept = NaN, dropped = NaN)
    sv = svdvals(J)
    rk = count(>(tol * sv[1]), sv)
    (rank = rk, n = n, kept = sv[rk] / sv[1], dropped = rk < n ? sv[rk + 1] / sv[1] : 0.0)
end

"""Best residual of fitting law `b` to law `a` at fixed θa (relative rate residuals at
`npts` random points, Levenberg–Marquardt on log θb, complex-step Jacobian, several starts;
the first start takes θa's values where the names agree)."""
function fit_into(a::PLaw, θa, b::PLaw, keq, rng; npts = 0, starts = 8, iters = 300)
    nb = length(b.names)
    mets = union(a.mets, b.mets)
    np_ = npts > 0 ? npts : 2 * max(nb, length(a.names)) + 12
    pts = [randconc(rng, mets) for _ in 1:np_]
    target = [real(prate(a, θa, keq, c)) for c in pts]
    sc = abs.(target) .+ 1e-14
    function rj(u)
        θ = exp.(clamp.(u, -50, 50))
        r = zeros(np_); J = zeros(np_, nb)
        for (i, c) in enumerate(pts)
            row, v = logjac_row(b, θ, keq, c)
            r[i] = (v - target[i]) / sc[i]
            J[i, :] .= row .* v ./ sc[i]
        end
        r, J
    end
    namea = Dict(x => j for (j, x) in enumerate(a.names))
    best = Inf
    for st in 1:starts
        u = st == 1 ? [haskey(namea, x) ? log(θa[namea[x]]) : randn(rng) for x in b.names] :
                      randn(rng, nb) .* 2.0
        nb == 0 && (best = maximum(abs, [(real(prate(b, Float64[], keq, c)) - t) / s
                                          for (c, t, s) in zip(pts, target, sc)]); break)
        r, J = rj(u)
        f = sum(abs2, r)
        λ = 1e-3
        for _ in 1:iters
            all(isfinite, J) || break
            A = J' * J
            δ = -((A + λ * Diagonal(diag(A) .+ 1e-14)) \ (J' * r))
            r2, J2 = rj(u .+ δ)
            f2 = sum(abs2, r2)
            if isfinite(f2) && f2 < f
                u .+= δ; r, J, f = r2, J2, f2; λ = max(λ / 3, 1e-14)
            else
                λ *= 5
            end
            (maximum(abs, r) < 1e-11 || λ > 1e14) && break
        end
        best = min(best, maximum(abs, r))
        best < 1e-9 && break
    end
    best
end

"""Same family both ways at two random θ each: (a in b, b in a) best residuals."""
function same_family(a::PLaw, b::PLaw, rng; ndraw = 2)
    keq = lu10(rng, -1, 1)
    ab = maximum(fit_into(a, [lu10(rng, -1, 1) for _ in a.names], b, keq, rng) for _ in 1:ndraw)
    ba = maximum(fit_into(b, [lu10(rng, -1, 1) for _ in b.names], a, keq, rng) for _ in 1:ndraw)
    ab, ba
end

"""Identity check when the two laws have the same fitted names: max relative difference at
random points with the same named values."""
function same_by_name(a::PLaw, b::PLaw, rng; npts = 20)
    Set(a.names) == Set(b.names) || return NaN
    θ = Dict(x => lu10(rng, -1, 1) for x in a.names)
    keq = lu10(rng, -1, 1)
    w = 0.0
    for _ in 1:npts
        c = randconc(rng, union(a.mets, b.mets))
        va = real(prate(a, [θ[x] for x in a.names], keq, c))
        vb = real(prate(b, [θ[x] for x in b.names], keq, c))
        w = max(w, abs(va - vb) / max(abs(va), abs(vb), 1e-300))
    end
    w
end

"""Growth scan: the largest log-slope d log v/d log t between t = 1e9 and 1e12 along rays:
all substrates together (products, modifiers 0) and each substrate alone (others at random
fixed values), and the same for products (forward = substrates)."""
function growth(pl, ms, rng; ndraw = 2)
    out = Dict{Symbol, Float64}(:fwd => -Inf, :rev => -Inf)
    for _ in 1:ndraw
        θ = [lu10(rng, -1, 1) for _ in pl.names]
        keq = lu10(rng, -1, 1)
        for (side, grow) in ((:fwd, ms.subs), (:rev, ms.prods))
            base = Dict(x => lu10(rng, -1, 1) for x in grow)
            for ray in [Any[grow]; Any[[x] for x in grow if length(grow) > 1]]
                v(t) = (c = Dict{Symbol, Float64}(x => (x in ray ? t : 1.0) * base[x] for x in grow);
                        abs(real(prate(pl, θ, keq, c))))
                s = log10(v(1e12) / v(1e9)) / 3
                isfinite(s) && (out[side] = max(out[side], s))
            end
        end
    end
    out
end

"""Best max relative residual of fitting law `b` (log θ, Levenberg–Marquardt, complex-step
Jacobian, `starts` random starts, first from `u0` when given) to the values `target` at the
concentration points `pts` with Keq fixed."""
function fit_values(b::PLaw, keq, pts, target, rng; starts = 12, iters = 400, u0 = nothing)
    nb = length(b.names)
    np_ = length(pts)
    sc = abs.(target) .+ 1e-14
    function rj(u)
        θ = exp.(clamp.(u, -50, 50))
        r = zeros(np_); J = zeros(np_, nb)
        for (i, c) in enumerate(pts)
            row, v = logjac_row(b, θ, keq, c)
            r[i] = (v - target[i]) / sc[i]
            J[i, :] .= row .* v ./ sc[i]
        end
        r, J
    end
    best = Inf
    for st in 1:starts
        u = st == 1 && u0 !== nothing ? copy(u0) : randn(rng, nb) .* 2.0
        r, J = rj(u)
        f = sum(abs2, r)
        isfinite(f) || continue
        λ = 1e-3
        for _ in 1:iters
            all(isfinite, J) || break
            A = J' * J
            δ = -((A + λ * Diagonal(diag(A) .+ 1e-14)) \ (J' * r))
            r2, J2 = rj(u .+ δ)
            f2 = sum(abs2, r2)
            if isfinite(f2) && f2 < f
                u .+= δ; r, J, f = r2, J2, f2; λ = max(λ / 3, 1e-14)
            else
                λ *= 5
            end
            (maximum(abs, r) < 1e-11 || λ > 1e14) && break
        end
        best = min(best, maximum(abs, r))
        best < 1e-9 && break
    end
    best
end
