# ABOUTME: The numeric growth flags of v4_n.jl rechecked: exact leading exponents of the law along
# ABOUTME: the flagged direction (both sides) and the law at t = 1e60, 1e120, 1e240.
include(joinpath(@__DIR__, "vlib.jl"))
setprecision(BigFloat, 1024)

"""Numerator and denominator of the slice law (plain, or n=1 MWC without sites)."""
function law_polys(L)
    L isa S.Law && return L.H * L.Q, L.D
    T = typeof(L.H).parameters[1]
    QI, DI = S.inactive_terms(L)
    Lp = S.term(T, S.mono(L.vars.index[:L] => 1))
    L.H * (L.fI * L.active.Q + Lp * L.fA * QI), L.fI * L.active.D + Lp * L.fA * DI
end
lead(P, gv, zv, w) = maximum((sum(w .* [Float64(S.mexp(mo, g)) for g in gv])
                              for mo in keys(P.terms) if all(z -> S.mexp(mo, z) == 0, zv));
                             init = -Inf)

lines = readlines(joinpath(VER, "logs", "v4_n.log"))
segs = Dict{Symbol, Vector{String}}()
cur = String[]
for l in lines
    m = match(r"^(\w+_rxn) \(", l)
    if m === nothing
        push!(cur, l)
    else
        segs[Symbol(m[1])] = cur
        global cur = String[]
    end
end
rng = Xoshiro(3)
for name in (:ter_bi_rxn, :ter_ter_rxn, :pyruvate_carboxylase_rxn)
    flags = [m for l in segs[name] for m in (match(r"law grows: \((\d+), (\d+)\) w=\[([^\]]*)\]", l),)
             if m !== nothing]
    isempty(flags) && continue
    data = load_run(name)
    byid = Dict(r[:id] => r for recs in data.rounds for r in recs)
    worst = -Inf
    plateau = 0
    for f in flags
        id = (parse(Int, f[1]), parse(Int, f[2]))
        w = parse.(Float64, split(f[3], ","))
        m = byid[id][:smech]
        L = S.rate_law(m)
        N, D = law_polys(L)
        idx(xs) = [L.vars.index[x] for x in xs]
        sides = [(m.substrates, m.products), (m.products, m.substrates)]
        sides = filter(s -> length(s[1]) == length(w), sides)
        for (grow, zero) in sides
            d = lead(N, idx(grow), idx(zero), w) - lead(D, idx(grow), idx(zero), w)
            worst = max(worst, d)
            vals = random_vals(m, rng)
            base = random_conc(m, rng)
            vs = map((big(10.0)^30, big(10.0)^60, big(10.0)^120, big(10.0)^240)) do t
                c = copy(base)
                foreach(z -> c[z] = big(0.0), zero)
                foreach(((x, e),) -> c[x] = base[x] * t^e, zip(grow, w))
                law_value(L, lawx(L, vals, c))
            end
            conv = abs(vs[4] - vs[3]) < 1e-3 * abs(vs[4])
            conv && (plateau += 1)
            @printf("   %s %s grow %s: exponent gap %.3g; v at t=1e30,1e60,1e120,1e240: %s %s %s %s\n",
                    name, id, grow, d, fmtb(vs[1]), fmtb(vs[2]), fmtb(vs[3]), fmtb(vs[4]))
        end
    end
    @printf("%s: %d numeric flags rechecked; largest exact exponent gap %.3g (positive = unbounded); %d side-scans level off by t = 1e240\n",
            name, length(flags), worst, plateau)
end
