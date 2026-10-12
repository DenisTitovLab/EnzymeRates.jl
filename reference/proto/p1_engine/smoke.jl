# ABOUTME: Smoke test of Engine.bare_rate_law on the plain specs against the package's raw polynomials.
# ABOUTME: Prints per-spec agreement and statistics; stdlib only.
const C = joinpath(@__DIR__, "..", "common")
include(joinpath(C, "Oracle.jl")); include(joinpath(C, "Engine.jl"))
topoly(pp) = Dict{Vector{Pair{Symbol, Int}}, Rational{Int}}(mono => c for (mono, c) in pp)
for m in Pop.load_population("specs")
    m.allosteric && continue
    try
        L = Engine.bare_rate_law(m)
        num, den, dfree = Engine.grouped(L, m)
        V = L.vars
        ok = (Poly.symbolic(num, V), Poly.symbolic(den, V), Poly.symbolic(dfree, V)) ==
             (topoly(m.raw_num), topoly(m.raw_den), topoly(m.raw_dfree))
        s = L.stats
        println(rpad(m.id, 10), " ok=", ok, " G=", s[:segments], " blocks=", s[:block_sizes],
                " cyc=", s[:cycles], " den=", s[:den_terms], "/", length(m.raw_den),
                " num=", s[:num_terms], " hill=", s[:hill_terms], " omega=", s[:omega_terms],
                " dead=", s[:dead_end], " reclose=", s[:re_closing], "  ", m.source)
    catch e
        println(rpad(m.id, 10), " ERROR ", sprint(showerror, e)[1:min(end, 300)], "  ", m.source)
    end
end
