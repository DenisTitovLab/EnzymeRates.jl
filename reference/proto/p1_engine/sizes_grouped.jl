# ABOUTME: Denominator sizes after kinetic-group renaming: expanded vs σ-factored vs σ-and-block-factored.
# ABOUTME: Shows what factored code generation would gain on the named (grouped) law.
const C = joinpath(@__DIR__, "..", "common")
include(joinpath(C, "Population.jl")); include(joinpath(C, "Engine.jl"))
using Printf, Statistics
const POPS = ["specs", "init_uniuni", "init_unibi", "init_bibi", "init_bibi_pp", "expand_bibi", "handwritten"]
tot = zeros(Int, 3); better_blk = 0; ratios = Float64[]
for pop in POPS
    t = zeros(Int, 3)
    for m in Pop.load_population(pop)
        m.allosteric && continue
        L = Engine.bare_rate_law(m; check = false)
        sub = Dict(v => Poly.mono(Poly.var!(L.vars, m.step_names[j][i]) => 1)
                   for (j, vs) in enumerate(L.stepvars) for (i, v) in enumerate(vs))
        g(p) = Poly.subst(p, v -> get(sub, v, nothing))
        den = g(L.den)
        sig = [g(s) for s in L.stats[:sigma]]
        F = L.stats[:factors]
        G = length(sig)
        # σ-only factored: Σ|σ_g| + Σ|D_g| with D_g = Π_b F_b(a(b, g)), all named.
        Dg = [Poly.expand(Poly.P{Int}, [g(F[b][L.stats[:attach][b][x]]) for b in eachindex(F)]) for x in 1:G]
        sig_only = sum(Poly.nterms, sig) + (G > 1 ? sum(Poly.nterms, Dg) : 0)
        sig_blk = sum(Poly.nterms, sig) + sum((Poly.nterms(g(f)) for Fb in F for f in values(Fb)); init = 0)
        t .+= (Poly.nterms(den), sig_only, sig_blk)
        global better_blk += sig_blk < sig_only
        push!(ratios, sig_blk / Poly.nterms(den))
    end
    tot .+= t
    @printf("%-13s named den terms: expanded %6d | σ-factored %6d | σ+block-factored %6d\n", pop, t...)
end
@printf("all           named den terms: expanded %6d | σ-factored %6d | σ+block-factored %6d; blocks help in %d mechanisms; σ+block/expanded median %.2f\n",
        tot..., better_blk, median(ratios))
