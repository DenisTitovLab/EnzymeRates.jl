# ABOUTME: Whether a rank-deficient child has the same law as its parent, tested both ways by
# ABOUTME: Levenberg–Marquardt on rate values (not coefficients): parent fits child, child fits parent.
include(joinpath(@__DIR__, "vlib.jl"))

const RX = (:ter_bi_rxn, :ter_ter_rxn, :pyruvate_carboxylase_rxn)

include(joinpath(@__DIR__, "vfit.jl"))

rng = Xoshiro(20261011 + 5)
for name in RX
    data = load_run(name)
    byid = Dict(r[:id] => r for recs in data.rounds for r in recs)
    fit = [r for recs in data.rounds for r in recs if r[:fitted]]
    Lph = [r for r in fit if haskey(get(r, :explain, Dict()), :L_phantom) &&
           r[:parent] != 0]
    other = [r for r in fit if :P in r[:labels] && !haskey(get(r, :explain, Dict()), :L_phantom) &&
             r[:parent] != 0]
    full = [r for r in fit if !(:P in r[:labels]) && !(:Ptoday in r[:labels]) &&
            get(r, :allo, false) && !(:D in r[:labels])]
    println("== $name: L-phantom children $(length(Lph)), other P children $(length(other)), full-rank non-D allosteric $(length(full))")
    for (lab, pool, k) in (("L_phantom", Lph, 6), ("other P", other, 3), ("full rank", full, 4))
        isempty(pool) && continue
        for r in pool[sort!(randperm(rng, length(pool))[1:min(k, length(pool))])]
            p = byid[r[:parent]]
            pc, cp = both_ways(r[:smech], p[:smech], rng)
            @printf("   %-9s %s (slice %d rank %d) vs parent %s (slice %d rank %d): parent fits child %.1e, child fits parent %.1e -> %s\n",
                    lab, r[:id], r[:n_slice], r[:rank], p[:id], p[:n_slice], p[:rank], pc, cp,
                    pc < 1e-8 && cp < 1e-8 ? "SAME LAW" :
                    cp < 1e-8 ? "parent's law inside child's" : "different")
        end
        flush(stdout)
    end
end
