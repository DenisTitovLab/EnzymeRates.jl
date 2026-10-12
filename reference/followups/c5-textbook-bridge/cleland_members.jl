# ABOUTME: C5 follow-up: fitted count and identifiable rank of the laws P4 found to be Cleland textbook
# ABOUTME: forms, plus the Theorell-Chance formula dimension; stdlib only.
include("cleland_dims.jl")
println("theorell-chance: 7 free; V2 from V1*Kp*Kiq/(V2*Kia*Kb) = Keq   dimension ",
        dim(theorell_chance, (x, Keq) -> (V1 = x[1], Ka = x[2], Kb = x[3], Kp = x[4], Kq = x[5], Kia = x[6],
            Kiq = x[7], V2 = x[1] * x[4] * x[7] / (Keq * x[6] * x[3]), Keq = Keq), 7, (:A, :B, :P, :Q), rng))
using Serialization
sk = Dict(r.id => r for r in deserialize("/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/proto/skeptic_P1-engine/data/rank_pop.jls"))
fam = Dict("ordered" => ["expand_bibi/349", "expand_bibi/358", "handwritten/7", "handwritten/8", "specs/5"],
           "pingpong" => ["expand_bibi/303", "expand_bibi/304", "expand_bibi/305", "handwritten/11",
                          "init_bibi_pp/261", "init_bibi_pp/262", "specs/7"],
           "re_random" => ["init_bibi/55", "init_bibi_pp/55", "specs/13"],
           "theorell_chance (first 8 of 19)" => ["expand_bibi/201", "expand_bibi/202", "handwritten/9",
               "handwritten/19", "init_bibi/4", "init_bibi/16", "init_bibi/36", "init_bibi/42"])
for (f, ids) in sort(collect(fam); by = first)
    println(rpad(f, 32), join(["$id $(sk[id].nθ)/$(sk[id].rank)" for id in ids], ", "))
end
