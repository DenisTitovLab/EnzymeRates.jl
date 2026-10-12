# ABOUTME: Dumps the census face fields of named records and the face flags of every :OnlyA-rejected
# ABOUTME: candidate, to compare the census's face check with the oracle on G1 candidates.
using Serialization
const CEN = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/census/small"
for (name, gids) in (("bi_bi_rxn.r2_100000", [17966, 17968, 14463]), ("bi_bi_pp_rxn", [3048]))
    d = deserialize(joinpath(CEN, "data", "$name.jls"))
    for g in gids
        r = d.recs[g]
        println("#$g $name faces=$(r[:faces]) kinds=$(r[:face_kind]) dvan=$(r[:dvan])")
        println("   face=$(r[:face])")
        println("   admit=$(r[:admit])")
    end
    adm = [r for r in d.recs if r[:status] === :ok && r[:fitted] && r[:admit] !== nothing]
    k(msg) = occursin("still join", msg) ? "G1_inside" : occursin("close a cycle", msg) ? "G1_cycle" : "island"
    for kind in ("G1_inside", "G1_cycle", "island")
        rs = [r for r in adm if k(r[:admit]) == kind]
        isempty(rs) && continue
        nf = count(r -> any(!=(:none), r[:face_kind]) || any(r[:dvan]), rs)
        println("$name $kind: $(length(rs)) rejected; with a census face flag: $nf")
    end
end
