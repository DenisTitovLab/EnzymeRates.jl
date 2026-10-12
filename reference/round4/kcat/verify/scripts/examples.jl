# ABOUTME: Prints the verifier's rows for the probe's worked kcat-zero examples (gaps, traps,
# ABOUTME: exits) so the report's examples can be checked one by one.
using Serialization
const VD = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/round4/kcat/verify"
rows = deserialize(joinpath(VD, "data", "zero_rows.jls"))
for (set, idx) in (("init_terter.part7", 4757), ("init_terter.part7", 5384), ("allo_bibi", 1023),
                   ("init_bibi", 57))
    i = findfirst(x -> x.set == set && x.index == idx, rows)
    i === nothing && (println(set, "/", idx, ": not in my rows"); continue)
    x = rows[i]
    println(set, "/", idx, ": class ", x.cls, " gaps full ", x.zf.gap, " core ", x.zc.gap,
            " ugap full ", x.zf.ugap, " core ", x.zc.ugap, " exact v>0 ", x.zf.pos && x.zf.num0 && x.zf.den0)
    for t in x.traps
        println("   trap X=", t.X, " core culprit ", t.core_culprit, " share ", round(t.share; digits = 6),
                " forms ", t.forms)
        foreach(e -> println("      ", e), t.exits)
    end
    x.cls === :finite && println("   kcat at 2 points ", [p.kc for p in x.pts], "; full law v(1e12 r) ", [p.s1[3] for p in x.pts])
end
