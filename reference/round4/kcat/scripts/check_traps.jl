# ABOUTME: Splits the traps of the kcat-zero plain laws by whether their substrate is a culprit of
# ABOUTME: the core law (exact gap < 0 there) and tallies the trap kinds and exits for each.
using Serialization, Printf
const KDATA = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/round4/kcat/data"
rows = [x for l in ("core", "expand_terter", "terter_sample")
        for x in deserialize(joinpath(KDATA, "rows_$l.jls")) if x.ok && x.cls === :zero]
println("kcat-zero laws: ", length(rows), " (plain ", count(x -> !x.allo, rows), ")")
subs = [:A, :B, :D]
tally = Dict{Tuple{Bool, String}, Int}()
withcore = 0
for x in rows
    x.allo && continue
    culprits = Set(subs[j] for j in 1:x.nsubs if x.zpc.gap[j] < 0)
    found = false
    for p in x.pts, t in p.z2.traps
        c = t.X in culprits
        found |= c
        k = !isempty(t.other) ? "other exits" :
            t.holdsP ? (isempty(t.blocked) ? "X + product/residual, no product-binding exit" :
                        "X + product/residual, product-binding exit") :
            (isempty(t.blocked) ? "X only, no product-binding exit" : "X only, product-binding exit")
        tally[(c, k)] = get(tally, (c, k), 0) + 1
    end
    global withcore += found
end
println("plain laws with a trap found for a core-law culprit substrate (at some point): ", withcore)
println("traps (laws x points x substrates) by (substrate is a core-law culprit, kind): ")
foreach(((k, v),) -> println("  ", k, " => ", v), sort!(collect(tally); by = string))
ex = [(x, t) for x in rows if !x.allo for t in x.pts[1].z2.traps
      if t.X in Set(subs[j] for j in 1:x.nsubs if x.zpc.gap[j] < 0) && isempty(t.blocked)]
println("core-culprit traps without a product-binding exit at point 1: ", length(ex))
for (x, t) in ex[1:min(end, 3)]
    println("  ", x.id, " X = ", t.X, " trap ", t.forms, " core gaps ", x.zpc.gap)
    foreach(l -> println("      ", l), x.lines)
end
