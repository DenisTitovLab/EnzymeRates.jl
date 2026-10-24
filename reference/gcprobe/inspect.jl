# ABOUTME: Prints the types held by the pipeline rows and by the population records that the
# ABOUTME: crashing workload deserializes, to see which objects survive between its phases.
include("/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/slice/skeptic/scripts/indep.jl")
rows = deserialize(joinpath(IDATA_IN, "pipe_terter1.jls"))
println(typeof(rows), " ", length(rows))
x = first(r for r in rows if r.ok)
for k in keys(x); println("  row.", k, " :: ", typeof(getfield(x, k))); end
rows2 = [x for k in 1:11 for x in deserialize(joinpath(IDATA_IN, "pipe_terter$k.jls")) if x.ok]
top = sort(rows2; by = x -> -x.nD)[1:40]
println("sets: ", unique(x.set for x in top))
open("/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/phase3/gcprobe/todo.tsv", "w") do io
    for x in top; println(io, x.id, '\t', x.set, '\t', x.index, '\t', x.nD); end
end
ps = deserialize(joinpath(IDATA_IN, "pop_$(top[1].set).jls"))
println(typeof(ps), " ", length(ps)); p = ps[top[1].index]
for k in fieldnames(typeof(p)); println("  pop.", k, " :: ", typeof(getfield(p, k))); end
recs = skeptic_records("expand_bibi")
println(typeof(recs), " ", length(recs)); r = recs[1]
println(typeof(r))
println("pop types: ", typeof(r.p)); for k in fieldnames(typeof(r.p)); println("  ", k, " :: ", typeof(getfield(r.p, k))); end
c = S.Converted(ps[top[1].index], top[1].index)
L = S.rate_law(c.mech); println(typeof(L))
