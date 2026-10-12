# ABOUTME: Writes hash(m) of every converted mechanism, run 1 from the saved data in order, run 2
# ABOUTME: rebuilt from the Pop records in reverse order after interning 20,000 extra symbols.
const PROTO = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/proto"
const DATA = joinpath(@__DIR__, "data")
const RUN = ARGS[1]
RUN == "2" && foreach(i -> Symbol("pad", i, "_", rand(UInt32)), 1:20_000)
include(joinpath(PROTO, "common", "Population.jl"))
include(joinpath(@__DIR__, "..", "src", "Slice.jl"))
include(joinpath(@__DIR__, "retag.jl"))
using Serialization

rows = String[]
sets = [split(l, '\t') for l in readlines(joinpath(DATA, "MANIFEST.tsv"))[2:end]]
RUN == "2" && reverse!(sets)
for r in sets
    conv = deserialize(joinpath(DATA, r[1] * ".jls"))
    pops = RUN == "2" ? deserialize(r[4]) : nothing
    order = RUN == "2" ? reverse(conv) : conv
    for c in order
        p = RUN == "2" ? (r[1] == "p5a_cases" ? retag(pops[c.index]) : pops[c.index]) : nothing
        m = RUN == "2" ? first(Slice.canonicalize(Slice.pop_raw(p))) : c.mech
        push!(rows, "$(r[1])/$(c.index)\t$(hash(m))\t$(hash(Slice.sig(m)))")
    end
end
# Control: a struct without an explicit hash, defined in a module built at run time.
ctl = Module(:Control)
Core.eval(ctl, :(struct C; x::Int; end))
push!(rows, "control\t$(hash(Base.invokelatest(ctl.C, 1)))\t0")
write(joinpath(@__DIR__, "logs", "hashes_run$RUN.tsv"), join(sort!(rows), "\n"), "\n")
println("run $RUN: ", length(rows) - 1, " mechanisms hashed")
