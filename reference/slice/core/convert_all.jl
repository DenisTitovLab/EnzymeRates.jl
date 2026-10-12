# ABOUTME: Converts every population and adversarial set to Slice mechanisms and saves them in
# ABOUTME: slice/core/data (Vector{Slice.Converted} per set), with a manifest and the failures.
const PROTO = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/proto"
const OUT = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/slice/core/data"
include(joinpath(PROTO, "common", "Population.jl"))
include(joinpath(@__DIR__, "..", "src", "Slice.jl"))
include(joinpath(@__DIR__, "hand_sets.jl"))
include(joinpath(@__DIR__, "retag.jl"))
using Serialization, Printf

common(n) = n => joinpath(PROTO, "common", "data", n * ".jls")
skp(dir, n, f = n) = n => joinpath(PROTO, dir, "data", f * ".jls")
mkpath(OUT)
serialize(joinpath(OUT, "pop_p1_direct.jls"), _testhelper_p1_direct())
serialize(joinpath(OUT, "pop_p2_hand.jls"), _testhelper_p2_hand())
const SETS = [
    common.(["specs", "handwritten", "init_uniuni", "init_unibi", "init_bibi",
             "init_bibi_pp", "expand_bibi", "allo_uniuni", "allo_bibi"]);
    [skp("skeptic_P5a-onlyA", n) for n in ("allo_uniuni_inh", "allo_unibi", "allo_bibi_pp",
                                           "allo_terbi")];
    [skp("skeptic_P5a-onlyA", "p5a_" * f, f) for f in ("adv_cases", "adv_cases2",
                                                      "adv_cases3", "adv_cases4", "cases")];
    [skp("skeptic_P1-engine", "p1_" * f, f) for f in ("adv_pkg", "adv_pkg2", "enum_inh")];
    "p1_direct" => joinpath(OUT, "pop_p1_direct.jls");
    [skp("skeptic_P2-solve", "p2_" * f, f) for f in ("adv_pkg", "adv_pkg2", "fuzz_dump",
                                                    "fuzz_dump_flip")];
    "p2_hand" => joinpath(OUT, "pop_p2_hand.jls")]

manifest = ["# set\tconverted\tfailed\tsource file"]
for (name, path) in SETS
    ps = deserialize(path)
    eltype(ps) <: Pop.PMech || (println(name, ": not a mechanism list: ", typeof(ps)); continue)
    out, failed = Slice.Converted[], String[]
    for (i, p) in enumerate(ps)
        try
            push!(out, Slice.Converted(name == "p5a_cases" ? retag(p) : p, i))
        catch e
            push!(failed, "$i\t$(p.id)\t$(p.source)\t" * first(sprint(showerror, e), 300))
        end
    end
    serialize(joinpath(OUT, name * ".jls"), out)
    write(joinpath(OUT, name * ".failed.txt"), join(failed, "\n"), isempty(failed) ? "" : "\n")
    push!(manifest, "$name\t$(length(out))\t$(length(failed))\t$path")
    @printf("%-20s %5d converted, %3d failed\n", name, length(out), length(failed))
    foreach(f -> println("    ", f), failed)
end
write(joinpath(OUT, "MANIFEST.tsv"), join(manifest, "\n"), "\n")
