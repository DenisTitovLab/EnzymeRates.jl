# ABOUTME: Inspects the formats of the skeptic and adversarial data files the core stage converts.
# ABOUTME: Prints each file's element type, count and a few ids; stdlib only.
include("/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/proto/common/Population.jl")
using Serialization
P = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/proto"
files = [joinpath(P, "skeptic_P5a-onlyA/data", f) for f in ("allo_uniuni_inh.jls", "allo_unibi.jls", "allo_bibi_pp.jls", "allo_terbi.jls", "adv_cases.jls", "adv_cases2.jls", "adv_cases3.jls", "adv_cases4.jls")]
append!(files, [joinpath(P, "skeptic_P1-engine/data", f) for f in ("adv_pkg.jls", "adv_pkg2.jls", "enum_inh.jls", "adv_check_rows.jls", "rank_pop.jls")])
append!(files, [joinpath(P, "skeptic_P2-solve/data", f) for f in ("adv_pkg.jls", "adv_pkg2.jls", "fuzz_dump.jls", "fuzz_dump_flip.jls")])
for f in files
    x = try deserialize(f) catch e; println(f, ": ERROR ", first(sprint(showerror, e), 200)); continue end
    print(replace(f, P => "<proto>"), ": ", typeof(x), " n=", x isa AbstractVector ? length(x) : "-")
    if x isa AbstractVector && !isempty(x) && first(x) isa Pop.PMech
        println("  ids: ", join([m.id * "(" * m.source[1:min(end, 40)] * ")" for m in x[1:min(end, 4)]], ", "),
                "  ok=", count(m -> m.ok, x), " allo=", count(m -> m.allosteric, x),
                " with step_names=", count(m -> !isempty(m.step_names) || !isempty(m.step_names_A), x))
    else
        println("  first: ", first(sprint(show, x isa AbstractVector ? first(x) : x), 300))
    end
end
