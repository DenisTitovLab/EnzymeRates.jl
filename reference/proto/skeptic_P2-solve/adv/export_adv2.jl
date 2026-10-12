# ABOUTME: Skeptic of probe P2: builds adversarial plain mechanisms with EnzymeRates and exports
# ABOUTME: them as Pop.PMech (package outputs included): second batch (Keq^(1/2), TC + iso).
using EnzymeRates, Random, Serialization, Printf, Statistics
const PROTO = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/proto"
const OUT = joinpath(PROTO, "skeptic_P2-solve", "data", "adv_pkg2.jls")
include(joinpath(PROTO, "common", "Population.jl"))
include(joinpath(PROTO, "common", "export", "export_lib.jl"))

mech(subs, prods, regs, body) = """
@enzyme_mechanism begin
    substrates: $subs
    products: $prods
    $(isempty(regs) ? "" : "regulators: " * regs)
    steps: begin
$body
    end
end"""

const DEFS = [
    # x = (E <--> F) taken twice per turnover: a Keq^(1/2) dependent in the probe's solve
    "half_keq" => mech("S", "P", "", """
        (E(S) <--> F(S), E <--> F)
        F(S) <--> E + P
        F + S <--> E(S)"""),
    # Theorell-Chance uni-uni closed by a free-enzyme isomerization
    "tc_iso_uni_uni" => mech("S", "P", "", """
        E + S <--> F + P
        F <--> E"""),
    # the fuzz/311 structure: a grouped TC pair and a free isomerization with zero flux
    "tc_pair_iso" => mech("S", "P", "", """
        (E + S <--> F + P, F + S <--> E + P)
        F <--> E
        E + S ⇌ F(S)"""),
]

out = Pop.PMech[]
for (i, (id, src)) in enumerate(DEFS)
    t = time()
    em = try
        Core.eval(Main, Meta.parse(src))
    catch e
        println(rpad(id, 40), "REJECTED by the package: ", first(sprint(showerror, e), 300))
        continue
    end
    pm = try
        to_pmech(ER._concrete(em), "adv2/$i", id; seed = 910_000 + i, do_kcat = false,
                 do_bench = false)
    catch e
        println(rpad(id, 40), "EXPORT FAILED: ", first(sprint(showerror, e), 300))
        continue
    end
    push!(out, pm)
    @printf("%-40s ok=%d fitted=%d ipc=%d %.1fs %s\n", id, pm.ok, length(pm.fitted_params),
            pm.independent_param_count, time() - t,
            isempty(pm.stage_errors) ? "" : string(pm.stage_errors))
    println("    fitted: ", pm.fitted_params)
    for d in pm.dependents; println("    dep ", d); end
    isempty(pm.pass2_rename) || println("    pass2: ", pm.pass2_rename)
    flush(stdout)
end
Pop.save_population("adv2", out; path = OUT)
println("wrote ", OUT, " (", length(out), " mechanisms)")
