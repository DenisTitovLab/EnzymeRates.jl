# ABOUTME: Smoke test of the core files: loads Slice, converts the specs population and prints
# ABOUTME: a few mechanisms with their orientation and names; stdlib only.
include("/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/proto/common/Population.jl")
include("/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/slice/src/Slice.jl")
ps = Pop.load_population("specs")
for (i, p) in enumerate(ps)
    c = Slice.Converted(p, i)
    m = c.mech
    i in (1, 7, 21, 22, 35) || continue
    println("== ", p.id, " ", p.source)
    fwd = Slice.orient(m)
    act = Slice.constant_names(m)
    ina = Slice.is_allosteric(m) ? Slice.constant_names(m; inactive = true) : nothing
    for (k, e) in enumerate(m.edges)
        println("  ", k, " g", m.group[k], " ", Slice.step_text(m, Slice.oriented(m, k, fwd)), "  ", act[k], ina === nothing ? "" : ina[k])
    end
    println("  sig length ", length(Slice.sig(m)), " roundtrip ", Slice.decode(Slice.sig(m)...) == m)
    Slice.check_boundary(m)
    println("  regulators ", Slice.regulator_names(m), " metabolites ", Slice.metabolites(m), " pkg ", p.metabolites)
end
