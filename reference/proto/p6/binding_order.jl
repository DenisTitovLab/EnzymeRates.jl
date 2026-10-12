# ABOUTME: Counts binding-like P6 edges whose take-up direction runs from the higher form index,
# ABOUTME: i.e. how often a pure a < b storage rule would store a binding as its release.
include(joinpath(@__DIR__, "harness.jl"))
nb = 0; nrev = 0; nplainrev = 0
for m in _testhelper_load_all()
    g = _testhelper_gm(_testhelper_toraw(m))
    for e in g.edges
        GM.is_binding(e) || continue
        global nb += 1
        e.a > e.b && (global nrev += 1; GM.is_plain_binding(g, e) && (global nplainrev += 1))
    end
end
println("binding-like edges $nb; stored from the higher form index $nrev (plain bindings among them: $nplainrev)")
