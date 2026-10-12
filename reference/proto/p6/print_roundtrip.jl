# ABOUTME: Prints every population mechanism as DSL text from its P6 GMech, parses and lowers the
# ABOUTME: text again, and checks the result equals the original (printer/lowering round trip).
include(joinpath(@__DIR__, "harness.jl"))
include(joinpath(@__DIR__, "DslLower.jl"))
using Printf

"""`m` as mechanism-macro body text."""
function dsl_text(m::GM.GMech)
    allo = !isempty(m.tags)
    out = ["substrates: " * join(m.substrates, ", "), "products: " * join(m.products, ", ")]
    inh = GM.inhibitors(m)
    isempty(inh) || push!(out, (allo ? "catalytic_inhibitors: " : "regulators: ") * join(inh, ", "))
    allo && push!(out, "catalytic_multiplicity: $(m.multiplicity)")
    push!(out, (allo ? "catalytic_steps" : "steps") * ": begin")
    for g in 1:GM.ngroups(m)
        ss = [GM.step_text(m, e) for (k, e) in enumerate(m.edges) if m.group[k] == g]
        body = length(ss) == 1 ? only(ss) : "(" * join(ss, ", ") * ")"
        push!(out, "    " * (allo ? (length(ss) == 1 ? "($body)" : body) * " :: $(m.tags[g])" : body))
    end
    push!(out, "end")
    if !isempty(m.sites)
        push!(out, "allosteric_regulators: " *
              join(("$l::$st" for s in m.sites for (l, st) in zip(s.ligands, s.states)), ", "))
        for s in m.sites
            push!(out, "regulatory_site(multiplicity = $(s.multiplicity)): begin\n    ligands: " *
                       join(s.ligands, ", ") * "\nend")
        end
    end
    "begin\n" * join(out, "\n") * "\nend"
end

ms = _testhelper_load_all()
bad = String[]
for m in ms
    g = _testhelper_gm(_testhelper_toraw(m))
    h = try _testhelper_gm(DslLower.lower(Meta.parse(dsl_text(g)))) catch e; push!(bad, m.id); continue end
    h == g || push!(bad, m.id)
end
@printf("printed, parsed and lowered back equal: %d of %d; failures: %s\n", length(ms) - length(bad),
        length(ms), join(first(bad, 8), ", "))
println(dsl_text(_testhelper_gm(_testhelper_toraw(only(filter(m -> m.id == "allo_uniuni/12", ms))))))
