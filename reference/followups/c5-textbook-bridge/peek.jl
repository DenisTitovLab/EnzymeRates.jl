# ABOUTME: C5 follow-up: prints the steps (package DSL style, stored direction), fitted count and rank of
# ABOUTME: chosen plain mechanisms, to describe the rank-deficient and Cleland-form laws.
include("/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/proto/common/Population.jl")
using Serialization
const SKDATA = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/proto/skeptic_P1-engine/data/rank_pop.jls"
sk = Dict(r.id => r for r in deserialize(SKDATA))
byid(id) = (p = split(id, "/"); Pop.load_population(String(p[1]))[parse(Int, p[2])])
met(x) = x.role === :Inh ? "$(x.name)::Inh" : string(x.name)
function form(sp)
    b = join(met.(sp.bound), ", ")
    res = isempty(sp.res_added) && isempty(sp.res_subtracted) ? "" :
          "; residual = " * join(string.(sp.res_added), " + ") * " - " * join(string.(sp.res_subtracted), " - ")
    isempty(b) && isempty(res) ? string(sp.conformation) : "$(sp.conformation)($b$res)"
end
side(m, f, ms) = join([form(m.species[f]); met.(ms)], " + ")
for id in ARGS
    m = byid(id); r = sk[id]
    println("$id  ($(m.source))  fitted $(r.nθ)  rank $(r.rank)  fitted names $(m.fitted_params)")
    for (g, js) in enumerate(m.groups)
        txt = [side(m, m.steps[j].from, m.steps[j].consumed) * (m.steps[j].re ? " ⇌ " : " <--> ") *
               side(m, m.steps[j].to, m.steps[j].released) for j in js]
        println("    ", length(txt) == 1 ? txt[1] : "(" * join(txt, ", ") * ")")
    end
end
