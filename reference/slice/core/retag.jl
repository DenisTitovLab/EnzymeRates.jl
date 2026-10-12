# ABOUTME: Restores the tags of the P5a grid cases the package rejected (exported from a
# ABOUTME: :NonequalAI twin), as P5a's own `retag` does; shared by the core-stage scripts.

"""P5a's `retag`: a rejected grid case (`uni:` source) was exported from a `:NonequalAI`
twin, so its tags come back from the source string; any other rejected case lost them."""
function retag(p)
    p.ok && return p
    startswith(p.source, "uni:") || error("a rejected case exported from a :NonequalAI " *
                                         "twin; its tags are not in the record")
    tS, tC, tP = Symbol.(split(split(p.source, ":")[3], ","))
    d = Dict(f => getfield(p, f) for f in fieldnames(Pop.PMech))
    d[:group_tags] = [any(x -> x.name === :S, p.steps[g[1]].consumed) ? tS :
                      any(x -> x.name === :P, p.steps[g[1]].consumed) ? tP : tC
                      for g in p.groups]
    Pop.PMech(; d...)
end
