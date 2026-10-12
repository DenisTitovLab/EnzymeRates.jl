# ABOUTME: P6 enumeration-style children on every population: RE/SS flip, group split and a
# ABOUTME: dead-end binding, each checked equal to the full constructor's result and timed.
include(joinpath(@__DIR__, "harness.jl"))
using Printf, Statistics

ms = _testhelper_load_all()
fails = Dict(k => String[] for k in (:roundtrip, :flip, :split, :deadend))
times = Dict(k => Float64[] for k in (:roundtrip, :flip, :split, :deadend))
napplied = Dict(k => 0 for k in keys(fails))
for (i, m) in enumerate(ms)
    raw0 = _testhelper_toraw(m)
    g, edge_of, flipped = _testhelper_build(raw0)
    gof = Pop.group_of_step(m)
    stepof = invperm(edge_of)
    # Lossless decode: raw(g) rebuilds g.
    napplied[:roundtrip] += 1
    _testhelper_gm(GM.raw(g)) == g || push!(fails[:roundtrip], m.id)
    push!(times[:roundtrip], _testhelper_time_ns(() -> _testhelper_gm(GM.raw(g)), 20))
    # RE/SS flip of group 1, against the full constructor on the exported input.
    pg = gof[stepof[findfirst(==(1), g.group)]]
    re = !m.steps[findfirst(==(pg), gof)].re
    s, p, steps, grp, kw = deepcopy(raw0)
    for j in eachindex(steps)
        grp[j] == pg && (steps[j] = GM.RawStep(steps[j].from, steps[j].to, steps[j].consumed, steps[j].released, re))
    end
    napplied[:flip] += 1
    GM.with_re(g, 1, re) == _testhelper_gm((s, p, steps, grp, kw)) || push!(fails[:flip], m.id)
    push!(times[:flip], _testhelper_time_ns(() -> GM.with_re(g, 1, re), 100))
    # Split: the last edge of the first multi-step group leaves it.
    G = GM.ngroups(g)
    multi = findfirst(l -> count(==(l), g.group) > 1, 1:G)
    if multi !== nothing
        k = findlast(==(multi), g.group)
        s, p, steps, grp, kw = deepcopy(raw0)
        grp[stepof[k]] = G + 1
        tags = isempty(kw.tags) ? kw.tags : [kw.tags; kw.tags[gof[stepof[k]]]]
        napplied[:split] += 1
        GM.with_split(g, [k]) == _testhelper_gm((s, p, steps, grp, merge(kw, (; tags)))) ||
            push!(fails[:split], m.id)
        push!(times[:split], _testhelper_time_ns(() -> GM.with_split(g, [k]), 100))
    end
    # Dead-end: a new inhibitor binds free E in a group of its own.
    free = findfirst(GM.is_free, g.forms)
    if free !== nothing
        conf = g.confs[g.forms[free].conf]
        st = GM.RawStep(GM.RawForm(conf, GM.Met[], Symbol[], Symbol[]),
                        GM.RawForm(conf, [GM.Met(:Inew, true)], Symbol[], Symbol[]),
                        [GM.Met(:Inew, true)], GM.Met[], true)
        child(r) = (r[1], r[2], [r[3]; st], [r[4]; maximum(r[4]) + 1],
                    merge(r[5], (; tags = isempty(r[5].tags) ? r[5].tags : [r[5].tags; :EqualAI])))
        napplied[:deadend] += 1
        _testhelper_gm(child(GM.raw(g))) == _testhelper_gm(child(raw0)) || push!(fails[:deadend], m.id)
        push!(times[:deadend], _testhelper_time_ns(() -> _testhelper_gm(child(GM.raw(g))), 20))
    end
end
println("operation   applied  mismatches  median ns  max ns")
for k in (:roundtrip, :flip, :split, :deadend)
    @printf("%-10s %8d %11d %10.1f %7.1f  %s\n", k, napplied[k], length(fails[k]),
            median(times[k]), maximum(times[k]), join(first(fails[k], 5), ", "))
end
