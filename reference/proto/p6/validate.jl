# ABOUTME: Checks the P6 graph type against the package's exported structure on every population:
# ABOUTME: lossless round trip, rendered names, orientation, inactive mask, metabolite lists, cycles.
include(joinpath(@__DIR__, "harness.jl"))
using Printf

ms = _testhelper_load_all()
bad = Dict{Symbol, Vector{String}}()
flag!(k, id) = push!(get!(bad, k, String[]), id)
ndiff_index = 0; nsteps = 0
tierhist = zeros(Int, 4); tier3_ids = String[]
recycles = 0
nplain = 0
for m in ms
    g, edge_of, flipped = _testhelper_build(_testhelper_toraw(m))
    global nsteps += length(m.steps)
    # Lossless: each package form maps to exactly one GM form, and each step's ends and lists match.
    formof = zeros(Int, length(m.species))
    ok = length(g.forms) == length(m.species) && length(g.edges) == length(m.steps)
    for (j, s) in enumerate(m.steps)
        e = g.edges[edge_of[j]]
        a, b, c, r = flipped[j] ? (e.b, e.a, e.rel, e.cons) : (e.a, e.b, e.cons, e.rel)
        for (i, f) in ((s.from, a), (s.to, b))
            formof[i] == 0 && (formof[i] = f)
            ok &= formof[i] == f
        end
        ok &= GM.mets_of(g, c) == sort(_testhelper_met.(s.consumed)) &&
              GM.mets_of(g, r) == sort(_testhelper_met.(s.released)) && e.re == s.re
    end
    ok &= allunique(formof) && all(>(0), formof)
    for (j, s) in enumerate(m.steps)                     # plain binding, read from the export
        f, t = m.species[s.from], m.species[s.to]
        plain = length(s.consumed) == 1 && isempty(s.released) &&
                sort([f.bound; s.consumed]; by = x -> (x.name, x.role)) == sort(t.bound; by = x -> (x.name, x.role)) &&
                sort(f.res_added) == sort(t.res_added) && sort(f.res_subtracted) == sort(t.res_subtracted)
        GM.is_plain_binding(g, g.edges[edge_of[j]]) == plain || flag!(:plain_binding, m.id)
        global nplain += plain
    end
    ok || flag!(:roundtrip, m.id)
    all(GM.form_name(g, formof[i]) == string(m.species[i].name) for i in eachindex(m.species)) ||
        flag!(:names, m.id)
    # Orientation: the package stores every step physical-forward.
    fwd = GM.physical_forward(g; tier3 = :name)
    all(fwd[edge_of[j]] == !flipped[j] for j in eachindex(m.steps)) || flag!(:orient_name, m.id)
    fwdi, tiers = GM.orient(g; tier3 = :index)
    global ndiff_index += count(fwdi .!= fwd)
    for t in tiers; tierhist[t + 1] += 1; end
    any(==(3), tiers) && push!(tier3_ids, m.id)
    # Groups: same partition of steps.
    gof = Pop.group_of_step(m)
    pairs_same = all((gof[i] == gof[j]) == (g.group[edge_of[i]] == g.group[edge_of[j]])
                     for i in eachindex(gof), j in eachindex(gof))
    pairs_same || flag!(:groups, m.id)
    if m.allosteric
        all(g.tags[g.group[edge_of[j]]] == m.group_tags[gof[j]] for j in eachindex(gof)) ||
            flag!(:tags, m.id)
        if !isempty(m.inactive_kept)
            mask = GM.inactive_mask(g)
            all(mask[edge_of[j]] == m.inactive_kept[j] for j in eachindex(gof)) ||
                flag!(:inactive_mask, m.id)
        end
    end
    GM.metabolites(g) == m.metabolites || flag!(:metabolites, m.id)
    GM.inhibitors(g) == sort(m.inhibitors) || flag!(:inhibitors, m.id)
    sort(GM.regulators(g)) == sort(m.regulators) || flag!(:regulators, m.id)
    # Analyses: E - V + C fundamental cycles, each closed and a multiple of the net reaction.
    a = try GM.analysis(g) catch e; flag!(:analysis_error, m.id); continue end
    root = collect(eachindex(g.forms))
    find(i) = root[i] == i ? i : (root[i] = find(root[i]))
    for e in g.edges; root[find(e.a)] = find(e.b); end
    ncomp = count(i -> find(i) == i, eachindex(g.forms))
    length(a.cycles) == length(g.edges) - length(g.forms) + ncomp || flag!(:cycle_count, m.id)
    for c in a.cycles
        div = zeros(Int, length(g.forms))
        for sk in c
            e = g.edges[abs(sk)]
            div[e.a] -= sign(sk); div[e.b] += sign(sk)
        end
        all(iszero, div) || flag!(:cycle_open, m.id)
    end
    GM.bottomless(g, a) === nothing || flag!(:bottomless, m.id)
    if !m.allosteric
        remask = BitVector([e.re for e in g.edges])
        nre = length(GM.analyze(g, remask).cycles)
        nre == Pop.re_cycle_closure(m)[1] || flag!(:re_cycles, m.id)
        global recycles += nre
    end
end
@printf("mechanisms %d, steps %d\n", length(ms), nsteps)
for k in (:roundtrip, :plain_binding, :names, :orient_name, :groups, :tags, :inactive_mask, :metabolites,
          :inhibitors, :regulators, :analysis_error, :cycle_count, :cycle_open, :bottomless,
          :re_cycles)
    v = get(bad, k, String[])
    @printf("%-15s failures %4d  %s\n", k, length(v), join(first(v, 8), ", "))
end
@printf("tier-3 by form index instead of rendered name: %d of %d steps change direction\n",
        ndiff_index, nsteps)
@printf("RE cycles across plain mechanisms: %d; plain bindings: %d\n", recycles, nplain)
@printf("steps decided by tier 0/1/2/3: %s; mechanisms with a tier-3 step: %d (%s)\n",
        tierhist, length(tier3_ids), join(first(tier3_ids, 12), ", "))
