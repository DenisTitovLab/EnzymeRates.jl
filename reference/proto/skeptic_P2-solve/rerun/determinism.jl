# ABOUTME: P2 determinism tests on every plain mechanism: (i) permuted storage gives identical solves;
# ABOUTME: (ii) regroupings that define the same model (split groups, merges along folds) give identical solves.
include(joinpath("/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/proto/common", "Thermo.jl"))
using Random
const T = Thermo

"""`m` with species, step order, group order, member order and step directions shuffled."""
function permuted(m::Pop.PMech, rng)
    ns = length(m.species); π = randperm(rng, ns); inv = invperm(π)   # new i holds old π[i]
    species = m.species[π]
    flip = rand(rng, Bool, length(m.steps))
    σ = randperm(rng, length(m.steps)); σinv = invperm(σ)               # new j holds old σ[j]
    steps = map(σ) do j
        s = m.steps[j]
        flip[j] ? Pop.PStep(inv[s.to], inv[s.from], s.released, s.consumed, s.re) :
                  Pop.PStep(inv[s.from], inv[s.to], s.consumed, s.released, s.re)
    end
    groups = [shuffle(rng, [σinv[j] for j in g]) for g in shuffle(rng, m.groups)]
    Pop.PMech(; id = m.id, source = m.source, substrates = m.substrates, products = m.products,
              inhibitors = m.inhibitors, metabolites = m.metabolites, species = species,
              steps = steps, groups = groups)
end

regrouped(m, groups) = Pop.PMech(; id = m.id, source = m.source, substrates = m.substrates,
    products = m.products, inhibitors = m.inhibitors, metabolites = m.metabolites,
    species = m.species, steps = m.steps, groups = groups)

signature(sol) = (T.independent_names(sol), T.dependent_table(sol),
                  sort!([T.colname(sol, d.col) => d.fold for d in sol.deps]; by = first))
labels(sol) = sort!([T.colname(sol, d.col) => d.label for d in sol.deps]; by = first)

rng = Xoshiro(20261009)
nm = 0; perm_runs = 0; perm_fail = 0; label_perm_fail = 0
split_tried = 0; split_same = 0; split_fail = 0; split_label_changes = 0
allsplit_same = 0; allsplit_fail = 0
merge_tried = 0; merge_fail = 0; merge_count_fail = 0
fails = String[]
for pop in Pop.populations(), m in Pop.load_population(pop)
    m.allosteric && continue
    global nm += 1
    sol = T.solve(m); sig = signature(sol); lab = labels(sol)
    for _ in 1:5
        global perm_runs += 1
        p = T.solve(permuted(m, rng))
        if signature(p) != sig
            global perm_fail += 1; push!(fails, "PERM $(m.id)")
        elseif labels(p) != lab
            global label_perm_fail += 1; push!(fails, "PERM-LABEL $(m.id)")
        end
    end
    n = length(sol.indep)
    # (ii-a) split one group of size >= 2 into singletons
    for (gi, g) in enumerate(m.groups)
        length(g) < 2 && continue
        global split_tried += 1
        groups = [[gg for (k, gg) in enumerate(m.groups) if k != gi]; [[j] for j in g]]
        s2 = T.solve(regrouped(m, groups))
        length(s2.indep) == n || continue
        global split_same += 1
        if signature(s2) != sig
            global split_fail += 1; push!(fails, "SPLIT $(m.id) group $gi")
        end
        labels(s2) != lab && (global split_label_changes += 1)
    end
    # (ii-b) split every group at once
    if any(g -> length(g) >= 2, m.groups)
        s3 = T.solve(regrouped(m, [[j] for g in m.groups for j in g]))
        if length(s3.indep) == n
            global allsplit_same += 1
            signature(s3) != sig && (global allsplit_fail += 1; push!(fails, "ALLSPLIT $(m.id)"))
        end
    end
    # (ii-c) merge the two groups of a Wegscheider fold (p = f across groups)
    gos = Pop.group_of_step(m)
    for d in sol.deps
        (d.fold && d.label === :wegscheider) || continue
        gp = gos[sol.cols[d.col].j]; gf = gos[sol.cols[d.terms[1].first].j]
        gp == gf && continue
        global merge_tried += 1
        groups = [[g for (k, g) in enumerate(m.groups) if k != gp && k != gf]; [[m.groups[gp]; m.groups[gf]]]]
        s4 = try
            T.solve(regrouped(m, groups))
        catch e
            push!(fails, "MERGE-ERR $(m.id): " * sprint(showerror, e)); nothing
        end
        s4 === nothing && continue
        length(s4.indep) == n || (global merge_count_fail += 1; push!(fails, "MERGE-COUNT $(m.id)"))
        signature(s4) != sig && (global merge_fail += 1; push!(fails, "MERGE $(m.id)"))
    end
end
println("plain mechanisms: $nm")
println("(i) permutations: $perm_runs runs, signature failures $perm_fail, label-only failures $label_perm_fail")
println("(ii-a) single-group splits: $split_tried tried, $split_same keep the count (thermodynamics ties the split back), failures $split_fail, label changes $split_label_changes")
println("(ii-b) split-all keeping the count: $allsplit_same, failures $allsplit_fail")
println("(ii-c) merges along Wegscheider folds: $merge_tried tried, count changes $merge_count_fail, failures $merge_fail")
foreach(println, fails[1:min(end, 40)])
