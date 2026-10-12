# ABOUTME: Totals of the P5b validation over every allosteric population, for the probe report.
# ABOUTME: Reads data/results_<pop>.jls and re-derives each mechanism for solve sizes; stdlib only.
using Serialization, Printf, Statistics
include(joinpath(@__DIR__, "Allo.jl"))

pops = ["specs", "handwritten", "allo_uniuni", "allo_bibi", "expand_bibi"]
R = Dict(p => deserialize(joinpath(@__DIR__, "data", "results_$(p).jls")) for p in pops)
M = Dict(p => Dict(m.id => m for m in Pop.load_population(p) if m.allosteric) for p in pops)
all_r = [r for p in pops for r in R[p]]
tally(v) = join(("$k=$(count(==(k), v))" for k in sort!(unique(v); by = string)), ", ")
println("mechanisms: ", length(all_r), " (", join(("$p $(length(R[p]))" for p in pops), ", "), ")")

# opposite folds split: a shared symbol vs twin (A/I copies)
for p in ("allo_bibi",)
    sh, tw = 0, 0
    for r in R[p]
        m = M[p][r[:id]]
        rA, rI = Dict(m.pass2_rename_A), Dict(m.pass2_rename_I)
        twin = Dict(m.step_names_A[j][k] => m.step_names_I_pkg[j][k] for j in eachindex(m.steps)
                    if m.inactive_kept[j] for k in 1:Allo.nslots(m.steps[j]))
        s = any(get(rI, y, nothing) == x for (x, y) in rA)
        t = any(haskey(twin, x) && haskey(twin, y) && twin[x] != x &&
                get(rI, twin[y], nothing) == twin[x] for (x, y) in rA)
        sh += s; tw += t && !s
    end
    println("$p opposite folds: shared symbol $sh, twin names only $tw")
end

# count differences against their causes, with rank deficits
causes = Dict{Any, Int}()
for r in all_r
    r[:count] == r[:pkg_count] && continue
    k = (pkg_minus_mine = r[:pkg_count] - r[:count],
         cause = !isempty(r[:undeclared_pkg]) && r[:pkg_one] ? :both :
                 !isempty(r[:undeclared_pkg]) ? :phantom_name : r[:pkg_one] ? :fold : :none,
         deficit_mine = r[:count] - r[:rank_mine][1], deficit_pkg = r[:pkg_count] - r[:rank_pkg][1])
    causes[k] = get(causes, k, 0) + 1
end
println("count differences (", sum(values(causes)), "):")
foreach(((k, n),) -> println("   $k => $n"), sort!(collect(causes); by = x -> string(x.first)))
same = filter(r -> r[:count] == r[:pkg_count], all_r)
println("equal counts: ", length(same), "; same model (package point satisfies mine) ",
        count(r -> r[:cons_err] < 1e-10, same), "; equal ranks ",
        count(r -> r[:rank_mine][1] == r[:rank_pkg][1], same), "; deficits ",
        tally([r[:count] - r[:rank_mine][1] for r in same]))
@printf("max errors: my v vs package sample %.1e; vs package law BigFloat %.1e; restriction v %.1e; n=1 oracle %.1e; detailed balance mine %.1e\n",
        maximum(r[:v_vs_sample] for r in same if r[:cons_err] < 1e-10),
        maximum(r[:v_vs_pkgbig] for r in same if r[:cons_err] < 1e-10),
        maximum(r[:restrict_v] for r in all_r), maximum(r[:oracle_n1] for r in all_r),
        maximum(r[:db_mine] for r in all_r))
println("phantom names: package ", count(r -> !isempty(r[:undeclared_pkg]), all_r), " mechanisms, ",
        sum(r -> length(r[:undeclared_pkg]), all_r), " names; mine ",
        count(r -> !isempty(r[:undeclared_mine]), all_r))
println("fold hazards: package opposite ", count(r -> r[:pkg_opp], all_r), ", one-state-only ",
        count(r -> r[:pkg_one] && !r[:pkg_opp], all_r), "; mine ", count(r -> r[:dep_in_law] > 0 || r[:opp_mine] > 0, all_r))
println("printed artifacts: package ", count(r -> r[:artifacts_pkg] > 0, all_r), " mechanisms (",
        sum(r -> r[:artifacts_pkg], all_r), "), mine ", count(r -> r[:artifacts_mine] > 0, all_r),
        "; tag-symmetric pairs mine ", sum(r -> r[:sym_same_text], all_r), "/", sum(r -> r[:sym_pairs], all_r),
        ", package same order ", sum(r -> r[:pkg_sym_same_order], all_r), "/", sum(r -> r[:pkg_sym_pairs], all_r))
println("normalization case: ", tally([r[:norm_pkg] for r in all_r]), "; mine differs ",
        count(r -> r[:norm_pkg] != r[:norm_mine], all_r))

# re-derive for the solve sizes, the RE-cycle substitutions and the weight kinds
rows, cols, maxint, nsub, dk = Int[], Int[], Int[], Ref(0), Dict{Any, Int}()
pp = Dict{Any, Int}()
for p in pops, m in values(M[p])
    md = Allo.derive(m)
    push!(rows, md.solve_stats.rows); push!(cols, length(md.cols)); push!(maxint, md.solve_stats.maxint)
    lawsyms = Set(md.T.names[v] for q in (md.A..., md.I...) for mo in keys(q.t) for (v, _) in mo)
    sub = 0
    for (st, keep) in ((:A, trues(length(m.steps))), (:I, md.keepI))
        re = [j for j in findall(keep) if m.steps[j].re]
        on = Set(re[i] for z in Allo.fundamental_cycles([(m.steps[j].from, m.steps[j].to) for j in re])
                 for i in findall(!iszero, z))
        sub += count(j -> (n = Allo.colname(m, (st, j, 1)); md.names.value[n] != (0, [n => 1]) &&
                           !haskey(md.names.folds, n)), on)
    end
    nsub[] += sub > 0
    conc = Set(md.T.index[x] for x in m.metabolites if haskey(md.T.index, x))
    kind(d) = !any(v in conc for mo in keys(d.t) for (v, _) in mo) ?
              (Allo.nterms(d) == 1 ? :mono : :poly) :
              any(mo -> !any(v in conc for (v, _) in mo), keys(d.t)) ? :met_pos : :met_zero
    k = (Allo.norm_case(md), kind(md.A[3]), kind(md.I[3]))
    dk[k] = get(dk, k, 0) + 1
    k[2] === :met_zero && (pp[:pingpong in m.family] = get(pp, :pingpong in m.family, 0) + 1)
end
println("joint solve: rows ≤ ", maximum(rows), ", columns ≤ ", maximum(cols), ", largest RREF entry ", maximum(maxint),
        "; mechanisms with a dependent RE constant on an RE cycle (written through independents) ", nsub[])
println("weights (case, d_A, d_I): ", join(("$k => $v" for (k, v) in sort!(collect(dk); by = string ∘ first)), "; "))
println("d_A without a concentration-free term, by ping-pong family: ", pp)
ts = [r[:derive_s] for r in all_r]
@printf("derive time per mechanism: median %.2f ms, p99 %.2f ms (first calls compile)\n",
        1e3 * median(ts), 1e3 * quantile(ts, 0.99))
