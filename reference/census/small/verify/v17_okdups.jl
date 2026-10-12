# ABOUTME: Estimates how many census-OK candidates (no W, D, P, N, Z) share their package-law family
# ABOUTME: with an earlier fit: each sampled OK candidate against its parent and earlier same-rank siblings.
include(joinpath(@__DIR__, "vlib.jl"))
rng = Xoshiro(171717)

function same_law(a, b)
    keq = lu10(rng, -1, 1)
    fit_into(a, [lu10(rng, -1, 1) for _ in a.names], b, keq, rng; starts = 4, iters = 200) < 1e-6 ||
        return false
    ab, ba = same_family(a, b, rng)
    if !(ab < 1e-8 && ba < 1e-8)       # one harder try before saying no
        keq = lu10(rng, -1, 1)
        ab = fit_into(a, [lu10(rng, -1, 1) for _ in a.names], b, keq, rng; starts = 30, iters = 600)
        ba = fit_into(b, [lu10(rng, -1, 1) for _ in b.names], a, keq, rng; starts = 30, iters = 600)
    end
    ab < 1e-8 && ba < 1e-8
end

for name in ["bi_bi_rxn", "bi_bi_pp_rxn"]
    d, mechs = loadrun(name)
    recs = d.recs
    fitted = [r for r in recs if r[:status] === :ok && r[:fitted]]
    first_key = Dict{String, Int}()
    isD = Set{Int}()
    for r in fitted
        haskey(first_key, r[:key]) ? push!(isD, r[:gid]) : (first_key[r[:key]] = r[:gid])
    end
    wflag(r) = get(r, :v0, false) || r[:re_turnover] !== nothing || any(!=(:none), r[:face_kind]) ||
               any(r[:dvan]) || r[:admit] !== nothing
    ok = [r[:gid] for r in fitted if r[:round] == 1 && !wflag(r) && !(r[:gid] in isD) &&
          r[:rank_pkg] == r[:np] && r[:vmax_fwd] && r[:vmax_rev] && recs[r[:parent]][:key] != r[:key]]
    pick = ok[sort!(randperm(rng, length(ok))[1:min(80, length(ok))])]
    bypar = Dict{Int, Vector{Int}}()
    for r in fitted
        r[:parent] == 0 && continue
        push!(get!(bypar, r[:parent], Int[]), r[:gid])
    end
    ndup = 0; nmoves = Dict{Symbol, Vector{Int}}()
    println("\n== $name: $(length(ok)) census-OK candidates in round 1; testing $(length(pick))")
    for g in pick
        r = recs[g]
        a = plaw(mechs[g])
        comps = [r[:parent]; [h for h in bypar[r[:parent]] if h < g && recs[h][:fitted] &&
                              recs[h][:rank_pkg] == r[:rank_pkg]]]
        hit = 0
        for h in comps
            recs[h][:rank_pkg] == r[:rank_pkg] || continue
            if same_law(a, plaw(mechs[h]))
                hit = h; break
            end
        end
        v = get!(nmoves, r[:move], [0, 0]); v[1] += 1
        if hit > 0
            ndup += 1; v[2] += 1
            @printf("  #%d (%s, np %d rank %d) same law as #%d (%s, np %d)\n      %s\n      %s\n", g, r[:move],
                    r[:np], r[:rank_pkg], hit, hit == r[:parent] ? "parent" : "earlier sibling",
                    recs[hit][:np], r[:text], recs[hit][:text])
            flush(stdout)
        end
    end
    println("  census-OK sampled: $(length(pick)); same package-law family as the parent or an earlier same-rank sibling: $ndup")
    println("  by move (tested, duplicates): ", join(("$k $(v[1]) $(v[2])" for (k, v) in nmoves), ", "))
end
