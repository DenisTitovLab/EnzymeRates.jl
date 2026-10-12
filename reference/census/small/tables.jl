# ABOUTME: Final census tables: per reaction and round the funnel, class counts with shares and
# ABOUTME: fit-cost weights, :OnlyA rejection kinds, and the structural cause of D duplicates.
# usage: tables.jl <data names, comma-separated>
include(joinpath(@__DIR__, "lib.jl"))

pc(a, b) = b == 0 ? "-" : @sprintf("%.1f", 100a / b)

function wlist(r)
    w = String[]
    get(r, :v0, false) && push!(w, "v=0")
    r[:re_turnover] === nothing || push!(w, "RE turnover")
    for (i, f) in enumerate(r[:faces])
        r[:face_kind][i] === :none || push!(w, string(r[:face_kind][i]))
        r[:dvan][i] && push!(w, "0/0")
    end
    r[:admit] === nothing || push!(w, "onlyA")
    unique!(w)
end
admit_kind(msg) = occursin("has none", msg) ? "no free form" :
                  occursin("flips at one free form", msg) ? "second free form" :
                  occursin("still join the two sides", msg) ? "G1: arc inside a component" :
                  occursin("close a cycle", msg) ? "G1: arcs close a cycle" :
                  occursin("reach the free enzyme only", msg) ? "island" : "other"

"""Whether two allosteric mechanisms differ only in the tags of groups that have no edge in
the inactive conformation (`inactive_mask`) of either."""
function tag_only_outside_inactive(a, b)
    (S.is_allosteric(a) && S.is_allosteric(b)) || return false
    a.edges == b.edges && a.forms == b.forms && a.group == b.group &&
        a.multiplicity == b.multiplicity && a.sites == b.sites || return false
    ma, mb = S.inactive_mask(a), S.inactive_mask(b)
    ma == mb || return false
    inside = Set(a.group[k] for k in eachindex(ma) if ma[k])
    all(g -> a.tags[g] == b.tags[g], inside)
end

function tables(name)
    d = deserialize(joinpath(CEN, "data", "$name.jls"))
    mechs = deserialize(joinpath(CEN, "data", "$name.mechs.jls"))
    recs = d.recs
    byg = Dict(r[:gid] => r for r in recs)
    fitted = [r for r in recs if r[:status] === :ok && r[:fitted]]
    first_of = Dict{String, Int}()
    for r in fitted
        k = r[:key]
        haskey(first_of, k) ? (r[:D] = first_of[k]) : (first_of[k] = r[:gid])
        r[:W] = wlist(r)
        r[:N] = !r[:vmax_fwd] || !r[:vmax_rev]
        p = r[:parent]
        r[:Z] = p != 0 && byg[p][:key] == r[:key]
        r[:P] = r[:rank_pkg] < r[:np]
        r[:bad] = !isempty(r[:W]) || haskey(r, :D) || r[:P] || r[:N] || r[:Z]
    end
    println("\n## $name")
    println("round | gen | seen | fresh | cx>337 | np>20 | err | eq_hash dup | FITTED | W | D | P | N | Z | any | OK | any: params, D terms (share of round)")
    for f in d.funnel
        rd = f[:round]
        fr = [r for r in fitted if r[:round] == rd]
        n = length(fr)
        c(k) = count(r -> k === :W ? !isempty(r[:W]) : k === :D ? haskey(r, :D) : r[k], fr)
        bad = [r for r in fr if r[:bad]]
        sp = sum(r[:np] for r in fr; init = 0)
        sd = sum(r[:dterms] for r in fr; init = 0)
        bp = sum(r[:np] for r in bad; init = 0)
        bd = sum(r[:dterms] for r in bad; init = 0)
        @printf("%d | %d | %d | %d | %d | %d | %d | %d | %d | %d (%s) | %d (%s) | %d (%s) | %d | %d (%s) | %d (%s) | %d (%s) | %s%%, %s%%\n",
                rd, f[:generated], f[:seen_skip], f[:fresh], f[:cx_skip], f[:param_skip], f[:error],
                f[:eqhash_dup], n, c(:W), pc(c(:W), n), c(:D), pc(c(:D), n), c(:P), pc(c(:P), n),
                c(:N), c(:Z), pc(c(:Z), n), length(bad), pc(length(bad), n), n - length(bad),
                pc(n - length(bad), n), pc(bp, sp), pc(bd, sd))
    end
    # :OnlyA rejection kinds
    ak = Dict{String, Int}()
    for r in fitted
        r[:admit] === nothing && continue
        k = admit_kind(r[:admit])
        ak[k] = get(ak, k, 0) + 1
    end
    isempty(ak) || println("  :OnlyA rejections by kind: ", join(("$k $v" for (k, v) in sort!(collect(ak))), ", "))
    # W rejections that are also rank-deficient / duplicates
    # D: structural cause
    Ds = [r for r in fitted if haskey(r, :D)]
    if !isempty(Ds)
        tagonly = 0
        bymove = Dict{Tuple{Symbol, Symbol}, Int}()
        for r in Ds
            a, b = to_slice(mechs[r[:gid]]), to_slice(mechs[r[:D]])
            tag_only_outside_inactive(a, b) && (tagonly += 1)
            k = (byg[r[:D]][:move], r[:move])
            bymove[k] = get(bymove, k, 0) + 1
        end
        println("  D: $(length(Ds)); differ from the earlier one only in tags of groups outside the inactive conformation: $tagonly")
        println("     (earlier move, this move): ", join(("$(k[1])→$(k[2]) $v" for (k, v) in sort!(collect(bymove); by = last, rev = true)), ", "))
    end
    # P by catalytic multiplicity and allosteric/plain
    Ps = [r for r in fitted if r[:P]]
    if !isempty(Ps)
        allo1 = count(r -> r[:allo] && r[:allo_n] == 1, Ps)
        allo2 = count(r -> r[:allo] && r[:allo_n] > 1, Ps)
        println("  P: $(length(Ps)); plain $(count(r -> !r[:allo], Ps)), allosteric n=1 $allo1, n>1 $allo2; ",
                "allosteric n=1 among all fitted: $(count(r -> r[:allo] && r[:allo_n] == 1, fitted)), of which P ",
                "$(count(r -> r[:allo] && r[:allo_n] == 1 && r[:P], fitted))")
        pm = Dict{Symbol, Int}()
        for r in Ps
            pm[r[:move]] = get(pm, r[:move], 0) + 1
        end
        println("     by move: ", join(("$k $v" for (k, v) in sort!(collect(pm); by = last, rev = true)), ", "))
        sl = count(r -> get(r, :rank_slice, 0) < get(r, :np_slice, 0), fitted)
        eqnp = count(r -> r[:np] == r[:np_slice], fitted)
        rk_mismatch_eqnp = count(r -> r[:np] == r[:np_slice] && r[:rank_pkg] != r[:rank_slice], fitted)
        println("     slice law rank-deficient: $sl of $(length(fitted)); package and slice counts equal on $eqnp, ",
                "and ranks differ on $rk_mismatch_eqnp of those")
    end
    fitted
end

foreach(tables, split(ARGS[1], ","))

# The seeds the seed screen keeps from a fit: which degeneracy condition each meets.
for name in (:bi_bi_pp_rxn,)
    rxn = getfield(Main, name)
    seeds, _ = identify_seeds(rxn)
    deg = filter(ER._degenerate, seeds)
    println("\n$name: $(length(deg)) degenerate seeds of $(length(seeds))")
    for m in deg
        gs = ER.steps(m)
        println("  V fwd $(!ER._has_vmax(gs, rxn, ER.Substrate)) rev $(!ER._has_vmax(gs, rxn, ER.Product)), ",
                "C $(ER._chemistry_equilibrates_both_sides(gs, rxn)): ", steps_text(m))
    end
end
