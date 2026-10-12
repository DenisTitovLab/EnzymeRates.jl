# ABOUTME: Recomputes the census class table from the saved records with independent code,
# ABOUTME: checks internal consistency (Z ⊂ D, np within D pairs, flags) and draws verification samples.
using Serialization, Printf, Random, Statistics

const CEN = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/census/small"
const VER = joinpath(CEN, "verify")
const RUNS = ["uni_uni_rxn", "uni_uni_with_reg", "uni_uni_allo", "uni_uni_allo_reg",
              "uni_uni_allo_2reg", "uni_uni_reg_and_inhibitor", "uni_bi_rxn", "bi_bi_rxn",
              "bi_bi_pp_rxn", "bi_bi_rxn.r2_100000", "bi_bi_pp_rxn.r2_1200"]

pc(a, b) = b == 0 ? "-" : @sprintf("%.1f%%", 100a / b)

"""The W reasons of a record, every flag the census stored (not only those tables.jl reads)."""
function wreasons(r)
    w = String[]
    haskey(r, :slice_err) && push!(w, "slice_err")
    get(r, :v0, false) && push!(w, "v0")
    get(r, :D0, false) && push!(w, "D0")
    get(r, :re_turnover, nothing) === nothing || push!(w, "re_turnover")
    for (i, k) in enumerate(get(r, :face_kind, Symbol[]))
        k === :none || push!(w, "face_" * string(k))
    end
    any(get(r, :dvan, Bool[])) && push!(w, "dvan")
    get(r, :admit, nothing) === nothing || push!(w, "admit")
    get(r, :role, nothing) === nothing || push!(w, "role")
    get(r, :names_err, nothing) === nothing || push!(w, "names")
    get(r, :unused, nothing) === nothing || push!(w, "unused")
    unique!(w)
end
admit_kind(msg) = occursin("has none", msg) ? "nofree" :
                  occursin("flips at one free form", msg) ? "secondfree" :
                  occursin("still join the two sides", msg) ? "G1_inside" :
                  occursin("close a cycle", msg) ? "G1_cycle" :
                  occursin("reach the free enzyme only", msg) ? "island" : "other"

samples = Dict{String, Any}()

for name in RUNS
    d = deserialize(joinpath(CEN, "data", "$name.jls"))
    recs = d.recs
    byg = Dict(r[:gid] => r for r in recs)
    @assert all(i -> recs[i][:gid] == i, eachindex(recs))
    println("\n## $name  (records $(length(recs)))")
    # funnel recomputed from records
    for f in d.funnel
        rd = f[:round]
        rr = [r for r in recs if r[:round] == rd]
        st(s) = count(r -> r[:status] === s, rr)
        nfit = count(r -> r[:status] === :ok && r[:fitted], rr)
        ninh = count(r -> r[:status] === :ok && !r[:fitted], rr)
        @printf("  round %d funnel: generated %d seen %d fresh %d (recs %d) cx %d (recs %d) np>20 %d err %d ok %d eqdup %d (recs %d) fitted %d (recs %d)\n",
                rd, f[:generated], f[:seen_skip], f[:fresh], length(rr), f[:cx_skip], st(:cx_skip),
                f[:param_skip], f[:error], f[:ok], f[:eqhash_dup], ninh, f[:fitted], nfit)
    end
    fitted = [r for r in recs if r[:status] === :ok && r[:fitted]]
    # eq_hash: every fitted eqh distinct; every inherited points to a fitted rep with same eqh
    eqs = [r[:eqh] for r in fitted]
    @printf("  fitted eq_hash distinct: %s; inherited reps consistent: %s\n", allunique(eqs),
            all(r -> byg[r[:rep]][:fitted] && byg[r[:rep]][:eqh] == r[:eqh],
                [r for r in recs if r[:status] === :ok && !r[:fitted]]))
    # classes, my own composition
    first_key = Dict{String, Int}()
    D = Dict{Int, Int}()
    for r in fitted           # gid order = processing order
        haskey(r, :key) || continue
        k = r[:key]
        haskey(first_key, k) ? (D[r[:gid]] = first_key[k]) : (first_key[k] = r[:gid])
    end
    W = Dict(r[:gid] => wreasons(r) for r in fitted)
    P = Dict(r[:gid] => (haskey(r, :rank_pkg) ? r[:rank_pkg] < r[:np] : missing) for r in fitted)
    Pneg = count(r -> get(r, :rank_pkg, 0) < 0, fitted)
    N = Dict(r[:gid] => (!r[:vmax_fwd] || !r[:vmax_rev]) for r in fitted)
    Z = Dict(r[:gid] => (r[:parent] != 0 && haskey(byg[r[:parent]], :key) &&
                         byg[r[:parent]][:key] == r[:key]) for r in fitted)
    nunranked = count(r -> !get(r, :ranked, false), fitted)
    println("  unranked fitted: $nunranked; rank_pkg < 0: $Pneg; slice_err: $(count(r -> haskey(r, :slice_err), fitted)); pkg_err: $(count(r -> haskey(r, :pkg_err), fitted))")
    rounds = sort!(unique(r[:round] for r in fitted))
    for rd in rounds
        fr = [r for r in fitted if r[:round] == rd]
        n = length(fr)
        g = [r[:gid] for r in fr]
        nW = count(x -> !isempty(W[x]), g)
        nD = count(x -> haskey(D, x), g)
        nP = count(x -> P[x] === true, g)
        nN = count(x -> N[x], g)
        nZ = count(x -> Z[x], g)
        any_ = [x for x in g if !isempty(W[x]) || haskey(D, x) || P[x] === true || N[x] || Z[x]]
        sp = sum(r[:np] for r in fr); sd = sum(r[:dterms] for r in fr)
        bp = sum(byg[x][:np] for x in any_; init = 0); bd = sum(byg[x][:dterms] for x in any_; init = 0)
        @printf("  round %d: fitted %d | W %d (%s) | D %d (%s) | P %d (%s) | N %d | Z %d (%s) | any %d (%s) | OK %d | any params %s, dterms %s\n",
                rd, n, nW, pc(nW, n), nD, pc(nD, n), nP, pc(nP, n), nN, nZ, pc(nZ, n),
                length(any_), pc(length(any_), n), n - length(any_), pc(bp, sp), pc(bd, sd))
        wk = Dict{String, Int}()
        for x in g, w in W[x]
            k = w == "admit" ? "admit:" * admit_kind(byg[x][:admit]) : w
            wk[k] = get(wk, k, 0) + 1
        end
        isempty(wk) || println("     W reasons: ", join(("$k $v" for (k, v) in sort!(collect(wk))), ", "))
        # Z ⊂ D?
        zn = count(x -> Z[x] && !haskey(D, x), g)
        zn > 0 && println("     Z but not D: $zn")
    end
    # D pairs at the package level: parameter counts
    dp = [(x, y) for (x, y) in D]
    if !isempty(dp)
        npdiff = count(((x, y),) -> byg[x][:np] != byg[y][:np], dp)
        mism = count(((x, y),) -> byg[x][:np] != byg[x][:np_slice] || byg[y][:np] != byg[y][:np_slice], dp)
        rkdiff = count(((x, y),) -> byg[x][:rank_pkg] != byg[y][:rank_pkg], dp)
        println("  D pairs $(length(dp)): package np differs $npdiff; a member with package np != slice np $mism; package rank differs $rkdiff")
        ex = [(x, y) for (x, y) in dp if byg[x][:rank_pkg] != byg[y][:rank_pkg]]
        for (x, y) in ex[1:min(end, 3)]
            println("     #$x np $(byg[x][:np]) rank $(byg[x][:rank_pkg]) (slice $(byg[x][:np_slice])/$(get(byg[x], :rank_slice, -9)))  vs  #$y np $(byg[y][:np]) rank $(byg[y][:rank_pkg]) (slice $(byg[y][:np_slice])/$(get(byg[y], :rank_slice, -9)))")
            println("       ", byg[x][:text], "\n       ", byg[y][:text])
        end
    end
    # P: package-only deficits and slice-only deficits
    pk = [r for r in fitted if P[r[:gid]] === true]
    if !isempty(pk)
        pkgonly = count(r -> get(r, :rank_slice, r[:np_slice]) == r[:np_slice], pk)
        sliceonly = count(r -> P[r[:gid]] === false && get(r, :rank_slice, 0) < r[:np_slice], fitted)
        println("  P $(length(pk)): slice law full rank (package-only deficit) $pkgonly; slice deficit but package full rank $sliceonly")
    end
    # N: numeric probe vs condition V
    nn = count(r -> get(r, :unb_fwd, false) || get(r, :unb_rev, false), fitted)
    println("  N condition V: $(count(r -> N[r[:gid]], fitted)); numeric probe unbounded: $nn; condition C: $(count(r -> r[:condC], fitted)); kcat err: $(count(r -> r[:kcat_err] !== nothing, fitted))")
    # samples
    rng = Xoshiro(hash(name) % 10^6)
    pick(xs, k) = isempty(xs) ? Int[] : xs[sort!(randperm(rng, length(xs))[1:min(k, length(xs))])]
    s = Dict{Symbol, Any}()
    gidsf = [r[:gid] for r in fitted]
    s[:W_island] = pick([x for x in gidsf if "admit" in W[x] && admit_kind(byg[x][:admit]) == "island"], 12)
    s[:W_G1] = pick([x for x in gidsf if "admit" in W[x] && startswith(admit_kind(byg[x][:admit]), "G1")], 14)
    s[:W_face] = [x for x in gidsf if any(startswith("face_"), W[x]) || "dvan" in W[x]]
    allo1 = [x for x in gidsf if isempty(W[x]) && byg[x][:allo] && get(byg[x], :allo_n, 1) == 1]
    s[:W_unflagged_allo] = pick(allo1, 12)
    s[:W_unflagged_plain] = pick([x for x in gidsf if isempty(W[x]) && !byg[x][:allo]], 12)
    s[:D_pairs] = [(x, D[x]) for x in pick(collect(keys(D)), 15)]
    s[:D_pairs_npdiff] = [(x, y) for (x, y) in dp if byg[x][:np] != byg[y][:np] || byg[x][:rank_pkg] != byg[y][:rank_pkg]]
    s[:D_unflagged] = pick([x for x in gidsf if !haskey(D, x) && byg[x][:round] > 0], 12)
    s[:P] = pick([x for x in gidsf if P[x] === true], 15)
    s[:P_unflagged] = pick([x for x in gidsf if P[x] === false], 12)
    s[:Z] = [(x, byg[x][:parent]) for x in pick([x for x in gidsf if Z[x]], 12)]
    s[:Z_unflagged] = [(x, byg[x][:parent]) for x in pick([x for x in gidsf if !Z[x] && byg[x][:parent] != 0], 12)]
    s[:N_probe] = [r[:gid] for r in fitted if get(r, :unb_fwd, false) || get(r, :unb_rev, false)]
    s[:N_sample] = pick(gidsf, 15)
    samples[name] = s
    println("  samples: ", join(("$k $(length(v))" for (k, v) in sort!(collect(s); by = first)), ", "))
end
serialize(joinpath(VER, "samples.jls"), samples)
println("\nwrote samples.jls")
