# ABOUTME: Aggregates saved census runs into the report tables: the funnel per round, the
# ABOUTME: classes of fitted candidates (W, D, P, N, Z, OK) with overlaps, costs and examples.
using Serialization, Printf, Statistics
const CEN = @__DIR__

pct(a, b) = b == 0 ? "-" : @sprintf("%.1f%%", 100a / b)
med(xs) = isempty(xs) ? NaN : median(xs)
p95(xs) = isempty(xs) ? NaN : quantile(xs, 0.95)

"""The well-posedness verdict of a record: the sub-classes it fails (empty: passes)."""
function wclasses(r)
    w = String[]
    haskey(r, :slice_err) && return ["no slice law: " * first(r[:slice_err], 60)]
    get(r, :v0, false) && push!(w, "v=0")
    r[:re_turnover] === nothing || push!(w, "RE turnover")
    fk = r[:face_kind]
    for (i, f) in enumerate(r[:faces])
        fk[i] === :bottomless && push!(w, "0/0 bottomless @" * f)
        r[:dvan][i] && push!(w, "0/0 D vanishes @" * f)
        fk[i] === :split && push!(w, "split @" * f)
        fk[i] === :trap && push!(w, "trap @" * f)
        fk[i] === :other && push!(w, "face other @" * f)
    end
    r[:admit] === nothing || push!(w, "no :OnlyA limit")
    r[:role] === nothing || push!(w, "role change")
    r[:names_err] === nothing || push!(w, "names")
    r[:unused] === nothing || push!(w, "unused reactant")
    unique!(w)
end
wshort(s) = replace(s, r" @.*" => "")

"""Never saturates: no maximal rate as the substrates grow or as the products grow
(condition V, structural), with the numeric probe for comparison."""
nclass(r) = haskey(r, :vmax_fwd) && (!r[:vmax_fwd] || !r[:vmax_rev])
nnum(r) = get(r, :unb_fwd, false) || get(r, :unb_rev, false)

function analyse(name)
    d = deserialize(joinpath(CEN, "data", "$name.jls"))
    recs = d.recs
    byg = Dict(r[:gid] => r for r in recs)
    println("\n", "="^100, "\n## ", d.label)
    # ── the funnel ─────────────────────────────────────────────────────────────────────
    println("\nFunnel per round (identify's order: expand → structural dedup → caps → eq_hash):")
    @printf("%-6s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s\n", "round", "parents", "raw", "gener.",
            "seen", "fresh", "cx>337", "np>20", "error", "eq_hash", "FITTED")
    for f in d.funnel
        @printf("%-6d %8s %8s %8d %8d %8d %8d %8d %8d %8d %8d\n", f[:round],
                string(get(f, :parents, "-")), string(get(f, :raw_children, "-")),
                f[:generated], f[:seen_skip], f[:fresh], f[:cx_skip], f[:param_skip],
                f[:error], f[:eqhash_dup], f[:fitted])
        f[:round] == 0 && println("       round 0 = $(f[:seed_source]): $(f[:seeds]) seeds, ",
            "$(f[:degenerate_seeds]) degenerate (seed screen; not fitted), replaced by ",
            "$(f[:flip_children_added]) non-degenerate flip children")
        f[:round] == 2 && println("       round 2 parents: $(f[:parents]) of $(f[:parent_pool]) round-1 frontier entries")
    end
    # ── classes of fitted candidates ───────────────────────────────────────────────────
    fitted = [r for r in recs if r[:status] === :ok && r[:fitted]]
    earlier = Dict{String, Int}()           # key => gid of the first fitted with it
    for r in fitted
        haskey(r, :key) || continue
        k = r[:key]
        if haskey(earlier, k)
            r[:D] = earlier[k]
        else
            earlier[k] = r[:gid]
        end
    end
    for r in fitted
        r[:W] = wclasses(r)
        r[:N] = nclass(r)
        p = byg[r[:gid]][:parent]
        r[:Z] = p != 0 && haskey(r, :key) && haskey(byg[p], :key) && byg[p][:key] == r[:key]
        r[:P] = get(r, :ranked, false) && haskey(r, :rank_pkg) && r[:rank_pkg] >= 0 &&
                r[:rank_pkg] < r[:np]
        r[:Pslice] = get(r, :ranked, false) && haskey(r, :rank_slice) &&
                     r[:rank_slice] < r[:np_slice]
    end
    cls(r) = (s = String[]; isempty(r[:W]) || push!(s, "W"); haskey(r, :D) && push!(s, "D");
              r[:P] && push!(s, "P"); r[:N] && push!(s, "N"); r[:Z] && push!(s, "Z");
              isempty(s) ? "OK" : join(s, "+"))
    for rd in 0:2
        fr = [r for r in fitted if r[:round] == rd]
        n = length(fr)
        n == 0 && continue
        ranked = count(r -> get(r, :ranked, false), fr)
        println("\nRound $rd: $n fitted (ranked $ranked)")
        nW = count(r -> !isempty(r[:W]), fr)
        nD = count(r -> haskey(r, :D), fr)
        nP = count(r -> r[:P], fr)
        nN = count(r -> r[:N], fr)
        nZ = count(r -> r[:Z], fr)
        nOK = count(r -> cls(r) == "OK", fr)
        @printf("  W %d (%s)  D %d (%s)  P %d (%s of ranked)  N %d (%s)  Z %d (%s)  OK %d (%s)\n",
                nW, pct(nW, n), nD, pct(nD, n), nP, pct(nP, ranked), nN, pct(nN, n), nZ,
                pct(nZ, n), nOK, pct(nOK, n))
        # W sub-classes
        ws = Dict{String, Int}()
        for r in fr, w in unique(wshort.(r[:W]))
            ws[w] = get(ws, w, 0) + 1
        end
        isempty(ws) || println("  W sub-classes: ", join(("$k $v" for (k, v) in sort!(collect(ws); by = last, rev = true)), "; "))
        wf = Dict{String, Int}()
        for r in fr, w in r[:W]
            occursin("@", w) || continue
            f = split(w, "@")[2]
            wf[f] = get(wf, f, 0) + 1
        end
        isempty(wf) || println("  W faces: ", join(("$k $v" for (k, v) in sort!(collect(wf))), "; "))
        # D detail
        if nD > 0
            src = Dict{Int, Int}()
            for r in fr
                haskey(r, :D) || continue
                src[byg[r[:D]][:round]] = get(src, byg[r[:D]][:round], 0) + 1
            end
            println("  D: first fitted in round ", join(("$k: $v" for (k, v) in sort!(collect(src))), ", "),
                    "; same move as that one: ", count(r -> haskey(r, :D) && byg[r[:D]][:move] == r[:move], fr))
        end
        # P detail
        if nP > 0
            parts = Dict{Symbol, Int}()
            prim = Dict{Symbol, Int}()
            pkgonly = 0
            defs = Dict{Int, Int}()
            for r in fr
                r[:P] || continue
                dpk = r[:np] - r[:rank_pkg]
                defs[dpk] = get(defs, dpk, 0) + 1
                if r[:np] - r[:rank_pkg] > get(r, :np_slice, r[:np]) - get(r, :rank_slice, r[:rank_pkg])
                    pkgonly += 1
                end
                dp = get(r, :deficit_parts, Dict{Symbol, Int}())
                for (k, v) in dp
                    parts[k] = get(parts, k, 0) + v
                end
                pk = isempty(dp) ? :package_only :
                     first(k for k in (:L_phantom_n1, :ss_noflux, :chain_flank, :nonequal_copy, :re_iso, :other) if haskey(dp, k))
                prim[pk] = get(prim, pk, 0) + 1
            end
            println("  P: deficit (np - rank) ", join(("$k: $v" for (k, v) in sort!(collect(defs))), ", "),
                    "; primary sub-class ", join(("$k $v" for (k, v) in sort!(collect(prim); by = last, rev = true)), ", "))
            println("     deficit units by cause: ", join(("$k $v" for (k, v) in sort!(collect(parts); by = last, rev = true)), ", "),
                    "; candidates whose package deficit exceeds the slice's: $pkgonly")
            nPs = count(r -> r[:Pslice], fr)
            println("     slice law deficit (rewrite): $nPs of $ranked ranked")
        end
        nNn = count(nnum, fr)
        nNf = count(r -> haskey(r, :vmax_fwd) && !r[:vmax_fwd], fr)
        println("  N: no Vmax forward $nNf, either way $nN; numeric probe unbounded $nNn; ",
                "agree $(count(r -> r[:N] == nnum(r), fr)) of $n")
        nC = count(r -> get(r, :condC, false), fr)
        nK = count(r -> get(r, :kcat_err, nothing) !== nothing, fr)
        println("  other: condition C (a reactant drops out) $nC; kcat precondition fails (fit ends as FitFailure) $nK")
        # overlaps
        combos = Dict{String, Int}()
        for r in fr
            c = cls(r)
            combos[c] = get(combos, c, 0) + 1
        end
        println("  combinations: ", join(("$k $v" for (k, v) in sort!(collect(combos); by = last, rev = true)), "; "))
        # by move
        mv = Dict{Symbol, Vector{Int}}()
        for r in fr
            v = get!(mv, r[:move], zeros(Int, 7))
            v[1] += 1
            isempty(r[:W]) || (v[2] += 1)
            haskey(r, :D) && (v[3] += 1)
            r[:P] && (v[4] += 1)
            r[:N] && (v[5] += 1)
            r[:Z] && (v[6] += 1)
            cls(r) == "OK" && (v[7] += 1)
        end
        println("  by move (fitted W D P N Z OK): ",
                join(("$k $(v[1]) $(v[2]) $(v[3]) $(v[4]) $(v[5]) $(v[6]) $(v[7])" for (k, v) in sort!(collect(mv))), "; "))
        # cost weights
        for (lab, sel) in (("all", fr), ("OK", [r for r in fr if cls(r) == "OK"]),
                           ("P", [r for r in fr if r[:P]]), ("D", [r for r in fr if haskey(r, :D)]),
                           ("W", [r for r in fr if !isempty(r[:W])]))
            isempty(sel) && continue
            np = [r[:np] for r in sel]
            dt = [r[:dterms] for r in sel if haskey(r, :dterms)]
            @printf("  size %-3s n=%-5d params median %.0f (sum %d); D terms median %.0f, p90 %.0f, max %d, sum %d\n",
                    lab, length(sel), med(np), sum(np), med(dt), isempty(dt) ? NaN : quantile(dt, 0.9),
                    isempty(dt) ? 0 : maximum(dt), isempty(dt) ? 0 : sum(dt))
        end
    end
    # reverse D: eq_hash duplicates with a different key from their representative
    inh = [r for r in recs if r[:status] === :ok && !r[:fitted]]
    rev = [r for r in inh if haskey(r, :key) && haskey(byg[r[:rep]], :key) && byg[r[:rep]][:key] != r[:key]]
    println("\nReverse: $(length(inh)) eq_hash duplicates in all; $(length(rev)) of them have another key than their representative")
    for r in rev[1:min(end, 2)]
        println("   #$(r[:gid]) vs rep #$(r[:rep]):\n     ", r[:text], "\n     ", byg[r[:rep]][:text])
    end
    # examples per class
    println("\nExamples:")
    function ex(lab, sel, k = 3)
        isempty(sel) && return
        println("  [$lab] $(length(sel)) candidates; e.g.")
        for r in sel[1:min(end, k)]
            extra = haskey(r, :D) ? " (same key as #$(r[:D]): $(byg[r[:D]][:text]))" : ""
            p = r[:parent] == 0 ? "" : " parent #$(r[:parent])"
            println("    #$(r[:gid]) r$(r[:round]) $(r[:move])$p np=$(r[:np]) rank=$(get(r, :rank_pkg, "-")) W=$(r[:W]) parts=$(get(r, :deficit_parts, "-"))\n      $(r[:text])$extra")
        end
    end
    allw = [r for r in fitted if !isempty(r[:W])]
    for w in unique(wshort(x) for r in allw for x in r[:W])
        ex("W " * w, [r for r in allw if any(x -> wshort(x) == w, r[:W])])
    end
    ex("D", [r for r in fitted if haskey(r, :D) && !r[:Z]])
    ex("Z", [r for r in fitted if r[:Z]])
    for pk in (:L_phantom_n1, :ss_noflux, :chain_flank, :nonequal_copy, :re_iso, :other)
        ex("P $pk", [r for r in fitted if r[:P] && haskey(get(r, :deficit_parts, Dict()), pk)])
    end
    ex("P package only", [r for r in fitted if r[:P] && isempty(get(r, :deficit_parts, Dict()))])
    ex("N", [r for r in fitted if r[:N]])
    ex("kcat precondition fails", [r for r in fitted if get(r, :kcat_err, nothing) !== nothing])
    # costs
    println("\nCost per candidate (seconds; median / p95 / max):")
    okr = [r for r in recs if r[:status] === :ok]
    for k in (:t_pkg, :t_convert, :t_law, :t_key, :t_support, :t_face, :t_admit, :t_vmax,
              :t_pkglaw, :t_unbounded, :t_rank_pkg, :t_rank_slice)
        xs = [r[k] for r in okr if haskey(r, k)]
        isempty(xs) && continue
        @printf("  %-13s n=%-5d %.2e / %.2e / %.2e\n", k, length(xs), med(xs), p95(xs), maximum(xs))
    end
    fitted
end

for name in split(ARGS[1], ",")
    analyse(name)
end
