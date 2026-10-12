# ABOUTME: Prints the census tables of saved runs: pipeline counts per round, the classes of
# ABOUTME: fitted candidates with overlaps and cost weights, check costs, and examples.

"""The faces messages of `r` whose text contains `word`."""
face_hits(r, word) = [f for f in get(r, :faces, Pair{String, String}[]) if occursin(word, f[2])]

"""The well-posedness sub-classes of a classified record (empty when it passes)."""
function wclasses(r)
    w = Symbol[]
    haskey(r, :conv_err) && push!(w, :W_ctor)
    get(r, :boundary, nothing) === nothing || push!(w, :W_boundary)
    get(r, :re_turn, nothing) === nothing || push!(w, :W_re_turnover)
    isempty(face_hits(r, "bottomless")) || push!(w, :W_bottomless)
    isempty(get(r, :dead, String[])) || push!(w, :W_D0_face)
    isempty(face_hits(r, "split")) || push!(w, :W_split)
    isempty(face_hits(r, "trap")) || push!(w, :W_trap)
    get(r, :v0, false) && push!(w, :W_v0)
    get(r, :admit, nothing) === nothing || push!(w, :W_admit)
    haskey(r, :law_err) && push!(w, :W_law_err)
    w
end

"""Whether the record's only face findings are on single-reactant faces."""
single_face_only(r) = !isempty(get(r, :faces, [])) &&
    all(f -> !(f[1] in ("P=0", "S=0")), r[:faces])

pct(a, b) = b == 0 ? "-" : @sprintf("%.1f%%", 100a / b)
q(v, p) = isempty(v) ? NaN : quantile(v, p)

"""Every class label of the fitted records of `rounds`, in run order (D needs the order)."""
function label!(rounds)
    seen_keys = Set{String}()
    for recs in rounds, r in recs
        r[:fitted] || continue
        haskey(r, :key) || continue
        k = r[:key]
        r[:dupkey] = k in seen_keys
        push!(seen_keys, k)
    end
    for recs in rounds, r in recs
        r[:fitted] || continue
        labels = wclasses(r)
        get(r, :dupkey, false) && push!(labels, :D)
        rk, ns = get(r, :rank, -1), get(r, :n_slice, -1)
        # P: the slice law (the rewrite's) has fewer identifiable directions than the slice
        # fits; Ptoday: the package law (what identify fits) has fewer than the package fits.
        rk >= 0 && ns > rk && push!(labels, :P)
        prk = get(r, :pkg_rank, rk >= 0 && r[:n_params] == ns && rk == ns ? rk : -1)
        r[:pkg_rank_used] = prk
        prk >= 0 && r[:n_params] > prk && push!(labels, :Ptoday)
        haskey(r, :pkg_rank) && r[:pkg_rank] != rk && push!(labels, :Xlaw)
        (get(r, :nosat_fwd, false) || get(r, :nosat_rev, false)) && push!(labels, :N)
        get(r, :same_key_as_parent, false) && push!(labels, :Z)
        r[:labels] = labels
    end
end

isW(l) = startswith(string(l), "W")

function print_tables(info, rounds; nexamples = 3)
    label!(rounds)
    # The key of each eq_hash's fitted representative, for eq_hash merges the key would split.
    rep_key = Dict{UInt64, String}()
    for recs in rounds, r in recs
        r[:fitted] && haskey(r, :key) && (rep_key[r[:eq_hash]] = r[:key])
    end
    println("="^100)
    @printf("%s: %d distinct init seeds; sample N0=%d (seed %d); full base tier: %s\n",
            info[:rxn], info[:n_seeds_all], info[:N0], info[:seed], info[:full])
    haskey(info, :full_memo) && println("full-run eq_hash memo: ", info[:full_memo])
    for (ri, recs) in enumerate(rounds)
        rnd = ri - 1
        println("-"^100)
        if rnd == 0
            r0, b = info[:r0], info[:r0_batch]
            @printf("ROUND 0: sampled seeds %d; seed screen: %d degenerate seeds removed, \
                    replaced by %d flip children (%d flip errors) -> %d base mechanisms\n",
                    r0.n_seeds, r0.n_degenerate,
                    r0.n_base - (r0.n_seeds - r0.n_degenerate), r0.n_base_fail, r0.n_base)
        else
            x, b = info[Symbol(:r, rnd)], info[Symbol(:r, rnd, :_batch)]
            @printf("ROUND %d: %d parents sampled from a pool of %d; %d parents threw in \
                    expand_mechanisms; %d children emitted, %d distinct (first parent \
                    wins); expansion %.1f s (max per parent %.2f s)\n", rnd, x.n_parents,
                    x.n_pool, x.n_fail_parents, x.n_emitted, x.n_children, x.t,
                    maximum(x.expand_times; init = 0.0))
            for (pid, msg) in x.fails
                println("    parent $(pid) threw: ", first(msg, 220))
            end
        end
        if rnd == 0 && haskey(info, :degenerate_slice) && !isempty(info[:degenerate_slice])
            dg = info[:degenerate_slice]
            @printf("  removed degenerate seeds under the slice: %d; no forward Vmax %d, no \
                    reverse Vmax %d, v = 0 %d, D = 0 on a face %d, face split/trap/bottomless \
                    %d, admit fails %d, law errors %d\n", length(dg),
                    count(d -> get(d, :nosat_fwd, false), dg),
                    count(d -> get(d, :nosat_rev, false), dg), count(d -> get(d, :v0, false), dg),
                    count(d -> !isempty(get(d, :dead, String[])), dg),
                    count(d -> !isempty(get(d, :faces, [])), dg),
                    count(d -> get(d, :admit, nothing) !== nothing, dg),
                    count(d -> haskey(d, :law_err) || haskey(d, :conv_err), dg))
            for d in dg[1:min(2, end)]
                println("    e.g. ", first(get(d, :text, ""), 400))
            end
        end
        @printf("  structural dedup: %d already seen in this run, %d produced by an \
                unsampled base mechanism of the full run (dropped)\n", b.n_seen, b.n_full)
        n = length(recs)
        cx = count(r -> r[:status] === :cx_skip, recs)
        ps = count(r -> r[:status] === :param_skip, recs)
        fl = count(r -> r[:status] === :fail, recs)
        ok = count(r -> r[:status] === :ok, recs)
        nf = count(r -> r[:fitted], recs)
        @printf("  fresh %d -> complexity skip (>%d) %d, param skip (>%d) %d, compile/render \
                errors %d -> compiled %d -> eq_hash already fit (memo) %d, eq_hash \
                duplicate in batch %d -> FITTED %d\n", n, CX_FILTER, cx, MAX_PARAMS, ps, fl,
                ok, b.n_memo, ok - b.n_memo - nf, nf)
        errs = Dict{String, Int}()
        for r in recs
            r[:status] === :fail || continue
            k = first(r[:err], 90)
            errs[k] = get(errs, k, 0) + 1
        end
        for (k, v) in sort!(collect(errs); by = last, rev = true)
            println("    error x$v: ", k)
        end
        dups = [r for r in recs if r[:status] === :ok && !r[:fitted] && haskey(r, :key) &&
                haskey(rep_key, r[:eq_hash])]
        @printf("  eq_hash duplicates (not fitted) %d, of which the key keeps %d apart from \
                their fitted representative\n", length(dups),
                count(r -> rep_key[r[:eq_hash]] != r[:key], dups))
        rnd > 0 && @printf("  fitted children whose eq_hash the full run's round 0 already \
                           holds: %d\n", count(r -> r[:fitted] && get(r, :in_full_memo, false),
                                                recs))
        fit = [r for r in recs if r[:fitted]]
        isempty(fit) && continue
        nfit = length(fit)
        sum_p = sum(r[:n_params] for r in fit)
        sum_d = sum(get(r, :nD, 0) for r in fit)
        @printf("  fitted: params median %.0f (max %d), law size (D terms) median %.0f \
                (p90 %.0f, max %d), V×τ median %.0f (max %d); sums: params %d, D terms %d\n",
                median([r[:n_params] for r in fit]), maximum(r[:n_params] for r in fit),
                median([get(r, :nD, 0) for r in fit]), q([get(r, :nD, 0) for r in fit], 0.9),
                maximum(get(r, :nD, 0) for r in fit), median([r[:cx] for r in fit]),
                maximum(r[:cx] for r in fit), sum_p, sum_d)
        nallo = count(r -> get(r, :allo, false), fit)
        @printf("  fitted allosteric %d of %d\n", nallo, nfit)
        cls = Dict{Symbol, Vector{Any}}()
        for r in fit, l in r[:labels]
            push!(get!(cls, l, Any[]), r)
        end
        anyW = [r for r in fit if any(isW, r[:labels])]
        okr = [r for r in fit if isempty(setdiff(r[:labels], [:Xlaw]))]
        row(lbl, rs) = @printf("    %-16s %5d  %6s of fitted | params %6s | D terms %6s\n",
                               lbl, length(rs), pct(length(rs), nfit),
                               pct(sum((r[:n_params] for r in rs); init = 0), sum_p),
                               pct(sum((get(r, :nD, 0) for r in rs); init = 0), sum_d))
        println("  classes among fitted (not exclusive; shares of count, of summed params, \
                of summed D terms):")
        row("W (any)", anyW)
        for l in sort!([k for k in keys(cls) if isW(k)])
            row("  " * string(l), cls[l])
        end
        for l in (:D, :P, :Ptoday, :N, :Z)
            row(string(l), get(cls, l, Any[]))
        end
        row("OK (none)", okr)
        row("Xlaw (note)", get(cls, :Xlaw, Any[]))
        # overlaps
        tags(r) = join(sort!(unique([isW(l) ? "W" : l in (:P, :Ptoday) ? "P" : string(l)
                                     for l in r[:labels] if l !== :Xlaw])), "+")
        combos = Dict{String, Int}()
        for r in fit
            t = tags(r)
            combos[isempty(t) ? "OK" : t] = get(combos, isempty(t) ? "OK" : t, 0) + 1
        end
        println("  exclusive combinations: ",
                join(("$k $v" for (k, v) in sort!(collect(combos); by = last, rev = true)), ", "))
        # face detail
        fsingle = count(r -> any(l -> l in (:W_split, :W_trap, :W_bottomless), r[:labels]) &&
                             single_face_only(r), fit)
        @printf("  face findings only on single-reactant faces: %d\n", fsingle)
        # P detail
        Pt = get(cls, :Ptoday, Any[])
        if !isempty(Pt)
            dt = [r[:n_params] - r[:pkg_rank_used] for r in Pt]
            @printf("  Ptoday: package deficiency (package count - package-law rank) total %d, \
                    median %.0f, max %d; of these also P (slice) %d\n", sum(dt), median(dt),
                    maximum(dt), count(r -> :P in r[:labels], Pt))
        end
        P = get(cls, :P, Any[])
        if !isempty(P)
            ex = Dict{Symbol, Int}()
            nP_pkgextra = count(r -> r[:n_params] > get(r, :n_slice, r[:n_params]), P)
            for r in P, (k, v) in get(r, :explain, Dict{Symbol, Int}())
                ex[k] = get(ex, k, 0) + v
            end
            defs = [r[:n_slice] - r[:rank] for r in P]
            @printf("  P: deficiency (slice count - slice-law rank) total %d, per candidate \
                    median %.0f max %d; candidates whose package count exceeds the slice count: \
                    %d (excess total %d); slice-level attribution: %s\n", sum(defs),
                    median(defs), maximum(defs), nP_pkgextra,
                    sum((max(0, r[:n_params] - get(r, :n_slice, r[:n_params])) for r in P);
                        init = 0), ex)
            subP = Dict{Symbol, Int}()
            for r in P
                for k in keys(get(r, :explain, Dict{Symbol, Int}()))
                    subP[k] = get(subP, k, 0) + 1
                end
                r[:n_params] > get(r, :n_slice, r[:n_params]) &&
                    (subP[:pkg_extra_names] = get(subP, :pkg_extra_names, 0) + 1)
            end
            println("  P candidates per sub-class: ", subP)
        end
        # never saturates detail
        nf_ = count(r -> get(r, :nosat_fwd, false), fit)
        nr_ = count(r -> get(r, :nosat_rev, false), fit)
        ns_ = count(r -> !isempty(get(r, :nosat_single, Symbol[])), fit)
        @printf("  N: no forward Vmax %d, no reverse Vmax %d; unbounded in a single reactant %d\n",
                nf_, nr_, ns_)
        # rank gaps
        kept = [r[:sv_kept] for r in fit if haskey(r, :sv_kept)]
        drop = [r[:sv_dropped] for r in fit if haskey(r, :sv_dropped) && r[:sv_dropped] > 0]
        @printf("  rank gaps: smallest kept singular ratio %.1e; largest dropped %.1e; \
                kept < 1e-6: %d; dropped > 1e-12: %d\n", minimum(kept; init = 1.0),
                maximum(drop; init = 0.0), count(<(1e-6), kept), count(>(1e-12), drop))
        # package rank validation
        pr = [r for r in fit if haskey(r, :pkg_rank)]
        if !isempty(pr)
            for why in (:needed, :sample)
                pw = [r for r in pr if r[:pkg_rank_why] === why]
                isempty(pw) && continue
                @printf("  package-law rank computed on %d fitted (%s): equals slice rank in \
                        %d; package count = slice count in %d; package law deficient in %d; \
                        slice law deficient in %d\n", length(pw), why,
                        count(r -> r[:pkg_rank] == r[:rank], pw),
                        count(r -> r[:n_params] == r[:n_slice], pw),
                        count(r -> r[:pkg_rank] < r[:pkg_rank_np], pw),
                        count(r -> r[:rank] < r[:n_slice], pw))
            end
            nshow = 0
            for r in pr
                r[:pkg_rank] == r[:rank] && continue
                (nshow += 1) > 4 && break
                @printf("    DISAGREE: pkg rank %d/%d (sv %.1e/%.1e) slice rank %d/%d: %s\n",
                        r[:pkg_rank], r[:pkg_rank_np], r[:pkg_rank_sv]..., r[:rank],
                        r[:n_slice], first(r[:text], 300))
            end
            perr = count(r -> haskey(r, :pkg_rank_err), fit)
            perr > 0 && println("  package rank errors: ", perr, " e.g. ",
                                first([r[:pkg_rank_err] for r in fit if haskey(r, :pkg_rank_err)]))
        end
        # examples
        println("  examples:")
        for l in [sort!([k for k in keys(cls) if isW(k)]); [:D, :P, :Ptoday, :N, :Z, :Xlaw]]
            rs = get(cls, l, Any[])
            isempty(rs) && continue
            for r in rs[1:min(nexamples, end)]
                extra = l in (:P, :Ptoday, :Xlaw) ? " package count $(r[:n_params]) rank $(r[:pkg_rank_used]); slice count $(get(r, :n_slice, -1)) rank $(r[:rank]) $(get(r, :explain, ""))" :
                        l === :N ? " single=$(get(r, :nosat_single, ""))" :
                        startswith(string(l), "W_split") || startswith(string(l), "W_trap") ||
                        startswith(string(l), "W_bottomless") ? " $(r[:faces][1:min(2, end)])" :
                        l === :W_D0_face ? " dead faces $(r[:dead])" :
                        l === :W_admit ? " $(first(r[:admit], 160))" :
                        l === :W_boundary ? " $(first(r[:boundary], 160))" : ""
                @printf("    [%s] id %s params %d D %d:%s\n      %s\n", l, r[:id], r[:n_params],
                        get(r, :nD, -1), extra, first(get(r, :text, ""), 900))
            end
        end
    end
end

"""Medians, p90 and maxima of the per-candidate check times over every round."""
function print_costs(rounds)
    println("-"^100)
    println("cost per candidate (ms): median / p90 / max over the candidates that ran the check")
    for k in (:t_cx, :t_dep, :t_render, :t_conv, :t_boundary, :t_admit, :t_returnover,
              :t_faces, :t_law, :t_key, :t_expand, :t_tests, :t_rank, :t_explain,
              :t_pkg_rank)
        v = [1e3 * r[k] for recs in rounds for r in recs if haskey(r, k) && r[k] > 0]
        isempty(v) && continue
        @printf("  %-14s n=%6d  %9.3f  %9.3f  %9.1f\n", k, length(v), median(v), q(v, 0.9),
                maximum(v))
    end
end
