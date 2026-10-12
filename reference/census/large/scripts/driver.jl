# ABOUTME: Runs the sampled identify rounds of one reaction (seeds, then expansions of sampled
# ABOUTME: parents), the slice classification of every compiled candidate, and saves records.

include(joinpath(@__DIR__, "census_lib.jl"))

"""Run the slice stage on every `:ok` record of `recs` in place; `key_of` maps a record id
to its key, for the parent comparison."""
function classify!(recs, rng, key_of)
    for r in recs
        r[:status] === :ok || continue
        out, sm = slice_stage(r[:mech], rng)
        merge!(r, out)
        r[:smech] = sm
        key_of[r[:id]] = get(r, :key, "")
        r[:parent] == 0 && continue
        pk = get(key_of, r[:parent], "")
        r[:same_key_as_parent] = !isempty(pk) && pk == get(r, :key, "-")
    end
end

"""The children of `parents` (records) as `_expand_parents` builds them: each parent
expanded alone, a throwing parent recorded, children deduplicated first parent wins."""
function expand_round(parents, rxn)
    children, parent_ids, fails = Any[], Tuple{Int, Int}[], Tuple{Tuple{Int, Int}, String}[]
    first_of = Dict{Any, Int}()
    n_emitted = 0
    times = Float64[]
    for p in parents
        t0 = time_ns()
        kids = try
            ER.expand_mechanisms(Union{ER.Mechanism, ER.AllostericMechanism}[p[:mech]], rxn)
        catch e
            push!(fails, (p[:id], errmsg(e)))
            push!(times, secs(t0))
            continue
        end
        push!(times, secs(t0))
        n_emitted += length(kids)
        for c in kids
            haskey(first_of, c) && continue
            first_of[c] = length(children) + 1
            push!(children, c)
            push!(parent_ids, p[:id])
        end
    end
    children, parent_ids, fails, n_emitted, times
end

function run_census(rxn_name::Symbol; N0 = 2000, N1 = 60, N2 = 0, seed = 20261011,
                    fullmemo = false, pkgrank_frac = 0.0, out = nothing)
    rxn = getfield(Main, rxn_name)
    rng = Xoshiro(seed)
    info = Dict{Symbol, Any}(:rxn => rxn_name, :N0 => N0, :N1 => N1, :N2 => N2, :seed => seed)
    T0 = time()
    t = @elapsed all_seeds = unique!(collect(ER.init_mechanisms(rxn)))
    n_all = length(all_seeds)
    info[:n_seeds_all] = n_all
    info[:t_init] = t
    @printf("%s: init_mechanisms %d distinct seeds in %.1f s; maxrss %.0f MB\n", rxn_name,
            n_all, t, Sys.maxrss() / 2^20)
    flush(stdout)
    idx0 = n_all > N0 ? sort!(randperm(Xoshiro(seed), n_all)[1:N0]) : collect(1:n_all)
    seeds = all_seeds[idx0]

    # The full run's base tier, as structure hashes (its `seen` after round 0), and with
    # `fullmemo` its eq_hash memo (the package stage on every base mechanism).
    full_seen = Set{UInt64}()
    n_deg_all = n_kids_all = n_fail_all = 0
    full_base = fullmemo ? Any[] : nothing
    t = @elapsed for m in all_seeds
        if ER._degenerate(m)
            n_deg_all += 1
            try
                for c in ER._expand_re_to_ss(m)
                    ER._assert_atom_conserving(c)
                    ER._degenerate(c) && continue
                    n_kids_all += 1
                    push!(full_seen, hash(c))
                    fullmemo && push!(full_base, c)
                end
            catch
                n_fail_all += 1
            end
        else
            push!(full_seen, hash(m))
            fullmemo && push!(full_base, m)
        end
    end
    info[:full] = (n_degenerate = n_deg_all, n_flip_children = n_kids_all,
                   n_flip_fail = n_fail_all, n_base_distinct = length(full_seen), t = t)
    @printf("full base tier: %d degenerate seeds, %d non-degenerate flip children, %d flip \
            failures, %d distinct base structures (%.1f s)\n", n_deg_all, n_kids_all,
            n_fail_all, length(full_seen), t)
    full_memo = Set{UInt64}()
    if fullmemo
        seen_f = Set{UInt64}()
        n_ok = 0
        t = @elapsed for m in full_base
            h = hash(m)
            h in seen_f && continue
            push!(seen_f, h)
            r = pkg_stage(m)
            r.status === :ok && (push!(full_memo, r.eq_hash); n_ok += 1)
        end
        info[:full_memo] = (n_ok = n_ok, n_eq_hash = length(full_memo), t = t)
        @printf("full base tier package stage: %d compiled under the caps, %d eq_hashes \
                (%.1f s)\n", n_ok, length(full_memo), t)
        full_base = nothing
    end
    all_seeds = nothing
    GC.gc()
    flush(stdout)

    # Round 0: the sampled seeds through `_base_tier`, then `_process_batch`.
    seen, memo = Set{UInt64}(), Set{UInt64}()
    key_of = Dict{Tuple{Int, Int}, String}()
    n_deg = count(ER._degenerate, seeds)
    t = @elapsed base, base_fail = ER._base_tier(seeds, rxn)
    info[:r0] = (n_seeds = length(seeds), n_degenerate = n_deg,
                 n_base = length(base), n_base_fail = length(base_fail), t_base = t)
    # The slice verdicts of the degenerate seeds the screen removes (never fitted today).
    degen = Dict{Symbol, Any}[]
    for m in seeds
        ER._degenerate(m) || continue
        dv, _ = slice_stage(m, rng)
        push!(degen, dv)
    end
    info[:degenerate_slice] = degen
    t = @elapsed recs0, c0 = process_batch!(base, seen, memo, nothing; rnd = 0)
    foreach(((i, r),) -> r[:id] = (0, i), enumerate(recs0))
    info[:r0_batch] = c0
    info[:t_pkg0] = t
    t = @elapsed classify!(recs0, rng, key_of)
    info[:t_slice0] = t
    @printf("round 0: %d seeds (%d degenerate) -> %d base; %d fitted; package %.0f s, \
            slice %.0f s; maxrss %.0f MB\n", length(seeds), n_deg, length(base),
            count(r -> r[:fitted], recs0), info[:t_pkg0], t, Sys.maxrss() / 2^20)
    flush(stdout)
    rounds = [recs0]
    for (rnd, NP, pseed) in ((1, N1, seed + 1), (2, N2, seed + 2))
        NP > 0 || break
        prev = rounds[end]
        pool = [r for r in prev if r[:status] === :ok && !get(r, :v0, false)]
        pick = length(pool) > NP ? sort!(randperm(Xoshiro(pseed), length(pool))[1:NP]) :
               collect(eachindex(pool))
        parents = pool[pick]
        t = @elapsed children, pids, fails, n_emit, etimes = expand_round(parents, rxn)
        info[Symbol(:r, rnd)] = (n_parents = length(parents), n_pool = length(pool),
                                 n_fail_parents = length(fails), fails = fails,
                                 n_emitted = n_emit, n_children = length(children), t = t,
                                 expand_times = etimes)
        t = @elapsed recs, cb = process_batch!(children, seen, memo, full_seen;
                                               parent = pids, rnd = rnd)
        foreach(((i, r),) -> r[:id] = (rnd, i), enumerate(recs))
        for r in recs
            r[:in_full_memo] = fullmemo && r[:status] === :ok && r[:eq_hash] in full_memo
        end
        info[Symbol(:r, rnd, :_batch)] = cb
        info[Symbol(:t_pkg, rnd)] = t
        t = @elapsed classify!(recs, rng, key_of)
        info[Symbol(:t_slice, rnd)] = t
        @printf("round %d: %d parents (%d threw) -> %d emitted, %d distinct children; \
                %d fresh; %d fitted; expand %.1f s, package %.0f s, slice %.0f s; \
                maxrss %.0f MB\n", rnd, length(parents), length(fails), n_emit,
                length(children), length(recs), count(r -> r[:fitted], recs),
                info[Symbol(:r, rnd)].t, info[Symbol(:t_pkg, rnd)], t, Sys.maxrss() / 2^20)
        flush(stdout)
        push!(rounds, recs)
    end

    # The package law's own rank: on every fitted candidate whose slice law is deficient or
    # whose package count differs from the slice count, and on a random `pkgrank_frac` of
    # the others (a check that their package law has the slice law's rank).
    if pkgrank_frac > 0
        prng = Xoshiro(seed + 7)
        for recs in rounds
            fit = [r for r in recs if r[:fitted] && get(r, :rank, -1) >= 0]
            for r in fit
                need = r[:rank] < r[:n_slice] || r[:n_params] != r[:n_slice]
                r[:pkg_rank_why] = need ? :needed : :sample
                need || rand(prng) < pkgrank_frac || continue
                t0 = time_ns()
                try
                    rk, np, kept, dropped = pkg_rank(r[:mech], prng)
                    r[:pkg_rank] = rk
                    r[:pkg_rank_np] = np
                    r[:pkg_rank_sv] = (kept, dropped)
                catch e
                    r[:pkg_rank_err] = errmsg(e)
                end
                r[:t_pkg_rank] = secs(t0)
            end
        end
    end
    info[:t_total] = time() - T0
    @printf("total %.0f s; maxrss %.0f MB\n", info[:t_total], Sys.maxrss() / 2^20)
    if out !== nothing
        serialize(out, (info = info, rounds = rounds))
        println("saved ", out)
    end
    info, rounds
end
