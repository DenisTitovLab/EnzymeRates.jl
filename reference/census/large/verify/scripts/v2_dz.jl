# ABOUTME: Classes D and Z: every flagged pair checked by law equality under the canonical
# ABOUTME: renaming, a sample by the oracle, and a renaming-blind fingerprint over all fitted.
include(joinpath(@__DIR__, "vlib.jl"))
setprecision(BigFloat, 256)

const RX = (:ter_bi_rxn, :ter_ter_rxn, :pyruvate_carboxylase_rxn)
rel(a, b) = Float64(abs(a - b) / max(abs(a), abs(b), big(1e-60)))

"""The fitted values of `y` mapped from those of `x` (`px`) through the canonical orders of
their laws (present constants position by position, absent ones drawn)."""
function renamed(cx, cy, px, rng)
    py = Dict{Symbol, BigFloat}()
    for (a, b) in zip(cx.order, cy.order)
        py[b] = px[a]
    end
    py
end

"""Values for law evaluation and for `all_parameters`: the fitted names of `m` from `fv`,
`Keq`, `E_total = 1`."""
pvals(m, fv, keq) = (; (n => fv[n] for n in S.fitted_params(m))..., Keq = keq,
                     E_total = big(1.0))

"""Law values of `mx` at `px` and of `my` at the renamed values, at `npts` concentrations;
with `orc` also the pruned oracle of each. Returns the largest relative differences."""
function pair_check(mx, my, rng; npts = 3, orc = false)
    Lx, Ly = S.rate_law(mx), S.rate_law(my)
    cx, cy = S.canonical(Lx), S.canonical(Ly)
    cx.key == cy.key || return (same_key = false, law = NaN, orc = NaN)
    fx = Dict{Symbol, BigFloat}(n => lu10(rng, -1, 1) for n in S.fitted_params(mx))
    fy = renamed(cx, cy, fx, rng)
    haskey(fx, :L) && (fy[:L] = fx[:L])      # L is fixed by name in the key, not renamed
    for n in S.fitted_params(my)
        haskey(fy, n) || (fy[n] = lu10(rng, -1, 1))
    end
    keq = lu10(rng, -1, 1)
    vx, vy = vals_from(mx, pvals(mx, fx, keq)), vals_from(my, pvals(my, fy, keq))
    wl, wo = 0.0, 0.0
    for _ in 1:npts
        conc = random_conc(mx, rng)
        lx, ly = law_value(Lx, lawx(Lx, vx, conc)), law_value(Ly, lawx(Ly, vy, conc))
        wl = max(wl, rel(lx, ly))
        if orc
            sp = random_speed(mx, rng)
            cxx = all_columns(mx, vx, rng)
            cyy = all_columns(my, vy, rng)
            ox = first(oracle(mx, cxx[1], cxx[2], vx, conc; Λ = big(10.0)^30, speed = sp))
            oy = first(oracle(my, cyy[1], cyy[2], vy, conc; Λ = big(10.0)^30,
                              speed = random_speed(my, rng)))
            wo = max(wo, rel(ox, oy), rel(ox, lx))
        end
    end
    (same_key = true, law = wl, orc = orc ? wo : NaN)
end

"""A fingerprint blind to renaming: the law's values with every fitted constant at 1.7
(then 0.3), `Keq` 2.3 and `L` 0.6, at three fixed concentration points, rounded to 40
digits, with the count of fitted constants a term holds."""
function fingerprint(m, pts)
    L = S.rate_law(m)
    consts = [n for (i, n) in enumerate(L.vars.names) if !L.vars.conc[i]]
    vals(cv) = Dict{Symbol, BigFloat}(n => (n === :Keq ? big(2.3) : n === :L ? big(0.6) : cv)
                                      for n in consts)
    fp = [round(law_value(L, lawx(L, vals(cv), conc)); sigdigits = 40)
          for cv in (big(1.7), big(0.3)) for conc in pts]
    (np = count(n -> !(n in (:Keq, :L)), consts), fp = fp)
end

"""How the D record `r` relates to the first record `o` with its key: rounds, parents,
allostery, and for two allosteric records with the same steps whether every group whose
tag differs lies outside the inactive conformation of both (no edge in `inactive_mask`)."""
function pair_kind(o, r)
    a, b = o[:smech], r[:smech]
    parts = [o[:id][1] == r[:id][1] ? "same round" : "earlier round",
             o[:id] == r[:parent] ? "the parent" :
             (r[:parent] != 0 && o[:parent] == r[:parent]) ? "sibling" : "other",
             string(S.is_allosteric(a) ? "allo" : "plain", "/", S.is_allosteric(b) ? "allo" : "plain")]
    if S.is_allosteric(a) && S.is_allosteric(b)
        same = a.forms == b.forms && a.edges == b.edges && a.group == b.group
        if same
            ia, ib = S.inactive_mask(a), S.inactive_mask(b)
            diff = findall(g -> a.tags[g] != b.tags[g], eachindex(a.tags))
            unreach = all(g -> all(k -> !ia[k] && !ib[k], S.members(a, g)), diff)
            push!(parts, "same steps", "tags differ in $(length(diff)) group(s)",
                  unreach ? "all outside both inactive conformations" : "SOME INSIDE")
        else
            push!(parts, "different steps")
        end
    end
    join(parts, ", ")
end

rng = Xoshiro(20261011 + 2)
for name in RX
    data = load_run(name)
    rounds = data.rounds
    byid = Dict(r[:id] => r for recs in rounds for r in recs)
    first_of = Dict{String, Any}()
    for recs in rounds, r in recs
        r[:fitted] && get!(first_of, r[:key], r)
    end
    fit = [r for recs in rounds for r in recs if r[:fitted]]
    println("="^100, "\n$name: $(length(fit)) fitted")

    # (a) every D: the law of the record equals the law of the first record with its key
    #     under the canonical renaming; a sample also by the oracle.
    D = [r for r in fit if :D in r[:labels]]
    worst, nbad, nsame = 0.0, 0, 0
    orc_pick = Set(D[sort!(randperm(rng, length(D))[1:min(12, length(D))])])
    worst_o = 0.0
    t0 = time()
    kinds = Dict{String, Int}()
    for r in D
        o = first_of[r[:key]]
        kinds[pair_kind(o, r)] = get(kinds, pair_kind(o, r), 0) + 1
        res = pair_check(o[:smech], r[:smech], rng; orc = r in orc_pick)
        res.same_key || (nbad += 1; continue)
        nsame += 1
        worst = max(worst, res.law)
        r in orc_pick && (worst_o = max(worst_o, res.orc))
        if r in orc_pick
            @printf("   oracle D %s of %s: law diff %.1e, oracle diff (both oracles and law) %.1e\n",
                    r[:id], o[:id], res.law, res.orc)
        end
    end
    @printf("(a) D flagged %d: key recomputed equal %d (%d not), worst law difference under the renaming %.1e; oracle sample %d worst %.1e (%.0f s)\n",
            length(D), nsame, nbad, worst, length(orc_pick), worst_o, time() - t0)
    println("    pair kinds: ", sort!(collect(kinds); by = last, rev = true))
    pk = Dict{String, Int}()
    for r in D
        o = first_of[r[:key]]
        k = string("package counts ", o[:n_params] == r[:n_params] ? "equal" : "differ",
                   ", package-law ranks ", get(o, :pkg_rank_used, -1) == get(r, :pkg_rank_used, -1) ? "equal" : "differ",
                   ", Xlaw ", (:Xlaw in o[:labels]) || (:Xlaw in r[:labels]))
        pk[k] = get(pk, k, 0) + 1
    end
    println("    today's package laws of the D pairs: ", sort!(collect(pk); by = last, rev = true))

    # (b) Z: each child against its parent, all by law and oracle
    Z = [r for r in fit if :Z in r[:labels]]
    for r in Z
        p = byid[r[:parent]]
        res = pair_check(p[:smech], r[:smech], rng; orc = true)
        @printf("   Z %s parent %s: same key %s, law diff %.1e, oracle diff %.1e; also D: %s\n",
                r[:id], p[:id], res.same_key, res.law, res.orc, :D in r[:labels])
    end

    # (c) fingerprint over every fitted record, in run order
    pts = [random_conc(fit[1][:smech], rng) for _ in 1:3]
    t0 = time()
    seen_fp = Dict{Any, Vector{Any}}()
    seen_key = Set{String}()
    tab = Dict{Tuple{Int, Bool, Bool}, Int}()      # (round, D by key, D by fingerprint)
    splits = Any[]
    for (ri, recs) in enumerate(rounds), r in recs
        r[:fitted] || continue
        f = fingerprint(r[:smech], pts)
        k = (f.np, f.fp)
        dk = r[:key] in seen_key
        df = haskey(seen_fp, k)
        df && !dk && push!(splits, (r, seen_fp[k]))
        dk && !df && println("   KEY DUP WITHOUT FINGERPRINT MATCH: ", r[:id])
        push!(get!(seen_fp, k, Any[]), r)
        push!(seen_key, r[:key])
        tab[(ri - 1, dk, df)] = get(tab, (ri - 1, dk, df), 0) + 1
    end
    @printf("(c) fingerprint over all fitted (%.0f s): (round, D by key, D by fingerprint) => %s\n",
            time() - t0, sort!(collect(tab)))
    println("   fingerprint matches an earlier record but the key does not: ", length(splits))
    for (r, others) in splits[1:min(6, end)]
        println("     ", r[:id], " (", r[:n_slice], " slice params, rank ", r[:rank], ") vs ",
                [(o[:id], o[:n_slice], o[:rank]) for o in others[1:min(3, end)]])
        println("       ", first(r[:text], 400))
        println("       ", first(others[1][:text], 400))
    end
    serialize(joinpath(VER, "data", "v2_splits_$(name).jls"),
              [(r[:id], [o[:id] for o in os]) for (r, os) in splits])
    flush(stdout)
end
