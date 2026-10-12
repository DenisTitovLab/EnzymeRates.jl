# ABOUTME: P8 driver over one or more populations (pooled): builds today's and the redesign's final
# ABOUTME: laws, keys them four ways, verifies every merge, and hunts missed merges by fingerprints.
include(joinpath(@__DIR__, "Laws.jl"))
include(joinpath(@__DIR__, "Factor.jl"))
using Random, Serialization, Printf, Statistics

const POPS = ARGS
const TAG = join(POPS, "+")
const OUT = joinpath(@__DIR__, "data")
mkpath(OUT)

logu(rng) = setprecision(() -> big(10)^(-2 + 4 * big(rand(rng))), BigFloat, Laws.PREC)
relerr(a, b) = Float64(abs(a - b) / max(abs(a), abs(b), big(1e-300)))
q(v) = isempty(v) ? "-" : @sprintf("median %.3g max %.3g", median(v), maximum(v))

"""Random values for the fitted names of `b`, `Keq`, `L` and the concentrations of `m`."""
function draw(rng, b::Laws.Built, m)
    v = Dict{Symbol, BigFloat}(p => logu(rng) for p in b.law.params)
    v[:Keq] = logu(rng); v[:L] = logu(rng)
    c = Dict{Symbol, BigFloat}(x => logu(rng) for x in unique([m.metabolites; m.regulators]))
    v, c
end

"""`vals` extended by the value of every dependent of `b.defs`."""
all_values(b::Laws.Built, vals) = setprecision(BigFloat, Laws.PREC) do
    out = copy(vals)
    for (d, rhs) in b.defs
        out[d] = prod((vals[s]^BigFloat(e) for (s, e) in rhs); init = big(1))
    end
    out
end

"""`v / E_total` at fitted values by the route that does not use the key: the oracle (plain),
the package's own exported law (MWC today) or P5b's evaluation (MWC redesign)."""
function independent_v(variant, b, m, vals, concs)
    variant === :today && b.kind === :plain && return Laws.oracle_today(m, vals, concs)
    variant === :today && return Laws.value(b, merge(vals, concs))
    Laws.reference_value(b, m, vals, concs)
end

"""Check of a built law against an independent evaluation: today's against the package's sample
rates, the redesign's against the oracle or P5b at one random point."""
function check_law(variant, b, m, rng)
    big_(ps) = Dict{Symbol, BigFloat}(k => big(x) for (k, x) in ps)
    variant === :today &&
        return maximum(relerr(Laws.value(b, merge(big_(s.params), big_(s.concs))), big(s.v))
                       for s in m.samples)
    vals, concs = draw(rng, b, m)
    relerr(Laws.value(b, merge(vals, concs)), Laws.reference_value(b, m, vals, concs))
end

"""SHA-256 of the normalized law with its own names (a key that does not see through renaming)."""
function exact_key(law::Canon.Law)
    lines = [string(i, "|", c, "|", join(("$s^$e" for (s, e) in mo), "*"))
             for (i, p) in enumerate(Canon.normalize(law)) for (mo, c) in p]
    Canon.bytes2hex(Canon.sha256(string(law.header, "|", join(sort(law.params), ","), "\n",
                                        join(sort!(lines), "\n"))))
end

"""Whether every definition of `from`, renamed by `mp`, holds among the definitions `to`: its
log-linear form reduces to zero once `to`'s dependents are written through their definitions."""
function defs_hold(from, mp, to)
    D = Dict(to)
    for (d, rhs) in from
        lin = Dict{Symbol, Rational{Int}}()
        add!(s, c) = haskey(D, s) ? (for (t, e) in D[s]; lin[t] = get(lin, t, 0) + c * e; end) :
                                    (lin[s] = get(lin, s, 0) + c)
        add!(get(mp, d, d), 1 // 1)
        for (s, e) in rhs
            add!(get(mp, s, s), -e)
        end
        all(iszero, values(lin)) || return false
    end
    true
end

kindof(s::Symbol) = (t = string(s); startswith(t, "Kd:") || startswith(t, "K:") ||
                     startswith(t, "K_") ? :K : startswith(t, "k:") || startswith(t, "k_") ? :k :
                     :other)

const ms = reduce(vcat, [Pop.load_population(p) for p in POPS])
const rng = Xoshiro(20261009)
const VARIANTS = (:today, :redesign)
const KEYS = (:k1, :k2, :k1c, :k2c)
recs = Dict{Any, Any}[]
built = Dict{Tuple{Symbol, Int}, Laws.Built}()
laws = Dict{Tuple{Symbol, Symbol, Int}, Canon.Law}()
canon = Dict{Tuple{Symbol, Symbol, Int}, Canon.Canonical}()
t_start = time()
for (i, m) in enumerate(ms)
    r = Dict{Any, Any}(:id => m.id, :eq_hash => m.eq_hash, :allosteric => m.allosteric,
                       :nfit => length(m.fitted_params), :family => m.family)
    for variant in VARIANTS
        builder = m.allosteric ? (variant === :today ? Laws.mwc_today : Laws.mwc_redesign) :
                  (variant === :today ? Laws.plain_today : Laws.plain_redesign)
        local b
        t0 = time_ns()
        try
            b = builder(m)
        catch e
            r[(variant, :error)] = sprint(showerror, e)
            continue
        end
        r[(variant, :build_ns)] = time_ns() - t0
        r[(variant, :check)] = check_law(variant, b, m, rng)
        r[(variant, :ndefs)] = length(b.defs)
        built[(variant, i)] = b
        for kk in KEYS
            local law
            tc = @elapsed law = kk === :k1 ? b.law : kk === :k2 ? b.alllaw :
                                first(Factor.cancel(kk === :k1c ? b.law : b.alllaw))
            if kk in (:k1c, :k2c)
                _, rem, unres = Factor.cancel(kk === :k1c ? b.law : b.alllaw)
                r[(variant, kk, :removed)] = rem
                r[(variant, kk, :unresolved)] = unres
            end
            c = Canon.canonical(law)
            tk = minimum(1:2) do _
                t1 = time_ns(); Canon.canonical(law); time_ns() - t1
            end
            # Invariance under renaming: two random renamings must give the same key, and the
            # induced map must carry the law onto its renamed copy exactly.
            inv_ok = true
            for _ in 1:2
                ren = Dict(p => Symbol("q", randstring(rng, 10)) for p in law.params)
                law2 = Canon.rename(law, ren)
                c2 = Canon.canonical(law2)
                inv_ok &= c2.key == c.key && Canon.same_under(Canon.name_map(c, c2), law, law2)
            end
            present = Set(s for comp in law.comps for mo in keys(comp) for (s, _) in mo)
            merge!(r, Dict((variant, kk, :key) => c.key, (variant, kk, :exact) => exact_key(law),
                           (variant, kk, :key_ns) => tk, (variant, kk, :cancel_s) => tc,
                           (variant, kk, :leaves) => c.leaves, (variant, kk, :hits) => c.hits,
                           (variant, kk, :invariant) => inv_ok,
                           (variant, kk, :common_factor) => Canon.common_factor_vars(law, rng),
                           (variant, kk, :absent) => count(p -> !(p in present), law.params),
                           (variant, kk, :nparams) => length(law.params),
                           (variant, kk, :nterms) => sum(length, law.comps),
                           (variant, kk, :enc_bytes) => sizeof(c.encoding)))
            laws[(variant, kk, i)] = law
            canon[(variant, kk, i)] = c
        end
    end
    push!(recs, r)
end
@printf("%s: %d mechanisms built and keyed in %.1f s\n", TAG, length(ms), time() - t_start)

# ─── Merge verification ────────────────────────────────────────────────

"""Verify member `i` against representative `r0` under the name map `mp` (representative's
names => member's): the keyed laws identical under the map, the definitions both ways (all-names
keys), and equal rates by the independent route at 3 random points, the member's fitted values
read off the representative's values of every name. Returns (exact ok, defs ok, max rel diff)."""
function verify(variant, kk, r0, i, mp)
    b0, b1 = built[(variant, r0)], built[(variant, i)]
    ok_exact = Canon.same_under(mp, laws[(variant, kk, r0)], laws[(variant, kk, i)])
    inv = Dict(v => k for (k, v) in mp)
    ok_defs = kk in (:k1, :k1c) || (defs_hold(b0.defs, mp, b1.defs) && defs_hold(b1.defs, inv, b0.defs))
    ok_exact, ok_defs, numeric(variant, r0, i, inv)
end

"""Largest relative difference of the two mechanisms' rates over 3 random points, the member's
fitted value of name `p` taken from the representative's value of `inv[p]`."""
function numeric(variant, r0, i, inv)
    b0, b1 = built[(variant, r0)], built[(variant, i)]
    worst = 0.0
    for _ in 1:3
        vals, concs = draw(rng, b0, ms[r0])
        allv = all_values(b0, vals)
        vals1 = Dict{Symbol, BigFloat}(p => allv[inv[p]] for p in b1.law.params)
        vals1[:Keq] = vals[:Keq]; vals1[:L] = vals[:L]
        worst = max(worst, relerr(independent_v(variant, b0, ms[r0], vals, concs),
                                  independent_v(variant, b1, ms[i], vals1, concs)))
    end
    worst
end

"""Classes of key `kk` of `variant` among `idx` (key => sorted member indices)."""
function classes_of(variant, kk, idx)
    cl = Dict{String, Vector{Int}}()
    for i in idx
        push!(get!(cl, recs[i][(variant, kk, :key)], Int[]), i)
    end
    foreach(sort!, values(cl))
    cl
end

"""Rename-invariant fingerprint of a fitted-name law: its rate at 3 points where every fitted
name takes one common value (concentrations, `Keq` and `L` fixed by name), 30 digits each."""
function fingerprint(b::Laws.Built, m)
    setprecision(BigFloat, 256) do
        map(1:3) do k
            vals = Dict{Symbol, BigFloat}(p => big(k + 1) / 3 for p in b.law.params)
            vals[:Keq] = big(k) / 2 + 1; vals[:L] = big(2k + 1) / 5
            for (j, x) in enumerate(sort(unique([m.metabolites; m.regulators])))
                vals[x] = big(j + k) / 7
            end
            @sprintf("%.30e", Laws.value(b, vals))
        end
    end
end

"""Each fitted name's response profile (rate with every name at 1/2 except this one at 3 and 1/5),
as text; a renaming permutes these profiles."""
function profiles(b::Laws.Built, m)
    setprecision(BigFloat, 256) do
        Dict(p => join(map((big(3), big(1) / 5)) do μ
                vals = Dict{Symbol, BigFloat}(s => big(1) / 2 for s in b.law.params)
                vals[p] = μ; vals[:Keq] = big(3) / 2; vals[:L] = big(7) / 5
                for (j, x) in enumerate(sort(unique([m.metabolites; m.regulators])))
                    vals[x] = big(j + 1) / 7
                end
                @sprintf("%.30e", Laws.value(b, vals))
            end, "|") for p in b.law.params)
    end
end

"""Bijections from `a`'s names to `b`'s that match equal profiles (all of them, up to `cap`)."""
function profile_maps(pa, pb; cap = 1000)
    sort(collect(values(pa))) == sort(collect(values(pb))) || return Dict{Symbol, Symbol}[]
    groups = Dict{String, Tuple{Vector{Symbol}, Vector{Symbol}}}()
    for (p, s) in pa; push!(get!(groups, s, (Symbol[], Symbol[]))[1], p); end
    for (p, s) in pb; push!(groups[s][2], p); end
    maps = [Dict{Symbol, Symbol}()]
    for (_, (xs, ys)) in groups
        perms = permutations_of(ys)
        maps = [merge(mp, Dict(zip(xs, ys2))) for mp in maps for ys2 in perms]
        length(maps) > cap && return maps[1:cap]
    end
    maps
end
permutations_of(v) = length(v) <= 1 ? [copy(v)] :
    [[v[k]; rest] for k in eachindex(v) for rest in permutations_of(v[setdiff(eachindex(v), k)])]

function summarize(variant, allo::Bool)
    idx = [i for i in eachindex(ms) if ms[i].allosteric == allo && haskey(built, (variant, i))]
    isempty(idx) && return nothing
    kind = allo ? "MWC" : "plain"
    nerr = count(i -> ms[i].allosteric == allo && !haskey(built, (variant, i)), eachindex(ms))
    println("── $TAG / $variant / $kind: $(length(idx)) laws, $nerr build errors")
    for i in eachindex(ms)
        e = get(recs[i], (variant, :error), nothing)
        e === nothing || ms[i].allosteric != allo || println("   error $(ms[i].id): ", first(e, 300))
    end
    chk = [recs[i][(variant, :check)] for i in idx]
    println("   law check (rel err vs ", variant === :today ? "package samples" : "oracle/P5b",
            "): ", q(chk), "; above 1e-12: ", count(>(1e-12), chk), "; build µs: ",
            q([recs[i][(variant, :build_ns)] / 1e3 for i in idx]))
    println("   classes: eq_hash ", length(unique(ms[i].eq_hash for i in idx)))
    out = Dict{Symbol, Any}()
    for kk in KEYS
        g(f) = [recs[i][(variant, kk, f)] for i in idx]
        lv, hv = g(:leaves), g(:hits)
        cf = [i for i in idx if !isempty(recs[i][(variant, kk, :common_factor)])]
        cl = classes_of(variant, kk, idx)
        sizes = sort!([length(v) for v in values(cl)]; rev = true)
        line = string("   [$kk] exact-name ", length(unique(g(:exact))), ", rename-proof ",
                      length(cl), " (largest ", first(sizes, 6), "); key µs ", q(g(:key_ns) ./ 1e3),
                      "; terms ", q(g(:nterms)), "; params ", q(g(:nparams)))
        kk in (:k1c, :k2c) && (line *= string("; laws with factors removed ",
                                              count(i -> !isempty(recs[i][(variant, kk, :removed)]), idx),
                                              ", unresolved ",
                                              count(i -> !isempty(recs[i][(variant, kk, :unresolved)]), idx),
                                              ", cancel µs ", q(g(:cancel_s) .* 1e6)))
        println(line)
        println("   [$kk] invariance failures ", count(!, g(:invariant)), "; search needed ",
                count(>(1), lv), " (max leaves ", maximum(lv; init = 0), "); symmetric laws ",
                count(>(1), hv), "; absent parameters in ", count(>(0), g(:absent)),
                "; common factor left in ", length(cf),
                isempty(cf) ? "" : string(" e.g. ", ms[cf[1]].id, " ",
                                          recs[cf[1]][(variant, kk, :common_factor)]))
        nmerge, bad_e, bad_d, bad_n, worst, kswap, ident, cross = 0, 0, 0, 0, 0.0, 0, 0, 0
        ctrl, ctrl_hit = 0, 0
        for v in values(cl)
            length(v) < 2 && continue
            length(unique(split(ms[i].id, "/")[1] for i in v)) > 1 && (cross += 1)
            for i in v[2:end]
                nmerge += 1
                mp = Canon.name_map(canon[(variant, kk, v[1])], canon[(variant, kk, i)])
                ok_e, ok_d, w = verify(variant, kk, v[1], i, mp)
                bad_e += !ok_e; bad_d += !ok_d; bad_n += w > 1e-20; worst = max(worst, w)
                any(kindof(a) != kindof(b) for (a, b) in mp) && (kswap += 1)
                all(a == b for (a, b) in mp) && (ident += 1)
                # negative control: swap the images of two fitted names; the rates must differ
                ps = built[(variant, v[1])].law.params
                if length(ps) >= 2 && ctrl < 60
                    a, b = ps[randperm(rng, length(ps))[1:2]]
                    bad = copy(mp); bad[a], bad[b] = mp[b], mp[a]
                    ctrl += 1
                    ctrl_hit += numeric(variant, v[1], i, Dict(y => x for (x, y) in bad)) > 1e-6
                end
            end
        end
        println("   [$kk] merges ", nmerge, " (classes spanning populations ", cross,
                "): exact failures ", bad_e, ", definition failures ", bad_d, ", numeric failures ",
                bad_n, @sprintf(" (max rel diff %.1e)", worst), "; swapped-map controls caught ",
                ctrl_hit, "/", ctrl, "; maps: identity ", ident, ", k<->K swaps ", kswap)
        eqc = Dict{UInt64, Vector{Int}}()
        for i in idx
            push!(get!(eqc, ms[i].eq_hash, Int[]), i)
        end
        multi = [v for v in values(eqc) if length(v) > 1]
        length(multi) > 0 && println("   [$kk] eq_hash classes with ≥ 2 members ", length(multi),
            ", split by the key ", count(v -> length(unique(recs[i][(variant, kk, :key)] for i in v)) > 1, multi))
        out[kk] = cl
    end
    # Missed merges: equal rename-invariant fingerprints but different keys (fully cancelled key).
    fp = Dict(i => fingerprint(built[(variant, i)], ms[i]) for i in idx)
    buckets = Dict{Vector{String}, Vector{Int}}()
    for i in idx
        push!(get!(buckets, fp[i], Int[]), i)
    end
    cand, renamings, examples = 0, 0, String[]
    for v in values(buckets)
        keys_ = unique(recs[i][(variant, :k1c, :key)] for i in v)
        length(keys_) < 2 && continue
        reps = [first(i for i in v if recs[i][(variant, :k1c, :key)] == k) for k in keys_]
        for x in 1:length(reps), y in x+1:length(reps)
            a, b = reps[x], reps[y]
            cand += 1
            pa, pb = profiles(built[(variant, a)], ms[a]), profiles(built[(variant, b)], ms[b])
            found = any(profile_maps(pa, pb)) do mp
                numeric(variant, a, b, Dict(y_ => x_ for (x_, y_) in mp)) < 1e-20
            end
            renamings += found
            length(examples) < 4 && push!(examples, string(ms[a].id, " vs ", ms[b].id,
                found ? " (a renaming)" : " (no renaming)", " params ",
                length(built[(variant, a)].law.params), "/", length(built[(variant, b)].law.params)))
        end
    end
    println("   fingerprint buckets ", length(buckets), "; key classes inside one bucket that differ ",
            cand, " pairs; of these, same function under a renaming (missed merges) ", renamings)
    for e in examples
        println("      ", e)
    end
    # the union of the fitted-name and all-names keys (cancelled)
    parent = Dict(i => i for i in idx)
    root(i) = parent[i] == i ? i : (parent[i] = root(parent[i]))
    for kk in (:k1c, :k2c), v in values(out[kk]), i in v[2:end]
        parent[root(i)] = root(v[1])
    end
    println("   union of k1c and k2c merges: ", length(unique(root(i) for i in idx)), " classes")
    out
end

cls = Dict((v, a) => summarize(v, a) for v in VARIANTS for a in (false, true))
serialize(joinpath(OUT, "$(TAG).jls"), (recs = recs, classes = cls))
println("saved ", joinpath(OUT, "$(TAG).jls"))
