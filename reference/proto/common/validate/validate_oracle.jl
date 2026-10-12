# ABOUTME: Validates Oracle.jl against thermodynamics and the package's stored outputs.
# ABOUTME: Checks detailed balance, Full-mode agreement and Reduced agreement on plain mechanisms.
include(joinpath(@__DIR__, "..", "Population.jl"))
include(joinpath(@__DIR__, "..", "Oracle.jl"))
using Random, Printf, Statistics

setprecision(BigFloat, Oracle.PREC)

"""BigFloat value of a dependent's monomial at `vals`."""
evalmono(mono, vals) = prod((vals[s]^BigFloat(e) for (s, e) in mono); init = big(1))

relerr(a, b) = abs(a - b) / max(abs(b), floatmin(Float64))
posneg(p) = (Pop.PPoly([t for t in p if t.second > 0]),
             Pop.PPoly([(t.first => -t.second) for t in p if t.second < 0]))

# The mechanisms: every plain spec and handwritten mechanism, and 200 plain enumerated ones.
pops = Dict(n => Pop.load_population(n) for n in Pop.populations())
plain(n) = haskey(pops, n) ? [m for m in pops[n] if !m.allosteric && m.ok] : Pop.PMech[]
pick(xs, n, seed) = length(xs) <= n ? xs :
    xs[sort!(randperm(Random.Xoshiro(seed), length(xs))[1:n])]
enum_small = [plain("init_uniuni"); plain("init_unibi")]
enum = [enum_small; pick(plain("init_bibi"), 60, 1); pick(plain("init_bibi_pp"), 60, 2)]
enum = [enum; pick(plain("expand_bibi"), 200 - length(enum), 3)]
cases = [("spec", plain("specs")); ("enumerated", enum); ("handwritten", plain("handwritten"))]
groups_of_cases = [(lab, m) for (lab, ms) in cases for m in ms]
println("cases: ", join(("$(lab) $(length(ms))" for (lab, ms) in cases), ", "))

rng = Random.Xoshiro(424242)
rows = NamedTuple[]
for (lab, m) in groups_of_cases
    t0 = time()
    # (i) detailed balance at Γ = Keq with bare consistent constants.
    db = BigFloat[]
    for _ in 1:2
        Keq = Oracle.logu(rng)
        k, _, _ = Oracle.consistent_constants(m, rng; Keq)
        c = Oracle.concs_at_ratio(m, rng, Keq)
        veq = Oracle.steady_state_v(m, k, c)
        c[m.products[1]] /= 10
        voff = Oracle.steady_state_v(m, k, c)
        push!(db, abs(veq) / abs(voff))
    end
    # (ii) Full mode: per-group random constants, as the package's Full samples drew them.
    nre, open_re = Pop.re_cycle_closure(m)
    hasre = nre > 0
    rename = Dict(m.pass2_rename)
    full_vs_pkg = Float64[]; full_vs_raw = Float64[]; pkg_fp = Float64[]
    for s in m.samples
        isnan(s.v_full) && continue
        vals = Dict{Symbol, BigFloat}(k => big(v) for (k, v) in s.full_params)
        merge!(vals, Dict{Symbol, BigFloat}(k => big(v) for (k, v) in s.concs))
        k = Oracle.constants_from_names(m, vals)
        vo = Oracle.steady_state_v(m, k, Dict(s.concs))
        vraw = Oracle.eval_poly(m.raw_num, vals) / Oracle.eval_poly(m.raw_den, vals)
        vpkg = Oracle.eval_poly(m.pkg_num, vals) / Oracle.eval_poly(m.pkg_den, vals)
        push!(full_vs_pkg, Float64(relerr(vo, big(s.v_full))))
        push!(full_vs_raw, Float64(relerr(vo, vraw)))
        push!(pkg_fp, Float64(relerr(big(s.v_full), vpkg)))
    end
    # (iii) Reduced: the package's independent values plus its dependents (re-evaluated in
    # BigFloat from their monomials), a Pass-2-folded name taking its target's value.
    red_vs_pkg = Float64[]; red_vs_poly = Float64[]; red_db = Float64[]; red_db_poly = Float64[]
    for s in m.samples
        isnan(s.v) && continue
        vals = Dict{Symbol, BigFloat}(k => big(v) for (k, v) in s.params)
        for (k, mono) in m.dependent_monomials
            vals[k] = evalmono(mono, vals)
        end
        for (a, b) in m.pass2_rename
            haskey(vals, a) || (vals[a] = vals[b])
        end
        k = Oracle.constants_from_names(m, vals)
        vo = Oracle.steady_state_v(m, k, Dict(s.concs))
        merge!(vals, Dict{Symbol, BigFloat}(k => big(v) for (k, v) in s.concs))
        vp = Oracle.eval_poly(m.pkg_num, vals) / Oracle.eval_poly(m.pkg_den, vals)
        push!(red_vs_pkg, Float64(relerr(vo, big(s.v))))
        push!(red_vs_poly, Float64(relerr(vo, vp)))
        # (iv) the package's constants at Γ = Keq: detailed balance of the oracle network and
        # of the package's own polynomials.
        c = Dict{Symbol, BigFloat}(k => big(v) for (k, v) in s.concs)
        Γ = prod(c[x] for x in m.products) / prod(c[x] for x in m.substrates)
        c[m.products[1]] *= vals[:Keq] / Γ
        veq = Oracle.steady_state_v(m, k, c)
        peq = Oracle.eval_poly(m.pkg_num, merge(vals, c))
        c[m.products[1]] /= 10
        voff = Oracle.steady_state_v(m, k, c)
        poff = Oracle.eval_poly(m.pkg_num, merge(vals, c))
        push!(red_db, Float64(abs(veq) / abs(voff)))
        push!(red_db_poly, Float64(abs(peq) / abs(poff)))
    end
    push!(rows, (; lab, id = m.id, source = m.source, nforms = length(m.species),
                 nsteps = length(m.steps), hasre, pass2 = !isempty(m.pass2_rename),
                 pass2_changes = m.pass2_changes_polys, open_re,
                 db = Float64(maximum(db)), full_vs_pkg, full_vs_raw, pkg_fp, red_vs_pkg,
                 red_vs_poly, red_db, red_db_poly, secs = time() - t0))
end

mx(v) = isempty(v) ? 0.0 : maximum(v)
println("\n(i) detailed balance |v(Γ=Keq)| / |v(Γ=Keq/10)|, bare consistent constants:")
for (lab, _) in cases
    rs = [r for r in rows if r.lab == lab]
    isempty(rs) && continue
    @printf("  %-11s n=%3d  max %.2e  median %.2e  fail(>1e-30) %d\n", lab, length(rs),
            maximum(r.db for r in rs), median([r.db for r in rs]), count(r -> r.db > 1e-30, rs))
end
println("\n(ii) Full mode at per-group random constants (4 samples per mechanism):")
for (lab, _) in cases
    rs = [r for r in rows if r.lab == lab]
    isempty(rs) && continue
    clean = [r for r in rs if !r.hasre && !r.pass2]
    @printf("  %-11s n=%3d  clean(no RE cycle, no Pass-2)=%3d  oracle vs exact raw polys max %.2e; \
            oracle vs package Float64 v_full: max %.2e, fail(>1e-12) %d\n", lab, length(rs),
            length(clean), mx(mx.(r.full_vs_raw for r in clean)),
            mx(mx.(r.full_vs_pkg for r in clean)), count(r -> mx(r.full_vs_pkg) > 1e-12, clean))
    for r in rs
        (r.hasre || r.pass2) || continue
        @printf("     excluded %-18s RE-cycle=%d (open %d) pass2=%d changes_polys=%d  oracle vs \
                v_full %.2e  oracle vs raw polys %.2e  %s\n", r.id, r.hasre, r.open_re, r.pass2,
                r.pass2_changes,
                mx(r.full_vs_pkg), mx(r.full_vs_raw), first(r.source, 60))
    end
    for r in clean
        mx(r.full_vs_pkg) > 1e-12 || continue
        @printf("     FAIL %-18s oracle vs v_full %.2e, package Float64 vs its exact polys %.2e, \
                oracle vs exact raw polys %.2e  %s\n", r.id, mx(r.full_vs_pkg), mx(r.pkg_fp),
                mx(r.full_vs_raw), first(r.source, 60))
    end
end
println("\n(iii) Reduced mode at the package's own dependent values (4 samples per mechanism):")
for (lab, _) in cases
    rs = [r for r in rows if r.lab == lab]
    isempty(rs) && continue
    @printf("  %-11s n=%3d  oracle vs exact package polys max %.2e; oracle vs package Float64 v: \
            max %.2e, fail(>1e-12) %d\n", lab, length(rs), mx(mx.(r.red_vs_poly for r in rs)),
            mx(mx.(r.red_vs_pkg for r in rs)), count(r -> mx(r.red_vs_pkg) > 1e-12, rs))
    for r in rs
        mx(r.red_vs_pkg) > 1e-12 || mx(r.red_vs_poly) > 1e-20 || continue
        @printf("     %-18s oracle vs v %.2e  oracle vs exact polys %.2e  %s\n", r.id,
                mx(r.red_vs_pkg), mx(r.red_vs_poly), first(r.source, 60))
    end
end
println("\n(iv) detailed balance at the package's own Reduced constants, Γ = Keq " *
        "(|v(Γ=Keq)| / |v(Γ=Keq/10)|, 4 samples per mechanism):")
for (lab, _) in cases
    rs = [r for r in rows if r.lab == lab]
    isempty(rs) && continue
    @printf("  %-11s n=%3d  oracle network max %.2e, fail(>1e-30) %d; package numerator \
            polynomial max %.2e, fail(>1e-30) %d\n", lab, length(rs), mx(mx.(r.red_db for r in rs)),
            count(r -> mx(r.red_db) > 1e-30, rs), mx(mx.(r.red_db_poly for r in rs)),
            count(r -> mx(r.red_db_poly) > 1e-30, rs))
end
@printf("\ntime: total %.1f s, max per mechanism %.2f s (%s, %d forms)\n", sum(r.secs for r in rows),
        maximum(r.secs for r in rows), rows[argmax([r.secs for r in rows])].id,
        rows[argmax([r.secs for r in rows])].nforms)
