# ABOUTME: Slice side of the well-posedness guard study: record loading, structural checks, the
# ABOUTME: law-based support test (bare and reduced), numeric checks and an RE-limit oracle.

const SP = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad"
const R = joinpath(SP, "round4/guard")
include(joinpath(SP, "proto/common/Population.jl"))
include(joinpath(SP, "slice/src/Slice.jl"))
const S = Slice
using Serialization, Random, Printf, LinearAlgebra

errmsg(e) = replace(first(sprint(showerror, e), 400), '\n' => ' ')
"""`nothing` when `f()` runs, else its error message."""
try_msg(f) = try
    f()
    nothing
catch e
    errmsg(e)
end

const CORE_SETS = ["specs", "handwritten", "init_uniuni", "init_unibi", "init_bibi",
    "init_bibi_pp", "expand_bibi", "allo_uniuni", "allo_bibi", "allo_uniuni_inh", "allo_unibi",
    "allo_bibi_pp", "allo_terbi", "p5a_adv_cases", "p5a_adv_cases2", "p5a_adv_cases3",
    "p5a_adv_cases4", "p5a_cases", "p1_adv_pkg", "p1_adv_pkg2", "p1_enum_inh", "p1_direct",
    "p2_adv_pkg", "p2_adv_pkg2", "p2_fuzz_dump", "p2_fuzz_dump_flip", "p2_hand"]
const PLAIN_ENUM = ["init_uniuni", "init_unibi", "init_bibi", "init_bibi_pp", "expand_bibi",
                    "p1_enum_inh"]
const ALLO_ENUM = ["allo_uniuni", "allo_bibi", "allo_uniuni_inh", "allo_unibi", "allo_bibi_pp",
                   "allo_terbi"]

# ─── plain-data specs for the package process ───────────────────────────────────────────

metspec(x::S.Met) = (x.name, x.modifier)
formspec(f::S.RawForm) = (f.conf, metspec.(f.bound), collect(f.added), collect(f.subtracted))
"""A mechanism's raw input as plain data (Base types only) for the package process."""
spec(subs, prods, steps, group, tags, mult, sites) =
    (substrates = collect(subs), products = collect(prods),
     steps = [(formspec(s.from), formspec(s.to), metspec.(s.consumed), metspec.(s.released),
               s.re) for s in steps],
     group = collect(Int, group), tags = collect(Symbol, tags), multiplicity = mult,
     sites = [(collect(s.ligands), collect(s.states), s.multiplicity) for s in sites])
function spec(m::S.Mechanism)
    (subs, prods, steps, group), kw = S.raw(m)
    spec(subs, prods, steps, group, kw.tags, kw.multiplicity, kw.sites)
end
spec(r::NamedTuple) = spec(r.substrates, r.products, r.steps, r.group, r.tags,
                           r.multiplicity, r.sites)

# ─── structural checks: the slice's and two the slice lacks ─────────────────────────────

"""The first rapid-equilibrium edge (in `mask`) that closes a cycle of RE edges whose net
metabolite change is not zero (turnover inside one RE segment: no finite rate), as text."""
function re_turnover(m, mask = trues(length(m.edges)))
    segment, parent, order = S.re_forest(m, mask)
    tree = Set(abs(p) for p in parent if p != 0)
    _, names, ex = S.segment_offsets(m, mask)
    row = [searchsortedfirst(names, x.name) for x in m.mets]
    for (k, e) in enumerate(m.edges)
        mask[k] && e.re && !(k in tree) || continue
        d = ex[:, e.b] .- ex[:, e.a]
        for j in eachindex(m.mets)
            d[row[j]] -= S.count_at(e.cons, j) - S.count_at(e.rel, j)
        end
        any(!iszero, d) && return S.step_text(m, e)
    end
    nothing
end

"""Declared reactants that no form or step holds as a reactant."""
function unused_reactants(m)
    used = Set(x.name for x in m.mets if !x.modifier)
    miss = [x for x in [m.substrates; m.products] if !(x in used)]
    isempty(miss) ? nothing : join(miss, ", ")
end

"""
The rapid-equilibrium segments of `mask` at the face where the names `zs` are zero, as the
fast dynamics see them: an RE edge's direction that takes up a `zs` name is dead, the other
direction is fast. A form survives when it holds its segment's fewest of every `zs` name
(`segment_offsets`), so its weight in the law stays positive on the face. The law is the
fast limit on the face when, per segment: some form survives (else bottomless: 0/0); the
surviving forms form one class over the RE edges between them; and no other closed class
of the fast graph (a trap of forms whose law weight is 0) is entered by an arc alive on the
face (a fast arc, or a steady-state direction that takes up no `zs` name). The first
violation as text, or `nothing`.
"""
function face_check(m, zs, mask = trues(length(m.edges)))
    segment, names, ex = S.segment_offsets(m, mask)
    zr = [r for (r, x) in enumerate(names) if x in zs]
    nf = length(m.forms)
    survive = [all(r -> ex[r, i] == 0, zr) for i in 1:nf]
    takes(p) = any(j -> S.count_at(p, j) > 0 && m.mets[j].name in zs, eachindex(m.mets))
    fast = [Int[] for _ in 1:nf]
    alive = Tuple{Int, Int}[]
    touched = falses(nf)
    for (k, e) in enumerate(m.edges)
        mask[k] || continue
        touched[e.a] = touched[e.b] = true
        takes(e.cons) || (push!(alive, (e.a, e.b)); e.re && push!(fast[e.a], e.b))
        takes(e.rel) || (push!(alive, (e.b, e.a)); e.re && push!(fast[e.b], e.a))
    end
    reach(i) = (seen = falses(nf); seen[i] = true; todo = [i];
                while !isempty(todo); u = pop!(todo); for v in fast[u]
                    seen[v] || (seen[v] = true; push!(todo, v)); end; end; seen)
    for s in 1:maximum(segment; init = 0)
        fs = [i for i in findall(==(s), segment) if touched[i]]
        length(fs) < 2 && continue
        sv = filter(i -> survive[i], fs)
        isempty(sv) && return "bottomless {" * join((S.form_text(m, i) for i in fs), ", ") * "}"
        r = reach(sv[1])
        all(i -> r[i], sv) || return "surviving forms " * S.form_text(m, sv[1]) * " and " *
            S.form_text(m, sv[findfirst(i -> !r[i], sv)]) * " split on the face"
        for i in fs
            survive[i] && continue
            ri = reach(i)
            any(j -> ri[j] && survive[j], fs) && continue
            trap = findall(ri)                      # closed: no fast arc to a survivor
            entered = any(((a, b),) -> b in trap && !(a in trap), alive)
            entered && return "trap {" * join((S.form_text(m, j) for j in trap), ", ") *
                               "} entered on the face"
        end
    end
    nothing
end

"""Every structural verdict of the slice on `m` (`nothing` = passes)."""
function slice_struct(m)
    bl = try_msg(() -> (fs = S.bottomless(m); fs === nothing ||
        error("bottomless {" * join((S.form_text(m, i) for i in fs), ", ") * "}")))
    role = (e = S.role_change(m); e === nothing ? nothing : "role change " * S.step_text(m, e))
    names = try_msg(() -> S.check_names(m))
    adm = try
        S.admit(m)
    catch e
        "admit error: " * errmsg(e)
    end
    fc(zs) = try
        f = face_check(m, zs)
        f === nothing && S.is_allosteric(m) && any(S.inactive_mask(m)) &&
            (f = face_check(m, zs, S.inactive_mask(m));
             f === nothing || (f = "inactive: " * f))
        f
    catch e
        "face_check error: " * errmsg(e)
    end
    (bottomless = bl, role = role, names = names, admit = adm,
     face_fwd = fc(m.products), face_rev = fc(m.substrates),
     re_cycle = (t = re_turnover(m); t === nothing ? nothing : "RE turnover cycle at " * t),
     no_turnover = try
         any(!iszero, S.cycles(m, S.orient(m))[2]) ? nothing : "no cycle turns the reaction over"
     catch e
         "cycles error: " * errmsg(e)
     end,
     unused = (u = unused_reactants(m); u === nothing ? nothing : "unused reactant " * u),
     solve = try_msg(() -> S.solve(m)))
end

# ─── the law-based support test ─────────────────────────────────────────────────────────

const Branch = Vector{Tuple{Any, Int}}      # a product of (polynomial, multiplicity)

minexp(p, v) = minimum(S.mexp(mo, v) for mo in keys(p.terms))
zero_branch(b) = any(((f, k),) -> k > 0 && iszero(f), b)
bmin(b, v) = sum((k * minexp(f, v) for (f, k) in b if k > 0); init = 0 // 1)

"""
The exact support test of `v ∝ N/D`, `N` and `D` sums of branches, each a product of
polynomials (positive coefficients for `D`'s factors; branches of one sum separated by a
variable, so no cancellation between them). After the joint shift to lowest terms in the
concentrations `concs`: `D ≡ 0` (`a`), `N ≡ 0` (`b`), and whether `D` (`dfwd`, `drev`) and
`N` (`nfwd`, `nrev`) vanish identically with every product (forward face) or every substrate
(reverse face) at zero. `mu` is the shift (the joint lowest exponents), all zero when the
law came in lowest terms.
"""
function support_test(N, D, subs, prods, concs)
    a, b = all(zero_branch, D), all(zero_branch, N)
    live = [x for x in [N; D] if !zero_branch(x)]
    mu = Dict(v => isempty(live) ? 0 // 1 : minimum(bmin(x, v) for x in live) for v in concs)
    attains(f, Z) = (lo = Dict(z => minexp(f, z) for z in Z);
                     any(mo -> all(z -> S.mexp(mo, z) == lo[z], Z), keys(f.terms)))
    survives(Bs, Z) = any(Bs) do x
        zero_branch(x) && return false
        all(z -> bmin(x, z) == mu[z], Z) && all(((f, k),) -> k == 0 || attains(f, Z), x)
    end
    (a = a, b = b, dfwd = !survives(D, prods), drev = !survives(D, subs),
     nfwd = !survives(N, prods), nrev = !survives(N, subs),
     shifted = any(!iszero, values(mu)))
end

conc_indices(m, vars) = ([vars.index[x] for x in m.substrates],
                         [vars.index[x] for x in m.products],
                         collect(1:length(S.metabolites(m))))

"""Branches of the bare law(s): plain `num/den`; allosteric the MWC combination of the
active (all edges) and inactive (`inactive_mask`) bare laws, flip weights their tree sums at
the flip form. Also each conformation's own branches."""
function bare_branches(m)
    bA = S.bare_law(m)
    nA, dA = Branch[[(bA.num, 1)]], Branch[[(bA.den, 1)]]
    S.is_allosteric(m) || return (N = nA, D = dA, vars = bA.vars, per = [(nA, dA)])
    imask = S.inactive_mask(m)
    T = Int
    one_ = one(S.Poly{T})
    bI = any(imask) ? S.bare_law(m, imask) : nothing
    fA, fI = bA.trees[1], bI === nothing ? one_ : bI.trees[1]
    numI, denI = bI === nothing ? (zero(S.Poly{T}), one_) : (bI.num, bI.den)
    Lp = S.term(T, S.mono(10^6 => 1))        # L: a variable no law uses
    W = [(one_ + S.term(T, S.mono(10^6 + i => 1)), s.multiplicity)
         for (i, s) in enumerate(m.sites)]     # a site's support: 1 + ligand terms
    n = m.multiplicity
    N = Branch[[(fI, n); (bA.num, 1); (bA.den, n - 1); W],
               [(Lp, 1); (fA, n); (numI, 1); (denI, n - 1); W]]
    D = Branch[[(fI, n); (bA.den, n); W], [(Lp, 1); (fA, n); (denI, n); W]]
    per = [(nA, dA)]
    bI === nothing || push!(per, (Branch[[(bI.num, 1)]], Branch[[(bI.den, 1)]]))
    (N = N, D = D, vars = bA.vars, per = per)
end

"""Branches of the reduced law (`S.rate_law`)."""
function reduced_branches(m, sol)
    L = S.rate_law(m, sol)
    if L isa S.Law
        return (N = Branch[[(L.H, 1), (L.Q, 1)]], D = Branch[[(L.D, 1)]], vars = L.vars, L = L)
    end
    A, T = L, eltype(values(L.H.terms))
    QI, DI = S.inactive_terms(A)
    Lp = S.term(T, S.mono(A.vars.index[:L] => 1))
    WA = [(s[1], s[3]) for s in A.sites]
    WI = [(s[2], s[3]) for s in A.sites]
    n = A.n
    N = Branch[[(A.H, 1); (A.fI, n); (A.active.Q, 1); (A.active.D, n - 1); WA],
               [(A.H, 1); (Lp, 1); (A.fA, n); (QI, 1); (DI, n - 1); WI]]
    D = Branch[[(A.fI, n); (A.active.D, n); WA], [(Lp, 1); (A.fA, n); (DI, n); WI]]
    (N = N, D = D, vars = A.vars, L = A)
end

"""The law-based verdicts of `m`: bare (always tried) and reduced (after the solve), with
times. `err` fields hold derivation errors."""
function law_test(m)
    out = Dict{Symbol, Any}()
    t0 = time_ns()
    bb = try
        bare_branches(m)
    catch e
        out[:bare_err] = errmsg(e)
        nothing
    end
    out[:t_bare_law] = (time_ns() - t0) / 1e9
    if bb !== nothing
        t1 = time_ns()
        s, p, c = conc_indices(m, bb.vars)
        out[:bare] = support_test(bb.N, bb.D, s, p, c)
        out[:bare_per] = [support_test(n, d, s, p, c) for (n, d) in bb.per]
        out[:t_bare_test] = (time_ns() - t1) / 1e9
    end
    t2 = time_ns()
    rb = try
        sol = S.solve(m)
        rb = reduced_branches(m, sol)
        out[:t_reduced_law] = (time_ns() - t2) / 1e9
        rb
    catch e
        out[:reduced_err] = errmsg(e)
        nothing
    end
    if rb !== nothing
        t3 = time_ns()
        s, p, c = conc_indices(m, rb.vars)
        out[:reduced] = support_test(rb.N, rb.D, s, p, c)
        out[:t_reduced_test] = (time_ns() - t3) / 1e9
        out[:law] = rb
    end
    out
end

# ─── numeric checks: the reduced law at random points, and the RE-limit oracle ──────────

logu(rng, lo = -1, hi = 1) = big(10)^(lo + (hi - lo) * big(rand(rng)))

"""Random parameter values for `m` (independents, `Keq`, `E_total = 1`) and every constant
by name (`all_parameters`)."""
function random_params(m, rng)
    names = S.fitted_params(m)
    p = (; (x => logu(rng) for x in names)..., Keq = logu(rng, -2, 2), E_total = big(1))
    p, Dict{Symbol, BigFloat}(pairs(S.all_parameters(m, p)))
end

"""Values of `vars` (concentrations from `conc`, the rest by name from `vals`)."""
varvals(vars, conc, vals) = [vars.conc[i] ? conc[x] : vals[x] for (i, x) in enumerate(vars.names)]

evalbranches(Bs, x) = sum((prod((S.evaluate(f, x)^k for (f, k) in b); init = big(1.0))
                           for b in Bs); init = big(0.0))

"""The MWC branches at multiplicity 1 with every regulator at zero (`W = 1`): the law the
two-conformation single-subunit network gives."""
function n1_branches(m, rb)
    rb.L isa S.Law && return rb.N, rb.D
    A = rb.L
    T = eltype(values(A.H.terms))
    QI, DI = S.inactive_terms(A)
    Lp = S.term(T, S.mono(A.vars.index[:L] => 1))
    (Branch[[(A.H, 1), (A.fI, 1), (A.active.Q, 1)], [(A.H, 1), (Lp, 1), (A.fA, 1), (QI, 1)]],
     Branch[[(A.fI, 1), (A.active.D, 1)], [(Lp, 1), (A.fA, 1), (DI, 1)]])
end

"""
The steady-state rate (net consumption of `m.substrates[1]` per enzyme) of the mass-action
network of `m` with every RE edge running at `Λ·speed[k]` (its ratio kept), in the active
constants; for an allosteric `m` the single-subunit two-conformation network: inactive copies
of the forms `inactive_mask` reaches, joined at the flip form by an RE step of ratio `L`.
"""
function oracle_v(m, vals, conc, Λ, speed)
    fwd = S.orient(m)
    allo = S.is_allosteric(m)
    nf = length(m.forms)
    imask = allo ? S.inactive_mask(m) : falses(length(m.edges))
    inI = falses(nf)
    if allo
        inI[1] = true
        foreach(k -> imask[k] && (inI[m.edges[k].a] = inI[m.edges[k].b] = true),
                eachindex(m.edges))
    end
    idx = zeros(Int, 2nf)
    n = 0
    for i in 1:nf; idx[i] = (n += 1); end
    for i in 1:nf; inI[i] && (idx[nf + i] = (n += 1)); end
    M = zeros(BigFloat, n, n)
    fl = Tuple{Int, Int, Int, BigFloat, BigFloat}[]
    x1 = m.substrates[1]
    cc(p) = prod((conc[x.name] for x in S.mets_of(m, p)); init = big(1.0))
    function add!(a, b, rf, rr, w)
        M[b, a] += rf; M[a, a] -= rf; M[a, b] += rr; M[b, b] -= rr
        push!(fl, (w, a, b, rf, rr))
    end
    for (state, names) in ((0, S.constant_names(m)),
                           (1, allo ? S.constant_names(m; inactive = true) : nothing))
        state == 1 && !allo && continue
        for (k, e) in enumerate(m.edges)
            state == 1 && !imask[k] && continue
            o = S.oriented(m, k, fwd)
            if o.re
                K = vals[names[k][1]]
                s = Λ * speed[k]
                kf, kr = S.is_binding(o) ? (s, s * K) : (s * K, s)
            else
                kf, kr = vals[names[k][1]], vals[names[k][2]]
            end
            w = count(x -> x.name === x1, S.mets_of(m, o.cons)) -
                count(x -> x.name === x1, S.mets_of(m, o.rel))
            add!(idx[state * nf + o.a], idx[state * nf + o.b], kf * cc(o.cons),
                 kr * cc(o.rel), w)
        end
    end
    allo && add!(idx[1], idx[nf + 1], Λ * speed[end] * vals[:L], Λ * speed[end], 0)
    M[n, :] .= 1
    rhs = zeros(BigFloat, n)
    rhs[n] = 1
    x = M \ rhs
    sum((w * (rf * x[a] - rr * x[b]) for (w, a, b, rf, rr) in fl if w != 0); init = big(0.0))
end

reldiff(x, y) = abs(x - y) / max(abs(x), abs(y), big(1e-30))

"""
Numeric checks of the reduced law at a random positive point and on the forward and reverse
faces: the law's value (NaN when `D = 0`), and the oracle at two RE speeds `Λ` and two random
speed patterns. Verdict per point: `:agree`, `:diverges` (the RE limit is infinite),
`:pattern` (the limit depends on how fast the RE steps are), `:mismatch` (the limit exists
and differs from the law), `:law_nan` with a unique oracle limit, or `:singular`.
"""
function numeric_checks(m, rb, rng)
    setprecision(BigFloat, 256) do
        p, vals = random_params(m, rng)
        vars = rb.vars
        conc = Dict{Symbol, BigFloat}(x => logu(rng, -1, 1) for x in S.metabolites(m))
        for l in S.regulators(m)
            conc[l] = big(0.0)                 # W = 1: the n = 1 view the oracle checks
        end
        merge!(vals, Dict(pairs(p)))
        N1, D1 = n1_branches(m, rb)
        nE = length(m.edges) + 1
        sp1 = [logu(rng) for _ in 1:nE]
        sp2 = [logu(rng) for _ in 1:nE]
        out = Dict{Symbol, Any}()
        for (face, zero_names) in ((:pos, Symbol[]), (:fwd, m.products), (:rev, m.substrates))
            c = copy(conc)
            foreach(x -> c[x] = big(0.0), zero_names)
            x = varvals(vars, c, vals)
            Nv, Dv = evalbranches(N1, x), evalbranches(D1, x)
            vlaw = iszero(Dv) ? big(NaN) : Nv / Dv
            # the full law's own N and D (multiplicity n) at the point, for finiteness
            Nf, Df = evalbranches(rb.N, x), evalbranches(rb.D, x)
            v = try
                o1 = oracle_v(m, vals, c, big(10)^10, sp1)
                o2 = oracle_v(m, vals, c, big(10)^20, sp1)
                o3 = oracle_v(m, vals, c, big(10)^20, sp2)
                scale = max(abs(o1), abs(o2), abs(o3), big(1e-12))
                verdict = abs(o2) > 1e4 * max(abs(o1), big(1e-12)) && abs(o2) > 1e3 ? :diverges :
                          abs(o2 - o3) > 1e-8 * scale + big(1e-14) ? :pattern :
                          isnan(vlaw) ? :law_nan :
                          abs(vlaw - o2) > 1e-8 * max(abs(vlaw), abs(o2)) + big(1e-14) ?
                          :mismatch : :agree
                (verdict = verdict, vlaw = Float64(vlaw), o = Float64.([o1, o2, o3]),
                 N = Float64(Nf), D = Float64(Df))
            catch e
                (verdict = :singular, vlaw = Float64(vlaw), o = Float64[], N = Float64(Nf),
                 D = Float64(Df), err = errmsg(e))
            end
            out[face] = v
        end
        out
    end
end

# ─── records ────────────────────────────────────────────────────────────────────────────

"""Every core record: `(key, set, mech or nothing, spec, ctor_err)`."""
function core_records()
    recs = NamedTuple[]
    for set in CORE_SETS
        conv = deserialize(joinpath(SP, "slice/core/data/$set.jls"))
        for c in conv
            push!(recs, (key = "$set/$(c.index)", set = set, prov = "core:$(c.id)",
                         mech = c.mech, spec = spec(c.mech), ctor_err = nothing))
        end
    end
    pops = deserialize(joinpath(SP, "slice/core/data/pop_p2_hand.jls"))
    done = Set(r.key for r in recs if r.set == "p2_hand")
    for (i, p) in enumerate(pops)
        "p2_hand/$i" in done && continue
        err = try_msg(() -> S.Converted(p, i))
        push!(recs, (key = "p2_hand/$i", set = "p2_hand", prov = "core:$(p.id)", mech = nothing,
                     spec = spec(S.pop_raw(p)), ctor_err = err))
    end
    recs
end

"""A random sample of the ter-ter Pop records (init part 1 and the expand set)."""
function terter_records(n1, n2, seed)
    rng = Xoshiro(seed)
    recs = NamedTuple[]
    for (file, n) in (("pop_init_terter.part1.jls", n1), ("pop_expand_terter.jls", n2))
        pops = deserialize(joinpath(SP, "slice/integrate/data", file))
        for i in sort(randperm(rng, length(pops))[1:n])
            p = pops[i]
            mech, err = try
                S.Converted(p, i).mech, nothing
            catch e
                nothing, errmsg(e)
            end
            push!(recs, (key = "terter/$file/$i", set = "terter_sample", prov = "terter:$(p.id)",
                         mech = mech, spec = spec(S.pop_raw(p)), ctor_err = err))
        end
    end
    recs
end
