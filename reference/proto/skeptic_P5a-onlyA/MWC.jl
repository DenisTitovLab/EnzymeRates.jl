# ABOUTME: Two-conformation (MWC formulation 1) oracle, joint A/I thermodynamic solve and
# ABOUTME: Newton-face limit engine for the P5a :OnlyA semantics probe; stdlib only.
const _COMMON = joinpath(@__DIR__, "..", "common")
isdefined(Main, :Pop) || include(joinpath(_COMMON, "Population.jl"))
isdefined(Main, :Poly) || include(joinpath(_COMMON, "Poly.jl"))
isdefined(Main, :Thermo) || include(joinpath(_COMMON, "Thermo.jl"))
isdefined(Main, :Engine) || include(joinpath(_COMMON, "Engine.jl"))

module MWC

import ..Pop, ..Thermo, ..Poly, ..Engine
using LinearAlgebra, Random

const PREC = 1024
const Q = Rational{Int64}
big1024(x) = setprecision(() -> BigFloat(x), BigFloat, PREC)

# ─── Tags and step kinds ─────────────────────────────────────────────

"""Allosteric tag of every flat step (`:EqualAI` for a plain mechanism)."""
step_tags(m::Pop.PMech) = (g = Pop.group_of_step(m);
    [isempty(m.group_tags) ? :EqualAI : m.group_tags[g[j]] for j in eachindex(m.steps)])

"""Whether each step belongs to a chemistry group: a kinetic group holding any step that is
not a plain binding (Thermo kind 2 or 3). Members share their constants, so a group has one
class; a plain binding grouped with a fused binding is chemistry, as the package's guard reads it."""
function is_chem(m::Pop.PMech; per_group = true)
    k = [c.kind != 1 for c in Thermo.canonical_steps(m)]
    per_group || return k
    g = Pop.group_of_step(m)
    [any(k[i] for i in m.groups[g[j]]) for j in eachindex(m.steps)]
end

"""RE flags of the inactive copy: an `:OnlyA` RE chemistry step becomes steady state,
so that both its rate constants can vanish (an RE step has no speed to scale)."""
function re_flags(m::Pop.PMech, st::Symbol; untie = false)
    st === :A && return [s.re for s in m.steps]
    tg, ch = step_tags(m), is_chem(m; per_group = !untie)
    [s.re && !(tg[j] === :OnlyA && ch[j]) for (j, s) in enumerate(m.steps)]
end

"""`m` with the RE flags `re` (structure only, for Thermo and Engine)."""
function with_re(m::Pop.PMech, re::AbstractVector{Bool})
    steps = [Pop.PStep(s.from, s.to, s.consumed, s.released, re[j])
             for (j, s) in enumerate(m.steps)]
    Pop.PMech(; id = m.id, source = m.source, substrates = m.substrates,
              products = m.products, inhibitors = m.inhibitors, regulators = m.regulators,
              metabolites = m.metabolites, species = m.species, steps = steps,
              groups = m.groups)
end

# ─── Joint A/I solve ─────────────────────────────────────────────────

"""A bare constant of one conformation: state, step, slot (Thermo convention), name and
Thermo priority."""
struct JCol
    st::Symbol
    j::Int
    r::Int
    name::String
    prio::NTuple{4, Int}
end

"""The joint log-linear solve over both conformations' bare constants. `cols` are in
elimination order; `dep[i] = (Keq exponent, [indep column => exponent])`."""
struct Joint
    m::Pop.PMech
    tags::Vector{Symbol}
    chem::Vector{Bool}
    re::Dict{Symbol, Vector{Bool}}
    cs::Dict{Symbol, Vector{Thermo.CStep}}
    cols::Vector{JCol}
    idx::Dict{Tuple{Symbol, Int, Int}, Int}
    C::Matrix{Int}
    rhs::Vector{Int}
    indep::Vector{Int}
    dep::Dict{Int, Tuple{Q, Vector{Pair{Int, Q}}}}
end

"""Column elimination order. `:pkg` eliminates the inactive copies first and, within a
state, follows Thermo's priority (larger key first, the larger name first); `:rev` reverses
the within-state order; `:afirst` eliminates the active copies first; `(:perm, seed)` uses
a random within-state order."""
function order_cols(cols::Vector{JCol}, mode)
    stkey(c) = mode === :afirst ? (c.st === :A ? 1 : 0) : (c.st === :I ? 1 : 0)
    if mode === :pkg || mode === :afirst
        return sort(cols; by = c -> (stkey(c), c.prio, c.name), rev = true)
    elseif mode === :rev
        return sort(cols; by = c -> (stkey(c), (.-c.prio)..., reverse(c.name)), rev = true)
    elseif mode isa Tuple && mode[1] === :perm
        rng = Random.Xoshiro(mode[2])
        r = Dict(c.name => rand(rng) for c in cols)
        return sort(cols; by = c -> (stkey(c), r[c.name]), rev = true)
    end
    error("unknown order $mode")
end

"""Steps of one state in canonical orientation with the state's RE flags."""
function state_cs(m::Pop.PMech, re::Vector{Bool})
    cs = Thermo.canonical_steps(m)
    fk = Thermo.formkey.(m.species)
    map(enumerate(cs)) do (j, c)
        c.re == re[j] && return c
        sa, sb = Thermo.sidekey(fk[c.a], c.cons), Thermo.sidekey(fk[c.b], c.rel)
        Thermo.CStep(c.a, c.b, c.cons, c.rel, false, c.flipped, c.kind, c.occ, c.fwd, false,
                     replace(c.key, "RE " => "SS "),
                     [string("k:", sa, " -> ", sb), string("k:", sb, " -> ", sa)])
    end
end

"""Solve the joint system of `m` (allosteric tags in `m.group_tags`): the cycles and
kinetic-group rows of both conformations, and an `A = I` tie for every constant of an
`:EqualAI` group. `order` picks the dependents (see `order_cols`). With `untie`, the
inactive copies of an `:OnlyA` group are not tied to each other (the group does not exist
in the inactive state) and each step keeps its own class (binding or chemistry)."""
function joint_solve(m::Pop.PMech; order = :pkg, untie = false)
    tags, chem = step_tags(m), is_chem(m; per_group = !untie)
    re = Dict(st => re_flags(m, st; untie) for st in (:A, :I))
    cs = Dict(st => state_cs(m, re[st]) for st in (:A, :I))
    raw = JCol[]
    for st in (:A, :I), c in Thermo.columns(cs[st])
        push!(raw, JCol(st, c.j, c.r, string(st, "|", c.name), c.prio))
    end
    cols = order_cols(raw, order)
    idx = Dict((c.st, c.j, c.r) => i for (i, c) in enumerate(cols))
    n = length(cols)
    rows, rhs = Vector{Int}[], Int[]
    for st in (:A, :I)
        loc = [(c.j, c.r) for c in cols if c.st === st]
        lidx = Dict(k => i for (i, k) in enumerate(loc))
        embed(r) = (g = zeros(Int, n);
                    for (k, i) in lidx; g[idx[(st, k...)]] += r[i]; end; g)
        cyc = Thermo.fundamental_cycles(length(m.species), cs[st])
        deltas = Thermo.step_delta(cs[st])
        for z in cyc
            push!(rows, embed(Thermo.cycle_row(cs[st], lidx, z)))
            push!(rhs, Thermo.cycle_multiple(m, deltas, z))
        end
        grps = untie && st === :I ? [g for g in m.groups if tags[g[1]] !== :OnlyA] : m.groups
        al = Thermo.group_alignment(m, cs[st], grps)
        for r in Thermo.group_rows(cs[st], lidx, al)
            push!(rows, embed(r)); push!(rhs, 0)
        end
    end
    for (j, t) in enumerate(tags)
        t === :EqualAI || continue
        for r in 1:(cs[:A][j].re ? 1 : 2)
            row = zeros(Int, n); row[idx[(:I, j, r)]] += 1; row[idx[(:A, j, r)]] -= 1
            push!(rows, row); push!(rhs, 0)
        end
    end
    C = isempty(rows) ? zeros(Int, 0, n) : permutedims(reduce(hcat, rows))
    M = zeros(Q, size(C, 1), n + 1)
    M[:, 1:n] .= C; M[:, n + 1] .= rhs
    piv, _ = Thermo.rref!(M)
    (n + 1) in piv && error("$(m.id): the joint constraints reduce to 0 = log Keq")
    indep = setdiff(1:n, piv)
    dep = Dict{Int, Tuple{Q, Vector{Pair{Int, Q}}}}()
    for (r, p) in enumerate(piv)
        dep[p] = (M[r, n + 1], [f => -M[r, f] for f in indep if M[r, f] != 0])
    end
    Joint(m, tags, chem, re, cs, cols, idx, C, rhs, indep, dep)
end

# ─── Scaling conventions ─────────────────────────────────────────────

"""Intended ε-exponent of every column under `conv` (`:i` or `:ii`): an inactive copy of an
`:OnlyA` RE binding K → ∞ (-1); an SS binding: under `:i` kon → 0 (+1) with koff fixed, under
`:ii` both → 0; chemistry (SS, or RE turned SS): both → 0. Every other column 0."""
function intended(J::Joint, conv::Symbol)
    a = zeros(Q, length(J.cols))
    for (i, c) in enumerate(J.cols)
        (c.st === :I && J.tags[c.j] === :OnlyA) || continue
        s = J.cs[:I][c.j]
        if J.chem[c.j]
            a[i] = 1
        elseif s.re
            a[i] = -1                          # Kd of a plain RE binding
        else                                   # SS plain binding; slot 1 takes it up
            a[i] = c.r == 1 ? 1 : (conv === :ii ? 1 : 0)
        end
    end
    a
end

"""Realized exponents along the path: the intended ones on the independents, each
dependent's from its row."""
function realized(J::Joint, aint::Vector{Q})
    a = zeros(Q, length(J.cols))
    for i in J.indep; a[i] = aint[i]; end
    for (p, (_, terms)) in J.dep
        a[p] = sum((e * aint[f] for (f, e) in terms); init = zero(Q))
    end
    a
end

"""Columns whose realized exponent breaks the convention: `partner` lists the columns that
are not inactive `:OnlyA` copies but move; `wrong` the inactive `:OnlyA` columns whose
realized sign differs from the intended one (a binding K not → ∞, a rate not → 0, a koff
that moves under `:i`)."""
function path_report(J::Joint, aint, a)
    partner, wrong = Int[], Int[]
    for (i, c) in enumerate(J.cols)
        onlyI = c.st === :I && J.tags[c.j] === :OnlyA
        if !onlyI
            a[i] != 0 && push!(partner, i)
        elseif sign(a[i]) != sign(aint[i])
            push!(wrong, i)
        end
    end
    partner, wrong
end

"""Sign every moved column must take along an :OnlyA path under `conv`: -1 for an RE
binding's Kd (→ ∞), +1 for a rate constant that vanishes, 0 for an SS binding's koff
under `:i`. Only inactive copies of :OnlyA groups move."""
function cone_signs(J::Joint, conv::Symbol)
    sig = Dict{Int, Int}()
    for (i, c) in enumerate(J.cols)
        (c.st === :I && J.tags[c.j] === :OnlyA) || continue
        s = J.cs[:I][c.j]
        sig[i] = J.chem[c.j] ? 1 : s.re ? -1 : (c.r == 1 || conv === :ii ? 1 : 0)
    end
    sig
end

"""Rational basis of the nullspace of an integer matrix."""
function nullbasis(A::Matrix{Int})
    nc = size(A, 2)
    M = zeros(Q, size(A, 1), nc + 1); M[:, 1:nc] .= A
    piv, _ = Thermo.rref!(M)
    map(setdiff(1:nc, piv)) do f
        v = zeros(Q, nc); v[f] = 1
        for (r, p) in enumerate(piv); v[p] = -M[r, f]; end
        v
    end
end

"""Feasible :OnlyA directions: exponents `a` that vanish on every column the tags do not
move, satisfy `C·a = 0` (every Haldane, Wegscheider, group and tie relation holds for all
ε), and give each moved column its required sign (`cone_signs`). Searches integer
combinations y ∈ {-R..R}^d of a nullspace basis; returns the distinct primitive
directions sorted by L1 norm, and `d`."""
function cone_directions(J::Joint, conv::Symbol; R = 3, maxgrid = 100_000)
    sig = cone_signs(J, conv)
    O = sort!(collect(keys(sig)))
    n = length(J.cols)
    isempty(O) && return [zeros(Q, n)], 0
    rows = [J.C[r, O] for r in axes(J.C, 1) if any(!iszero, J.C[r, O])]
    for (t, i) in enumerate(O)
        sig[i] == 0 && push!(rows, [u == t ? 1 : 0 for u in eachindex(O)])
    end
    A = isempty(rows) ? zeros(Int, 0, length(O)) : permutedims(reduce(hcat, rows))
    B = nullbasis(A)
    d = length(B)
    d == 0 && return Vector{Q}[], 0
    (2R + 1)^d > maxgrid && (R = max(1, floor(Int, (maxgrid^(1 / d) - 1) / 2)))
    found = Set{Vector{Int}}()
    for y in Iterators.product(ntuple(_ -> -R:R, d)...)
        (2R + 1)^d > 4maxgrid && break             # too many dimensions: sample only
        all(iszero, y) && continue
        aO = sum(y[k] * B[k] for k in 1:d)
        all(t -> (s = sig[O[t]]; s == 0 ? aO[t] == 0 : sign(aO[t]) == s), eachindex(O)) ||
            continue
        den = foldl(lcm, (denominator(x) for x in aO); init = 1)
        v = [Int(x * den) for x in aO]
        push!(found, v .÷ foldl(gcd, v; init = 0))
    end
    # random interior points by alternating projections (the nullspace and the sign box),
    # rounded to exact rational directions and checked exactly
    if d > 0
        N = Float64.(reduce(hcat, B))
        G = N' * N
        rng = Random.Xoshiro(17)
        box(x) = [(s = sig[O[t]]; s == 0 ? 0.0 : s > 0 ? max(x[t], 1.0) : min(x[t], -1.0))
                  for t in eachindex(O)]
        for _ in 1:12
            x = [(s = sig[O[t]]; s * (1 + 4rand(rng))) for t in eachindex(O)]
            y = zeros(d)
            for _ in 1:400
                y = G \ (N' * box(x)); x = N * y
            end
            for K in (1, 2, 4, 8, 16, 64)
                yq = [Q(round(Int, K * v), K) for v in y]
                aO = sum(yq[k] * B[k] for k in 1:d)
                all(t -> (s = sig[O[t]]; s == 0 ? aO[t] == 0 : sign(aO[t]) == s), eachindex(O)) || continue
                den = foldl(lcm, (denominator(x) for x in aO); init = 1)
                v = [Int(x * den) for x in aO]
                push!(found, v .÷ foldl(gcd, v; init = 0))
                break
            end
        end
    end
    dirs = sort!(collect(found); by = v -> (sum(abs, v), v))
    [(a = zeros(Q, n); for (t, i) in enumerate(O); a[i] = v[t]; end; a) for v in dirs], d
end

# ─── Values ──────────────────────────────────────────────────────────

"""Every bare column value from independent values and Keq."""
function all_values(J::Joint, x_indep::AbstractDict{Int}, keq)
    setprecision(BigFloat, PREC) do
        x = Dict{Int, BigFloat}(i => BigFloat(x_indep[i]) for i in J.indep)
        for (p, (k, terms)) in J.dep
            v = BigFloat(keq)^BigFloat(k)
            for (f, e) in terms; v *= x[f]^BigFloat(e); end
            x[p] = v
        end
        x
    end
end

"""Values along the path at `eps`: independents `hat * eps^a`, dependents from the rows."""
path_values(J::Joint, hat::AbstractDict{Int}, a, eps, keq) =
    all_values(J, Dict(i => setprecision(() -> BigFloat(hat[i]) * BigFloat(eps)^BigFloat(a[i]),
                                         BigFloat, PREC) for i in J.indep), keq)

"""Per-step constants of state `st` in the stored orientation, Oracle convention (`[K]`
RE, `[kf, kr]` SS)."""
function stored(J::Joint, x::AbstractDict{Int}, st::Symbol)
    map(eachindex(J.m.steps)) do j
        c = J.cs[st][j]
        if c.re
            K = x[J.idx[(st, j, 1)]]
            [c.flipped && !c.bindK ? 1 / K : K]
        else
            k1, k2 = x[J.idx[(st, j, 1)]], x[J.idx[(st, j, 2)]]
            c.flipped ? [k2, k1] : [k1, k2]
        end
    end
end

"""Random log-uniform independents in [1e-2, 1e2] and a Keq."""
random_hat(J::Joint, rng) =
    Dict(i => big1024(10.0^(4 * rand(rng) - 2)) for i in J.indep)

# ─── Networks ────────────────────────────────────────────────────────

"""Steady state of a network on `n` nodes with edges `(a, b, rf, rr, ω)` (`rf` drives
a → b, `rr` b → a, ω the copies of the first substrate one a → b turnover consumes):
populations summing to 1 and the net consumption of the first substrate."""
function solve_net(n::Int, edges)
    setprecision(BigFloat, PREC) do
        M = zeros(BigFloat, n, n)
        for (a, b, rf, rr, _) in edges
            M[b, a] += rf; M[a, a] -= rf; M[a, b] += rr; M[b, b] -= rr
        end
        M[n, :] .= 1
        rhs = zeros(BigFloat, n); rhs[n] = 1
        x = M \ rhs
        x, sum((ω * (rf * x[a] - rr * x[b]) for (a, b, rf, rr, ω) in edges if ω != 0);
               init = big(0))
    end
end

"""Catalytic edges of one conformation over `nodeof(form)` (or `nothing` to skip):
`k[j]` per step in stored orientation, `re[j]` its kind, `active[j]` whether it exists.
RE steps are fast pairs at `Λ`."""
function cat_edges(m::Pop.PMech, k, re, active, concs, Λ, nodeof)
    x1 = m.substrates[1]
    c(x) = BigFloat(concs[x.name])
    out = Tuple{Int, Int, BigFloat, BigFloat, Int}[]
    for (j, s) in enumerate(m.steps)
        active[j] || continue
        a, b = nodeof(s.from), nodeof(s.to)
        (a === nothing || b === nothing) && continue
        f, r = re[j] ? (Pop.is_binding(s) ? (Λ, Λ * k[j][1]) : (Λ, Λ / k[j][1])) :
               (BigFloat(k[j][1]), BigFloat(k[j][2]))
        rf = f * prod((c(y) for y in s.consumed); init = big(1))
        rr = r * prod((c(y) for y in s.released); init = big(1))
        ω = count(y -> y.name === x1, s.consumed) - count(y -> y.name === x1, s.released)
        push!(out, (a, b, rf, rr, ω))
    end
    out
end

"""The free form (no bound metabolite, no residual) that carries the flip."""
free_form(m::Pop.PMech) = findfirst(Pop.is_free, m.species)

"""Forms reachable from every free form over the `active` steps (the package's pruning)."""
function reachable(m::Pop.PMech, active)
    reach = Set(i for (i, sp) in enumerate(m.species) if Pop.is_free(sp))
    changed = true
    while changed
        changed = false
        for (j, s) in enumerate(m.steps)
            active[j] || continue
            if s.from in reach && !(s.to in reach)
                push!(reach, s.to); changed = true
            elseif s.to in reach && !(s.from in reach)
                push!(reach, s.from); changed = true
            end
        end
    end
    reach
end

"""Single-protomer steady state of one conformation: the free fraction and the flux per
protomer. Forms outside `keep` hold no mass."""
function conformation_state(m, k, re, active, concs, Λ; keep = 1:length(m.species))
    forms = sort!(collect(keep))
    loc = Dict(f => i for (i, f) in enumerate(forms))
    edges = cat_edges(m, k, re, active, concs, Λ, f -> get(loc, f, nothing))
    length(forms) == 1 && return big1024(1), big1024(0)
    x, J = solve_net(length(forms), edges)
    x[loc[free_form(m)]], J
end

"""Regulatory-site factor of one conformation: Π_sites (1 + Σ [lig]/K)^m over the ligands
present in that conformation (`regK[(site, lig)]` is `nothing` when absent)."""
function W_factor(m::Pop.PMech, regK, concs)
    setprecision(BigFloat, PREC) do
        prod((begin
                  z = big(1)
                  for lig in site.ligands
                      K = regK[(si, lig)]
                      K === nothing || (z += BigFloat(concs[lig]) / K)
                  end
                  z^site.multiplicity
              end for (si, site) in enumerate(m.sites)); init = big(1))
    end
end

"""Formulation-1 rate per site at catalytic multiplicity `n`, from each conformation's own
single-protomer network: an oligomer's conformation is a product of independent protomer
chains, and only the all-free oligomer flips, so conformation c holds weight
`L_c · W_c / π_c(free)^n` and turns over at its protomer flux. `kI`/`reI`/`activeI`/`keepI`
describe the inactive network (limit or deletion)."""
function v_factorized(m, kA, kI, reI, activeI, keepI, L, regKA, regKI, concs; n, Λ)
    setprecision(BigFloat, PREC) do
        reA = [s.re for s in m.steps]
        πA, JA = conformation_state(m, kA, reA, trues(length(m.steps)), concs, Λ)
        πI, JI = conformation_state(m, kI, reI, activeI, concs, Λ; keep = keepI)
        wA = W_factor(m, regKA, concs) / πA^n
        wI = BigFloat(L) * W_factor(m, regKI, concs) / πI^n
        (wA * JA + wI * JI) / (wA + wI)
    end
end

"""Regulator occupancy states: one entry per site copy, 0 (empty) or a ligand index."""
reg_states(m::Pop.PMech) = vec(collect(Iterators.product(
    [0:length(site.ligands) for site in m.sites for _ in 1:site.multiplicity]...)))

"""Explicit formulation-1 network for `nprot` protomers (1 or 2): states are
(conformation, protomer forms, regulator occupancy); every protomer runs its own catalytic
steps, regulators bind their site copies at rapid equilibrium, and only the all-free,
regulator-free state flips (A → I at rate L, back at 1). Returns the rate per site."""
function v_explicit(m, kA, kI, reI, activeI, keepI, L, regKA, regKI, concs; nprot, Λ)
    setprecision(BigFloat, PREC) do
        nf = length(m.species)
        regs = isempty(m.sites) ? [()] : reg_states(m)
        site_of = [si for (si, site) in enumerate(m.sites) for _ in 1:site.multiplicity]
        keep = Dict(:A => sort!(collect(1:nf)), :I => sort!(collect(keepI)))
        forms = Dict(st => vec([Tuple(t) for t in Iterators.product(ntuple(_ -> keep[st], nprot)...)])
                     for st in (:A, :I))
        # a conformation holds only the occupancies of the ligands it binds
        present(st, r) = all(ci -> r[ci] == 0 ||
            (st === :A ? regKA : regKI)[(site_of[ci], m.sites[site_of[ci]].ligands[r[ci]])] !== nothing,
            eachindex(r))
        rs = Dict(st => [r for r in regs if present(st, r)] for st in (:A, :I))
        nodes = Dict{Any, Int}()
        for st in (:A, :I), x in forms[st], r in rs[st]
            nodes[(st, x, r)] = length(nodes) + 1
        end
        edges = Tuple{Int, Int, BigFloat, BigFloat, Int}[]
        reA = [s.re for s in m.steps]
        x1 = m.substrates[1]
        cc(y) = BigFloat(concs[y.name])
        for st in (:A, :I)
            k, re, active = st === :A ? (kA, reA, trues(length(m.steps))) : (kI, reI, activeI)
            regK = st === :A ? regKA : regKI
            # every protomer runs its own catalytic steps; each edge is listed once, from
            # the joint state whose protomer p sits on the step's `from` form
            for x in forms[st], r in rs[st], p in 1:nprot, (j, s) in enumerate(m.steps)
                (active[j] && s.from == x[p] && s.to in keep[st]) || continue
                y = ntuple(q -> q == p ? s.to : x[q], nprot)
                f, rv = re[j] ? (Pop.is_binding(s) ? (Λ, Λ * k[j][1]) : (Λ, Λ / k[j][1])) :
                        (BigFloat(k[j][1]), BigFloat(k[j][2]))
                rf = f * prod((cc(z) for z in s.consumed); init = big(1))
                rr = rv * prod((cc(z) for z in s.released); init = big(1))
                ω = count(z -> z.name === x1, s.consumed) - count(z -> z.name === x1, s.released)
                push!(edges, (nodes[(st, x, r)], nodes[(st, y, r)], rf, rr, ω))
            end
            for x in forms[st], r in rs[st], (ci, si) in enumerate(site_of)
                r[ci] == 0 || continue
                for (li, lig) in enumerate(m.sites[si].ligands)
                    K = regK[(si, lig)]
                    K === nothing && continue
                    r2 = ntuple(q -> q == ci ? li : r[q], length(r))
                    push!(edges, (nodes[(st, x, r)], nodes[(st, x, r2)],
                                  Λ * BigFloat(concs[lig]), Λ * K, 0))
                end
            end
        end
        ff = free_form(m)
        x0 = ntuple(_ -> ff, nprot); r0 = isempty(m.sites) ? () : ntuple(_ -> 0, length(site_of))
        push!(edges, (nodes[(:A, x0, r0)], nodes[(:I, x0, r0)], BigFloat(L), big(1), 0))
        _, J = solve_net(length(nodes), edges)
        J / nprot
    end
end

# ─── Symbolic laws and leading faces ─────────────────────────────────

"""Shared symbol table of a combined law: independents `A|…`/`I|…`, `Keq`, `L`, the
metabolites and the regulator constants."""
mutable struct Table
    T::Poly.Vars
end

"""The bare law of state `st` (or a deletion graph) as `(N, Q, D)` over `tab`, with every
bare constant replaced by its monomial in the joint solve's independents and Keq."""
function state_law(J::Joint, st::Symbol, tab::Table)
    m = J.m
    pm = with_re(m, J.re[st])
    law = Engine.bare_rate_law(pm; check = false)
    T = tab.T
    keq = Poly.var!(T, :Keq)
    colmono(i) = haskey(J.dep, i) ?
        (begin
             k, terms = J.dep[i]
             all(isinteger, [k; last.(terms)]) || error("$(m.id): non-integer exponent")
             Poly.mono((keq => Int(k)),
                       (Poly.var!(T, Symbol(J.cols[f].name)) => Int(e) for (f, e) in terms)...)
         end) : Poly.mono(Poly.var!(T, Symbol(J.cols[i].name)) => 1)
    sub = Dict{Int, Poly.Mono}()
    for name in m.metabolites
        haskey(law.vars.index, name) && (sub[law.vars.index[name]] = Poly.mono(Poly.var!(T, name) => 1))
    end
    for (j, vs) in enumerate(law.stepvars)
        c = J.cs[st][j]
        if c.re
            mn = colmono(J.idx[(st, j, 1)])
            sub[vs[1]] = c.flipped && !c.bindK ? Poly.mpow(mn, -1) : mn
        else
            m1, m2 = colmono(J.idx[(st, j, 1)]), colmono(J.idx[(st, j, 2)])
            sub[vs[1]], sub[vs[2]] = c.flipped ? (m2, m1) : (m1, m2)
        end
    end
    f(v) = sub[v]
    Poly.subst(law.num, f), Poly.subst(law.den, f), Poly.subst(law.dfree, f)
end

"""Regulatory factor polynomial Π_sites (1 + Σ x_lig · K^-1)^m of one conformation, with
`regname[(site, lig)]` the constant's symbol or `nothing` when the ligand is absent."""
function W_poly(m::Pop.PMech, regname, tab::Table)
    T = tab.T
    W = one(Poly.P{Int})
    for (si, site) in enumerate(m.sites)
        z = one(Poly.P{Int})
        for lig in site.ligands
            K = regname[(si, lig)]
            K === nothing && continue
            z += Poly.term(Poly.P{Int}, Poly.mono(Poly.var!(T, lig) => 1, Poly.var!(T, K) => -1))
        end
        for _ in 1:site.multiplicity; W *= z; end
    end
    W
end

# An expression over polynomials: (:leaf, P), (:prod, [e...]), (:sum, [e...]).
leaf(p) = (:leaf, p)
eprod(xs...) = (:prod, collect(Any, xs))
esum(xs...) = (:sum, collect(Any, xs))
epow(x, n) = (:prod, Any[x for _ in 1:n])

"""Full expansion of an expression."""
expand_expr(e) = e[1] === :leaf ? e[2] :
    e[1] === :prod ? foldl(*, (expand_expr(x) for x in e[2]); init = one(Poly.P{Int})) :
    foldl(+, (expand_expr(x) for x in e[2]); init = zero(Poly.P{Int}))

"""Order vector of a monomial under the lexicographic directions `dirs` (each a
Dict var => weight)."""
mord(mn, dirs) = Tuple(sum((get(d, v, zero(Q)) * e for (v, e) in mn); init = zero(Q)) for d in dirs)

"""Initial form of a polynomial: the terms of lexicographically least order."""
function initial(p::Poly.P{Int}, dirs)
    iszero(p) && return nothing, p
    best = minimum(mord(mn, dirs) for mn in keys(p.t))
    best, Poly.P{Int}(Dict(mn => c for (mn, c) in p.t if mord(mn, dirs) == best))
end

"""Initial form of an expression: products multiply initial forms; a sum keeps the
children of least order and, if their initial forms cancel, expands the sum.
Returns (order, polynomial, cancelled::Bool)."""
function initial_expr(e, dirs)
    if e[1] === :leaf
        o, p = initial(e[2], dirs)
        return o, p, false
    elseif e[1] === :prod
        o, p, c = ntuple(_ -> zero(Q), length(dirs)), one(Poly.P{Int}), false
        for x in e[2]
            ox, px, cx = initial_expr(x, dirs)
            ox === nothing && return nothing, zero(Poly.P{Int}), c || cx
            o = o .+ ox
            p *= px; c |= cx
        end
        return o, p, c
    else
        parts = [initial_expr(x, dirs) for x in e[2]]
        parts = [q for q in parts if q[1] !== nothing]
        isempty(parts) && return nothing, zero(Poly.P{Int}), false
        best = minimum(q[1] for q in parts)
        p = foldl(+, (q[2] for q in parts if q[1] == best); init = zero(Poly.P{Int}))
        c = any(q[3] for q in parts)
        if iszero(p)
            o2, p2 = initial(expand_expr(e), dirs)
            return o2, p2, true
        end
        return best, p, c
    end
end

"""The formulation-1 combine at multiplicity n in polynomial (cross-weighted) form:
num = N_A Q_A^(n-1) D_I^n W_A + L N_I Q_I^(n-1) D_A^n W_I, den likewise with Q^n."""
function combine(lawA, lawI, WA, WI, Lp, n)
    NA, QA, DA = lawA; NI, QI, DI = lawI
    num = esum(eprod(leaf(NA), epow(leaf(QA), n - 1), epow(leaf(DI), n), leaf(WA)),
               eprod(leaf(Lp), leaf(NI), epow(leaf(QI), n - 1), epow(leaf(DA), n), leaf(WI)))
    den = esum(eprod(epow(leaf(QA), n), epow(leaf(DI), n), leaf(WA)),
               eprod(leaf(Lp), epow(leaf(QI), n), epow(leaf(DA), n), leaf(WI)))
    num, den
end

"""Leading-order value of num/den under `dirs` at `vals` (indexed by table variable):
the ratio of initial forms when the orders agree, 0 when the numerator is of higher
order, ±Inf when lower. Also returns both orders and whether a cancellation occurred."""
function lead_value(num, den, dirs, vals)
    on, pn, cn = initial_expr(num, dirs)
    od, pd, cd = initial_expr(den, dirs)
    v = setprecision(BigFloat, PREC) do
        on === nothing ? big(0) :
        on == od ? Poly.evaluate(pn, vals) / Poly.evaluate(pd, vals) :
        on > od ? big(0) : sign(Poly.evaluate(pn, vals)) * BigFloat(Inf)
    end
    (v = v, on = on, od = od, cancel = cn || cd, pn = pn, pd = pd)
end

end # module
