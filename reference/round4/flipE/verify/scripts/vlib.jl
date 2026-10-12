# ABOUTME: Independent two-conformation mass-action oracle: maps parameter names to steps by
# ABOUTME: the naming convention, solves Haldane itself, flips only the all-free-E oligomer.
using EnzymeRates, LinearAlgebra, Random, Printf
const ER = EnzymeRates
setprecision(BigFloat, 512)

"""One step as written: `from + cons -> to + rel`, rapid equilibrium or not, with its tag.
Forms are package species names (conformation followed by the bound metabolites)."""
struct VStep
    from::String
    cons::Vector{String}
    to::String
    rel::Vector{String}
    re::Bool
    tag::Symbol
end

"""DSL text of a form: `ES` is written `E(S)` through the `dslform` map, others as is."""
dslside(f, ms, dslform) = join([get(dslform, f, f); ms], " + ")
dslline(s::VStep, dslform) = dslside(s.from, s.cons, dslform) * (s.re ? " ⇌ " : " <--> ") *
                             dslside(s.to, s.rel, dslform) * " :: " * String(s.tag)

function dsltext(subs, prods, steps, dslform; n = 1)
    "@allosteric_mechanism begin\n    substrates: $(join(subs, ", "))\n" *
    "    products: $(join(prods, ", "))\n    catalytic_multiplicity: $n\n" *
    "    catalytic_steps: begin\n" *
    join(["        " * dslline(s, dslform) for s in steps], "\n") * "\n    end\nend"
end

build(src) = try
    Core.eval(Main, Meta.parse(src))
catch err
    err isa LoadError && (err = err.error)
    sprint(showerror, err)
end

tagpre(t, σ) = t === :EqualAI ? "" : (σ === :A ? "A_" : "I_")
instate(s::VStep, σ) = σ === :A || s.tag !== :OnlyA
isbinding(s::VStep) = length(s.cons) + length(s.rel) == 1

"""Names, in conformation σ, of the step's constants and its forward equilibrium ratio
[to]·Πrel/([from]·Πcons) as a function of a name→value lookup `v`. Every candidate name
the convention allows is tried; `have(name)` says whether the package knows a name."""
function stepconsts(s::VStep, σ, have)
    t = tagpre(s.tag, σ)
    sd(f, ms) = join([f; ms], "_")
    if !s.re
        kf = "k_" * t * sd(s.from, s.cons) * "_to_" * sd(s.to, s.rel)
        kr = "k_" * t * sd(s.to, s.rel) * "_to_" * sd(s.from, s.cons)
        return (names = [kf, kr], ratio = v -> v(kf) / v(kr), rates = v -> (v(kf), v(kr)))
    end
    if isbinding(s)
        if length(s.cons) == 1        # F + M ⇌ B: dissociation constant K_B_to_F_M
            K = "K_" * t * s.to * "_to_" * s.from * "_" * s.cons[1]
            return (names = [K], ratio = v -> 1 / v(K), rates = nothing)
        else                          # B ⇌ F + M
            K = "K_" * t * s.from * "_to_" * s.to * "_" * s.rel[1]
            return (names = [K], ratio = v -> v(K), rates = nothing)
        end
    end
    K1 = "K_" * t * s.from * "_to_" * s.to
    K2 = "K_" * t * s.to * "_to_" * s.from
    have(K1) && return (names = [K1], ratio = v -> v(K1), rates = nothing)
    have(K2) && return (names = [K2], ratio = v -> 1 / v(K2), rates = nothing)
    error("no RE constant named $K1 or $K2")
end

"""Steps of the unique cycle (the 2-core of the form graph) in traversal order with
orientation ±1, starting at free E."""
function cycle_of(steps)
    deg = Dict{String, Int}()
    alive = trues(length(steps))
    for s in steps; deg[s.from] = get(deg, s.from, 0) + 1; deg[s.to] = get(deg, s.to, 0) + 1; end
    changed = true
    while changed
        changed = false
        for (j, s) in enumerate(steps)
            alive[j] || continue
            if deg[s.from] == 1 || deg[s.to] == 1
                alive[j] = false; deg[s.from] -= 1; deg[s.to] -= 1; changed = true
            end
        end
    end
    cyc = findall(alive)
    isempty(cyc) && return Tuple{Int, Int}[]
    start = any(j -> "E" in (steps[j].from, steps[j].to), cyc) ? "E" : steps[cyc[1]].from
    out, cur, used = Tuple{Int, Int}[], start, Set{Int}()
    while length(out) < length(cyc)
        j = findfirst(j -> !(j in used) && cur in (steps[j].from, steps[j].to), cyc)
        jj = cyc[j]; push!(used, j)
        if steps[jj].from == cur
            push!(out, (jj, 1)); cur = steps[jj].to
        else
            push!(out, (jj, -1)); cur = steps[jj].from
        end
    end
    out
end

"""All constant values: the package's fitted ones plus each conformation's dependent one,
solved from that conformation's Haldane relation over the cycle (Π oriented ratios = Keq
along the substrate-to-product direction). Asserts every non-cycle constant is fitted."""
function solve_values(steps, fitted::Dict{String, BigFloat}, Keq, subs, prods)
    vals = copy(fitted)
    cyc = cycle_of(steps)
    for σ in (:A, :I)
        active = [(j, o) for (j, o) in cyc if instate(steps[j], σ)]
        length(active) == length(cyc) || continue         # the cycle is broken in σ
        # Orientation: does the traversal consume the substrates?
        net = Dict{String, Int}()
        for (j, o) in active, (ms, sg) in ((steps[j].cons, -o), (steps[j].rel, o)), m in ms
            net[m] = get(net, m, 0) + sg
        end
        dir = all(get(net, x, 0) == -1 for x in subs) ? 1 :
              all(get(net, x, 0) == 1 for x in subs) ? -1 : error("bad cycle orientation")
        target = Keq^dir
        unknown = String[]
        for (j, o) in active
            c = stepconsts(steps[j], σ, x -> haskey(vals, x) || haskey(fitted, x))
            for x in c.names
                haskey(vals, x) || push!(unknown, x)
            end
        end
        unknown = unique(unknown)
        if isempty(unknown)
            # Haldane must hold already.
            prodv = prod(begin
                c = stepconsts(steps[j], σ, x -> haskey(vals, x))
                c.ratio(x -> vals[x])^o
            end for (j, o) in active)
            abs(prodv / target - 1) < big(10)^-60 ||
                error("Haldane violated in $σ with no free constant: $(Float64(prodv / target))")
            continue
        end
        length(unknown) == 1 || error("state $σ: unknowns $unknown")
        x = only(unknown)
        # The oriented product is linear in x or in 1/x: probe with x = 1 and x = 2.
        f(xv) = (vals[x] = xv; prod(begin
            c = stepconsts(steps[j], σ, y -> haskey(vals, y))
            c.ratio(y -> vals[y])^o
        end for (j, o) in active))
        p1, p2 = f(big(1)), f(big(2))
        e = round(Int, log2(p2 / p1))
        abs(e) == 1 || error("exponent $e for $x")
        vals[x] = (target / p1)^(1 // e)
        abs(f(vals[x]) / target - 1) < big(10)^-100 || error("Haldane solve failed")
    end
    # Every constant of every step present in some state must now be known.
    for s in steps, σ in (:A, :I)
        instate(s, σ) || continue
        c = stepconsts(s, σ, y -> haskey(vals, y))
        for x in c.names
            haskey(vals, x) || error("constant $x of step $(s.from)->$(s.to) in $σ unknown")
        end
    end
    vals
end

"""Per-protomer steady-state rate (net uptake of `x1`) of the concerted network: `n`
protomers independent within a conformation; only the all-`flip` oligomer interconverts
between conformations, with [flip…flip]_I/[flip…flip]_A = L (a bridge edge, so its speed
is irrelevant). RE steps run at Λ·(ratio·Πcons, Πrel); `v∞ ≈ 2v(2Λ) − v(Λ)`. Inactive
forms that the flip form cannot reach along σ = I steps are left out."""
function oracle_rate(steps, vals, concs, L, x1; n = 1, flip = "E", Λ = big(10)^40,
                     blocked = Dict{Int, Tuple{BigFloat, BigFloat}}())
    function once(Λ)
        c(ms) = prod((concs[m] for m in ms); init = big(1))
        trans = Dict(σ => Tuple{String, String, BigFloat, BigFloat, Int}[] for σ in (:A, :I))
        # `blocked[j] = (kf, kr)` keeps :OnlyA step j in the inactive conformation as a
        # steady-state step with those (already slowed) rates instead of deleting it.
        for (j, (kf, kr)) in blocked
            s = steps[j]
            push!(trans[:I], (s.from, s.to, kf * c(s.cons), kr * c(s.rel),
                              count(==(x1), s.cons) - count(==(x1), s.rel)))
        end
        for σ in (:A, :I), s in steps
            instate(s, σ) || continue
            k = stepconsts(s, σ, y -> haskey(vals, y))
            look(y) = vals[y]
            rf, rr = s.re ? (Λ * k.ratio(look) * c(s.cons), Λ * c(s.rel)) :
                            (k.rates(look)[1] * c(s.cons), k.rates(look)[2] * c(s.rel))
            w = count(==(x1), s.cons) - count(==(x1), s.rel)
            push!(trans[σ], (s.from, s.to, rf, rr, w))
        end
        # Inactive forms reachable from the flip form along σ = I steps.
        reach = Set([flip]); ch = true
        while ch
            ch = false
            for (a, b, _, _, _) in trans[:I]
                (a in reach) != (b in reach) && (push!(reach, a, b); ch = true)
            end
        end
        filter!(t -> t[1] in reach && t[2] in reach, trans[:I])
        forms = Dict(:A => sort!(unique([f for t in trans[:A] for f in (t[1], t[2])])),
                     :I => sort!(collect(reach)))
        nodes = Tuple{Symbol, Vector{String}}[]
        for σ in (:A, :I)
            for st in Iterators.product(ntuple(_ -> forms[σ], n)...)
                push!(nodes, (σ, collect(st)))
            end
        end
        idx = Dict(nd => i for (i, nd) in enumerate(nodes))
        N = length(nodes)
        M = zeros(BigFloat, N, N)
        flux = Tuple{Int, Int, BigFloat, BigFloat, Int}[]
        for (σ, st) in nodes, i in 1:n, (a, b, rf, rr, w) in trans[σ]
            st[i] == a || continue
            st2 = copy(st); st2[i] = b
            p, q = idx[(σ, st)], idx[(σ, st2)]
            M[q, p] += rf; M[p, p] -= rf; M[p, q] += rr; M[q, q] -= rr
            w == 0 || push!(flux, (p, q, rf, rr, w))
        end
        allf = fill(flip, n)
        if haskey(idx, (:I, allf)) && haskey(idx, (:A, allf))
            p, q = idx[(:A, allf)], idx[(:I, allf)]
            M[q, p] += L; M[p, p] -= L; M[p, q] += 1; M[q, q] -= 1
        end
        M[N, :] .= 1
        rhs = zeros(BigFloat, N); rhs[N] = 1
        x = M \ rhs
        v = sum((w * (rf * x[p] - rr * x[q]) for (p, q, rf, rr, w) in flux); init = big(0))
        mI = sum((x[i] for (i, nd) in enumerate(nodes) if nd[1] === :I); init = big(0))
        v / n, mI
    end
    v1, m1 = once(Λ)
    v2, _ = once(2Λ)
    2v2 - v1, m1
end

"""The package's rate for an `AllostericEnzymeMechanism` `m` (E_total = 1)."""
function pkg_rate(m, fitted::Dict{String, BigFloat}, Keq, concs)
    fp = ER.fitted_params(m)
    prm = NamedTuple{(fp..., :Keq, :E_total)}((Tuple(fitted[String(p)] for p in fp)..., Keq,
                                              big(1)))
    mets = Tuple(ER.metabolites(m))
    BigFloat(ER.rate_equation(m, NamedTuple{mets}(Tuple(concs[String(x)] for x in mets)), prm))
end

relerr(a, b) = abs(a - b) / max(abs(b), big(10)^-300)
lu(rng, lo, hi) = exp(BigFloat(log(lo) + (log(hi) - log(lo)) * rand(rng)))
