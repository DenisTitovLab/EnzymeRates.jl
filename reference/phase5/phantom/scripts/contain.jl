# ABOUTME: Numeric family containment both ways between a one-subunit variant and a comparator law
# ABOUTME: (its L → 0 law or its enumeration parent's): Float64 LM fits with restarts, BigFloat polish.
# usage: contain.jl <data dir> <selection tsv: rxn gid comparator> <out tsv>
include(joinpath(@__DIR__, "plib.jl"))
include(joinpath(SP, "phase3", "census", "census_common.jl"))

const DATA, SEL, OUTF = ARGS[1], ARGS[2], ARGS[3]
const T0 = time()
const RECS = Dict{String, Dict{Int, Vector{String}}}()
function recs(rxn)
    get!(RECS, rxn) do
        rows = read_rows(joinpath(DATA, "recs_$rxn.tsv"))
        Dict(parse(Int, r[1]) => r for r in rows[2:end])
    end
end
roles(rxn) = startswith(rxn, "uniuni") ? ([:S], [:P]) : ([:A, :B], [:P, :Q])

"""A law to fit: its branches, the names it fits and the values it holds fixed."""
struct Target
    B::BLaw
    free::Vector{Symbol}
end

"""`log v` of `T.B` and its gradient over `log t` of the free names at each point; `z`
the free names' log bases, `fixed` the other values (bases by name), `X` the
concentration points."""
function logv(Tg::Target, z::Vector{F}, fixed::Dict{Symbol, F}, X, Q) where {F}
    B = Tg.B
    cols = colsof(B, Tg.free)
    t = Vector{F}(undef, length(B.vars.names))
    for (v, x) in enumerate(B.vars.names)
        B.vars.conc[v] || (t[v] = get(fixed, x, zero(F)))
    end
    for (j, v) in enumerate(cols)
        t[v] = exp(z[j])
    end
    y = Vector{F}(undef, length(X))
    J = Matrix{F}(undef, length(X), length(cols))
    for (r, xr) in enumerate(X)
        for (v, x) in enumerate(B.vars.names)
            B.vars.conc[v] && (t[v] = xr[x])
        end
        N, D, Ng, Dg = point(B, t, cols, Q)
        y[r] = log(abs(N / D))     # H sets the sign, shared by both laws
        J[r, :] .= Ng ./ N .- Dg ./ D
    end
    y, J
end

"""Levenberg–Marquardt on `Σ (log v - target)^2` from `z`; returns `(z, cost)`."""
function lm(Tg, z, fixed, X, Q, target; iters = 300)
    F = eltype(z)
    y, J = logv(Tg, z, fixed, X, Q)
    r = y .- target
    cost = sum(abs2, r)
    λ = F(1e-3)
    for _ in 1:iters
        A = J' * J
        g = J' * r
        δ = try
            -((A + λ * Diagonal(diag(A) .+ F(1e-30))) \ g)
        catch
            λ *= 10
            continue
        end
        z2 = z .+ δ
        ok = all(isfinite, z2) && all(abs.(z2) .< 300)
        if ok
            y2, J2 = try
                logv(Tg, z2, fixed, X, Q)
            catch
                (fill(F(NaN), length(X)), J)
            end
            r2 = y2 .- target
            c2 = sum(abs2, r2)
            if isfinite(c2) && c2 < cost
                z, y, J, r, cost = z2, y2, J2, r2, c2
                λ = max(λ / 3, F(1e-30))
                cost < F(1e-60) && break
                continue
            end
        end
        λ *= 4
        λ > F(1e20) && break
    end
    z, cost
end

"""Whether `f`'s law at random constants (`Lfix` values fixed) lies in `g`'s family: the
largest relative error of the best fit of `g` (`gfixed` held fixed), Float64 LM with 8
restarts and a 256-bit polish."""
function contained(f::Target, g::Target, gfixed_names::Vector{Symbol}, rng;
                   restarts::Int = 6, fixval = nothing)
    Qf, Qg = denlcm(f.B), denlcm(g.B)
    concs = [x for (v, x) in enumerate(g.B.vars.names) if g.B.vars.conc[v]]
    concf = [x for (v, x) in enumerate(f.B.vars.names) if f.B.vars.conc[v]]
    allc = unique!([concs; concf])
    npt = 3 * length(g.free) + 10
    X64 = [Dict(x => 10.0^(4 * rand(rng) - 2) for x in allc) for _ in 1:npt]
    fvals = Dict{Symbol, Float64}()
    for (v, x) in enumerate(f.B.vars.names)
        f.B.vars.conc[v] || (fvals[x] = 10.0^(2 * rand(rng) - 1))
    end
    shared = Dict{Symbol, Float64}(:Keq => fvals[:Keq])
    for x in gfixed_names
        shared[x] = fixval === nothing ? 10.0^(2 * rand(rng) - 1) : fixval
    end
    setprecision(BigFloat, 256) do
        XB = [Dict(x => big(v) for (x, v) in xr) for xr in X64]
        fbig = Dict(x => big(v)^(1 / big(Qf)) for (x, v) in fvals)
        tf, _ = logv(Target(f.B, Symbol[]), BigFloat[], fbig, XB, Qf)
        gfix64 = Dict(x => v^(1 / Qg) for (x, v) in shared)
        gfixB = Dict(x => big(v)^(1 / big(Qg)) for (x, v) in shared)
        best = (Inf, Float64[])
        for _ in 1:restarts
            z0 = [log(10.0^(2 * rand(rng) - 1)) / Qg for _ in g.free]
            z, c = lm(g, z0, gfix64, X64, Qg, Float64.(tf); iters = 200)
            c < best[1] && (best = (c, z))
        end
        zB, cB = lm(g, big.(best[2]), gfixB, XB, Qg, tf; iters = 60)
        Float64(sqrt(cB / length(XB)))
    end
end

function main()
    sel = [split(l, '\t') for l in eachline(SEL) if !isempty(l)]
    done = Set{String}()
    isfile(OUTF) && foreach(l -> push!(done, join(split(l, '\t')[1:3], '\t')), eachline(OUTF))
    open(OUTF, "a") do io
        for (rxn, gs, cmp) in sel
            id = join((rxn, gs, cmp), '\t')
            id in done && continue
            gid = parse(Int, gs)
            R = recs(rxn)
            r = R[gid]
            subs, prods = roles(rxn)
            rng = Xoshiro(gid + 7)
            t = time()
            row = try
                sm = parse_neutral(r[10], subs, prods)
                A = S.rate_law(sm)
                f = Target(blaw(A), filter(!=(:L), collect(S.fitted_params(sm))))
                g = if cmp == "LP"
                    pm = parse_neutral(r[10], subs, prods; plain = true)
                    Target(blaw(S.law(pm)), collect(S.fitted_params(pm)))
                else
                    pr = R[parse(Int, r[3])]
                    pm = parse_neutral(pr[10], subs, prods)
                    Target(blaw(S.rate_law(pm)), collect(S.fitted_params(pm)))
                end
                g_has_L = :L in g.free
                # child in comparator: the variant at random (θ, L), the comparator fitted
                ab = contained(f, g, Symbol[], rng)
                # comparator in child: L fixed at a random value, then L free
                fL = Target(f.B, f.free)
                ba_fix = contained(g, fL, [:L], rng)
                ba_free = contained(g, Target(f.B, [f.free; :L]), Symbol[], rng)
                [rxn, gs, cmp, @sprintf("%.1e", ab), @sprintf("%.1e", ba_fix),
                 @sprintf("%.1e", ba_free), string(Int(g_has_L)),
                 @sprintf("%.1f", time() - t)]
            catch e
                [rxn, gs, cmp, "error: " * replace(first(sprint(showerror, e), 150), r"[\t\n]" => " ")]
            end
            println(io, join(row, '\t'))
            flush(io)
            time() - T0 > 520 && (println("stopping for time"); break)
        end
    end
    @printf("finished in %.0f s, maxrss %.0f MB\n", time() - T0, Sys.maxrss() / 2^20)
end

basename(PROGRAM_FILE) == "contain.jl" && main()
