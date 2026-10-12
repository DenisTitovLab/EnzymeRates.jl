# ABOUTME: Class P by independent ranks: BigFloat central differences of the slice law at new
# ABOUTME: points, of the mass-action oracle, and of the package's own rate_equation.
include(joinpath(@__DIR__, "vlib.jl"))

const RX = (:ter_bi_rxn, :ter_ter_rxn, :pyruvate_carboxylase_rxn)

"""The pivoted Gram–Schmidt norms of `dv/dlog θ` (rows scaled to unit max) of `f(θ, conc)`
over the names `names` at one random point `θ0` (`10^[-1.5,1.5]`), `2n+6` random
concentration points, central differences of step `h` in log θ."""
function fd_norms(f, names, mkconc, rng; h)
    np = length(names)
    θ0 = Dict{Symbol, BigFloat}(n => lu10(rng, -1.5, 1.5) for n in names)
    J = zeros(BigFloat, 2np + 6, np)
    for r in axes(J, 1)
        conc = mkconc()
        for (j, n) in enumerate(names)
            θp, θm = copy(θ0), copy(θ0)
            θp[n] *= exp(h)
            θm[n] *= exp(-h)
            J[r, j] = (f(θp, conc) - f(θm, conc)) / (2h)
        end
        mx = maximum(abs, view(J, r, :))
        mx > 0 && (J[r, :] ./= mx)
    end
    pivoted_norms(J)
end

rank_of(norms, tol) = count(>(tol), norms)

"""The slice law's rank at a new point (512 bits, h = 1e-40, tolerance 1e-30)."""
function slice_rank(m, rng)
    L = S.rate_law(m)
    names = collect(S.fitted_params(m))
    keq = lu10(rng, -1, 1)
    f(θ, conc) = law_value(L, lawx(L, merge(θ, Dict(:Keq => keq)), conc))
    nr = fd_norms(f, names, () -> random_conc(m, rng; lo = -1.5, hi = 1.5), rng;
                  h = big(10.0)^-40)
    rank_of(nr, 1e-30), length(names), nr
end

"""The oracle's rank (pruned n=1 network, Λ = 1e30, h = 1e-22, tolerance 1e-18)."""
function oracle_rank(m, rng)
    names = collect(S.fitted_params(m))
    keq = lu10(rng, -1, 1)
    sp = random_speed(m, rng)
    function f(θ, conc)
        p = (; (n => θ[n] for n in names)..., Keq = keq, E_total = big(1.0))
        vals = vals_from(m, p)
        cols, x, _, _ = all_columns(m, vals, rng)
        first(oracle(m, cols, x, vals, conc; Λ = big(10.0)^30, speed = sp))
    end
    nr = fd_norms(f, names, () -> random_conc(m, rng; lo = -1.5, hi = 1.5), rng;
                  h = big(10.0)^-22)
    rank_of(nr, 1e-18), length(names), nr
end

"""The package law's rank: `rate_equation` of the compiled mechanism in BigFloat, its
own fitted parameters (512 bits, h = 1e-40, tolerance 1e-30)."""
function package_rank(mp, rng)
    em = ER.compile_mechanism(mp)
    names = collect(ER.fitted_params(em))
    mets = collect(ER.metabolites(em))
    keq = lu10(rng, -1, 1)
    function f(θ, conc)
        p = (; (n => θ[n] for n in names)..., Keq = keq, E_total = big(1.0))
        c = (; (x => conc[x] for x in mets)...)
        Base.invokelatest(ER.rate_equation, em, c, p)
    end
    nr = fd_norms(f, names, () -> Dict{Symbol, BigFloat}(x => lu10(rng, -1.5, 1.5) for x in mets),
                  rng; h = big(10.0)^-40)
    rank_of(nr, 1e-30), length(names), nr
end

short(nr, rk) = @sprintf("kept min %.1e, dropped max %.1e", Float64(nr[rk]),
                         rk < length(nr) ? Float64(nr[rk + 1]) : 0.0)

rng = Xoshiro(20261011 + 3)
tally = Dict{Tuple{String, Bool, Bool}, Int}()   # (kind, census P, my slice-law P)
for name in RX
    data = load_run(name)
    for (ri, recs) in enumerate(data.rounds)
        fit = [r for r in recs if r[:fitted] && get(r, :rank, -1) >= 0]
        P = [r for r in fit if :P in r[:labels]]
        notP = [r for r in fit if !(:P in r[:labels])]
        Ponly = [r for r in fit if :Ptoday in r[:labels] && !(:P in r[:labels])]
        border = [r for r in fit if get(r, :sv_kept, 1.0) < 1e-6]
        bykind = Dict{String, Vector{Any}}()
        for r in P
            k = join(sort!(string.(collect(keys(get(r, :explain, Dict()))))), "+")
            push!(get!(bykind, k, Any[]), r)
        end
        pickP = Any[]
        for (k, rs) in bykind
            append!(pickP, rs[sort!(randperm(rng, length(rs))[1:min(4, length(rs))])])
        end
        pickP = pickP[1:min(10, end)]
        pickN = notP[sort!(randperm(rng, length(notP))[1:min(10, length(notP))])]
        pickT = Ponly[sort!(randperm(rng, length(Ponly))[1:min(4, length(Ponly))])]
        println("-- $name round $(ri - 1): P $(length(P)) $(Dict(k => length(v) for (k, v) in bykind)), not P $(length(notP)), package-only $(length(Ponly)), borderline $(length(border))")
        for (lab, rs) in (("P", pickP), ("notP", pickN), ("pkg-only", pickT), ("borderline", border))
            for (i, r) in enumerate(rs)
                m = r[:smech]
                rk1, n1, nr1 = slice_rank(m, rng)
                rk2, _, nr2 = slice_rank(m, rng)
                line = @sprintf("   %-10s %s census slice %d/%d pkg %d/%d | new points %d/%d (%s), %d/%d",
                                lab, r[:id], r[:rank], r[:n_slice], r[:pkg_rank_used],
                                r[:n_params], rk1, n1, short(nr1, rk1), rk2, n1)
                if lab in ("P", "notP") && i <= 2 || lab == "pkg-only" && i <= 2
                    rko, _, nro = oracle_rank(m, rng)
                    line *= @sprintf(" | oracle %d/%d (%s)", rko, n1, short(nro, rko))
                end
                if lab == "pkg-only" || (lab == "P" && i <= 3)
                    rkp, npk, nrp = try
                        package_rank(r[:mech], rng)
                    catch e
                        -1, -1, BigFloat[]
                    end
                    line *= rkp < 0 ? " | package error" :
                            @sprintf(" | package %d/%d (%s)", rkp, npk, short(nrp, rkp))
                end
                println(line, "\n      ", string(get(r, :explain, "")))
                k = (lab, :P in r[:labels], rk1 < n1 || rk2 < n1)
                tally[k] = get(tally, k, 0) + 1
            end
            flush(stdout)
        end
    end
end
println("tally (sample, census P, new slice-law rank deficient) => count: ", sort!(collect(tally)))
