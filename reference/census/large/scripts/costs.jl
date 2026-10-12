# ABOUTME: Warm per-check timings on a sample of saved fitted candidates (each check run three
# ABOUTME: times, minimum kept), and the seed screen's degenerate seeds split into V and C.

include(joinpath(@__DIR__, "census_lib.jl"))
include(joinpath(@__DIR__, "summarize.jl"))

"""The minimum time in µs of three runs of `f`, after one warm-up run."""
function tmin(f)
    f()
    minimum(@elapsed(f()) for _ in 1:3) * 1e6
end

"""The support tests alone: `N ≡ 0` and the faces where `D` vanishes, in lowest terms."""
function support_only(m, vars, N, D)
    sh = S.conc_shift((N, D), vars)
    N2, D2 = N * sh, D * sh
    faces = [[vars.index[x] for x in zs] for (_, zs) in face_sets(m)]
    survives(P, Z) = any(mo -> all(z -> S.mexp(mo, z) == 0, Z), keys(P.terms))
    iszero(N2), count(Z -> !survives(D2, Z), faces)
end

"""The joint saturation test alone (no forward or reverse Vmax)."""
function saturation_only(m, vars, N, D)
    sh = S.conc_shift((N, D), vars)
    N2, D2 = N * sh, D * sh
    si = [vars.index[x] for x in m.substrates]
    pi_ = [vars.index[x] for x in m.products]
    function top(P, zero, grow)
        best, found = -1 // 1, false
        for mo in keys(P.terms)
            all(z -> S.mexp(mo, z) == 0, zero) || continue
            d = sum((S.mexp(mo, x) for x in grow); init = 0 // 1)
            best = found ? max(best, d) : d
            found = true
        end
        found, best
    end
    unb(g, z) = (a = top(N2, z, g); b = top(D2, z, g); a[1] && b[1] && a[2] > b[2])
    unb(si, pi_), unb(pi_, si)
end

rng = Xoshiro(17)
const NPICK = parse(Int, ARGS[1])
for name in Symbol.(ARGS[2:end])
    data = deserialize(joinpath(@__DIR__, "..", "data", "$(name).jls"))
    fit = [r for recs in data.rounds for r in recs if r[:fitted] && get(r, :rank, -1) >= 0]
    pick = fit[sort!(randperm(rng, length(fit))[1:min(NPICK, length(fit))])]
    T = Dict{String, Vector{Float64}}()
    add(k, v) = push!(get!(T, k, Float64[]), v)
    for r in pick
        mp, m = r[:mech], r[:smech]
        add("package: _eq_complexity", tmin(() -> ER._eq_complexity(mp)))
        em = ER.compile_mechanism(mp)
        add("package: fitted count (solve)", tmin(() -> ER._dependent_param_exprs(typeof(em))))
        add("package: rate_equation_string", tmin(() -> ER.rate_equation_string(em)))
        add("convert to slice", tmin(() -> to_slice(mp)))
        add("check_boundary (bottomless P=0/S=0, role, names)", tmin(() -> try_msg(() -> S.check_boundary(m))))
        add("admit (:OnlyA island rule)", tmin(() -> S.admit(m)))
        add("RE turnover", tmin(() -> re_turnover(m)))
        imask = S.is_allosteric(m) ? S.inactive_mask(m) : falses(length(m.edges))
        add("face check, all $(length(face_sets(m))) faces", tmin(() -> for (_, zs) in face_sets(m)
            face_check(m, zs); any(imask) && face_check(m, zs, imask); end))
        add("edge flux (structural)", tmin(() -> S.edge_flux(m)))
        local sol, L
        add("slice solve + law", tmin(() -> (sol = S.solve(m); L = S.rate_law(m, sol))))
        add("rename-proof key", tmin(() -> S.canonical(L)))
        local vars, N, D, Nb
        add("expand law (N, D)", tmin(() -> ((vars, N, D, Nb) = expand_law_checked(L))))
        add("support tests (N≡0, D=0 on faces)", tmin(() -> support_only(m, vars, N, D)))
        add("saturation (no Vmax)", tmin(() -> saturation_only(m, vars, N, D)))
        fn = slice_fitted(m, sol)
        add("rank (Jacobian + SVD)", tmin(() -> law_jacobian(vars, N, D, fn, Xoshiro(1))))
        add("package-law rank (re-derive + complex step)", tmin(() -> pkg_rank(mp, Xoshiro(1))))
    end
    println("$name: warm cost per candidate in µs (median / p90 / max) on $(length(pick)) fitted")
    for (k, v) in sort!(collect(T); by = x -> median(x[2]))
        @printf("  %-52s %10.1f %10.1f %10.1f\n", k, median(v), quantile(v, 0.9), maximum(v))
    end
    flush(stdout)
end
