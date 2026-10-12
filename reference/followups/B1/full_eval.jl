# ABOUTME: Follow-up B1: re-evaluates P5a's :OnlyA cone limit over ALL enumerated directions (P5a used 16),
# ABOUTME: for (ii), (iu), untied R_D with free ratios (iiu) and with bridge kinetic ratios kept (Dk).
# Usage: julia full_eval.jl <population> [first last]   → data/full_<population>[_first].jls
include(joinpath(@__DIR__, "b1_lib.jl"))
include(joinpath(S, "proto/p5a/limits_lib.jl"))

"""P5a's `MWC.cone_directions` with `extra` equality rows (pairs of column indices that must take equal
exponents); same grid and interior sampler, same seeds."""
function cone_dirs(J, conv::Symbol; extra = Tuple{Int, Int}[], R = 3, maxgrid = 100_000)
    sig = MWC.cone_signs(J, conv)
    O = sort!(collect(keys(sig)))
    n = length(J.cols)
    isempty(O) && return [zeros(MWC.Q, n)], 0
    pos = Dict(i => t for (t, i) in enumerate(O))
    rows = [J.C[r, O] for r in axes(J.C, 1) if any(!iszero, J.C[r, O])]
    for (t, i) in enumerate(O)
        sig[i] == 0 && push!(rows, [u == t ? 1 : 0 for u in eachindex(O)])
    end
    for (i1, i2) in extra
        row = zeros(Int, length(O)); row[pos[i1]] += 1; row[pos[i2]] -= 1
        push!(rows, row)
    end
    A = isempty(rows) ? zeros(Int, 0, length(O)) : permutedims(reduce(hcat, rows))
    B = MWC.nullbasis(A)
    d = length(B)
    d == 0 && return Vector{MWC.Q}[], 0
    (2R + 1)^d > maxgrid && (R = max(1, floor(Int, (maxgrid^(1 / d) - 1) / 2)))
    found = Set{Vector{Int}}()
    okdir(aO) = all(t -> (s = sig[O[t]]; s == 0 ? aO[t] == 0 : sign(aO[t]) == s), eachindex(O))
    prim(aO) = (den = foldl(lcm, (denominator(x) for x in aO); init = 1);
                v = [Int(x * den) for x in aO]; v .÷ foldl(gcd, v; init = 0))
    for y in Iterators.product(ntuple(_ -> -R:R, d)...)
        (2R + 1)^d > 4maxgrid && break
        all(iszero, y) && continue
        aO = sum(y[k] * B[k] for k in 1:d)
        okdir(aO) && push!(found, prim(aO))
    end
    N = Float64.(reduce(hcat, B)); G = N' * N
    rng = Random.Xoshiro(17)
    box(x) = [(s = sig[O[t]]; s == 0 ? 0.0 : s > 0 ? max(x[t], 1.0) : min(x[t], -1.0)) for t in eachindex(O)]
    for _ in 1:12
        x = [(s = sig[O[t]]; s * (1 + 4rand(rng))) for t in eachindex(O)]
        y = zeros(d)
        for _ in 1:400; y = G \ (N' * box(x)); x = N * y; end
        for K in (1, 2, 4, 8, 16, 64)
            aO = sum(MWC.Q(round(Int, K * y[k]), K) * B[k] for k in 1:d)
            okdir(aO) || continue
            push!(found, prim(aO)); break
        end
    end
    dirs = sort!(collect(found); by = v -> (sum(abs, v), v))
    [(a = zeros(MWC.Q, n); for (t, i) in enumerate(O); a[i] = v[t]; end; a) for v in dirs], d
end

"""Kinetic doors of `m` under the untied R_D typing that are bridges of the full step graph."""
function bridge_kinetic_doors(m)
    a = analyze(m, :D)
    [j for j in a.doors if !a.aff[j] && j in a.bridges]
end

"""Full-direction verdict for one mechanism and convention (:ii, :iu, :iiu, :Dk)."""
function full_eval(m, conv; maxdirs = parse(Int, get(ENV, "B1_MAXDIRS", "3000")))
    rng = Random.Xoshiro(11)                                 # P5a light(): same draws
    untie = conv in (:iu, :iiu, :Dk)
    c = conv === :iu ? :i : :ii
    J = MWC.joint_solve(m; order = :pkg, untie)
    hat = MWC.random_hat(J, rng); keq, L = logu(rng), logu(rng)
    rA, rI = draw_reg(m, rng)
    concs = Dict(x => logu(rng) for x in m.metabolites)
    extra = Tuple{Int, Int}[]
    if conv === :Dk
        for j in bridge_kinetic_doors(m)
            J.cs[:I][j].re && error("bridge kinetic door $j is RE in I")
            push!(extra, (J.idx[(:I, j, 1)], J.idx[(:I, j, 2)]))
        end
    end
    r = Dict{Symbol, Any}(:id => m.id, :conv => conv, :nextra => length(extra))
    t = @elapsed begin
        dirs, d = cone_dirs(J, c; extra)
        r[:cone_dim], r[:cone_n] = d, length(dirs)
        r[:v_del] = v_deletion(J, hat, keq, L, rA, rI, concs)
        if !isempty(dirs)
            Sx = sym_setup(J, hat, keq, L, rA, rI, concs)
            res = [sym_eval(Sx, a) for a in dirs[1:min(end, maxdirs)]]
            vs = [s.v for s in res]
            r[:cone_v] = vs[1]
            r[:cone_spread] = maximum(reldiff(v, vs[1], 0) for v in vs)
            r[:cone_seq] = any(s -> s.ordI[2] !== nothing && s.ordI[3] !== nothing && s.ordI[2] < s.ordI[3], res)
            r[:cone_vs_del] = maximum(reldiff(v, r[:v_del], 0) for v in vs)
            r[:n_eval] = length(vs)
            r[:n_distinct] = length(unique(round.(Float64.(vs); sigdigits = 12)))
            r[:any_zero] = any(iszero, vs)
            # P5a's 16-direction pick, for comparison
            pick = unique([1; round.(Int, range(1, length(dirs); length = min(16, length(dirs))))])
            pick = filter(<=(length(vs)), pick)
            r[:spread16] = maximum(reldiff(vs[k], vs[1], 0) for k in pick)
        end
    end
    r[:time] = t
    r
end

if abspath(PROGRAM_FILE) == @__FILE__
    name = ARGS[1]
    pop = filter(m -> :OnlyA in m.group_tags, load_allo(name))
    length(ARGS) > 1 && (pop = pop[parse(Int, ARGS[2]):min(end, parse(Int, ARGS[3]))])
    if length(ARGS) > 3
        want = Set(ARGS[4:end]); pop = filter(m -> m.id in want, pop)
    end
    recs = Dict{Symbol, Any}[]
    t0 = time()
    for (k, m) in enumerate(pop), conv in (:ii, :iu, :iiu, :Dk)
        r = try
            full_eval(m, conv)
        catch e
            Dict{Symbol, Any}(:id => m.id, :conv => conv, :error => first(sprint(showerror, e), 200))
        end
        push!(recs, r)
        conv === :Dk && k % 50 == 0 && (@printf("%s %d/%d %.0fs\n", name, k, length(pop), time() - t0); flush(stdout))
    end
    mkpath(joinpath(@__DIR__, "data"))
    tag = length(ARGS) > 3 ? "_ids" : length(ARGS) > 1 ? "_$(ARGS[2])" : ""
    serialize(joinpath(@__DIR__, "data", "full_$(name)$(tag).jls"), recs)
    for conv in (:ii, :iu, :iiu, :Dk)
        rs = filter(r -> r[:conv] === conv, recs)
        vs = [verdict(r) for r in rs]
        cnt = Dict{String, Int}(); for v in vs; cnt[v] = get(cnt, v, 0) + 1; end
        @printf("%s (%s): %s | max dirs evaluated %d, mean %.0f; time %.1f s\n", name, conv,
                join(("$k=$v" for (k, v) in sort!(collect(cnt); by = first)), " "),
                maximum((get(r, :n_eval, 0) for r in rs); init = 0), sum(get(r, :n_eval, 0) for r in rs) / max(1, length(rs)),
                sum(get(r, :time, 0.0) for r in rs))
        for r in filter(r -> haskey(r, :error), rs)[1:min(end, 3)]; println("   error ", r[:id], ": ", r[:error]); end
    end
    @printf("total %.0f s\n", time() - t0)
end
