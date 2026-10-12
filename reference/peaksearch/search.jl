# ABOUTME: Runs the adversarial peak search over every mechanism variant, product setting and Keq,
# ABOUTME: then re-verifies each best point in BigFloat and by an inner sup over concentrations.

include(joinpath(@__DIR__, "search_lib.jl"))
using Serialization

const MODE = isempty(ARGS) ? "pilot" : ARGS[1]
const BUDGET = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 20_000
const KEYS = length(ARGS) >= 3 ? split(ARGS[3], ",") : ["T", "g", "a", "b", "c", "d1", "d2", "f"]
const SETTINGS = length(ARGS) >= 4 ? Symbol.(split(ARGS[4], ",")) : [:zero, :Q, :P, :PQ]
const LKS = length(ARGS) >= 5 ? parse.(Float64, split(ARGS[5], ",")) : [-3.0, 0.0, 3.0]
const OUT = joinpath(@__DIR__, "results_$(MODE).jls")

setprecision(BigFloat, 256)

function verify(m, setting, lK, x)
    nθ = ntheta(m)
    rb = evalR(m, big.(x), setting, lK; T = BigFloat)
    kf, kr, _ = rates(m, big.(x[1:nθ]), lK; T = BigFloat)
    # the saturating limit, pushed further out
    vs80 = vsat_dirs(m, kf, kr; dirs = ((80, 80), (80, 40), (40, 80)))
    Rin, lcin = inner_sup(m, x[1:nθ], lK, setting)
    (Rbig = rb.R, vsdirs = Float64.(rb.vsdirs), vs80 = Float64.(vs80), Rinner = Rin,
     lcinner = lcin)
end

function main()
    m = MECHS["T"]
    x = vcat(zeros(ntheta(m)), [1.0, 1.0])
    evalR(m, x, :zero, 0.0)
    t = @elapsed for _ in 1:20_000
        evalR(m, x, :zero, 0.0)
    end
    @printf("timing: %.2f µs per objective evaluation (textbook, 4 GTH solves)\n", t / 20_000 * 1e6)
    results = Dict{Any, Any}()
    for key in KEYS, setting in SETTINGS, lK in LKS
        m = MECHS[key]
        t = @elapsed res = search(m, setting, lK; seed = hash((key, setting, lK)) % UInt32,
                                  budget = BUDGET)
        nθ = ntheta(m)
        if isempty(res.x)
            @printf("%-3s %-5s log10Keq=%5.1f  no feasible point found\n", key, setting, lK)
            continue
        end
        lc = res.x[nθ+1:end]
        ver = verify(m, setting, lK, res.x)
        results[(key, setting, lK)] = (; res..., ver...)
        @printf("%-3s %-5s log10Keq=%5.1f  R=%.12f  Rbig-1=%+.3e  Rinner=%.12f  at log10(A,B[,P][,Q])=%s  evals=%d  %.1fs\n",
                key, setting, lK, res.R, Float64(ver.Rbig - 1), ver.Rinner,
                string(round.(lc; digits = 3)), res.evals, t)
        ds = ver.vsdirs
        spread = (maximum(ds) - minimum(ds)) / maximum(ds)
        @printf("      v_sat: diag %.6e, A>>B %.6e, B>>A %.6e (rel spread %.1e; at 1e80: %.6e); R vs diagonal v_sat = %.6g\n",
                ds[1], ds[2], ds[3], spread, maximum(ver.vs80), res.R * maximum(ds) / ds[1])
        flush(stdout)
        serialize(OUT, results)
    end
    # Box restricted to A, B <= 1e2 so the optimizer cannot sit on the saturating corner
    if MODE == "full"
        println("\n-- zero products, A, B restricted to 1e-8..1e2 (no saturating corner) --")
        for key in ("T", "c", "b", "a", "g")
            m = MECHS[key]
            res = search(m, :zero, 0.0; seed = 99, budget = BUDGET, chi = 2.0)
            lc = res.x[ntheta(m)+1:end]
            rb = evalR(m, big.(res.x), :zero, 0.0; T = BigFloat, chi = 2.0)
            @printf("%-3s zero  restricted  R=%.12f  Rbig-1=%+.3e at log10(A,B)=%s\n",
                    key, res.R, Float64(rb.R - 1), string(round.(lc; digits = 3)))
            results[(key, :zero_restricted, 0.0)] = res
            flush(stdout)
        end
        serialize(OUT, results)
    end
    println("\n-- best points --")
    for (k, r) in sort(collect(results); by = p -> string(p[1]))
        k[2] === :zero_restricted && continue
        m = MECHS[k[1]]
        @printf("%s  R=%.6g\n", string(k), r.R)
        print(describe(m, r.x[1:ntheta(m)], k[3]))
    end
end

main()
