# ABOUTME: Targeted re-fit diagnostic for issue (2) within-hash loss spread.
# ABOUTME: Tests whether CMA-ES maxtime=60 was the root cause; refits worst buckets at maxtime=600.

using Pkg
Pkg.activate(temp=true)
Pkg.develop(path=joinpath(@__DIR__, ".."))
Pkg.add(["CSV", "DataFrames", "Random", "OptimizationPyCMA"])

using EnzymeRates, CSV, DataFrames, Random, OptimizationPyCMA
using EnzymeRates: _canonical_rate_eq_hash, FittingProblem, fit_rate_equation

# --------------------------------------------------------------------------
# Configuration
# --------------------------------------------------------------------------
const SPREAD_THRESHOLD = 0.001   # 0.1% — fitter-noise level
const TOP_K = 5                   # worst-spread buckets per n
const N_RESTARTS = 10
const MAXTIME = 600.0
const LDH_KEQ = 20000.0
const RESULT_CSV = joinpath(@__DIR__, "refit_results.csv")

# --------------------------------------------------------------------------
# Data loading
# --------------------------------------------------------------------------
function load_ldh_data()
    raw = CSV.read(joinpath(@__DIR__, "LDH_data.csv"), DataFrame)
    raw.group = string.(raw.Article, "/", raw.Fig)
    select(raw, [:group, :Rate, :Pyruvate, :NADH, :Lactate, :NAD])
end

# --------------------------------------------------------------------------
# Bad-bucket identification
# --------------------------------------------------------------------------
function find_bad_buckets(n::Int)
    csv_path = joinpath(@__DIR__, "params_estimate_$(n).csv")
    df = CSV.read(csv_path, DataFrame)
    df = filter(:loss => l -> l < 1.0, df)   # drop sentinel losses

    mechs = [eval(Meta.parse(row.mechanism_type))() for row in eachrow(df)]
    df.canonical_hash = [string(_canonical_rate_eq_hash(m), base=16, pad=16) for m in mechs]
    df.mech = mechs

    bad = NamedTuple[]
    for g in groupby(df, :canonical_hash)
        nrow(g) < 2 && continue
        lo, hi = extrema(g.loss)
        spread = (hi - lo) / lo
        spread > SPREAD_THRESHOLD || continue
        push!(bad, (hash=g.canonical_hash[1], spread=spread, lo=lo, hi=hi, df=DataFrame(g)))
    end
    sort!(bad; by=b -> b.spread, rev=true)
    return first(bad, min(TOP_K, length(bad)))
end

# --------------------------------------------------------------------------
# Incremental CSV save / resume
# --------------------------------------------------------------------------
function init_or_resume_results()
    if isfile(RESULT_CSV)
        return CSV.read(RESULT_CSV, DataFrame)
    end
    DataFrame(
        n_params=Int[], hash=String[], mech_type=String[],
        old_loss=Float64[], new_loss=Float64[],
        fit_time_s=Float64[], status=String[],
    )
end

function append_result!(results_df, row)
    push!(results_df, row)
    CSV.write(RESULT_CSV, results_df)
end

function already_fit(results_df, n::Int, hash::String, mech_type::String)
    isempty(results_df) && return false
    any(
        (results_df.n_params .== n) .&
        (results_df.hash .== hash) .&
        (results_df.mech_type .== mech_type),
    )
end

# --------------------------------------------------------------------------
# Fitting
# --------------------------------------------------------------------------
function refit_mech(m, data, _old_loss)
    fp = FittingProblem(m, data; Keq=LDH_KEQ)
    t = @elapsed fit = fit_rate_equation(fp, PyCMAOpt(); n_restarts=N_RESTARTS, maxtime=MAXTIME)
    return fit.loss, t
end

# --------------------------------------------------------------------------
# Summary printer
# --------------------------------------------------------------------------
function print_summary(results)
    println("\n############ SUMMARY ############")
    for n in (5, 6, 7, 8)
        sub = filter(r -> r.n_params == n, results)
        isempty(sub) && continue
        println("\nn=$n:")
        for hg in groupby(sub, :hash)
            okrows = filter(r -> r.status == "ok", hg)
            if nrow(okrows) >= 2
                old_lo, old_hi = extrema(okrows.old_loss)
                new_lo, new_hi = extrema(okrows.new_loss)
                old_spread = (old_hi - old_lo) / old_lo
                new_spread = (new_hi - new_lo) / new_lo
                verdict = if new_spread < SPREAD_THRESHOLD
                    "COLLAPSED"
                elseif new_spread < old_spread / 2
                    "MOSTLY"
                elseif new_spread < old_spread
                    "PARTIAL"
                else
                    "PERSISTS"
                end
                println(
                    "  $(hg.hash[1]): " *
                    "old=$(round(old_spread*100, digits=2))% → " *
                    "new=$(round(new_spread*100, digits=2))% " *
                    "[$verdict] (n=$(nrow(okrows)))",
                )
            else
                println(
                    "  $(hg.hash[1]): insufficient finite refits " *
                    "($(nrow(okrows))/$(nrow(hg)))",
                )
            end
        end
    end
end

# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------
function main()
    data = load_ldh_data()
    println("LDH data: $(nrow(data)) rows, $(length(unique(data.group))) groups")

    results = init_or_resume_results()
    println("Resuming with $(nrow(results)) prior results")

    for n in (5, 6, 7, 8)
        println("\n############ n=$n ############")
        bad = find_bad_buckets(n)
        println("Top $(length(bad)) bad buckets at n=$n:")
        for (i, b) in enumerate(bad)
            println(
                "  $i: hash=$(b.hash) " *
                "spread=$(round(b.spread*100, digits=2))% " *
                "mechs=$(nrow(b.df))",
            )
        end

        for (bi, b) in enumerate(bad)
            println("\n--- bucket $bi/$(length(bad)) hash=$(b.hash) ---")
            for (mi, row) in enumerate(eachrow(b.df))
                mt = row.mechanism_type
                if already_fit(results, n, b.hash, mt)
                    println("  mech $mi/$(nrow(b.df)): SKIP (already fit)")
                    continue
                end
                println(
                    "  mech $mi/$(nrow(b.df)): " *
                    "old_loss=$(round(row.loss, sigdigits=5)) refitting...",
                )
                t0 = time()
                status, new_loss, fit_time = try
                    nl, ft = refit_mech(row.mech, data, row.loss)
                    ("ok", nl, ft)
                catch e
                    msg = sprint(showerror, e)
                    ("error: $(msg[1:min(80, length(msg))])", NaN, time() - t0)
                end
                println(
                    "    -> new_loss=$(round(new_loss, sigdigits=5)) " *
                    "status=$(status) ($(round(fit_time, digits=1))s)",
                )
                append_result!(results, (
                    n_params=n, hash=b.hash,
                    mech_type=mt, old_loss=row.loss,
                    new_loss=new_loss, fit_time_s=fit_time,
                    status=status,
                ))
            end
        end
    end

    print_summary(results)
end

main()
