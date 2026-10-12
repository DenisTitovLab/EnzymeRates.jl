# ABOUTME: Fit-cost calibration in one Julia process: the candidate selection (select.jl), then the
# ABOUTME: timed fits (fit.jl). ARGS: fitting budget in seconds, maxtime per restart in seconds.
include("/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/census/fitcost/select.jl")
let t = time()
    main_select()
    @printf("selection phase %.0f s\n", time() - t)
end
GC.gc()
include(joinpath(FCD, "fit.jl"))
main_fit(; budget = length(ARGS) >= 1 ? parse(Float64, ARGS[1]) : 3600.0,
         maxtime = length(ARGS) >= 2 ? parse(Float64, ARGS[2]) : 600.0,
         fit_cap = length(ARGS) >= 3 ? parse(Float64, ARGS[3]) : 900.0)
