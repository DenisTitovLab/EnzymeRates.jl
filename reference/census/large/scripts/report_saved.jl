# ABOUTME: Prints the census tables and check costs of saved runs; usage:
# ABOUTME: report_saved.jl <rxn> [<rxn> ...] (reads data/<rxn>.jls).
include(joinpath(@__DIR__, "census_lib.jl"))
include(joinpath(@__DIR__, "summarize.jl"))
for name in ARGS
    data = deserialize(joinpath(@__DIR__, "..", "data", "$(name).jls"))
    print_tables(data.info, data.rounds)
    print_costs(data.rounds)
end
