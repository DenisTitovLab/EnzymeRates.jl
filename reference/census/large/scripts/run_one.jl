# ABOUTME: Runs the census of one reaction from the command line and prints its tables;
# ABOUTME: usage: run_one.jl <rxn> <N0> <N1> <N2> <fullmemo:0|1> <pkgrank_frac>.
include(joinpath(@__DIR__, "driver.jl"))
include(joinpath(@__DIR__, "summarize.jl"))
rxn = Symbol(ARGS[1])
N0, N1, N2 = parse.(Int, ARGS[2:4])
fullmemo = ARGS[5] == "1"
frac = parse(Float64, ARGS[6])
out = joinpath(@__DIR__, "..", "data", "$(rxn).jls")
info, rounds = run_census(rxn; N0, N1, N2, fullmemo, pkgrank_frac = frac, out)
print_tables(info, rounds)
print_costs(rounds)
