# ABOUTME: Prints the slowest pipeline rows of a saved run with their stage times; usage:
# ABOUTME: peek_rows.jl label [n]
include(joinpath(@__DIR__, "lib.jl"))
rows = filter(x -> x.ok, deserialize(joinpath(IDATA, "pipe_$(ARGS[1]).jls")))
n = length(ARGS) > 1 ? parse(Int, ARGS[2]) : 12
for x in first(sort(rows; by = x -> -x.t_derive), n)
    @printf("%-14s %-18s derive %7.1f ms: solve %.2f law %.2f key %.2f print %.2f code %.2f kcat %.2f; nD %d trees %s\n",
            x.set, x.id, 1e3x.t_derive, 1e3x.t_solve, 1e3x.t_law, 1e3x.t_key, 1e3x.t_print,
            1e3x.t_code, 1e3x.t_kcat, x.nD, x.trees)
end
