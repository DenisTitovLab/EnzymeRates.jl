# ABOUTME: Item 8: wall time of a second-process check run as test_compile_budget.jl runs its
# ABOUTME: subprocesses: package load alone, and keys plus generated bodies of the 42 specs and 6 gate laws.
using Printf
const K = @__DIR__
julia = Base.julia_cmd().exec[1]
load = `$julia --startup-file=no --project=$K/env -e 'using EnzymeRates'`
work = `$julia --startup-file=no -e "include(\"$K/klib.jl\"); rows = tsv(\"$K/data/x5_small.tsv\")[1:48]; for r in rows; m = parse_neutral(r[end], names_of(r[4]), names_of(r[5])); sol = S.solve(m); L = S.rate_law(m, sol); println(canonical_named(L)[1].key[1:16], \" \", bytes2hex(sha256(string(S.body(S.program(m, sol, L)))))[1:16]); end"`
for (name, cmd) in (("load EnzymeRates", load), ("load EnzymeRates", load),
                    ("keys and bodies of 48 laws (slice)", work))
    out = IOBuffer()
    t = @elapsed run(pipeline(cmd; stdout = out))
    lines = split(String(take!(out)), '\n'; keepempty = false)
    @printf("%-36s %.1f s, %d output lines\n", name, t, length(lines))
    name == "keys and bodies of 48 laws (slice)" && open(joinpath(K, "data", "subproc_keys.txt"), "w") do io
        foreach(l -> println(io, l), lines)
    end
end
