# ABOUTME: Unit-rule census driver: runs `analyze` on mechanisms i0:i1 of one population until a
# ABOUTME: time budget runs out, writing data/<pop>_<i0>_<last>.jls; args: <pop> <i0> <i1> [budget_s].
include(joinpath(@__DIR__, "ulib.jl"))
include(joinpath(@__DIR__, "pops.jl"))

pop, i0, i1 = ARGS[1], parse(Int, ARGS[2]), parse(Int, ARGS[3])
budget = length(ARGS) >= 4 ? parse(Float64, ARGS[4]) : Inf
P = population(pop)
i1 = min(i1, length(P))
@printf("%s: %d mechanisms, running %d:%d within %.0f s\n", pop, length(P), i0, i1, budget)
rng = Xoshiro(20261014 + i0)
res = []
t0 = time()
stop = i0 - 1
for i in i0:i1
    time() - t0 > budget && break
    global stop = i
    id, text, subs, prods = P[i]
    t = time()
    row = try
        m = parse_neutral(text, subs, prods)
        merge((; id, i, ok = true, err = ""), analyze(m, rng))
    catch err
        (; id, i, ok = false, err = first(sprint(showerror, err), 300))
    end
    push!(res, merge(row, (; secs = time() - t)))
    i % 200 == 0 && @printf("  %d done, %.0f s, maxrss %.0f MB\n", i, time() - t0, Sys.maxrss() / 2^20)
end
out = joinpath(R5, "data", "$(pop)_$(i0)_$(stop).jls")
serialize(out, res)
println("NEXT ", stop + 1)
@printf("wrote %s: %d rows (%d errors) in %.0f s; maxrss %.0f MB\n", out, length(res),
        count(r -> !r.ok, res), time() - t0, Sys.maxrss() / 2^20)
