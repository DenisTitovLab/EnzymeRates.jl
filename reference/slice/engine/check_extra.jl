# ABOUTME: Engine-stage extras: complexity against the package's _eq_complexity, the inactive
# ABOUTME: law of every allosteric record (smoke), and the spanning-tree cap's error.
include(joinpath(@__DIR__, "lib.jl"))
recs = load_records()
pk = [r for r in recs if r.p.eq_complexity >= 0]
same = count(r -> S.complexity(r.c.mech) == r.p.eq_complexity, pk)
@printf("complexity(m) equals the package's _eq_complexity in %d of %d records \
        (plain %d, allosteric %d)\n", same, length(pk), count(r -> !S.is_allosteric(r.c.mech), pk),
        count(r -> S.is_allosteric(r.c.mech), pk))
for r in pk
    S.complexity(r.c.mech) == r.p.eq_complexity ||
        println("  ", r.set, " ", r.p.id, ": mine ", S.complexity(r.c.mech), " package ",
                r.p.eq_complexity)
end
allo = [r for r in recs if S.is_allosteric(r.c.mech)]
nok, nempty, nfail, npos, nzero, t, nact = 0, 0, String[], 0, 0, 0.0, 0
for r in allo
    m = r.c.mech
    sol = try S.solve(m) catch; continue end
    mask = S.inactive_mask(m)
    try
        La = S.law(m, sol)
        all(>(0), values(La.Q.terms)) && !iszero(La.Q) && (global nact += 1)
    catch err
        push!(nfail, "$(r.set) $(r.p.id) active: $(sprint(showerror, err))")
    end
    any(mask) || (global nempty += 1; continue)
    try
        global t += @elapsed L = S.law(m, sol, mask; inactive = true)
        global nok += 1
        all(>(0), values(L.Q.terms)) && (global npos += 1)
        iszero(L.Q) && (global nzero += 1)
        S.complexity(m, mask)
    catch err
        push!(nfail, "$(r.set) $(r.p.id): $(sprint(showerror, err))")
    end
end
@printf("allosteric records %d: active law with a positive nonzero Q in %d; inactive law \
        derived %d (Q positive %d, of which Q = 0 %d), empty inactive mask %d, failures %d; \
        %.2f s in all\n", length(allo), nact, nok, npos, nzero, nempty, length(nfail), t)
foreach(println, first(nfail, 10))
for n in 4:7
    pairs = [(i, j) for i in 1:n for j in i + 1:n]
    msg = try
        "$(length(S.spanning_trees(n, pairs))) trees"
    catch err
        sprint(showerror, err)
    end
    println("complete graph K$n ($(n^(n - 2)) trees x $n nodes): ", msg)
end
