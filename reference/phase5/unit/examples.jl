# ABOUTME: Printed examples for the unit-rule report: coefficient forms under the spec rule and
# ABOUTME: under Cleland's N where they differ, and today's spec law lengths against the new ones.
include(joinpath(@__DIR__, "ulib.jl"))
include(joinpath(@__DIR__, "pops.jl"))

"""The coefficient form's `v` line and `N` line of plain law `L` over normalizer `nz`."""
function vlines(m, L, sol, nz)
    _, coef = plain_texts(m, L, sol, nz)
    filter(l -> startswith(l, "v =") || startswith(l, "N =") || startswith(l, "d_0") ||
                startswith(l, "n_0"), split(coef, "\n"))
end

want = Dict("spec" => ["spec:Segel Bi Uni Uni Uni PP Ter Bi", "spec:Segel Ping Pong Bi Bi"],
            "pingpong" => ["6e46f55b0b68e5d9"])
for (pop, ids) in want, (id, text, subs, prods) in population(pop)
    id in ids || continue
    m = parse_neutral(text, subs, prods)
    sol = S.solve(m)
    L = S.rate_law(m, sol)
    println("\n=== ", id, "\n", text)
    for (label, nz) in (("spec rule", spec_norm(m, L)), ("Cleland", cleland_norm(m, L)),
                        ("order sqfree", spec_norm(m, L; order = :sqfree)))
        nz === nothing && continue
        println("--- ", label, ": ", nz.kind, " λ = ", lamtext(nz, L.vars.names))
        foreach(println, vlines(m, L, sol, nz))
    end
end

println("\n=== Lengths of today's printed spec laws and the new full law (spec rule)")
today = today_laws()
rng = Xoshiro(5)
for (id, text, subs, prods) in population("spec")
    m = parse_neutral(text, subs, prods)
    r = analyze(m, rng)
    @printf("%-60s today %6d  new %6d  coefficient form %6d\n", first(id, 60),
            length(today[id]), r.len_full, r.len_coef)
end
