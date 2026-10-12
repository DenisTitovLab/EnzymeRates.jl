# ABOUTME: Digests every stage output of every converted record under one copy of the slice
# ABOUTME: sources; usage: regress.jl <src dir> <out.jls> [sets…]. Two runs show the edits.
const SCRATCH = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad"
include(joinpath(SCRATCH, "proto", "common", "Population.jl"))
include(joinpath(ARGS[1], "Slice.jl"))
using Serialization, SHA, Printf
const S = Slice
const DATA = joinpath(SCRATCH, "slice", "core", "data")

"""A polynomial as sorted text."""
ptext(p) = join(sort!(["$(c)*$(mo)" for (mo, c) in p.terms]), " + ")
"""A law (plain or MWC) as text."""
ltext(L::S.Law) = join([ptext(L.H), ptext(L.Q), ptext(L.D), ptext.(L.sigma)..., ptext.(L.trees)...,
                        string(L.shift), string(L.vars.names)], "\n")
ltext(A) = join([ptext(A.H), ltext(A.active), A.inactive === nothing ? "-" : ltext(A.inactive),
                 ptext(A.fA), ptext(A.fI), [ptext(s[1]) * "|" * ptext(s[2]) * "|$(s[3])"
                                            for s in A.sites]..., string(A.n),
                 string(A.vars.names)], "\n")
"""`f()` as text, or the error's message."""
attempt(f) = try
    string(f())
catch err
    "ERROR " * sprint(showerror, err)
end
dig(s) = bytes2hex(sha256(s))[1:16]
mtext(cs) = join(sort!([string(S.sig(c)) for c in cs]), "\n")

out = Dict{Tuple{String, Int}, Vector{Pair{Symbol, String}}}()
t0 = time()
for r in [split(l, '\t') for l in readlines(joinpath(DATA, "MANIFEST.tsv"))[2:end]]
    set = String(r[1])
    length(ARGS) > 2 && !(set in ARGS[3:end]) && continue
    for c in deserialize(joinpath(DATA, set * ".jls"))
        m = c.mech
        length(ARGS) > 2 && (println(set, " ", c.index, " ", c.id); flush(stdout))
        parts = Pair{Symbol, String}[]
        add(k, f) = push!(parts, k => dig(attempt(f)))
        add(:admit, () -> S.admit(m))
        add(:boundary, () -> (S.check_boundary(m); "ok"))
        add(:complexity, () -> S.complexity(m))
        add(:names, () -> (S.constant_names(m), S.constant_names(m; inactive = true),
                           S.regulator_names(m)))
        sol = try S.solve(m) catch; nothing end
        add(:solve, () -> (sol.cols, sol.indep, [(d.col, d.terms, d.keq, d.label)
                                                  for d in sol.deps], sol.fold))
        add(:fitted, () -> S.fitted_params(m))
        L = try S.is_allosteric(m) ? S.allo_law(m, sol) : S.law(m, sol) catch; nothing end
        add(:law, () -> ltext(L))
        add(:printed, () -> S.rate_equation_string(m, L, sol))
        add(:coefficients, () -> S.coefficient_form(m, L))
        add(:kcat, () -> begin
            vars, ts = S.is_allosteric(m) ? S.kcat_terms(m, S.allo_law(m, sol, S.kcat_mask(m))) :
                       S.kcat_terms(m, sol)
            join([ptext(n) * " / " * ptext(d) for (n, d) in ts], "\n") * string(vars.names)
        end)
        add(:body, () -> S.body(S.program(m)))
        add(:bottomless, () -> S.bottomless(m))
        # the dead-end move enumerates every subset of the bound modifiers (2^20 placements
        # for p1_adv_pkg's adv_pkg/20, as the package's move does): moves up to 6 modifiers
        if length(S.inhibitors(m)) <= 6
            add(:flip, () -> mtext(S.flip_children(m)))
            add(:split, () -> mtext(S.split_children(m)))
            add(:dead_end, () -> mtext(S.dead_end_children(m, [:I])))
        end
        add(:flux, () -> (S.edge_flux(m), S.catalytic_core(m)))
        out[(set, c.index)] = parts
    end
    @printf("%-18s done at %.0f s\n", set, time() - t0)
    flush(stdout)
end
serialize(ARGS[2], out)
println("records: ", length(out))
