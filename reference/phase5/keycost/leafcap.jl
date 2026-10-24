# ABOUTME: Item 12: what the key's search does at and near the leaf cap, on toy laws whose fitted
# ABOUTME: constants are symmetric (k interchangeable constants: k! leaves; m tied pairs: 2^m).
include(joinpath(@__DIR__, "klib.jl"))

"""A plain law over concentration S and constants `ks`, with `H = S`, `Q = 1` and
`D = 1 + Σ_j S^d(j)·Π ks[j]`: `groups[j]` lists the constants of term-group `j`."""
function toy(k::Int; pairs::Bool = false)
    vars = S.Vars()
    s = S.var!(vars, :S; conc = true)
    ks = [S.var!(vars, Symbol("k", i)) for i in 1:k]
    S.var!(vars, :Keq)
    one_ = S.term(Int, S.Mono())
    D = copy(one_.terms)
    D = S.Poly{Int}(D)
    for (i, v) in enumerate(ks)
        # pairs: constants 2j-1 and 2j share the concentration power j, so only the
        # pair members tie; otherwise every constant multiplies S alike
        d = pairs ? cld(i, 2) : 1
        S.addto!(D, S.mono(s => d, v => 1), 1)
    end
    S.Law(vars, S.term(Int, S.mono(s => 1)), one_, D, [D], [D], S.Mono())
end

for (label, L) in [("k=6 symmetric", toy(6)), ("k=7 symmetric", toy(7)),
                   ("k=8 symmetric", toy(8)), ("k=9 symmetric", toy(9)),
                   ("m=12 tied pairs", toy(24; pairs = true)),
                   ("m=16 tied pairs", toy(32; pairs = true)),
                   ("m=17 tied pairs", toy(34; pairs = true))]
    t = @elapsed r = try
        canonical_named(L)
    catch e
        e
    end
    if r isa Exception
        @printf("%-16s %6d terms: %s after %.2f s\n", label, key_terms(L),
                sprint(showerror, r), t)
    else
        @printf("%-16s %6d terms: %d leaves, %d nodes, %.3f s, key %s\n", label, key_terms(L),
                r[2], r[3], t, r[1].key[1:16])
    end
    flush(stdout)
end
