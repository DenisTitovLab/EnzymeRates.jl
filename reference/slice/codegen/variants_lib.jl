# ABOUTME: Experiment helpers: codegen variant builders B (dependents as locals) and C
# ABOUTME: (expanded denominator) over the slice, with P7 loaded for comparison.
include(joinpath(@__DIR__, "lib.jl"))
include(joinpath(PROTO, "p7", "Encode.jl"))

"""Variant B: the bare law's segment sums and tree sums with each dependent (not a fold)
as a local computed once, `Q` and `H` from the law."""
function program_b(m, sol, L::S.Law{T}) where {T}
    b = S.bare_law(m, trues(length(m.edges)), T)
    nc = length(S.metabolites(m))
    vars, defs = deepcopy(L.vars), Pair{Int, S.Poly{T}}[]
    names = S.column_names(m, sol)
    vidx(i) = vars.index[names[i]]
    deps = Dict(d.col => d for d in sol.deps)
    dloc = Dict{Int, Int}()
    dvar(col) = get!(dloc, col) do
        d = deps[col]
        S.local!(vars, defs, S.term(T, S.mono((vidx(f) => e for (f, e) in d.terms)...,
                                                vars.index[:Keq] => d.keq)))
    end
    function sub(p)
        out = S.Poly{T}()
        for (mo, c) in p.terms
            acc = S.Mono()
            for (v, e) in mo
                if v <= nc
                    acc = S.mmul(acc, S.mono(v => e))
                    continue
                end
                col = sol.index[S.Col(b.cols[v - nc]..., false)]
                col = get(sol.fold, col, col)
                acc = S.mmul(acc, S.mono((haskey(deps, col) ? dvar(col) : vidx(col)) => e))
            end
            S.addto!(out, acc, c)
        end
        out
    end
    sig, trs = map(sub, b.sigma), map(sub, b.trees)
    f(p) = S.factor!(vars, defs, p)
    den = sum(S.term(T, L.shift) * f(s) * f(t) for (s, t) in zip(sig, trs))
    S.Program(vars, defs, f(L.H) * f(L.Q), den)
end

"""Variant A: the law's own segment and tree sums, dependents substituted."""
function program_a(L::S.Law{T}) where {T}
    vars, defs = deepcopy(L.vars), Pair{Int, S.Poly{T}}[]
    f(p) = S.factor!(vars, defs, p)
    den = sum(S.term(T, L.shift) * f(s) * f(t) for (s, t) in zip(L.sigma, L.trees))
    S.Program(vars, defs, f(L.H) * f(L.Q), den)
end

"""Variant C: the expanded denominator `D` by Horner, `H` and `Q` as locals."""
function program_c(L::S.Law{T}) where {T}
    vars, defs = deepcopy(L.vars), Pair{Int, S.Poly{T}}[]
    f(p) = S.factor!(vars, defs, p)
    S.Program(vars, defs, f(L.H) * f(L.Q), L.D)
end

