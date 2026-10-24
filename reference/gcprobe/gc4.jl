# ABOUTME: Re-runs the slice's `law` stage by stage, with the GC allowed to run only in the
# ABOUTME: stages listed in GCON (letters A–H), to find the stage whose code the crash needs.
include("/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/slice/src/Slice.jl")
const S = Slice
include("/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/phase3/gcprobe/mechs_src.jl")
const NM = parse(Int, get(ENV, "NM", "6"))
const ROUNDS = parse(Int, get(ENV, "ROUNDS", "4"))
const GCON = get(ENV, "GCON", "ABCDEFGH")

stage(f, c) = c in GCON ? f() : (GC.enable(false); x = f(); GC.enable(true); x)

function mylaw(m)
    T = Int
    mask = trues(length(m.edges)); inactive = false
    sol = stage(() -> S.solve(m), 'A')
    b = stage(() -> S.bare_law(m, mask, T), 'B')
    vars, names, value, keq, nc = stage('C') do
        vars, names = S.conc_vars(m), S.column_names(m, sol)
        value = Dict(i => S.mono(S.var!(vars, names[i]) => 1) for i in sol.indep)
        keq, nc = S.var!(vars, :Keq), length(S.metabolites(m))
        for d in sol.deps
            value[d.col] = foldl(S.mmul, (S.mpow(value[f], e) for (f, e) in d.terms);
                                 init = S.mono(keq => d.keq))
        end
        vars, names, value, keq, nc
    end
    sub(p) = S.subst(p, v -> v <= nc ? nothing :
                           value[sol.index[S.Col(b.cols[v - nc]..., inactive)]])
    H = S.term(T, S.conc_mono(vars, m.substrates)) -
        S.term(T, S.mmul(S.conc_mono(vars, m.products), S.mono(keq => -1)))
    Q = stage(() -> S.divexact(sub(b.num), H), 'D')
    D = stage(() -> sub(b.den), 'E')
    shift = stage(() -> S.conc_shift((Q, D), vars), 'F')
    Q2, D2 = stage(() -> (Q * shift, D * shift), 'G')
    sg, tr = stage(() -> (map(sub, b.sigma), map(sub, b.trees)), 'H')
    S.Law(vars, H, Q2, D2, sg, tr, shift)
end

function main()
    ms = [S.Mechanism(a...; kw...) for (a, kw) in RAWS[1:NM]]
    println("built ", length(ms), " GCON=", GCON); flush(stdout)
    for r in 1:ROUNDS, (k, m) in enumerate(ms)
        L = mylaw(m)
        println(r, " ", k, " ", S.nterms(L.D)); flush(stdout)
        get(ENV, "FULLGC", "0") == "1" && GC.gc(true)
    end
    println("finished")
end
main()
