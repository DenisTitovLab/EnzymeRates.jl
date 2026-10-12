# ABOUTME: Card C4 exploration: the plain laws without a constant term, their concentration support,
# ABOUTME: affine rank, K-form relation classes (full law and products at zero), and timing.
include(joinpath(@__DIR__, "kform_lib.jl"))

t0 = time()
n = 0
for pop in Pop.populations(), m in Pop.load_population(pop)
    m.allosteric && continue
    global n += 1
    L = P4.law(m)
    D = P4.cmap(L, L.den)
    haskey(D, Mono()) && continue
    S = P4.cmono(L, m.substrates)
    mons = sort(collect(keys(D)); by = mo -> monokey(L, mo))
    k = kform(L, D)
    prods = Set(L.conc[x] for x in m.products)
    D0 = Dict(mo => c for (mo, c) in D if !any(((v, _),) -> v in prods, mo))
    k0 = kform(L, D0)
    piSX = all(x -> haskey(D, mmul(S, mpow(mono(L.conc[x] => 1), -1))), m.substrates)
    println(rpad(m.id, 18), " ", rpad(m.source[1:min(end, 40)], 40), " fam=", join(m.family, ","))
    println("   support (", length(mons), "): ", join([P4.fmt_mono(mo, L.vars) for mo in mons], " "))
    println("   rank ", k.rank, " of ", 1 + length(k.mets), " | N determined ", k.n_determined,
            " | relations ", length(k.rels), " ", join(string.(k.cls), ","),
            " | ΠS/X all present ", piSX, " | c(ΠS) terms ", haskey(D, S) ? nterms(D[S]) : 0)
    println("   P=Q=0: support ", join([P4.fmt_mono(mo, L.vars) for mo in sort(collect(keys(D0)); by = mo -> monokey(L, mo))], " "),
            " | rank ", k0.rank, " | N determined ", k0.n_determined, " | relations ", join(string.(k0.cls), ","))
end
@printf("plain %d, %.1f s\n", n, time() - t0)
