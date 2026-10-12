# ABOUTME: Card C4 verification: units of the no-constant-term normalizers (Cleland and graded-basis N)
# ABOUTME: over all 49 laws, and the shapes of K_X = d_0/d_X when the constant coefficient d_0 is a sum.
include(joinpath(@__DIR__, "kform_lib.jl"))

shape(p) = nterms(p) == 1 ? "mono" : "sum"
kshape = Dict{Tuple{String, String}, Int}()
ncle = 0; ncle_ok = 0; ngr = 0; ngr_ok = 0
for pop in Pop.populations(), m in Pop.load_population(pop)
    m.allosteric && continue
    L = P4.law(m)
    D = P4.cmap(L, L.den)
    uD = only(Set(P4.unit(L, t) for t in keys(L.den.t)))       # unit of every term of D
    ucoef(mo) = uD .- (P4.degree(mo), 0)                       # unit of the coefficient d_m
    if haskey(D, Mono())
        d0 = D[Mono()]
        nterms(d0) > 1 || continue
        for mo in keys(D)
            P4.degree(mo) == 1 || continue
            a, b = P4.reduce_ratio(d0, D[mo])                  # K_X = a/b
            k = (shape(a), shape(b)); kshape[k] = get(kshape, k, 0) + 1
        end
        continue
    end
    S = P4.cmono(L, m.substrates)
    if all(x -> haskey(D, mmul(S, mpow(mono(L.conc[x] => 1), -1))), m.substrates)
        global ncle += 1
        n = length(m.substrates)
        u = reduce((p, q) -> p .+ q, (ucoef(mmul(S, mpow(mono(L.conc[x] => 1), -1))) for x in m.substrates)) .-
            (n - 1) .* ucoef(S)
        u == uD && (global ncle_ok += 1)
    end
    k = kform(L, D)
    if k.N_K !== nothing
        global ngr += 1
        u = reduce((p, q) -> p .+ q, (Int(numerator(e)) .* ucoef(b) for (b, e) in zip(k.basis, k.N_K[1, :]) if e != 0);
                   init = (0, 0))
        uK = [reduce((p, q) -> p .+ q, (Int(numerator(e)) .* ucoef(b) for (b, e) in zip(k.basis, k.N_K[i, :]) if e != 0);
                     init = (0, 0)) for i in 2:size(k.N_K, 1)]
        (u == uD && all(==((1, 0)), uK)) && (global ngr_ok += 1)
    end
end
println("no constant term: Cleland N = Π_X c(ΠS/X)/c(ΠS)^(n-1) defined in ", ncle, "; has constant-term units in ", ncle_ok)
println("no constant term: graded-basis N defined in ", ngr, "; N has constant-term units and every K_X units of concentration in ", ngr_ok)
println("d_0 a sum: shapes of K_X = d_0/d_X (numerator, denominator) after cancelling content and linear factors: ", kshape,
        "; total ", sum(values(kshape)))
