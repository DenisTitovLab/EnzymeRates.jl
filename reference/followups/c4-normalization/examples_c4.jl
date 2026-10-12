# ABOUTME: Card C4 examples: denominators divided by the constant-term coefficient (each term as c/K),
# ABOUTME: and K-forms D/N = Σ γ_m Π (X/K_X)^m_X of laws without a constant term, in kinetic constants.
include(joinpath(@__DIR__, "kform_lib.jl"))
include(joinpath(@__DIR__, "..", "..", "proto", "p4", "Build.jl"))

rt(L, a, b) = b == one(PI) ? P4.fmt_factored(a, L.vars) :
              P4.fmt_factored(a, L.vars) * " / (" * P4.fmt_factored(b, L.vars) * ")"
cm(L, mo) = isempty(mo) ? "1" : P4.fmt_mono(mo, L.vars)

"""D/d_0 term by term: each concentration monomial c^m with K_m = d_0/d_m (reduced), so the term
reads c^m/K_m; the numerator likewise as n_m/d_0."""
function show_const(label, m)
    L = P4.law(m)
    D, Nm = P4.cmap(L, L.den), P4.cmap(L, L.num)
    d0 = D[Mono()]
    println("\n=== ", label, " [", m.id, "]")
    for d in L.sol.deps
        println("  dependent (", d.label, "): ", L.sol.cols[d.col].name)
    end
    println("  d_0 = ", P4.fmt_poly(d0, L.vars), "   (", nterms(d0), " term(s))")
    for mo in sort(collect(keys(D)); by = mo -> monokey(L, mo))
        isempty(mo) && continue
        a, b = P4.reduce_ratio(d0, D[mo])
        println("  D: ", cm(L, mo), " / K,  K = ", rt(L, a, b))
    end
    for mo in sort(collect(keys(Nm)); by = mo -> monokey(L, mo))
        a, b = P4.reduce_ratio(Nm[mo], d0)
        println("  num: ", cm(L, mo), " × ", rt(L, a, b))
    end
    k = kform(L, D)
    println("  K-form relations: ", isempty(k.rels) ? "none" :
            join([cm(L, mo) * "=" * string(c) for ((mo, _), c) in zip(k.rels, k.cls)], ", "))
    for ((mo, λ), c) in zip(k.rels, k.cls)
        c === :ratio || continue
        ab = cratio(D, k.basis, -λ; extra = mo)
        println("    γ(", cm(L, mo), ") = ", ab === nothing ? "(fractional)" : rt(L, ab...))
    end
end

"""K-form of `D` (or of its products-free part): N and K_X over the basis coefficients and in kinetic
constants, units, and each γ_m."""
function show_kform(label, m; products_zero = false)
    L = P4.law(m)
    D = P4.cmap(L, L.den)
    prods = Set(L.conc[x] for x in m.products)
    products_zero && (D = Dict(mo => c for (mo, c) in D if !any(((v, _),) -> v in prods, mo)))
    k = kform(L, D)
    println("\n=== ", label, " [", m.id, "]", products_zero ? " (products at zero)" : "")
    for d in L.sol.deps
        println("  dependent (", d.label, "): ", L.sol.cols[d.col].name)
    end
    println("  support: ", join([cm(L, mo) for mo in sort(collect(keys(D)); by = mo -> monokey(L, mo))], " "))
    println("  basis: ", join([cm(L, b) for b in k.basis], " "), "; rank ", k.rank, " of ", 1 + length(k.mets),
            "; N determined ", k.n_determined)
    uD = only(Set(P4.unit(L, t) for (mo, c) in D for t in keys((c * mo).t)))
    ucoef(mo) = uD .- (P4.degree(mo), 0)
    if k.N_K !== nothing
        names = ["N"; ["K_" * string(x) for x in k.mets]]
        for (i, nm) in enumerate(names)
            ex = k.N_K[i, :]
            u = reduce((a, b) -> a .+ b, (Int(numerator(e)) .* ucoef(bb) for (bb, e) in zip(k.basis, ex) if e != 0);
                       init = (0, 0))
            ab = cratio(D, k.basis, ex)
            println("  ", nm, " = ", cprod(L, k.basis, ex), " = ", ab === nothing ? "(fractional)" : rt(L, ab...),
                    "   units ", u, i == 1 ? " (constant-term units " * string(uD) * ")" : "")
        end
    end
    println("  relations: ", isempty(k.rels) ? "none" :
            join([cm(L, mo) * "=" * string(c) for ((mo, _), c) in zip(k.rels, k.cls)], ", "))
    for ((mo, λ), c) in zip(k.rels, k.cls)
        ab = cratio(D, k.basis, -λ; extra = mo)
        println("    γ(", cm(L, mo), ") = ", cprod(L, vcat([mo], k.basis), vcat([one(QQ)], -λ)), " = ",
                ab === nothing ? "(fractional)" : rt(L, ab...))
    end
end

show_const("SS uni-uni", byid("specs/1"))
show_const("SS ordered bi-bi (Segel)", byid("specs/5"))
show_const("RE random bi-bi, enumerated (grouped K)", byid("init_bibi/55"))
re_random = Build.mech("re_random_bibi", [:A, :B], [:P, :Q],
    ["E + A ⇌ E(A)", "E + B ⇌ E(B)", "E(A) + B ⇌ E(A,B)", "E(B) + A ⇌ E(A,B)",
     "E(A,B) <--> E(P,Q)", "E(P,Q) ⇌ E(Q) + P", "E(P,Q) ⇌ E(P) + Q", "E(P) ⇌ E + P",
     "E(Q) ⇌ E + Q"])
show_const("RE random bi-bi, textbook (independent K)", re_random)

show_kform("classic ping-pong bi-bi (Segel)", byid("specs/7"))
show_kform("classic ping-pong bi-bi (Segel)", byid("specs/7"); products_zero = true)
show_kform("residual ping-pong, SS chemistry", byid("handwritten/11"))
show_kform("residual ping-pong, RE, dead-end", byid("init_bibi_pp/61"))
show_kform("residual ping-pong, RE", byid("init_bibi_pp/62"))
show_kform("RE Theorell-Chance step", byid("handwritten/21"))
show_kform("hexa uni ping-pong (Segel)", byid("specs/17"); products_zero = true)
