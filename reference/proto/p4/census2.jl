# ABOUTME: P4 census 2: rule-A reference monomiality, graded-order ties, and which Cleland textbook forms
# ABOUTME: are identities of each plain bi-bi law (Schwartz-Zippel test, exact arithmetic mod 2^61-1).
include(joinpath(@__DIR__, "P4.jl"))
include(joinpath(@__DIR__, "bridge_text.jl"))
using Random
using .Poly: Mono, mono, P, term, nterms
const PI = P{Int}
include(joinpath(@__DIR__, "Fp.jl"))

"""Cleland bi-bi constants at a point: the coefficient ratios of `cleland_bibi`, as field elements."""
function cleland_at(L, x, m)
    N, D = P4.cmap(L, L.num), P4.cmap(L, L.den)
    c(xs...) = (p = get(D, P4.cmono(L, xs), PI()); iszero(p) ? nothing : evalp(p, x))
    n(xs...) = (p = get(N, P4.cmono(L, xs), PI()); iszero(p) ? nothing : evalp(p, x))
    r(a, b) = b === nothing ? nothing : a === nothing ? Fp(0) : a / b
    neg(a) = a === nothing ? nothing : -a
    ki(pairs...) = (for (t, tx) in pairs
                        (c(t...) === nothing || c(tx...) === nothing) || return r(c(t...), c(tx...))
                    end; nothing)
    A, B, Pr, Q = m.substrates[1], m.substrates[2], m.products[1], m.products[2]
    (V1 = r(n(A, B), c(A, B)), V2 = r(neg(n(Pr, Q)), c(Pr, Q)), Ka = r(c(B), c(A, B)),
     Kb = r(c(A), c(A, B)), Kp = r(c(Q), c(Pr, Q)), Kq = r(c(Pr), c(Pr, Q)),
     Kia = ki(((), (A,)), ((Pr,), (A, Pr))), Kib = ki(((Pr, Q), (B, Pr, Q)), ((), (B,))),
     Kip = ki(((A, B), (A, B, Pr)), ((), (Pr,))), Kiq = ki(((), (Q,)), ((B,), (B, Q))),
     Keq = r(n(A, B), neg(n(Pr, Q))))
end

rng = Xoshiro(3)
mono_A = 0; mono_B = 0; nB = 0; n = 0; ties = Int[]; nnoconst = 0; nbibi = 0
forms = Dict{String, Vector{String}}()
for pop in Pop.populations(), m in Pop.load_population(pop)
    m.allosteric && continue
    global n += 1
    L = P4.law(m)
    D = P4.cmap(L, L.den)
    S = P4.cmono(L, m.substrates)
    ref = haskey(D, Mono()) ? Mono() : S
    nterms(D[ref]) == 1 && (global mono_A += 1)
    if haskey(D, S)
        global nB += 1
        nterms(D[S]) == 1 && (global mono_B += 1)
    end
    if !haskey(D, Mono())
        global nnoconst += 1
        lo = minimum(P4.degree, keys(D))
        push!(ties, count(mo -> P4.degree(mo) == lo, keys(D)))
    end
    length(m.substrates) == 2 && length(m.products) == 2 || continue
    global nbibi += 1
    La = P4.law(m; cycles = :re)
    hits = Set{String}()
    for t in (ordered, theorell_chance, pingpong, re_random)
        ok = true
        for _ in 1:3
            x = [Fp(rand(rng, 1:Int(PR - 1))) for _ in La.vars.names]
            k = cleland_at(La, x, m)
            xc = (A = x[La.conc[m.substrates[1]]], B = x[La.conc[m.substrates[2]]],
                  P = x[La.conc[m.products[1]]], Q = x[La.conc[m.products[2]]])
            ok &= try t(k, xc) == evalp(La.num, x) / evalp(La.den, x) catch; false end
            ok || break
        end
        ok && push!(hits, string(nameof(t)))
    end
    push!(get!(forms, isempty(hits) ? "none" : join(sort(collect(hits)), "+"), String[]), m.id)
end
println("plain mechanisms: ", n)
println("rule A (constant term, else ΠS): reference coefficient a single monomial in ", mono_A, " of ", n)
println("rule B (always ΠS): ΠS present in ", nB, "; its coefficient a single monomial in ", mono_B)
println("no constant term: ", nnoconst, "; monomials at the lowest degree (graded-order ties): ",
        sort(collect(Dict(t => count(==(t), ties) for t in unique(ties)))))
println("\nbi-bi mechanisms: ", nbibi)
for (f, ids) in sort(collect(forms); by = first)
    println("  textbook identity ", rpad(f, 30), length(ids), "  e.g. ", first(ids, 8))
end
