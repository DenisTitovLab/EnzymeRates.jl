# ABOUTME: Card C4 census over the 806 plain laws: units after dividing by the constant term, shape of
# ABOUTME: the reference coefficient and its groups, and the K-form classes (with products, and at P=Q=0).
include(joinpath(@__DIR__, "kform_lib.jl"))

cnt = Dict{String, Int}()
bump!(k) = (cnt[k] = get(cnt, k, 0) + 1)
ids = Dict{String, Vector{String}}()
note!(k, id) = push!(get!(ids, k, String[]), id)
d0terms = Int[]
shapes_lin = Dict{Tuple{String, String}, Int}()
shapes_all = Dict{Tuple{String, String}, Int}()
shape(p) = nterms(p) == 1 ? "mono" : "sum"
rng = Xoshiro(11)
t0 = time()
for pop in Pop.populations(), m in Pop.load_population(pop)
    m.allosteric && continue
    bump!("plain")
    L = P4.law(m)
    D, Nm = P4.cmap(L, L.den), P4.cmap(L, L.num)
    uD = Set(P4.unit(L, t) for t in keys(L.den.t))
    uN = Set(P4.unit(L, t) for t in keys(L.num.t))
    homog = length(uD) == 1 && length(uN) == 1 && only(uN) == only(uD) .+ (0, -1)
    homog && bump!("homogeneous")
    prods = Set(L.conc[x] for x in m.products)
    D0 = Dict(mo => c for (mo, c) in D if !any(((v, _),) -> v in prods, mo))
    if haskey(D, Mono())
        bump!("const")
        d0 = D[Mono()]
        u0 = Set(P4.unit(L, t) for t in keys(d0.t))
        # every term of D/d_0 dimensionless and every term of num/d_0 of unit 1/time
        ok = length(u0) == 1 && all(P4.unit(L, t) == only(u0) for t in keys(L.den.t)) &&
             all(P4.unit(L, t) == only(u0) .+ (0, -1) for t in keys(L.num.t))
        ok && bump!("const: D/d0 dimensionless termwise, num/d0 1/time termwise")
        push!(d0terms, nterms(d0))
        nterms(d0) == 1 ? bump!("const: d0 single monomial") : bump!("const: d0 sum")
        if nterms(d0) > 1
            for (mo, c) in D
                isempty(mo) && continue
                a, b = P4.reduce_ratio(c, d0)
                k = (shape(a), shape(b))
                shapes_all[k] = get(shapes_all, k, 0) + 1
                P4.degree(mo) == 1 && (shapes_lin[k] = get(shapes_lin, k, 0) + 1)
            end
        end
        k = kform(L, D; rng)
        cls = Set(k.cls)
        lab = isempty(k.rels) ? "no relation" : :ratio in cls ? "ratio" : :number in cls ? "number" : "one"
        bump!("const: K-form " * lab); note!("const: K-form " * lab, m.id)
    else
        bump!("noconst")
        S = P4.cmono(L, m.substrates)
        nterms(D[S]) == 1 ? bump!("noconst: c(ΠS) single monomial") : bump!("noconst: c(ΠS) sum")
        cle = all(x -> haskey(D, mmul(S, mpow(mono(L.conc[x] => 1), -1))), m.substrates)
        cle ? bump!("noconst: Cleland N defined") : (bump!("noconst: Cleland N undefined");
                                                     note!("noconst: Cleland N undefined", m.id))
        for (tag, DD) in (("full", D), ("P=Q=0", D0))
            k = kform(L, DD; rng)
            cls = Set(k.cls)
            lab = !k.n_determined ? "N undetermined" : isempty(k.rels) ? "one (no relation)" :
                  :ratio in cls ? "ratio" : :number in cls ? "number" : "one"
            bump!("noconst $tag: " * lab); note!("noconst $tag: " * lab, m.id)
            if tag == "full" && k.n_determined && cle
                # does the graded-basis N equal the Cleland N = Π_X c(ΠS/X) / c(ΠS)^(|S|-1)?
                x = randpoint(L, rng)
                ex = k.N_K === nothing ? nothing : k.N_K[1, :]
                if ex !== nothing && any(e -> denominator(e) != 1, ex)
                    bump!("noconst: graded N needs a fractional power")
                elseif ex !== nothing
                    Ng = prod(evalq(DD[b], x)^Int(numerator(e)) for (b, e) in zip(k.basis, ex) if e != 0)
                    Nc = prod(evalq(D[mmul(S, mpow(mono(L.conc[y] => 1), -1))], x) for y in m.substrates) /
                         evalq(D[S], x)^(length(m.substrates) - 1)
                    Ng == Nc ? bump!("noconst: graded N == Cleland N") :
                               (bump!("noconst: graded N != Cleland N"); note!("graded N != Cleland N", m.id))
                end
            end
        end
    end
end
@printf("time %.1f s\n", time() - t0)
for (k, v) in sort(collect(cnt); by = first)
    println(rpad(k, 70), v)
end
println("d0 term counts: ", sort(collect(Dict(t => count(==(t), d0terms) for t in unique(d0terms)))))
println("d0 a sum: group shapes (numerator, denominator) of d_m/d_0 after cancelling content and linear factors")
println("  linear terms: ", shapes_lin)
println("  all terms:    ", shapes_all)
for (k, v) in sort(collect(ids); by = first)
    println(k, " (", length(v), "): ", first(v, 12))
end
