# ABOUTME: Independent check of EnzymeRates.rate_equation (Full and Reduced modes) against this
# ABOUTME: verifier's exact solver at the product-activation optimum and the hand-built example.

using EnzymeRates, Printf, Serialization
isdefined(@__MODULE__, :VL) || include(joinpath(@__DIR__, "vlib.jl"))
using .VL: QQ, toQ, qpow10, rate_exact, decode, MECHS

const em = @enzyme_mechanism begin
    substrates: A, B
    products: P, Q
    steps: begin
        E + A <--> E(A)
        E + B <--> E(B)
        E(A) + B <--> E(A, B)
        E(B) + A <--> E(A, B)
        E(A, B) <--> E(P, Q)
        E(P, Q) <--> E(Q) + P
        E(P, Q) <--> E(P) + Q
        E(Q) <--> E + Q
        E(P) <--> E + P
    end
end
const mT = MECHS["T"]

"""A side of a reaction as (bound ligands, free ligands), parsed from a name like "EA_B"."""
function side(s::AbstractString)
    toks = split(s, "_")
    form = toks[1]
    startswith(form, "E") || error("unexpected form $form")
    (sort(collect(form[2:end])), sort([only(t) for t in toks[2:end]]))
end
myside(m, state, ligs) = (sort(collect(String(m.states[state])[2:end])),
                          sort([only(String(VL.LIG[l])) for l in ligs]))

"""Map every package rate-constant name to one of my steps' kf or kr."""
function mapping(kf, kr)
    vals = Dict{Symbol,Float64}()
    for p in parameters(em, Full)
        p === :E_total && continue
        s = String(p)
        startswith(s, "k_") || error("unexpected parameter $s")
        a, b = split(s[3:end], "_to_")
        sa, sb = side(a), side(b)
        hits = Float64[]
        for (i, e) in enumerate(mT.edges)
            f = myside(mT, e.u, e.cons)
            t = myside(mT, e.w, e.rel)
            (sa, sb) == (f, t) && push!(hits, Float64(kf[i]))
            (sa, sb) == (t, f) && push!(hits, Float64(kr[i]))
        end
        length(hits) == 1 || error("$s matched $(length(hits)) steps")
        vals[p] = only(hits)
    end
    vals
end

function compare(label, kf, kr, lK, points)
    vals = mapping(kf, kr)
    full = NamedTuple{parameters(em, Full)}(Tuple(p === :E_total ? 1.0 : vals[p]
                                                  for p in parameters(em, Full)))
    fp = EnzymeRates.fitted_params(em)
    red = NamedTuple{(fp..., :Keq, :E_total)}((Tuple(vals[p] for p in fp)..., exp10(lK), 1.0))
    mets = metabolites(em)
    kfQ, krQ = toQ.(kf), toQ.(kr)
    println(label, "  (Full params: ", length(parameters(em, Full)) - 1, ", fitted: ", length(fp), ")")
    for c in points
        conc = NamedTuple{mets}(Tuple(Float64(c[findfirst(==(x), VL.LIG)]) for x in mets))
        vF = rate_equation(em, conc, full, Full)
        vR = rate_equation(em, conc, red, Reduced)
        vx = Float64(rate_exact(mT, kfQ, krQ, Tuple(toQ(x) for x in c))[1])
        @printf("   (A,B,P,Q) = %-36s  Full %.12e  Reduced %.12e  exact %.12e  max rel diff %.1e\n",
                string(Float64.(c)), vF, vR, vx, max(abs(vF - vx), abs(vR - vx)) / abs(vx))
    end
end

let kf = [1.0 for _ in 1:9], kr = [1.0 for _ in 1:9]
    kf[9] = kr[9] = 1e-3
    pts = [(1e20, 1e20, 0.0, 0.0), (1e8, 1e8, 0.0, 1e-2), (1e8, 1e8, 0.0, 1.0),
           (1e8, 1e8, 0.0, 1e4), (1e8, 1e8, 0.0, 1e8), (1.0, 1.0, 0.0, 10.0)]
    compare("T example: all constants 1, E(P) <-> E + P at 1e-3 both ways, Keq = 1", kf, kr, 0.0, pts)
end
let RES = deserialize(joinpath(@__DIR__, "..", "results_full.jls")), r = RES[("T", :Q, 0.0)]
    kf, kr = decode(mT, r.x[1:15], 0.0)
    lc = r.x[16:end]
    pts = [(1e20, 1e20, 0.0, 0.0), (exp10(lc[1]), exp10(lc[2]), 0.0, exp10(lc[3])),
           (1e3, 1e5, 0.0, 1e2), (1e8, 1e8, 1e-3, 1e4)]
    compare("T :Q Keq = 1 search optimum", Float64.(kf), Float64.(kr), 0.0, pts)
end

# Sanity of the cycle-product check: it must see three independent cycles in the textbook and
# flag a 1% change of one chord's reverse constant.
let RES = deserialize(joinpath(@__DIR__, "..", "results_full.jls")), r = RES[("T", :Q, 0.0)]
    kf, kr = decode(mT, r.x[1:15], 0.0)
    dev0, rows0 = VL.thermo_check(mT, kf, kr, 0.0)
    println("thermo check at the T :Q optimum: ", length(rows0), " cycles (step, n, log10 prod, dev) = ", rows0)
    for i in (4, 8, 9, 1)
        kr2 = copy(kr)
        kr2[i] *= big(1.01)
        dev, _ = VL.thermo_check(mT, kf, kr2, 0.0)
        @printf("   kr[%d] x 1.01 -> max dev %.6f (expected log10(1.01) = %.6f)\n", i, Float64(dev),
                log10(1.01))
    end
end
