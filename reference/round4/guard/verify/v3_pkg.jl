# ABOUTME: Package-side checks, part 2: RE turnover and zero-numerator laws in the package, how
# ABOUTME: loss! treats v ≡ 0 and NaN, every plain flagged child against the oracle, (c) samples.
include(joinpath(@__DIR__, "vpkglib.jl"))

"""Max relative deviation of the exact-face RE limit from the package law on each face
(products zero; substrates zero), for `npar` parameter draws; NaN law reported as such."""
function face_dev(em, ref, rng; npar = 2)
    mech = ER.Mechanism(em)
    names = ER.fitted_params(em)
    mets = collect(ER.metabolites(em))
    rxn = ER.reaction(mech)
    subs = [ER.name(x) for x in ER.substrates(rxn)]
    prods = [ER.name(x) for x in ER.products(rxn)]
    out = Dict{Symbol, Vector{String}}(:fwd => String[], :rev => String[])
    for _ in 1:npar
        params = (; (x => 10.0^(2rand(rng) - 1) for x in names)..., Keq = 10.0^(2rand(rng) - 1),
                  E_total = 1.0)
        steps = mass_action(mech, params)
        base = Dict{Symbol, BigFloat}(x => big(10.0^(2rand(rng) - 1)) for x in mets)
        nt(c) = (; (x => Float64(c[x]) for x in mets)...)
        for (f, zs) in ((:fwd, prods), (:rev, subs))
            c0 = copy(base); foreach(z -> c0[z] = big(0), zs)
            law = ER.rate_equation(em, nt(c0), params)
            fl = face_limits(steps, base, zs, ref, rng; npat = 10)
            push!(out[f], isnan(law) ? "lawNaN" : any(isnan, fl.vals) ? "nonunique" :
                  fmt(maximum(abs.(fl.vals .- law)) / abs(law)))
        end
    end
    out
end

const B_re_turnover = :(@enzyme_mechanism begin
    substrates: A, B
    products: P, Q
    steps: begin
        E + A ⇌ E(A)
        E + P ⇌ E(P)
        (E + Q ⇌ E(Q), E(P) + Q ⇌ E(P, Q))
        E(A) + B ⇌ E(P, Q)
        E(Q) + P <--> E(P, Q)
    end
end)
const N0_del_cat = @enzyme_mechanism begin
    substrates: A, B
    products: P, Q
    steps: begin
        E + B ⇌ E(B)
        E + P ⇌ E(P)
        E(B) + A ⇌ E(A, B)
        E(A, B) <--> E(P, Q)
    end
end
const G_dead_end = @enzyme_mechanism begin
    substrates: S
    products: P, Q
    steps: begin
        E + S ⇌ E(S)
        E(S) <--> E(P, Q)
        E(Q) + P ⇌ E(P, Q)
        E + Q <--> E(Q)
        E(P) + Q <--> E(P, Q)
    end
end

function main()
    rng = Xoshiro(99)

    println("## RE turnover beside an SS route, in the package (mut/expand_bibi/453/flip2)")
    try
        em = eval(B_re_turnover)
        println("  built; rate_equation_string: ", first(ER.rate_equation_string(em), 120))
    catch e
        println("  package error: ", errline(e))
    end

    println("\n## v ≡ 0 law (mut/init_bibi/26/del_cat) in the package")
    s = ER.rate_equation_string(N0_del_cat)
    println("  ", strip(last(split(s, '\n'; keepempty = false))))
    println("  _degenerate: ", ER._degenerate(ER.Mechanism(N0_del_cat)))
    names = ER.fitted_params(N0_del_cat)
    data = (group = [1, 1, 1, 2, 2], Rate = [0.5, 1.0, 1.5, 0.2, 0.4],
            A = [0.1, 1.0, 10.0, 1.0, 1.0], B = [1.0, 1.0, 1.0, 0.1, 10.0],
            P = zeros(5), Q = zeros(5))
    fp = ER.FittingProblem(N0_del_cat, data; Keq = 10.0)
    ls = [ER.loss!(randn(rng, length(names)) .* 3, fp) for _ in 1:5]
    println("  loss! at 5 random parameter points (relative mode): ", ls)
    fpa = ER.FittingProblem(N0_del_cat, data; Keq = 10.0, scale_k_to_kcat = nothing)
    println("  loss! (absolute mode): ", [ER.loss!(randn(rng, length(names)) .* 3, fpa) for _ in 1:3])
    p = (; (x => 1.0 for x in names)..., Keq = 10.0, E_total = 1.0)
    try
        ER.rescale_parameter_values(N0_del_cat, p; scale_k_to_kcat = 1.0)
        println("  rescale_parameter_values: no error")
    catch e
        println("  rescale_parameter_values error: ", errline(e))
    end

    println("\n## D = 0 on the forward face (hand ss_dead_end_at_face) in the package")
    names = ER.fitted_params(G_dead_end)
    face = (group = [1, 1, 1], Rate = [0.5, 1.0, 1.5], S = [0.1, 1.0, 10.0], P = zeros(3), Q = zeros(3))
    off = (group = [1, 1, 1], Rate = [0.5, 1.0, 1.5], S = [0.1, 1.0, 10.0], P = fill(0.01, 3),
           Q = fill(0.01, 3))
    fpf = ER.FittingProblem(G_dead_end, face; Keq = 10.0)
    fpo = ER.FittingProblem(G_dead_end, off; Keq = 10.0)
    xs = [randn(rng, length(names)) .* 3 for _ in 1:5]
    println("  loss! on initial-rate data (P = Q = 0): ", [ER.loss!(x, fpf) for x in xs])
    println("  loss! on the same points with P = Q = 0.01: ", [ER.loss!(x, fpo) for x in xs])
    println("  _degenerate: ", ER._degenerate(ER.Mechanism(G_dead_end)))

    println("\n## every plain flagged enumerator child: package law vs exact-face RE limit")
    fl = deserialize(joinpath(R, "data", "enum_flagged.jls"))
    for (i, r) in enumerate(fl)
        isempty(r.child.tags) || (println("  $i allosteric, skipped"); continue)
        rxn = dummy_rxn(r.child)
        em = ER.compile_mechanism(mybuild(r.child, rxn))
        d = face_dev(em, first(r.child.substrates), rng)
        @printf("  %2d %-45s max|face-law|/|law| fwd %s  rev %s\n", i, first(r.face, 45),
                join(d[:fwd], ","), join(d[:rev], ","))
    end

    println("\n## the probe's (c) records (D = 0 on a face): package law and the oracle")
    keys_ = deserialize(joinpath(R, "data", "face_keys.jls"))
    specs = Dict(r.key => r.spec for r in deserialize(joinpath(R, "data", "specs_full.jls")))
    n = 0
    for (mode, key, f) in keys_
        mode == "full" && haskey(specs, key) || continue
        sp = specs[key]
        isempty(sp.tags) || continue
        m = try mybuild(sp, dummy_rxn(sp)) catch e; println("  $key ctor: ", errline(e)); continue end
        n += 1; n > 4 && break
        em = ER.compile_mechanism(m)
        mech = ER.Mechanism(em)
        names = ER.fitted_params(em); mets = collect(ER.metabolites(em))
        params = (; (x => 10.0^(2rand(rng) - 1) for x in names)..., Keq = 10.0^(2rand(rng) - 1), E_total = 1.0)
        steps = mass_action(mech, params)
        base = Dict{Symbol, BigFloat}(x => big(10.0^(2rand(rng) - 1)) for x in mets)
        zs = f === :fwd ? sp.products : sp.substrates
        c0 = copy(base); foreach(z -> c0[z] = big(0), zs)
        law = ER.rate_equation(em, (; (x => Float64(c0[x]) for x in mets)...), params)
        il = interior_limits(steps, base, zs, first(sp.substrates), rng)
        flv = face_limits(steps, base, zs, first(sp.substrates), rng; npat = 5)
        @printf("  %-32s %s face: law %s; interior limits %s; exact face %s (nterm %d)\n", key, f,
                string(law), join(fmt.([r[end] for r in il]), "/"), join(unique(fmt.(flv.vals)), "/"), flv.nterm)
    end
end

main()
