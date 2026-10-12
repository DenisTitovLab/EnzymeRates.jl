# ABOUTME: Class W against the mass-action oracle: the :OnlyA island rule by vanishing limits
# ABOUTME: under several potentials, and the face classes by exact-face values at two speeds.
include(joinpath(@__DIR__, "vlib.jl"))

const RX = (:ter_bi_rxn, :ter_ter_rxn, :pyruvate_carboxylase_rxn)
const NS = 12
rel(a, b) = Float64(abs(a - b) / max(abs(a), abs(b), big(1e-60)))

"""The n=1 network with the `:OnlyA` steps vanishing, under four potentials: `A` (every
component off the flip form's arc tree above it), `A2` (the same with other free constants,
increments and speeds), `B` (sources off the tree level with it), `C` (below it); each at
t = 30 and 45. Returns the law and the pattern values at each of `npts` points."""
function vanish_test(m, L, rng; npts = 2)
    st = onlya_structure(m)
    fwd = S.orient(m)
    pts = []
    for _ in 1:npts
        vals = random_vals(m, rng)
        conc = random_conc(m, rng)
        sp = random_speed(m, rng)
        vl = law_value(L, lawx(L, vals, conc))
        c1 = all_columns(m, vals, rng)
        c2 = all_columns(m, vals, rng)
        out = Dict{String, Tuple{BigFloat, BigFloat}}()
        for (lab, off, cx) in (("A", 3, c1), ("A2", 3, c2), ("B", 0, c1), ("C", -3, c1))
            phi = potentials(st, off, rng)
            sig = map(eachindex(m.edges)) do k
                o = S.oriented(m, k, fwd)
                abs(phi[st.comp[o.a]] - phi[st.comp[o.b]]) / 2 + 2 + big(rand(rng))
            end
            v = [first(oracle(m, cx[1], cx[2], vals, conc; Λ = big(10.0)^12, speed = sp,
                              mode = :vanish, phi = phi, t = big(tt), sigma = sig))
                 for tt in (30, 45)]
            out[lab] = (v[1], v[2])
        end
        push!(pts, (law = vl, out = out, resid = c1[3], nfree = c1[4]))
    end
    st, pts
end

"""Verdict of the vanishing test: `:unique_law` when every pattern converges to the law,
`:pattern_dependent` when the islands-up patterns give the law and another pattern does
not, else a description."""
function vanish_verdict(pts)
    conv = all(p -> all(((a, b),) -> rel(a, b) < 1e-8, values(p.out)), pts)
    upok = all(p -> rel(p.out["A"][2], p.law) < 1e-8 && rel(p.out["A2"][2], p.law) < 1e-8, pts)
    alleq = all(p -> all(((a, b),) -> rel(b, p.law) < 1e-8, values(p.out)), pts)
    alleq && return :unique_law
    upok && any(p -> rel(p.out["B"][2], p.law) > 1e-6 || rel(p.out["C"][2], p.law) > 1e-6,
                pts) && return conv ? :pattern_dependent : :pattern_dependent_noconv
    Symbol("other_conv=$(conv)_up=$(upok)")
end

"""Face values of `m`: the law and the oracle on the exact face at Λ = 1e14 and 1e20 with
one speed pattern and at 1e20 with another (speeds in 10^[-3,3]), per face; the verdict
`:agree`, `:pattern`, `:mismatch`, `:diverges`, `:nonunique` or `:law_nan`."""
function face_test(m, L, rng)
    vals = random_vals(m, rng)
    conc = random_conc(m, rng)
    cols, x, _, _ = all_columns(m, vals, rng)
    sp1, sp2 = random_speed(m, rng; dec = 3), random_speed(m, rng; dec = 3)
    out = []
    for (fname, zs) in [("pos", Symbol[]); vfaces(m)]
        c = copy(conc)
        foreach(z -> c[z] = big(0.0), zs)
        vl = law_value(L, lawx(L, vals, c))
        o1, n1 = oracle(m, cols, x, vals, c; Λ = big(10.0)^14, speed = sp1)
        o2, n2 = oracle(m, cols, x, vals, c; Λ = big(10.0)^20, speed = sp1)
        o3, n3 = oracle(m, cols, x, vals, c; Λ = big(10.0)^20, speed = sp2)
        sc = max(abs(o1), abs(o2), abs(o3), big(1e-40))
        v = (n1 != 1 || n2 != 1 || n3 != 1) ? :nonunique :
            abs(o2) > 1e4 * max(abs(o1), big(1e-40)) ? :diverges :
            abs(o2 - o3) > 1e-8 * sc ? :pattern :
            isnan(vl) ? :law_nan :
            rel(vl, o2) > 1e-8 ? :mismatch : :agree
        push!(out, (face = fname, verdict = v, law = vl, o = (o1, o2, o3),
                    gap = isnan(vl) ? NaN : max(rel(vl, o2), rel(vl, o3))))
    end
    out
end

census_face(r) = [f for f in get(r, :faces, []) if any(w -> occursin(w, f[2]),
                                                        ("split", "trap", "bottomless"))]

rng = Xoshiro(20261011 + 1)
println("="^100, "\nPART 1: the island rule (W_admit) against vanishing limits")
summary1 = Dict{Tuple{Symbol, Bool, Symbol}, Int}()
for name in RX
    data = load_run(name)
    for (ri, recs) in enumerate(data.rounds)
        ri == 1 && continue
        allo = [r for r in recs if r[:fitted] && get(r, :allo, false)]
        flagged = [r for r in allo if :W_admit in r[:labels]]
        clean = [r for r in allo if r[:admit] === nothing]
        for (lab, pool) in (("flagged", flagged), ("unflagged", clean))
            isempty(pool) && continue
            pick = pool[sort!(randperm(rng, length(pool))[1:min(NS, length(pool))])]
            println("-- $name round $(ri - 1) $lab: $(length(pick)) of $(length(pool))")
            for r in pick
                m = r[:smech]
                L = S.rate_law(m)
                st, pts = vanish_test(m, L, rng)
                v = vanish_verdict(pts)
                k = (name, lab == "flagged", v)
                summary1[k] = get(summary1, k, 0) + 1
                gaps = join((@sprintf("B %.1e C %.1e A2 %.1e", rel(p.out["B"][2], p.law),
                                      rel(p.out["C"][2], p.law), rel(p.out["A2"][2], p.law))
                             for p in pts), " | ")
                @printf("   %s %-22s lost %d free %d resid %.0e | %s\n", r[:id], v,
                        length(st.lost), pts[1].nfree, Float64(pts[1].resid), gaps)
            end
            flush(stdout)
        end
    end
end
println("island-rule summary (reaction, census-flagged, oracle verdict) => count:")
foreach(kv -> println("   ", kv), sort!(collect(summary1); by = string))

println("="^100, "\nPART 2: face classes against exact-face oracle values")
summary2 = Dict{Tuple{Bool, Symbol}, Int}()
for name in RX
    data = load_run(name)
    fit = [r for recs in data.rounds for r in recs if r[:fitted]]
    flagged = [r for r in fit if !isempty(census_face(r))]
    plain = [r for r in fit if isempty(census_face(r)) && !get(r, :allo, false)]
    allo = [r for r in fit if isempty(census_face(r)) && get(r, :allo, false) &&
            r[:admit] === nothing]
    sets = [("flagged", flagged),
            ("unflagged plain", plain[sort!(randperm(rng, length(plain))[1:NS])]),
            ("unflagged allosteric", allo[sort!(randperm(rng, length(allo))[1:NS])])]
    for (lab, rs) in sets
        println("-- $name $lab: $(length(rs))")
        for r in rs
            m = r[:smech]
            L = S.rate_law(m)
            ft = face_test(m, L, rng)
            bad = [f for f in ft if f.verdict !== :agree]
            k = (lab == "flagged", isempty(bad) ? :all_agree : :some_face_fails)
            summary2[k] = get(summary2, k, 0) + 1
            cf = join(("$(f[1]): $(f[2])" for f in census_face(r)), "; ")
            oc = join((@sprintf("%s %s (gap %.1e)", f.face, f.verdict, f.gap) for f in bad),
                      "; ")
            println("   $(r[:id]) census [$(cf)] oracle [$(isempty(oc) ? "all faces agree" : oc)]")
        end
        flush(stdout)
    end
end
println("face summary (census-flagged, oracle) => count: ", sort!(collect(summary2); by = string))
