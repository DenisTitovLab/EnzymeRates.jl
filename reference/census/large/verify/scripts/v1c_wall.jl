# ABOUTME: W recounted over every fitted candidate: the island rule rebuilt from components and
# ABOUTME: arcs against the census's admit, and every face against the exact-face oracle.
include(joinpath(@__DIR__, "vlib.jl"))
setprecision(BigFloat, 256)

const RX = (:ter_bi_rxn, :ter_ter_rxn, :pyruvate_carboxylase_rxn)
rel(a, b) = Float64(abs(a - b) / max(abs(a), abs(b), big(1e-60)))

"""Every face of `m` against the oracle at one point and one ordering (Λ = 1e20); a face
that disagrees with the law is redone at two more orderings and Λ = 1e26 and classed
`:pattern`, `:mismatch` or `:nonunique` (or `:law_nan`)."""
function faces_all(m, L, rng)
    vals = random_vals(m, rng)
    conc = random_conc(m, rng)
    cols, x, _, _ = all_columns(m, vals, rng)
    sp = random_speed(m, rng; dec = 3)
    bad = Tuple{String, Symbol, Float64}[]
    for (fname, zs) in vfaces(m)
        c = copy(conc)
        foreach(z -> c[z] = big(0.0), zs)
        vl = law_value(L, lawx(L, vals, c))
        o, n = oracle(m, cols, x, vals, c; Λ = big(10.0)^20, speed = sp)
        n == 1 && !isnan(vl) && rel(vl, o) < 1e-8 && continue
        os = [oracle(m, cols, x, vals, c; Λ = big(10.0)^26, speed = random_speed(m, rng; dec = 3))
              for _ in 1:2]
        any(t -> t[2] != 1, os) || n != 1 ? push!(bad, (fname, :nonunique, NaN)) :
        rel(os[1][1], os[2][1]) > 1e-8 ? push!(bad, (fname, :pattern, rel(os[1][1], os[2][1]))) :
        isnan(vl) ? push!(bad, (fname, :law_nan, NaN)) :
        rel(vl, os[1][1]) > 1e-8 ? push!(bad, (fname, :mismatch, rel(vl, os[1][1]))) : nothing
    end
    bad
end

census_face(r) = [f for f in get(r, :faces, []) if any(w -> occursin(w, f[2]),
                                                        ("split", "trap", "bottomless"))]

rng = Xoshiro(20261011 + 12)
for name in RX
    data = load_run(name)
    println("="^100, "\n$name")
    # (a) the island rule rebuilt, against admit, every fitted allosteric record
    msgs = Dict{String, Int}()
    tab = Dict{Tuple{Int, Bool, Bool}, Int}()
    for (ri, recs) in enumerate(data.rounds), r in recs
        r[:fitted] && get(r, :allo, false) || continue
        a = r[:admit]
        a === nothing || (k = first(a, 60); msgs[k] = get(msgs, k, 0) + 1)
        st = onlya_structure(r[:smech])
        key = (ri - 1, a !== nothing, !isempty(st.lost))
        tab[key] = get(tab, key, 0) + 1
    end
    println("(a) (round, admit rejects, rebuilt island rule finds a lost component) => count: ",
            sort!(collect(tab)))
    println("    admit messages (first 60 characters): ", msgs)
    # (b) every fitted record's faces by the oracle
    t0 = time()
    ftab = Dict{Tuple{Int, Bool, Bool}, Int}()
    kinds = Dict{Symbol, Int}()
    news = Any[]
    for (ri, recs) in enumerate(data.rounds), r in recs
        r[:fitted] || continue
        m = r[:smech]
        bad = faces_all(m, S.rate_law(m), rng)
        cf = !isempty(census_face(r))
        key = (ri - 1, cf, !isempty(bad))
        ftab[key] = get(ftab, key, 0) + 1
        foreach(b -> kinds[b[2]] = get(kinds, b[2], 0) + 1, bad)
        if !isempty(bad) && !cf
            push!(news, (r[:id], bad, r[:text]))
        end
        if cf || !isempty(bad)
            println("   $(r[:id]) census $(census_face(r)) oracle $(bad)")
        end
    end
    @printf("(b) faces (%.0f s): (round, census face class, oracle finds a face) => %s; oracle face kinds %s\n",
            time() - t0, sort!(collect(ftab)), kinds)
    for (id, bad, txt) in news[1:min(6, end)]
        println("   ORACLE-ONLY $id $bad\n      ", first(txt, 600))
    end
    flush(stdout)
end
