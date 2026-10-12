# ABOUTME: Harder overshoot searches where the probe's 3-point search found none or few: unflagged
# ABOUTME: plain laws with SS abortive cycles, MWC n=1 without them, MWC n=2 without tag effects.
include(joinpath(@__DIR__, "vlib.jl"))
include(joinpath(@__DIR__, "overshoot_lib.jl"))

t0 = time()
rng = Xoshiro(2024)
rowsP = [x for l in ("core", "expand_terter", "terter_sample")
         for x in deserialize(joinpath(KDATA, "rows_$l.jls")) if x.ok && x.cls === :finite]
flag = Set((x.set, x.index) for x in rowsP if x.ratio > 1.01)
recs = Dict{Tuple{String, Int}, Any}()
for set in CORE
    foreach(r -> recs[(r.set, r.index)] = r.m, core_records(set))
end
foreach(r -> recs[(r.set, r.index)] = r.m, pop_records("expand_terter"))
picks = deserialize(joinpath(KDATA, "terter_sample.jls"))
for set in unique(first.(picks))
    foreach(r -> recs[(r.set, r.index)] = r.m, pop_records(set, [i for (s, i, _) in picks if s == set]))
end

"""The probe's conformation-dependence flag: a substrate binding in a :NonequalAI or :OnlyA
group, or an SS substrate binding with chemistry in such a group."""
function conf_dependent(m)
    fwd = S.orient(m)
    os = [S.oriented(m, k, fwd) for k in eachindex(m.edges)]
    nm(p) = S.mets_of(m, p)
    hp(i) = m.forms[i].res != 0 || any(x -> x.name in m.products, nm(m.forms[i].bound))
    takes_sub(o) = any(x -> x.name in m.substrates && !x.modifier, nm(o.cons))
    chem(o) = !any(x -> x.name in m.products, nm(o.cons)) &&
              (any(x -> x.name in m.products, nm(o.rel)) || (hp(o.b) && !hp(o.a)))
    tag(k) = m.tags[m.group[k]]
    any(k -> takes_sub(os[k]) && tag(k) in (:NonequalAI, :OnlyA), eachindex(os)) ||
        (any(k -> takes_sub(os[k]) && !os[k].re, eachindex(os)) &&
         any(k -> chem(os[k]) && tag(k) in (:NonequalAI, :OnlyA), eachindex(os)))
end

function run_group(label, keys, npts)
    hits = Any[]
    for k in keys
        r = test_law(recs[k], npts, rng; plo = -3.0, phi = 3.0)
        r.ok && r.ratio > 1.01 && push!(hits, (k, r))
    end
    @printf("%s: %d laws searched (%d points each, constants 1e-3..1e3, box 1e-8..1e8): sup > 1.01 kcat in %d (%.0f s)\n",
            label, length(keys), npts, length(hits), time() - t0)
    for (k, r) in hits[1:min(3, end)]
        @printf("   %s/%d ratio %.4g (kcat %.4g, sup %.4g at log10 %s) oracle %.1e\n", k[1], k[2], r.ratio,
                r.kc, r.sup, string(round.(r.at; digits = 2)), r.orc)
        foreach(l -> println("      ", l), steps_text(recs[k]))
    end
    flush(stdout)
    hits
end

fin = [(x.set, x.index) for x in rowsP]
plainTTT = [k for k in fin if !S.is_allosteric(recs[k]) && overshoot_class(recs[k])[3] && !(k in flag)]
run_group("plain, SS step at an abortive form on a non-turnover cycle, not flagged by the probe",
          plainTTT[randperm(rng, length(plainTTT))[1:min(300, end)]], 6)
allo1 = [k for k in fin if S.is_allosteric(recs[k]) && recs[k].multiplicity == 1 &&
         !overshoot_class(recs[k])[3]]
run_group("MWC n = 1 without an SS step at an abortive form on a non-turnover cycle",
          allo1[randperm(rng, length(allo1))[1:min(300, end)]], 6)
allo2 = [k for k in fin if S.is_allosteric(recs[k]) && recs[k].multiplicity == 2 &&
         !overshoot_class(recs[k])[3] && !conf_dependent(recs[k])]
run_group("MWC n = 2 without abortive SS cycles or conformation-dependent substrate binding",
          allo2[randperm(rng, length(allo2))[1:min(300, end)]], 6)
@printf("done in %.0f s, maxrss %.0f MB\n", time() - t0, Sys.maxrss() / 2^20)
