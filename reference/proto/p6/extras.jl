# ABOUTME: P6 side checks: GMech as a compile-time type parameter (Sig) and back, serialization
# ABOUTME: round trip, a rendered-name collision, tied-group frequency and inactive-state reorientation.
include(joinpath(@__DIR__, "harness.jl"))
using Printf, Serialization, Statistics

ms = _testhelper_load_all()
gms = [_testhelper_gm(_testhelper_toraw(m)) for m in ms]

# ── GMech as a type parameter: encode the canonical fields as nested tuples, decode by collect.
sig(m::GM.GMech) = (Tuple(m.substrates), Tuple(m.products), Tuple((x.name, x.inh) for x in m.mets),
    Tuple(m.confs), Tuple((f.bound, f.res, f.conf) for f in m.forms),
    Tuple((e.a, e.b, e.cons, e.rel, e.re) for e in m.edges), Tuple(m.group), Tuple(m.tags),
    m.multiplicity, Tuple((Tuple(s.ligands), Tuple(s.states), s.multiplicity) for s in m.sites))
function unsig(@nospecialize(s::Tuple))
    GM.GMech(collect(Symbol, s[1]), collect(Symbol, s[2]), [GM.Met(x...) for x in s[3]],
             collect(Symbol, s[4]), [GM.Form(x...) for x in s[5]], [GM.Edge(x...) for x in s[6]],
             collect(Int32, s[7]), collect(Symbol, s[8]), s[9],
             [GM.Site(collect(Symbol, a), collect(Symbol, b), c) for (a, b, c) in s[10]],
             Ref{Union{Nothing, GM.Analysis}}(nothing))
end
struct _testhelper_Mech{S} end
leaves(x) = x isa Tuple ? sum(leaves, x; init = 0) : 1
"""Mirror of the package's `_sig_of` step encoding (atoms and reaction tuple left out)."""
function pkgsig(m::Pop.PMech)
    met(x) = (x.role === :S ? :Substrate : x.role === :P ? :Product : :CompetitiveInhibitor, x.name)
    sp(i) = (Tuple(met.(m.species[i].bound)), m.species[i].conformation,
             (Tuple((:Substrate, x) for x in m.species[i].res_added),
              Tuple((:Product, x) for x in m.species[i].res_subtracted)))
    Tuple(Tuple((sp(m.steps[j].from), sp(m.steps[j].to), Tuple(met.(m.steps[j].consumed)),
                 Tuple(met.(m.steps[j].released)), m.steps[j].re) for j in g) for g in m.groups)
end
roundtrip_bad = 0; types_bad = 0; t_unsig = Float64[]; nleaf = Int[]; nleaf_pkg = Int[]
for (i, g) in enumerate(gms)
    s = sig(g)
    T = _testhelper_Mech{s}
    S2 = T.parameters[1]
    unsig(S2) == g || (global roundtrip_bad += 1)
    g2 = _testhelper_gm(_testhelper_toraw(ms[i]; rng = Xoshiro(i), flip = trues(length(ms[i].steps))))
    _testhelper_Mech{sig(g2)} === T || (global types_bad += 1)
    push!(t_unsig, _testhelper_time_ns(() -> unsig(S2), 50))
    push!(nleaf, leaves(s)); push!(nleaf_pkg, leaves(pkgsig(ms[i])))
end
@printf("Sig: decode(encode(m)) != m: %d; equal variants giving a different type: %d of %d\n",
        roundtrip_bad, types_bad, length(gms))
@printf("Sig decode (no re-canonicalization): median %.0f ns, max %.0f ns\n", median(t_unsig), maximum(t_unsig))
@printf("Sig leaves: GMech median %d (max %d) vs package step encoding median %d (max %d)\n",
        median(nleaf), maximum(nleaf), median(nleaf_pkg), maximum(nleaf_pkg))

# ── Serialization round trip drops the cache and keeps ==.
ser_bad = 0
for g in gms
    GM.analysis(g)
    io = IOBuffer(); serialize(io, g); seekstart(io)
    h = deserialize(io)
    (h == g && hash(h) == hash(g) && h.cache[] === nothing) || (global ser_bad += 1)
end
@printf("serialize/deserialize: not equal or cache kept: %d of %d\n", ser_bad, length(gms))

# ── Two distinct forms rendering one package name: content keeps them apart, names collide.
F(ms...) = GM.RawForm(:E, collect(GM.Met, ms), Symbol[], Symbol[])
nadp = GM.GMech([:NAD, :NADP], [:P, :Q],
    [GM.RawStep(F(), F(GM.Met(:NAD, false)), [GM.Met(:NAD, false)], GM.Met[], true),
     GM.RawStep(F(GM.Met(:NAD, false)), F(GM.Met(:NAD, false), GM.Met(:P, false)), [GM.Met(:P, false)], GM.Met[], true),
     GM.RawStep(F(), F(GM.Met(:NADP, false)), [GM.Met(:NADP, false)], GM.Met[], true)],
    [1, 2, 3])
println("forms of the NAD/P/NADP example: ", [GM.form_name(nadp, i) for i in eachindex(nadp.forms)],
        " — distinct forms: ", length(nadp.forms))

# ── How often the tie rules decide anything.
tied_groups = 0; tier23 = 0
for g in gms
    _, tier = GM.orient(g)
    global tier23 += count(>(1), tier)
    for l in 1:GM.ngroups(g)
        count(k -> g.group[k] == l && tier[k] > 1, eachindex(g.edges)) >= 2 && (global tied_groups += 1)
    end
end
@printf("steps decided by tier 2 or 3: %d; kinetic groups with two or more such steps: %d\n", tier23, tied_groups)

# ── Today's inactive state rebuild re-orients steps relative to the full graph.
allo = [m for m in ms if m.allosteric && !isempty(m.inactive_reversed)]
nrev = sum(count(m.inactive_reversed) for m in allo)
nmech = count(m -> any(m.inactive_reversed), allo)
@printf("package inactive rebuild: %d steps stored reversed vs the full graph, in %d of %d allosteric mechanisms\n",
        nrev, nmech, length(allo))
nnames = count(m -> any(j -> m.inactive_kept[j] && !isempty(m.step_names_I_pkg[j]) &&
                        m.step_names_I_pkg[j] != m.step_names_I[j], eachindex(m.steps)), allo)
@printf("allosteric mechanisms whose inactive-state names differ from the full-graph rendering: %d\n", nnames)
