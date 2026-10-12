# ABOUTME: Classifies the mechanisms whose leading-order kcat differs from the package (P5a):
# ABOUTME: dead-end leaf present, size of the saturated face, and the per-corner values.
include(joinpath(@__DIR__, "pkgvals.jl"))
ids = split(read(joinpath(@__DIR__, "logs", "kcat.log"), String), '\n')
diff = [strip(split(split(l, "differs:")[2])[1]) for l in ids if occursin("differs:", l)]
err = [strip(split(split(l, "):")[2])[1]) for l in ids if occursin("package error", l)]
pops = Dict{String, Vector{Pop.PMech}}()
function get_m(id)
    p, i = split(id, "/")
    get!(() -> Pop.load_population(p), pops, p)[parse(Int, i)]
end
for id in [diff; err]
    m = get_m(id)
    T = Poly.Vars()
    num = Dict{Poly.Mono,Int}(); den = Dict{Poly.Mono,Int}()
    if !m.allosteric
        for (p, d) in ((m.pkg_num, num), (m.pkg_den, den)), (mono, c) in p
            Poly.addto!(d, Poly.mono((Poly.var!(T, s) => e for (s, e) in mono)...), Int(c))
        end
        prods = Set(T.index[p] for p in m.products if haskey(T.index, p))
        oth = Set(T.index[p] for p in setdiff(m.metabolites, [m.substrates; m.products]) if haskey(T.index, p))
        keep(d) = Poly.P{Int}(Dict(mn => c for (mn, c) in d if !any(v -> v.first in prods || v.first in oth, mn)))
        ds = [Dict(T.index[s] => MWC.Q(-1) for s in m.substrates if haskey(T.index, s))]
        on, pn = MWC.initial(keep(num), ds); od, pd = MWC.initial(keep(den), ds)
        cm(mn) = [p for p in mn if T.names[p.first] in m.substrates]
        nface = length(unique(cm(mn) for mn in keys(pd.t)))
        println(rpad(id, 18), " dead_end_leaf=", :dead_end_leaf in m.family, " inhibitor=", :inhibitor in m.family,
                " ord num/den (corner 0) ", on, "/", od, " den-face patterns ", nface)
    else
        println(rpad(id, 18), " allosteric dead_end_leaf=", :dead_end_leaf in m.family, " tags=", m.group_tags)
    end
end
