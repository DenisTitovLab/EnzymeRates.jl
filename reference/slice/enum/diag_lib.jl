# ABOUTME: Test helpers of the split diagnosis: export to Pop.PMech, the split's units, its
# ABOUTME: partner rule, and a P5b-style count (inactive cycles of the pruned graph only).

"""A slice mechanism as a Pop.PMech (structure only), so the package can rebuild it."""
function _testhelper_to_pmech(m::Slice.Mechanism, id)
    role(x) = x.modifier ? :Inh : x.name in m.substrates ? :S : :P
    pm(p) = [Pop.PMet(x.name, role(x)) for x in Slice.mets_of(m, p)]
    species = map(eachindex(m.forms)) do i
        f = m.forms[i]
        res = [x.name for x in Slice.mets_of(m, f.res)]
        Pop.PSpecies(pm(f.bound), m.confs[f.conf], filter(in(m.substrates), res),
                     filter(in(m.products), res), Symbol(Slice.form_name(m, i)))
    end
    Pop.PMech(; id, source = id, m.substrates, m.products, species,
              steps = [Pop.PStep(e.a, e.b, pm(e.cons), pm(e.rel), e.re) for e in m.edges],
              groups = [Slice.members(m, g) for g in 1:Slice.ngroups(m)],
              allosteric = Slice.is_allosteric(m), cat_multiplicity = m.multiplicity,
              group_tags = copy(m.tags),
              sites = [Pop.PSite(s.ligands, s.multiplicity, s.states) for s in m.sites])
end

"""The split move's units of `m` and the child of a unit set, as `split_children` builds them."""
function _testhelper_units(m)
    fwd, flux, G = Slice.orient(m), Slice.edge_flux(m), Slice.ngroups(m)
    units = [(g, p..., [k for q in p if !m.edges[q[1]].re && !any(flux[q]) for k in q])
             for g in 1:G for p in Slice.bipartitions(m, fwd, g)]
    function child(sel)
        labels, tags, re = Int.(m.group), copy(m.tags), [e.re for e in m.edges]
        for (i, (g, _, part, rev)) in enumerate(units[sel])
            labels[part] .= G + i
            Slice.is_allosteric(m) && push!(tags, m.tags[g])
            re[rev] .= true
        end
        Slice.regroup(m, labels, tags, re)
    end
    units, child
end

"""The split move's partner rule for `units` of `m`, as `split_children` applies it."""
function _testhelper_partners(m, units)
    segment = first(Slice.re_forest(m))
    near = [m.edges[Slice.members(m, g)[1]].re ?
            Set(segment[m.edges[k].a] for k in Slice.members(m, g)) : Set{Int}()
            for g in 1:Slice.ngroups(m)]
    function (sel)
        isempty(sel) && return eachindex(units)
        used = [units[u][1] for u in sel]
        touched = union((near[g] for g in used)...)
        [u for (u, (g,)) in enumerate(units) if !(g in used) && !isdisjoint(near[g], touched)]
    end
end

"""Fundamental cycles of the subgraph of the edges in `mask`, as signed edge counts in the
edges' storage directions."""
function _testhelper_cycles(m, mask)
    _, parent, order = Slice.re_forest(Slice.with_re(m, mask))
    up(i) = (k = parent[i]; e = m.edges[abs(k)]; k > 0 ? Int(e.a) : Int(e.b))
    depth = zeros(Int, length(m.forms))
    foreach(i -> parent[i] == 0 || (depth[i] = depth[up(i)] + 1), order)
    tree = falses(length(m.edges))
    foreach(p -> p == 0 || (tree[abs(p)] = true), parent)
    map(findall(k -> mask[k] && !tree[k], eachindex(m.edges))) do k
        z = zeros(Int, length(m.edges))
        z[k] = 1
        x, y = Int(m.edges[k].b), Int(m.edges[k].a)
        while x != y
            if depth[x] >= depth[y]
                p = parent[x]; z[abs(p)] -= sign(p); x = up(x)
            else
                p = parent[y]; z[abs(p)] += sign(p); y = up(y)
            end
        end
        z
    end
end

"""P5b's count: one class per group slot in the active state; in the inactive state one class
per `:NonequalAI` group that keeps an edge, `:EqualAI` edges on their active class and no
`:OnlyA` columns; rows are the cycles of the full graph and of the kept graph only."""
function _testhelper_p5b_count(m)
    fwd, G = Slice.orient(m), Slice.ngroups(m)
    keep = Slice.inactive_mask(m)
    slots(g) = m.edges[Slice.members(m, g)[1]].re ? 1 : 2
    col(k, inactive) = inactive && m.tags[m.group[k]] === :NonequalAI ? G + m.group[k] :
                       Int(m.group[k])
    rows = Vector{Int}[]
    states = Slice.is_allosteric(m) ? ((false, trues(length(m.edges))), (true, keep)) :
             ((false, trues(length(m.edges))),)
    for (inactive, mask) in states, z in _testhelper_cycles(m, mask)
        r = zeros(Int, 2G)
        foreach(k -> z[k] == 0 || (r[col(k, inactive)] += z[k] * (fwd[k] ? 1 : -1)),
                eachindex(z))
        push!(rows, r)
    end
    I = [g for g in 1:G if Slice.is_allosteric(m) && m.tags[g] === :NonequalAI &&
         any(keep[Slice.members(m, g)])]
    M = isempty(rows) ? zeros(Rational{Int}, 0, 2G) :
        Rational{Int}.(permutedims(reduce(hcat, rows)))
    sum(slots, 1:G) + sum(slots, I; init = 0) - length(Slice.rref!(M, 2G))
end
