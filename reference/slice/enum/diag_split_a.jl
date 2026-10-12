# ABOUTME: Split diagnosis, step 1: for every allosteric parent whose split children differ from
# ABOUTME: the package's, builds each unit set's child and exports them for the package's count.
include(joinpath(@__DIR__, "load.jl"))
include(joinpath(@__DIR__, "diag_lib.jl"))

out, mechs = Pop.PMech[], Slice.Mechanism[]
meta = NamedTuple[]
for set in ("expand_bibi", "allo_bibi")
    rows = deserialize(joinpath(S, "slice/enum/data/cmp_$set.jls"))
    for (ri, r) in enumerate(rows)
        x = r.split
        isempty(x.missing) && isempty(x.extra) && continue
        m = r.mech
        units, child = _testhelper_units(m)
        groups = [u[1] for u in units]
        # every unit set of at most three units from distinct groups
        sets = [[u] for u in eachindex(units)]
        for a in eachindex(units), b in a+1:length(units)
            groups[a] != groups[b] || continue
            push!(sets, [a, b])
            for c in b+1:length(units)
                groups[c] != groups[a] && groups[c] != groups[b] && push!(sets, [a, b, c])
            end
        end
        id0 = "$set/$(r.index)"
        push!(out, _testhelper_to_pmech(m, id0))
        push!(mechs, m)
        cands = [child(s) for s in sets]
        for (s, c) in zip(sets, cands)
            push!(out, _testhelper_to_pmech(c, id0 * "/" * join(s, ",")))
            push!(mechs, c)
        end
        where(cs) = [findfirst(==(c), cands) for c in cs]
        push!(meta, (; set, row = ri, index = r.index, id = id0, nunits = length(units),
                     groups, sets, pkg_sets = where(x.pkg), mine_sets = where(x.mine)))
    end
end
@printf("%d parents, %d mechanisms exported\n", length(meta), length(out))
@printf("package children not found among sets of <= 3 units: %d; mine: %d\n",
        count(isnothing, (i for mt in meta for i in mt.pkg_sets)),
        count(isnothing, (i for mt in meta for i in mt.mine_sets)))
@printf("round trip of the export: %d of %d\n",
        count(((p, m),) -> Slice.Converted(p, 0).mech == m, zip(out, mechs)), length(out))
serialize(joinpath(S, "slice/enum/data/diag_split_export.jls"), out)
serialize(joinpath(S, "slice/enum/data/diag_split_meta.jls"), meta)
