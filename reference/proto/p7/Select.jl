# ABOUTME: The P7 measurement set: all 42 specs and the 30 largest-denominator mechanisms of
# ABOUTME: init_bibi and expand_bibi (plain: package den terms; allosteric: both states' den terms).
isdefined(Main, :Pop) || include(joinpath(@__DIR__, "..", "common", "Population.jl"))

den_size(m) = m.allosteric ? length(m.state_polys_A[2]) + length(m.state_polys_I[2]) : length(m.pkg_den)

"""`(population, mechanism)` pairs in measurement order."""
function selected()
    out = Tuple{String, Pop.PMech}[]
    for m in Pop.load_population("specs"); push!(out, ("specs", m)); end
    for pop in ("init_bibi", "expand_bibi")
        ms = Pop.load_population(pop)
        for m in sort(ms; by = den_size, rev = true, alg = MergeSort)[1:30]   # stable: id order on ties
            push!(out, (pop, m))
        end
    end
    out
end
