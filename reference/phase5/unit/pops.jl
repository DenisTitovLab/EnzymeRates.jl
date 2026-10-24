# ABOUTME: Population lists for the unit-rule census, read from data/<pop>.tsv as export_pops.jl
# ABOUTME: wrote them: (id, neutral text, substrates, products), plus today's law for the specs.

const R5 = @__DIR__

"""A TSV field with census_common.jl's escapes undone."""
function unesc(s::AbstractString)
    io, i, cs = IOBuffer(), 1, collect(s)
    while i <= length(cs)
        c = cs[i]
        if c == '\\' && i < length(cs)
            n = cs[i + 1]
            print(io, n == 't' ? '\t' : n == 'n' ? '\n' : n == 'r' ? '\r' : n)
            i += 2
        else
            print(io, c)
            i += 1
        end
    end
    String(take!(io))
end

rows(path) = [[unesc(f) for f in split(l, '\t')] for l in eachline(path) if !isempty(l)]

names_of(s) = isempty(s) ? Symbol[] : Symbol.(split(s, ","))

"""
    population(name) -> Vector{(id, text, substrates, products)}

The rows of `data/<name>.tsv` (`spec`, `bibi`, `pingpong`, `terter`), in file order.
"""
population(name) = [(r[1], r[4], names_of(r[2]), names_of(r[3]))
                    for r in rows(joinpath(R5, "data", name * ".tsv"))]

"""Today's printed law of each spec, by id."""
today_laws() = Dict(r[1] => r[5] for r in rows(joinpath(R5, "data", "spec.tsv")))
