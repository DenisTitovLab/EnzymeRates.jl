# ABOUTME: Builds plain Pop.PMech mechanisms from step text such as "E(A) + B ⇌ E(A,B)" (module Build).
# ABOUTME: Releases are stored as the bindings they reverse, as the package stores them; P4 probe.
isdefined(Main, :Pop) || include(joinpath(@__DIR__, "..", "common", "Population.jl"))

module Build

import ..Pop

"""A form token `E`, `E(A)` or `F(A,B)` as (conformation, bound names)."""
function form(tok)
    m = match(r"^([A-Za-z]+)(?:\((.*)\))?$", strip(tok))
    m === nothing && error("bad form $tok")
    m[1], m[2] === nothing ? String[] : String.(strip.(split(m[2], ",")))
end

"""The mechanism of `steps` (strings `side ⇌ side` for RE, `side <--> side` for SS; a side is a
form followed by `+ metabolite`s). `groups` lists flat step indices per kinetic group (default
singletons). Metabolites outside `substrates` and `products` are inhibitors."""
function mech(id, substrates, products, steps; groups = nothing)
    role(x) = Symbol(x) in substrates ? :S : Symbol(x) in products ? :P : :Inh
    keys = Tuple{String, Vector{String}}[]
    fidx(f) = (k = findfirst(==(f), keys); k === nothing ? (push!(keys, f); length(keys)) : k)
    parsed = map(steps) do s
        re = occursin("⇌", s)
        l, r = strip.(split(s, re ? "⇌" : "<-->"))
        sl, sr = strip.(split(l, "+")), strip.(split(r, "+"))
        (fidx(form(sl[1])), String.(sl[2:end]), fidx(form(sr[1])), String.(sr[2:end]), re)
    end
    order = sortperm(keys; by = k -> (length(k[2]), k[1] * join(sort(k[2]))))
    pos = invperm(order)
    fname(k) = k[1] * join(sort(k[2]))
    species = [Pop.PSpecies(Pop.PMet[Pop.PMet(Symbol(b), role(b)) for b in sort(k[2])],
                            Symbol(k[1]), Symbol[], Symbol[], Symbol(fname(k))) for k in keys[order]]
    pm(xs) = Pop.PMet[Pop.PMet(Symbol(x), role(x)) for x in xs]
    side(f, xs) = join([fname(keys[f]); xs], "_")
    psteps, names = Pop.PStep[], Vector{Symbol}[]
    for (a, ca, b, rb, re) in parsed
        isempty(ca) && !isempty(rb) && ((a, ca, b, rb) = (b, rb, a, ca))   # store a release as its binding
        push!(psteps, Pop.PStep(pos[a], pos[b], pm(ca), pm(rb), re))
        push!(names, re ? [Symbol("K_", side(b, rb), "_to_", side(a, ca))] :
                          [Symbol("k_", side(a, ca), "_to_", side(b, rb)), Symbol("k_", side(b, rb), "_to_", side(a, ca))])
    end
    mets = Symbol.(unique([String.(substrates); String.(products);
                           [x for p in parsed for x in [p[2]; p[4]]]]))
    Pop.PMech(; id, source = "built:" * id, substrates = collect(Symbol, substrates),
              products = collect(Symbol, products), metabolites = mets, species, steps = psteps,
              groups = groups === nothing ? [[j] for j in eachindex(psteps)] : groups,
              step_names = names, family = [:built])
end

end # module
