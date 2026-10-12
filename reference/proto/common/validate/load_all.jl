# ABOUTME: Loads every population and runs the README example, as a smoke test of the files.
# ABOUTME: Stdlib only.
include(joinpath(@__DIR__, "readme_example.jl"))
n = 0
for p in Pop.populations()
    global n += length(Pop.load_population(p))
end
println("populations: ", Pop.populations(), "; mechanisms: ", n)
