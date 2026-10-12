# ABOUTME: Evaluates the package's rate_equation on the face (products or substrates at zero)
# ABOUTME: of accepted mechanisms whose slice law is 0/0 there; counts NaN results.
using EnzymeRates, Serialization, Random, Printf
const ER = EnzymeRates
const R = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/round4/guard"
include(joinpath(R, "pkg_build.jl"))
faces = deserialize(joinpath(R, "data", "face_keys.jls"))
specs = Dict{String, Any}()
for mode in ("full", "hand"), r in deserialize(joinpath(R, "data", "specs_$mode.jls"))
    specs[r.key] = r.spec
end
rng = Xoshiro(4)
nan = 0; fin = 0; rej = 0
for (mode, key, face) in faces
    sp = specs[key]
    m = try build(sp) catch e; global rej += 1; println(key, " rejected: ", first(sprint(showerror, e), 120)); continue end
    em = ER.compile_mechanism(m)
    names = ER.fitted_params(em)
    params = (; (x => 10.0^(2rand(rng) - 1) for x in names)..., Keq = 10.0^(4rand(rng) - 2), E_total = 1.0)
    mets = ER.metabolites(em)
    zs = face === :fwd ? sp.products : sp.substrates
    concs = (; (x => (x in zs ? 0.0 : 10.0^(2rand(rng) - 1)) for x in mets)...)
    pos = (; (x => 10.0^(2rand(rng) - 1) for x in mets)...)
    v0 = ER.rate_equation(em, concs, params)
    v1 = ER.rate_equation(em, pos, params)
    isnan(v0) ? (global nan += 1) : (global fin += 1)
    @printf("%-40s %s face: v = %s   (at a positive point %.4g)\n", key, face, string(v0), v1)
end
println("faces: $(length(faces)); package rate_equation NaN on the face: $nan; finite: $fin; rejected by the package: $rej")
