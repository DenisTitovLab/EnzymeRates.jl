# ABOUTME: Prints, for fitted candidates whose package law and slice law differ in dimension,
# ABOUTME: the package's and the slice's fitted names; usage: xlaw.jl <rxn> <n>.
include(joinpath(@__DIR__, "census_lib.jl"))
include(joinpath(@__DIR__, "summarize.jl"))
name, n = Symbol(ARGS[1]), parse(Int, ARGS[2])
data = deserialize(joinpath(@__DIR__, "..", "data", "$(name).jls"))
label!(data.rounds)
xs = [r for recs in data.rounds for r in recs if r[:fitted] && :Xlaw in get(r, :labels, Symbol[])]
println("$name: $(length(xs)) fitted candidates whose package law has another rank than the slice law")
pkgextra = Dict{Int, Int}()
for r in xs
    d = r[:n_params] - r[:n_slice]
    pkgextra[d] = get(pkgextra, d, 0) + 1
end
println("package count - slice count among them: ", pkgextra)
for r in xs[1:min(n, end)]
    em = ER.compile_mechanism(r[:mech])
    _, indep = ER._dependent_param_exprs(typeof(em))
    sl = slice_fitted(r[:smech], S.solve(r[:smech]))
    norm(x) = replace(string(x), "_res_+" => "_res_plus", "_-" => "_minus")
    pk = Set(norm.(indep)); sk = Set(string.(sl))
    println("-- id $(r[:id]) package $(length(indep)) rank $(r[:pkg_rank_used]); slice $(length(sl)) rank $(r[:rank])")
    println("   ", first(r[:text], 600))
    println("   package only: ", sort!(collect(setdiff(pk, sk))))
    println("   slice only:   ", sort!(collect(setdiff(sk, pk))))
end
