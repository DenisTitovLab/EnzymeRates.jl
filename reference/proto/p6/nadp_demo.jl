# ABOUTME: Shows two distinct P6 forms that render one package form name (E(NAD, P) and E(NADP)).
# ABOUTME: Content identity keeps them apart; the rendered parameter names would collide.
include(joinpath(@__DIR__, "harness.jl"))
src = read(joinpath(@__DIR__, "extras.jl"), String)
i, j = findfirst("F(ms...)", src), findfirst("println(\"forms of the NAD/P/NADP", src)
include_string(Main, src[first(i):first(j)-1])
println("forms: ", [GM.form_name(nadp, k) for k in eachindex(nadp.forms)], " — ", length(nadp.forms), " distinct forms")
