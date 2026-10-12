# ABOUTME: Skeptic P5a: exact check that the untied (iu) cone of AD1/AD11 is empty: prints the
# ABOUTME: constraint rows restricted to the moved columns, their signs and the nullspace basis.
include(joinpath(@__DIR__, "MWC.jl"))
using Serialization
adv = deserialize(joinpath(@__DIR__, "data", "adv_cases.jls"))
for id in ("AD1_mixed_square_SS", "AD11_mixed_two_plain"), untie in (false, true)
    m = only(filter(x -> x.id == id, adv))
    J = MWC.joint_solve(m; order = :pkg, untie)
    sig = MWC.cone_signs(J, :i)
    O = sort!(collect(keys(sig)))
    println("── $id  untie=$untie  moved columns (sign: -1 K→∞, +1 rate→0, 0 fixed):")
    for i in O; println("     ", rpad(J.cols[i].name, 44), " sign ", sig[i]); end
    rows = [J.C[r, O] for r in axes(J.C, 1) if any(!iszero, J.C[r, O])]
    println("   constraint rows restricted to the moved columns:")
    for r in unique(rows); println("     ", r); end
    for (t, i) in enumerate(O); sig[i] == 0 && push!(rows, [u == t ? 1 : 0 for u in eachindex(O)]); end
    B = MWC.nullbasis(permutedims(reduce(hcat, rows)))
    println("   nullspace basis (after fixing the sign-0 columns): ", isempty(B) ? "{0}" : B)
    dirs, d = MWC.cone_directions(J, :i)
    println("   cone_directions: d=$d, found ", length(dirs))
end
