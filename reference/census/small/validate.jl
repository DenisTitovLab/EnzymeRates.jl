# ABOUTME: Checks that the census classifiers fire on known cases: the guard probe's hand-built
# ABOUTME: cases, its 34 flagged enumerator children, named mutants, and a non-saturating law.
include(joinpath(@__DIR__, "lib.jl"))
include(joinpath(SCR, "round4", "guard", "pkg_build.jl"))
const GD = joinpath(SCR, "round4", "guard", "data")

wsum(r) = begin
    w = String[]
    get(r, :v0, false) && push!(w, "v=0")
    get(r, :re_turnover, nothing) === nothing || push!(w, "RE turnover")
    for (i, f) in enumerate(get(r, :faces, String[]))
        r[:face_kind][i] === :none || push!(w, string(r[:face_kind][i], "@", f))
        r[:dvan][i] && push!(w, "Dvanish@" * f)
    end
    get(r, :admit, nothing) === nothing || push!(w, "admit")
    w
end

function check(label, m; rank = true)
    rec = Dict{Symbol, Any}()
    out = String[]
    try
        sm, sol, L, br = slice_facts!(rec, m)
        push!(out, "W=$(wsum(rec))")
        if rank
            r = slice_rank(L, br, Xoshiro(2))
            push!(out, "slice rank=$(r.rank)/$(r.n)")
        end
    catch e
        push!(out, "slice ERROR " * first(errmsg(e), 120))
    end
    try
        pkg_facts!(rec, m, Xoshiro(1); rank = false)
        push!(out, "vmax=$(rec[:vmax_fwd])/$(rec[:vmax_rev]) unb=$(rec[:unb_fwd])/$(rec[:unb_rev]) " *
              "condC=$(rec[:condC]) kcat_err=$(rec[:kcat_err] === nothing ? "-" : first(rec[:kcat_err], 40))")
    catch e
        push!(out, "package ERROR " * first(errmsg(e), 120))
    end
    println(rpad(label, 46), " ", join(out, " | "))
    flush(stdout)
end

println("── hand-built cases (guard probe) ──")
for r in deserialize(joinpath(GD, "specs_hand.jls"))
    m = try
        build(r.spec)
    catch e
        println(rpad(r.key, 46), " package build fails: ", first(errmsg(e), 100)); continue
    end
    check(r.key, m)
end

println("\n── the probe's flagged enumerator children (face split/trap expected) ──")
fl = deserialize(joinpath(GD, "enum_flagged.jls"))
nflag = 0
for (i, r) in enumerate(fl)
    m = build(r.child)
    rec = Dict{Symbol, Any}()
    try
        slice_facts!(rec, m)
        w = wsum(rec)
        isempty(w) || (global nflag += 1)
        i <= 6 && println("  ", r.level, " ", r.move, " probe: ", first(r.face, 60), "\n     census: ", w)
    catch e
        println("  error ", errmsg(e))
    end
end
println("flagged by the census: $nflag of $(length(fl))")

println("\n── named mutants of the probe ──")
specs = deserialize(joinpath(GD, "specs_full.jls"))
want = Set(["mut/init_bibi/26/del_cat", "mut/expand_bibi/453/flip2", "mut/init_bibi_pp/224/del_cat",
            "mut/allo_terbi/90/to_ss", "mut/p1_enum_inh/196/del_cat", "mut/p1_enum_inh/710/del_cat"])
for r in specs
    r.key in want || continue
    m = try
        build(r.spec)
    catch e
        println(rpad(r.key, 46), " package build fails: ", first(errmsg(e), 100)); continue
    end
    check(r.key, m)
end

println("\n── a law that never saturates (verify report) ──")
nonsat = @enzyme_mechanism begin
    substrates: S
    products: P
    steps: begin
        E + S <--> E(S)
        E(S) ⇌ E(P)
        E(P) ⇌ E + P
    end
end
check("E + S <--> E(S), E(S) ⇌ E(P), E(P) ⇌ E + P", ER._concrete(nonsat))
