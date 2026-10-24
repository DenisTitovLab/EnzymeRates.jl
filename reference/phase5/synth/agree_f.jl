# ABOUTME: Law-type variant F against the slice (U) on the law study's census sample: equal law parts,
# ABOUTME: printed law, coefficient form, key and code body; law time interleaved per mechanism; live bytes.
include("/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/phase45/research5/law/scripts/lib.jl")
const U = load_variant("U")
include("/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/phase45/research5/synth/variantF/SliceF.jl")
const F = Base.invokelatest(getfield, Main, :SliceF)
t0 = time()

"""The law's parts as BigInt dictionaries; the inactive conformation as `(QI, DI)`."""
function parts(S, L)
    hasproperty(L, :active) || return (:plain, bdict(L.H), bdict(L.Q), bdict(L.D))
    QI, DI = S.inactive_terms(L)
    (:allo, bdict(L.H), bdict(L.active.Q), bdict(L.active.D), bdict(QI), bdict(DI), bdict(L.fA),
     bdict(L.fI), [(bdict(a), bdict(b), k) for (a, b, k) in L.sites], L.n)
end
function outputs(S, m)
    sol = S.solve(m)
    L = S.rate_law(m, sol)
    sol, L, (parts(S, L), S.rate_equation_string(m, L, sol), S.coefficient_form(m, L),
             S.canonical(L).key, S.body(S.program(m, sol, L)))
end

rows = NamedTuple[]
nbuild = 0
for r in census_rows()
    try
        push!(rows, (set = r.set, id = r.id, mu = slice_mech(U, r.net), mf = slice_mech(F, r.net)))
    catch
        global nbuild += 1
    end
end
println("julia ", VERSION, "; mechanisms: ", length(rows), "; rows neither builds: ", nbuild); flush(stdout)

diffs, errs, nallo, nflip = Dict{String, Int}(), 0, 0, 0
ratios = Float64[]
for (k, r) in enumerate(rows)
    ou = try outputs(U, r.mu) catch e; sprint(showerror, e) end
    of = try outputs(F, r.mf) catch e; sprint(showerror, e) end
    if ou isa String || of isa String
        ou isa String && of isa String && first(ou, 40) == first(of, 40) ? (global errs += 1) :
            (diffs["error in one only"] = get(diffs, "error in one only", 0) + 1)
        continue
    end
    for (j, part) in enumerate(("law", "string", "coefficients", "key", "code"))
        ou[3][j] == of[3][j] || (diffs[part] = get(diffs, part, 0) + 1)
    end
    if of[2] isa F.AlloLaw
        global nallo += 1
        global nflip += F.flip_only(of[2].inactive)
    end
    order = isodd(k) ? ((U, r.mu, ou), (F, r.mf, of)) : ((F, r.mf, of), (U, r.mu, ou))
    t = Dict{Module, Float64}()
    for (S, m, o) in order
        t[S] = best(() -> S.rate_law(m, o[1]))
    end
    push!(ratios, t[F] / t[U])
    k % 500 == 0 && (println("  ", k, " done, ", round(Int, time() - t0), " s"); flush(stdout))
end
println("errors in both: ", errs, "; differences F vs U: ", isempty(diffs) ? "none" : diffs)
println("allosteric laws: ", nallo, ", of them with a flip-only inactive conformation: ", nflip)
@printf("law time F/U over %d laws: median %.3f p10 %.3f p90 %.3f\n", length(ratios),
        q(ratios, 0.5), q(ratios, 0.1), q(ratios, 0.9))
for (label, S, field) in (("U", U, :mu), ("F", F, :mf))
    foreach(r -> try S.rate_law(getfield(r, field)) catch end, rows[1:50])
    GC.gc(true); GC.gc(true)
    lb0 = Base.gc_live_bytes()
    keep = Any[]
    for r in rows
        try
            push!(keep, S.rate_law(getfield(r, field)))
        catch
        end
    end
    GC.gc(true); GC.gc(true)
    @printf("live bytes %s: %.2f KiB per law held (%d laws)\n", label,
            (Base.gc_live_bytes() - lb0) / 1024 / length(keep), length(keep))
    keep = nothing
end
println("finished in ", round(Int, time() - t0), " s")
