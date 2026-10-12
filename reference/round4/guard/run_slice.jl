# ABOUTME: Runs the slice side of the guard study on every record (core sets, skeptic cases,
# ABOUTME: mutants): structural checks, law-based support test, numeric checks; saves results.
# Usage: run_slice.jl <mode> where mode is smoke (a few records) or full.
include(joinpath(@__DIR__, "lib.jl"))
include(joinpath(@__DIR__, "advbuild.jl"))
include(joinpath(@__DIR__, "mutants.jl"))

const MODE = isempty(ARGS) ? "smoke" : ARGS[1]
const OUT = joinpath(R, "data")
mkpath(OUT)

function evaluate(r, rng)
    m = r.mech
    base = (key = r.key, set = r.set, prov = r.prov, ctor_err = r.ctor_err)
    m === nothing && return base
    st = slice_struct(m)
    lt = law_test(m)
    num = nothing
    if haskey(lt, :law)
        num = try
            numeric_checks(m, lt[:law], rng)
        catch e
            Dict(:err => errmsg(e))
        end
        delete!(lt, :law)
    end
    merge(base, (allo = S.is_allosteric(m), nforms = length(m.forms),
                 nedges = length(m.edges), nre = count(e -> e.re, m.edges),
                 struct_ = st, law = lt, num = num))
end

function main()
    t0 = time()
    recs = [core_records(); adv_records(); MODE == "smoke" ? NamedTuple[] : terter_records(40, 20, 7)]
    muts = mutant_records(MODE == "smoke" ? 15 : 1100, 20261010)
    @printf("records: %d core+adv, %d mutants (%.1f s)\n", length(recs), length(muts),
            time() - t0)
    flush(stdout)
    allr = [recs; muts]
    if MODE == "smoke"
        rng = Xoshiro(1)
        allr = allr[sort!(unique([1:40; randperm(rng, length(allr))[1:200]]))]
    end
    serialize(joinpath(OUT, "specs_$MODE.jls"),
              [(key = r.key, set = r.set, spec = r.spec) for r in allr])
    rng = Xoshiro(20261011)
    out = Vector{Any}(undef, length(allr))
    for (i, r) in enumerate(allr)
        out[i] = try
            evaluate(r, rng)
        catch e
            (key = r.key, set = r.set, prov = r.prov, ctor_err = r.ctor_err,
             fatal = errmsg(e))
        end
        i % 500 == 0 && (@printf("%d/%d %.0f s\n", i, length(allr), time() - t0); flush(stdout))
    end
    serialize(joinpath(OUT, "slice_$MODE.jls"), out)
    @printf("done %d records in %.0f s\n", length(out), time() - t0)
end

main()
