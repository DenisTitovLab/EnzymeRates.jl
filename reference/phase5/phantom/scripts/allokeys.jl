# ABOUTME: Writes the dedup key of every one-subunit allosteric candidate's MWC law, so the study
# ABOUTME: can count the fits a phantom rule saves net of the key's own deduplication.
include(joinpath(@__DIR__, "plib.jl"))
include(joinpath(SP, "phase3", "census", "census_common.jl"))
const DATA = ARGS[1]
for rxn in ARGS[2:end]
    subs, prods = startswith(rxn, "uniuni") ? ([:S], [:P]) : ([:A, :B], [:P, :Q])
    rows = read_rows(joinpath(DATA, "recs_$rxn.tsv"))
    open(joinpath(DATA, "allokeys_$rxn.tsv"), "w") do io
        for r in rows[2:end]
            r[5] == "1" && r[6] == "1" && r[9] == "ok" || continue
            k = try
                canonical_named(S.allo_law(parse_neutral(r[10], subs, prods))).key
            catch e
                "error"
            end
            println(io, r[1], '\t', k)
        end
    end
    println(rxn, " done, maxrss ", round(Sys.maxrss() / 2^20), " MB")
end
