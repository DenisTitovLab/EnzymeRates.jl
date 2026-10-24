# ABOUTME: Copy of the skeptic's bbR4.jl (the F6 GC-crash workload) with switches that each
# ABOUTME: change one input: how rows load, how expand_bibi records are built, which sets run.
include("/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/slice/skeptic/scripts/indep.jl")

const NTOP = parse(Int, get(ENV, "NTOP", "40"))
const ROWS = get(ENV, "ROWS", "full")      # full | none (todo from todo.tsv)
const RECS = get(ENV, "RECS", "deser")     # deser | rebuild (Converted built from Pop records) | skip
const NEXP = parse(Int, get(ENV, "NEXP", "3"))
const TODO = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/phase3/gcprobe/todo.tsv"

function load_todo()
    if ROWS == "full"
        rows = [x for k in 1:11
                for x in deserialize(joinpath(IDATA_IN, "pipe_terter$k.jls")) if x.ok]
        top = sort(rows; by = x -> -x.nD)[1:NTOP]
        return [(x.id, x.set, x.index) for x in top], top
    end
    t = [split(l, '\t') for l in readlines(TODO)][1:NTOP]
    [(String(a), String(b), parse(Int, c)) for (a, b, c) in t], nothing
end

function expand_records()
    RECS == "deser" && return Dict(r.c.index => r for r in skeptic_records("expand_bibi"))
    row = only(r for r in (split(l, '\t') for l in readlines(joinpath(CORE, "MANIFEST.tsv")))
               if r[1] == "expand_bibi")
    pops = deserialize(row[4])
    Dict(i => (set = "expand_bibi", p = p, c = S.Converted(p, i)) for (i, p) in enumerate(pops))
end

function main()
    picked, top = load_todo()
    todo = [[("expand_bibi", "expand_bibi", i) for i in (626, 642, 658)[1:NEXP]]; picked]
    RECS == "skip" && (todo = picked)
    GC.gc()
    for set in unique(t[2] for t in todo)
        ps = set == "expand_bibi" ? nothing : deserialize(joinpath(IDATA_IN, "pop_$set.jls"))
        recs = set == "expand_bibi" ? expand_records() : nothing
        for (label, _, index) in filter(t -> t[2] == set, todo)
            c = ps === nothing ? recs[index].c : S.Converted(ps[index], index)
            L = S.rate_law(c.mech)
            nD = L isa S.Law ? S.nterms(L.D) : S.nterms(L.active.D)
            println(label, " ", nD)
            flush(stdout)
        end
        ps = recs = nothing
        GC.gc()
    end
    println("finished ", top === nothing ? 0 : length(top))
end

main()
