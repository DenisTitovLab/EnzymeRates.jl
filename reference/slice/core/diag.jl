# ABOUTME: Diagnostics for the core stage: the method instances a fresh Sig length compiles,
# ABOUTME: constructor cost by phase, and the tag mutation that compared equal.
const PROTO = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/proto"
include(joinpath(PROTO, "common", "Population.jl"))
include(joinpath(@__DIR__, "..", "src", "Slice.jl"))
using Serialization, Printf, Random, Statistics
const S = Slice
const DATA = joinpath(@__DIR__, "data")

# 1. Which method instances a fresh Sig length compiles in the generated call.
function snapshot()
    out = Set{Core.MethodInstance}()
    Base.visit(Core.methodtable) do m
        for mi in Base.specializations(m)
            push!(out, mi)
        end
    end
    out
end
@generated function nedges(::S.EnzymeMechanism{T}) where {T}
    length(S.decode(T...).edges)
end
ms = [c.mech for c in deserialize(joinpath(DATA, "init_bibi.jls"))]
nedges(S.compile(ms[1]))
lens = Set([length(S.sig(ms[1]))])
for m in ms[2:end]
    length(S.sig(m)) in lens && continue
    push!(lens, length(S.sig(m)))
    em = S.compile(m)
    before = snapshot()
    nedges(em)
    new = setdiff(snapshot(), before)
    println("fresh Sig length ", length(S.sig(m)), ": ", length(new), " new instances")
    for mi in first(sort!(collect(new); by = mi -> string(mi.specTypes)), 80)
        println("    ", first(string(mi.specTypes), 150))
    end
    length(lens) >= 3 && break
end

# 2. Constructor cost by phase.
function time_ns_per(f, k; reps = 11)
    f()
    median(map(_ -> (t0 = time_ns(); for _ in 1:k; Base.donotdelete(f()); end;
                     (time_ns() - t0) / k), 1:reps))
end
pops = Pop.load_population("init_bibi")
raws = [S.pop_raw(p) for p in pops]
t_all = [time_ns_per(() -> S.canonicalize(r), 50) for r in raws[1:5:end]]
t_orient = [time_ns_per(() -> S.orient(m), 50) for m in ms[1:5:end]]
t_check = [time_ns_per(() -> S.check_groups(m, S.orient(m)), 50) for m in ms[1:5:end]]
t_conn = [time_ns_per(() -> S.connected(m), 50) for m in ms[1:5:end]]
@printf("init_bibi constructor %.0f ns; orient %.0f; check_groups+orient %.0f; connected %.0f\n",
        median(t_all), median(t_orient), median(t_check), median(t_conn))

# 3. The tag mutation that compared equal (seeded as in check_core is not reproducible here,
# so try every tag change of every allosteric record and report any equal or same-hash one).
for set in ("specs", "handwritten", "allo_uniuni", "allo_bibi", "expand_bibi")
    for (i, p) in enumerate(Pop.load_population(set))
        p.allosteric || continue
        r = S.pop_raw(p)
        m = first(S.canonicalize(r))
        for g in eachindex(r.tags), t in S.TAGS
            t === r.tags[g] && continue
            tags = copy(r.tags)
            tags[g] = t
            m2 = try
                first(S.canonicalize(r.substrates, r.products, r.steps, r.group; tags,
                                     r.multiplicity, r.sites))
            catch
                continue
            end
            (m2 == m || hash(m2) == hash(m)) &&
                println("tag mutation equal: $(p.id) group $g $(r.tags[g]) -> $t: ",
                        m2 == m, " ", hash(m2) == hash(m), " tags ", m.tags, " vs ", m2.tags)
        end
    end
end
