# ABOUTME: Helpers shared by the solve-stage check scripts: loading the converted sets with
# ABOUTME: their Pop records, maps between Pop names and solve columns, exact ranks, values.
const SCRATCH = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad"
const PROTO = joinpath(SCRATCH, "proto")
const DATA = joinpath(SCRATCH, "slice", "core", "data")
isdefined(Main, :Oracle) || include(joinpath(PROTO, "common", "Oracle.jl"))
isdefined(Main, :Slice) || include(joinpath(SCRATCH, "slice", "src", "Slice.jl"))
include(joinpath(SCRATCH, "slice", "core", "retag.jl"))
using Serialization, Printf, Random, Statistics, LinearAlgebra
const S = Slice
setprecision(BigFloat, 512)

const HANDBUILT = ("p1_direct", "p2_hand", "p2_fuzz_dump", "p2_fuzz_dump_flip")

"""Every converted record as `(set, p, c)`; `sets` restricts to the named sets."""
function load_records(sets = nothing)
    out = []
    for r in [split(l, '\t') for l in readlines(joinpath(DATA, "MANIFEST.tsv"))[2:end]]
        name = String(r[1])
        sets === nothing || name in sets || continue
        pops = deserialize(r[4])
        for c in deserialize(joinpath(DATA, name * ".jls"))
            p = name == "p5a_cases" ? retag(pops[c.index]) : pops[c.index]
            push!(out, (set = name, p = p, c = c))
        end
    end
    out
end

"""Exact rank of an integer or rational matrix."""
function rank_exact(A::AbstractMatrix)
    size(A, 1) == 0 && return 0
    M = Rational{BigInt}.(A)
    length(S.rref!(M, size(M, 2)))
end

"""
For Pop step `j`'s stored constant `r` (Oracle convention), the column slot of its edge and
the sign of its log in that slot's log: `pop_flip` read backwards. `names` selects the
step's state (`inactive` reads the inactive column; `reversed` says the Pop names of that
state run against the Pop step, as the package's pruned inactive graph may store it).
"""
function pop_slot(c, fwd, j, r; reversed = false)
    k = c.edge_of[j]
    e = c.mech.edges[k]
    same = (fwd[k] != c.flipped[j]) != reversed
    e.re ? (1, same || S.is_binding(e) ? 1 : -1) : (same ? r : 3 - r, 1)
end

"""Per-edge constant values (`constant_names` order) of the active state, or of the inactive
state (an `:EqualAI` edge reads its active name), read by name from `allp`; an inactive
edge outside the law (`mask`) gets placeholder ones, which the oracle never reads."""
function edge_values(m, allp; inactive = false, mask = trues(length(m.edges)))
    act = S.constant_names(m)
    ina = inactive ? S.constant_names(m; inactive = true) : act
    map(eachindex(m.edges)) do k
        mask[k] || return [big(1) for _ in act[k]]
        names = inactive && m.tags[m.group[k]] !== :EqualAI ? ina[k] : act[k]
        [BigFloat(getfield(allp, n)) for n in names]
    end
end

"""Random BigFloat parameter values (log-uniform in [1e-2, 1e2]) for `parameters(m)`."""
random_params(m, rng) = (; (n => Oracle.logu(rng) for n in S.parameters(m))...)
