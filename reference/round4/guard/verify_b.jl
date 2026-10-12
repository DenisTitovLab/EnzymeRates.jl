# ABOUTME: Checks that law test (b) N≡0 coincides with "no turnover cycle or RE turnover", and
# ABOUTME: that law test (c) coincides with bottomless or an SS-level dead end.
using Serialization
const R = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/round4/guard"
c = Dict{String, Int}(); bump(k) = c[k] = get(c, k, 0) + 1
for s in deserialize(joinpath(R, "data", "slice_full.jls"))
    haskey(s, :law) && haskey(s.law, :reduced) || continue
    t, st = s.law[:reduced], s.struct_
    nt, rt, un = st.no_turnover !== nothing, st.re_cycle !== nothing, st.unused !== nothing
    bump("b=$(t.b) noturn=$(nt) reTurn=$(rt) unused=$(un)")
end
foreach(println, sort!(collect(c)))
