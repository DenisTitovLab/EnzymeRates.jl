# ABOUTME: Per-conformation bare-law verdicts of allosteric records against the combined law's.
using Serialization
const R = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/round4/guard"
sl = deserialize(joinpath(R, "data", "slice_full.jls"))
c = Dict{String, Int}(); bump(k) = c[k] = get(c, k, 0) + 1
ex = Dict{String, String}()
for s in sl
    get(s, :allo, false) && haskey(s, :law) && haskey(s.law, :bare_per) || continue
    per, comb = s.law[:bare_per], s.law[:bare]
    bump("allosteric with bare laws")
    length(per) == 2 || (bump("inactive state = flip form alone"); continue)
    a, i = per
    i.b && !comb.b && (bump("inactive N≡0, combined N≢0"); get!(ex, "iN", s.key))
    a.b && !comb.b && (bump("active N≡0, combined N≢0"); get!(ex, "aN", s.key))
    (i.dfwd || i.drev) && !(comb.dfwd || comb.drev) && (bump("inactive D=0 on a face, combined not"); get!(ex, "iD", s.key))
    (a.dfwd || a.drev) && !(comb.dfwd || comb.drev) && (bump("active D=0 on a face, combined not"); get!(ex, "aD", s.key))
end
println(sort!(collect(c))); println(ex)
