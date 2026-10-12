# ABOUTME: Debug: face-check verdicts against the oracle's face verdicts in a run.
using Serialization
const R = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/round4/guard"
out = deserialize(joinpath(R, "data", "slice_$(ARGS[1]).jls"))
for r in out
    haskey(r, :struct_) || continue
    f, v = r.struct_.face_fwd, r.struct_.face_rev
    o = haskey(r, :num) && r.num isa Dict && !haskey(r.num, :err) ? (r.num[:fwd].verdict, r.num[:rev].verdict) : (:none, :none)
    (f !== nothing || v !== nothing || o != (:agree, :agree)) && println(r.key, " fwd=", f, " rev=", v, " oracle=", o)
end
