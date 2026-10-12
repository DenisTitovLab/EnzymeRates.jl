# ABOUTME: Mutant follow-up: what kcat = 0 and "never saturates" mean among RE/SS-toggled
# ABOUTME: mutants (no turnover, no forward rate, or substrate inhibition), with examples.
using Serialization, Printf, Statistics
const VD = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/round4/kcat/verify"
rows = deserialize(joinpath(VD, "data", "mutant_rows.jls"))
cnt(f, xs) = (d = Dict{Any, Int}(); foreach(x -> (k = f(x); d[k] = get(d, k, 0) + 1), xs);
              sort!(collect(d); by = string))
L = filter(x -> x.law, rows)
kind(x) = x.Qzero ? :no_turnover_Q0 : !x.zf.den0 ? :zero_over_zero_at_P0 :
          !x.zf.num0 ? :no_forward_rate_at_P0 : :forward_rate
z = filter(x -> x.cls === :zero, L)
println("kcat-zero mutants: ", length(z), " by (allo, verdict, kind): ", cnt(x -> (x.allo, x.verdict, kind(x)), z))
zf = filter(x -> kind(x) === :forward_rate, z)
println("  with a forward rate: ", length(zf), "; exact v>0 at P=0 (all coefficients positive): ",
        count(x -> x.zf.pos && x.zf.num0 && x.zf.den0, zf),
        "; some substrate alone drives v to 0: ", count(x -> any(g -> g !== nothing && g < 0, x.zf.gap), zf),
        "; uniform gap < 0: ", count(x -> x.zf.ugap !== nothing && x.zf.ugap < 0, zf))
inf = filter(x -> x.cls === :infinite, L)
println("never-saturating mutants: ", length(inf), " by (allo, verdict, kind): ", cnt(x -> (x.allo, x.verdict, kind(x)), inf))
println("  of derived mutants with a forward rate: ", count(x -> kind(x) === :forward_rate, L))
println("  growing substrate gaps (full law, X alone): ", cnt(x -> sort([g === nothing ? 99 : Int(g) for g in x.zf.gap]), inf))
for x in inf[1:min(3, end)]
    println("  example ", x.parent, " ", x.how, " allo ", x.allo, " verdict ", x.verdict, " ugap ", x.zf.ugap)
    foreach(l -> println("      ", l), x.lines)
end
# which parents: how many distinct parents give never-saturating children
println("  distinct parents: ", length(unique(x.parent for x in inf)), " of ", length(unique(x.parent for x in L)))
# enumerated parents' own classes are in zero_rows: compare
zr = deserialize(joinpath(VD, "data", "zero_rows.jls"))
par = Dict("$(x.set)/$(x.index)" => x for x in zr if x.ok)
println("  parents' own kcat class: ", cnt(x -> haskey(par, x.parent) ? par[x.parent].cls : :missing, inf))
