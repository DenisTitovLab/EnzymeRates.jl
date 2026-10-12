# ABOUTME: Sanity counts over saved census runs: errors per stage, package vs slice parameter
# ABOUTME: counts and ranks, Float64/BigFloat rank agreement; prints a few disagreements.
using Serialization, Printf
const CEN = @__DIR__
for name in split(ARGS[1], ",")
    d = deserialize(joinpath(CEN, "data", "$name.jls"))
    recs = d.recs
    ok = [r for r in recs if r[:status] === :ok]
    fit = [r for r in ok if r[:fitted]]
    ranked = [r for r in fit if get(r, :ranked, false)]
    c(f, xs) = count(f, xs)
    println("== $(d.label): fresh $(length(recs)), ok $(length(ok)), fitted $(length(fit)), ranked $(length(ranked))")
    println("   slice_err $(c(r -> haskey(r, :slice_err), ok)), pkg_err $(c(r -> haskey(r, :pkg_err), ok)), rank_err $(c(r -> haskey(r, :rank_err), ranked))")
    for r in ok
        haskey(r, :slice_err) && (println("     slice_err #$(r[:gid]): ", first(r[:slice_err], 160)); break)
    end
    for r in ok
        haskey(r, :pkg_err) && (println("     pkg_err #$(r[:gid]): ", first(r[:pkg_err], 160)); break)
    end
    npm = [r for r in ok if haskey(r, :np_slice) && r[:np] != r[:np_slice]]
    println("   np pkg != slice: $(length(npm)) of $(length(ok)) ok")
    rkm = [r for r in ranked if haskey(r, :rank_pkg) && haskey(r, :rank_slice) && r[:rank_pkg] != r[:rank_slice]]
    println("   rank pkg != slice: $(length(rkm)) of $(length(ranked)) ranked")
    for r in rkm[1:min(end, 3)]
        println("     #$(r[:gid]) np $(r[:np])/$(r[:np_slice]) rank $(r[:rank_pkg])/$(r[:rank_slice]) gaps $(r[:rank_pkg_gap]) $(r[:rank_slice_gap])\n       $(r[:text])")
    end
    bigr = c(r -> haskey(r, :rank_pkg_gap) && r[:rank_pkg_gap][3], ranked)
    neg = c(r -> get(r, :rank_pkg, 0) < 0, ranked)
    disagree = c(r -> haskey(r, :f64_rank_agrees) && !r[:f64_rank_agrees], ranked)
    println("   pkg rank redone in BigFloat: $bigr; pkg rank non-finite: $neg; slice Float64 rank disagrees: $disagree")
    gk = [r[:rank_slice_gap][1] for r in ranked if haskey(r, :rank_slice_gap)]
    gd = [r[:rank_slice_gap][2] for r in ranked if haskey(r, :rank_slice_gap)]
    isempty(gk) || @printf("   slice rank: smallest kept pivot %.1e, largest dropped %.1e\n", minimum(gk), maximum(gd))
    gk = [r[:rank_pkg_gap][1] for r in ranked if haskey(r, :rank_pkg_gap)]
    gd = [r[:rank_pkg_gap][2] for r in ranked if haskey(r, :rank_pkg_gap)]
    isempty(gk) || @printf("   pkg rank: smallest kept %.1e, largest dropped %.1e\n", minimum(gk), maximum(gd))
end
