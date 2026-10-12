# ABOUTME: Largest residual of my constraints at the package's sample points, equal-count vs differing.
# ABOUTME: Stdlib only (probe P5b).
using Serialization
rs = [r for p in ("specs", "handwritten", "allo_uniuni", "allo_bibi", "expand_bibi")
      for r in deserialize(joinpath(@__DIR__, "data", "results_$(p).jls"))]
eq = filter(r -> r[:count] == r[:pkg_count], rs); df = filter(r -> r[:count] != r[:pkg_count], rs)
println("equal counts: max residual ", maximum(r[:cons_err] for r in eq), " over ", length(eq))
println("differing counts: min residual ", minimum(r[:cons_err] for r in df), " median ", sort([r[:cons_err] for r in df])[end ÷ 2], " over ", length(df))
