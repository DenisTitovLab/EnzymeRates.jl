# ABOUTME: Follow-up verifier run in the package environment: d1 same-ray rates and a textbook slice
# ABOUTME: counterexample (stdlib only), then the EnzymeRates rate_equation cross-check.

for f in ("d1ray.jl", "pkgcheck.jl")
    println("\n######## ", f)
    t = time()
    try
        include(joinpath(@__DIR__, f))
    catch err
        println("FAILED: ", sprint(showerror, err, catch_backtrace()))
    end
    println("(", round(time() - t; digits = 1), " s)")
    flush(stdout)
end
