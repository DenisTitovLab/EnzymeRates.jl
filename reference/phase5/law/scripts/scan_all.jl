# ABOUTME: Runs scan.jl for several law-type variants in one Julia process (VARIANTS, in order), each
# ABOUTME: in a fresh module, with markers on stdout and stderr so split_scan.py can scan each segment.
for v in split(get(ENV, "VARIANTS", "C,B,A,U"), ',')
    println(stderr, "=== VARIANT ", v); flush(stderr)
    println("=== VARIANT ", v); flush(stdout)
    ENV["VARIANT"] = v
    path = joinpath(@__DIR__, "scan.jl")
    # a `module` expression, so the module has its own `include`
    Core.eval(Main, Expr(:module, true, Symbol("Run", v), Expr(:block, :(include($path)))))
    println(stderr, "=== END VARIANT ", v); flush(stderr)
end
