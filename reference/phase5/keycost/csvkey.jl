# ABOUTME: Item 6: how CSV.jl reads back a hex key column of 16 or 64 digits written by CSV.write,
# ABOUTME: in a one-row file (best_equation.csv) and a many-row file; measured on random digests.
using CSV, DataFrames, SHA, Random, Printf
println("CSV ", pkgversion(CSV), "; Julia ", VERSION)

"""How a one-row CSV with `key = s` reads back: `:string` when the value survives, else the
parsed type."""
function readback(s)
    io = IOBuffer()
    CSV.write(io, DataFrame(key = [s], loss = [0.5]))
    df = CSV.read(IOBuffer(take!(io)), DataFrame)
    v = df.key[1]
    v isa AbstractString && v == s ? :string : Symbol(typeof(v))
end

rng = Xoshiro(20261014)
for len in (16, 64)
    n = 200_000
    counts = Dict{Symbol, Int}()
    examples = Dict{Symbol, Tuple{String, Any}}()
    for i in 1:n
        s = bytes2hex(sha256(string(i, ":", rand(rng, UInt64))))[1:len]
        r = readback(s)
        counts[r] = get(counts, r, 0) + 1
        r === :string || get!(examples, r, (s, nothing))
    end
    @printf("len %d: %d one-row files; %s\n", len, n,
            join(("$k $v ($(@sprintf("%.2e", v / n)))" for (k, v) in sort(collect(counts))), ", "))
    for (k, (s, _)) in examples
        io = IOBuffer()
        CSV.write(io, DataFrame(key = [s]))
        println("  example ", k, ": ", s, " reads back as ",
                repr(CSV.read(IOBuffer(take!(io)), DataFrame).key[1]))
    end
end
# Classes of 16-hex strings a number parser could take (exact fractions of 16^16)
# (Parsers.jl takes 'e' and 'f' as exponent markers: the Float64 example above)
p_digits = (10 / 16)^16
p_ef = 2 * 14 * 10.0^15 / 16.0^16   # digits with one interior 'e' or 'f'
@printf("16 hex: all digits %.2e, digits with one interior e or f %.2e, sum %.2e\n", p_digits,
        p_ef, p_digits + p_ef)
@printf("64 hex: all digits %.2e, digits with one interior e or f %.2e\n", (10 / 16)^64,
        2 * 62 * big(10)^63 / big(16)^64)
# A many-row file: one String row in the column keeps the whole column String
io = IOBuffer()
CSV.write(io, DataFrame(key = ["0123456789012345", "abcdef0123456789"]))
println("two rows (numeric + hex): ", repr(CSV.read(IOBuffer(take!(io)), DataFrame).key))
# Birthday bound for false merges at n keys
for n in (10^4, 10^6, 10^8)
    @printf("n = %.0e keys: P(two share a 64-bit prefix) ≈ %.1e; 256-bit ≈ %.1e\n", n,
            Float64(n)^2 / 2.0^65, big(n)^2 / big(2)^257)
end
