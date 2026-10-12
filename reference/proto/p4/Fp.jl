# ABOUTME: Arithmetic in the prime field mod 2^61 - 1 (type Fp) and polynomial evaluation in it, for
# ABOUTME: Schwartz-Zippel identity tests in the P4 census scripts.
const PR = Int128(2)^61 - 1

"""An element of the prime field mod 2^61 - 1."""
struct Fp
    v::Int128
end
Fp(x::Integer) = Fp(mod(Int128(x), PR))
Base.:+(a::Fp, b::Fp) = Fp(mod(a.v + b.v, PR))
Base.:-(a::Fp, b::Fp) = Fp(mod(a.v - b.v, PR))
Base.:-(a::Fp) = Fp(mod(-a.v, PR))
Base.:*(a::Fp, b::Fp) = Fp(mod(a.v * b.v, PR))
Base.inv(a::Fp) = (a.v == 0 && throw(DivideError()); Fp(powermod(a.v, PR - 2, PR)))
Base.:/(a::Fp, b::Fp) = a * inv(b)
Base.:^(a::Fp, k::Integer) = k >= 0 ? Fp(powermod(a.v, k, PR)) : inv(a)^(-k)
for op in (:+, :-, :*, :/)
    @eval Base.$op(a::Fp, b::Integer) = $op(a, Fp(b))
    @eval Base.$op(a::Integer, b::Fp) = $op(Fp(a), b)
end
evalp(p::PI, x) = foldl(+, (Fp(c) * foldl(*, (x[v]^e for (v, e) in mo); init = Fp(1)) for (mo, c) in p.t);
                        init = Fp(0))

