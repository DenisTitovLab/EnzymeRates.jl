# ABOUTME: Reproducer: law type C's Poly layout (an immutable struct holding one Dict) merged with
# ABOUTME: another type in a caller loop; ARGS[1] = vec | nothing | alias (Poly as a Dict alias).
const Mono = Vector{Pair{Int, Rational{Int}}}
struct PolyW{T}
    terms::Dict{Mono, T}
end
const PolyA{T} = Dict{Mono, T}
const MODE = ARGS[1]
const P = MODE == "alias" ? PolyA{Int} : PolyW{Int}
mk(n) = MODE == "alias" ? PolyA{Int}(Mono([i => 1 // 1]) => i for i in 1:n) :
                          PolyW{Int}(Dict{Mono, Int}(Mono([i => 1 // 1]) => i for i in 1:n))
nterms(p::PolyW) = length(p.terms)
nterms(p::PolyA) = length(p)
struct Law
    H::P
    D::P
end
struct AlloLaw
    H::P
    sites::Vector{Tuple{P, P, Int}}
end
@noinline make(allo::Bool, n) = allo ? AlloLaw(mk(n), [(mk(n), mk(n), 1)]) : Law(mk(n), mk(n))

function main(rounds)
    for r in 1:rounds
        L = make(false, 1000)                       # Union{Law, AlloLaw}, as rate_law under C
        y = MODE == "nothing" ? (isodd(r) ? L.H : nothing) :
                                (L isa Law ? L.H : L.sites)
        println(r, " ", y === nothing ? 0 : y isa Vector ? length(y) : nterms(y))
        GC.gc(true)
        println("full gc ", r, " ok")
    end
end
main(3)
