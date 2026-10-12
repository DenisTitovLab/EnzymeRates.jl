# ABOUTME: P7 codegen (module Gen): a positional rate law encoded in a type (RateLaw) and its
# ABOUTME: @generated Float64 evaluator, plus a per-mechanism NamedTuple wrapper; stdlib only.
module Gen

"""
A rate law as one homogeneous `NTuple{N, Int}` type parameter, every variable by position:
parameters `1:np` (fitted…, Keq, E_total), dependents `np+1:np+nd` (power products of
parameters), then the concentrations. Layout of `CODE`:
`np, nc, nd, nn, num, den, scale, recip, horner, salt`, then per dependent `len, (v, e)…`, then
per node `nterms`, and per term `coef, nmono, (v, e)…, nfacs, (node, pow)…`. A node may use
earlier nodes; `v = x[scale] · node[num] / node[den]`. `recip` precomputes reciprocals instead
of dividing per term, `horner` emits every node by greedy Horner, and `salt` keeps otherwise
equal laws apart (to measure without sharing).
"""
struct RateLaw{CODE} end

"""The decoded law (concrete types, so the generator compiles once for every law)."""
struct IR
    np::Int
    nc::Int
    deps::Vector{Vector{Pair{Int, Int}}}
    nodes::Vector{Vector{Tuple{Int, Vector{Pair{Int, Int}}, Vector{Pair{Int, Int}}}}}
    out::NTuple{3, Int}
    recip::Bool
    horner::Bool
end

function decode(@nospecialize(code))
    c = Int[code...]                   # splatting compiles once; indexing would box per call
    k = Ref(10)
    next() = c[k[] += 1]
    pairs(n) = [next() => next() for _ in 1:n]
    deps = [pairs(next()) for _ in 1:c[3]]
    nodes = [[(next(), pairs(next()), pairs(next())) for _ in 1:next()] for _ in 1:c[4]]
    IR(c[1], c[2], deps, nodes, (c[5], c[6], c[7]), c[8] == 1, c[9] == 1)
end

xs(i) = Symbol(:x, i)
rs(i) = Symbol(:r, i)
ns(k) = Symbol(:n, k)

"""Balanced binary tree of `op` over `xs`: every call has two operands (no varargs path)."""
nest(op, xs) = length(xs) == 1 ? xs[1] :
    Expr(:call, op, nest(op, xs[1:end÷2]), nest(op, xs[end÷2+1:end]))

"""A term's factors as `atom => power`: a positive exponent as `x`, a negative one as the
reciprocal `r` (or as a divisor when `recip` is false), a node reference as `n`."""
function atoms(t, recip)
    up, down = Pair{Symbol, Int}[], Pair{Symbol, Int}[]
    for (v, e) in t[2]
        e > 0 ? push!(up, xs(v) => e) : recip ? push!(up, rs(v) => -e) : push!(down, xs(v) => -e)
    end
    append!(up, [ns(f) => p for (f, p) in t[3]])
    up, down
end

expand(fs) = Any[a for (a, e) in fs for _ in 1:e]

"""A node as a flat sum of products: positive terms minus negative terms."""
function flatsum(node, recip)
    pos, neg = Any[], Any[]
    for t in node
        up, down = atoms(t, recip)
        f = expand(up)
        abs(t[1]) == 1 || pushfirst!(f, Float64(abs(t[1])))
        p = isempty(f) ? 1.0 : nest(:*, f)
        isempty(down) || (p = :($p / $(nest(:*, expand(down)))))
        push!(t[1] > 0 ? pos : neg, p)
    end
    isempty(neg) && return nest(:+, pos)
    isempty(pos) ? :(-$(nest(:+, neg))) : :($(nest(:+, pos)) - $(nest(:+, neg)))
end

"""A node by greedy Horner: take out the factor that most terms share, recurse on both parts."""
function horner(ts::Vector{Tuple{Int, Dict{Symbol, Int}}})
    c0 = sum((c for (c, a) in ts if isempty(a)); init = 0)
    rest = [(c, a) for (c, a) in ts if !isempty(a)]
    isempty(rest) && return Float64(c0)
    cnt = Dict{Symbol, Int}()
    for (_, a) in rest, s in keys(a); cnt[s] = get(cnt, s, 0) + 1; end
    best = first(sort!(collect(cnt); by = p -> (-p.second, string(p.first))))[1]
    dec(a) = (b = copy(a); b[best] == 1 ? delete!(b, best) : (b[best] -= 1); b)
    inner = horner([(c, dec(a)) for (c, a) in rest if haskey(a, best)])
    ex = inner == 1.0 ? best : :($best * $inner)
    others = [(c, a) for (c, a) in rest if !haskey(a, best)]
    c0 == 0 || push!(others, (c0, Dict{Symbol, Int}()))
    isempty(others) ? ex : :($ex + $(horner(others)))
end
hornersum(node) = horner([(t[1], Dict(atoms(t, true)[1])) for t in node])

"""The body of the evaluator: positional loads, dependents, reciprocals, nodes, `v`."""
function body(ir::IR)
    recip = ir.recip || ir.horner
    nd = length(ir.deps)
    ex = Any[Expr(:meta, :noinline)]
    append!(ex, [:($(xs(i)) = p[$i]) for i in 1:ir.np])
    append!(ex, [:($(xs(ir.np + nd + i)) = c[$i]) for i in 1:ir.nc])
    for (k, d) in enumerate(ir.deps)
        up = Any[xs(v) for (v, e) in d if e > 0 for _ in 1:e]
        dn = Any[xs(v) for (v, e) in d if e < 0 for _ in 1:-e]
        u = isempty(up) ? 1.0 : nest(:*, up)
        push!(ex, :($(xs(ir.np + k)) = $(isempty(dn) ? u : :($u / $(nest(:*, dn))))))
    end
    if recip
        neg = sort!(unique(v for node in ir.nodes for t in node for (v, e) in t[2] if e < 0))
        append!(ex, [:($(rs(v)) = 1.0 / $(xs(v))) for v in neg])
    end
    for (k, node) in enumerate(ir.nodes)
        push!(ex, :($(ns(k)) = $(ir.horner ? hornersum(node) : flatsum(node, recip))))
    end
    push!(ex, :(return $(xs(ir.out[3])) * $(ns(ir.out[1])) / $(ns(ir.out[2]))))
    Expr(:block, ex...)
end

@generated evaluate(::RateLaw{CODE}, c::NTuple{NC, Float64}, p::NTuple{NP, Float64}) where
    {CODE, NC, NP} = body(decode(CODE))

"""A mechanism singleton: its NamedTuple `rate_equation` method is defined by `define!`."""
struct Mech{ID} end

function rate_equation end

"""Define `rate_equation(::Mech{id}, concs, params)`: read the concentrations and parameters by
name into tuples and call the shared positional kernel of `law`."""
function define!(id::Symbol, law::RateLaw, cnames, pnames)
    c = Expr(:tuple, (:(concs.$n) for n in cnames)...)
    p = Expr(:tuple, (:(params.$n) for n in pnames)...)
    @eval rate_equation(::Mech{$(QuoteNode(id))}, concs::NamedTuple, params::NamedTuple) =
        evaluate($law, $c, $p)
    Mech{id}()
end

end # module
