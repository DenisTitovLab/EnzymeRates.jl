# per-node op counts of flat and factored laws
include(joinpath(@__DIR__, "Encode.jl"))
countops(b) = (n = Ref(0); w(x) = (x isa Expr && (x.head === :call && (n[] += 1); foreach(w, x.args))); w(b); n[])
ms = [Pop.load_population("expand_bibi"); Pop.load_population("init_bibi")]
for id in ARGS
    m = only(filter(x -> x.id == id, ms))
    for (name, law) in (("flat", P7.flat_law(m)), ("fact", P7.factored_law(m)))
        ir = Gen.decode(Tuple(code(law, (true, false, 0))))
        println(id, " ", name, ": deps ", length(ir.deps), " nodes ", length(ir.nodes), " num=", ir.out[1], " den=", ir.out[2])
        for (k, node) in enumerate(ir.nodes)
            ops = countops(Gen.flatsum(node, true))
            deg = sum(sum(abs(e) for (v, e) in t[2]; init = 0) + sum(p for (f, p) in t[3]; init = 0) for t in node)
            println("   node $k: terms ", length(node), " ops ", ops, " total degree ", deg, " facs ", sum(length(t[3]) for t in node))
        end
        println("   whole body ops ", countops(Gen.body(ir)))
    end
end
