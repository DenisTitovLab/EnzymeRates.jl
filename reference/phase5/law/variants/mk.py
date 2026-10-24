# ABOUTME: Builds the law-type variants of the slice (U as is, A widened to BigInt, B a Dict
# ABOUTME: alias, C checked Int only) as module directories; each patch must match exactly once.
import os, sys
SP = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad"
SRC = f"{SP}/slice/src"
OUT = f"{SP}/phase45/research5/law/variants"
FILES = ["types.jl", "names.jl", "poly.jl", "convert.jl", "thermo.jl", "api.jl", "engine.jl",
         "normalize.jl", "print.jl", "allosteric.jl", "admission.jl", "print_allo.jl",
         "codegen.jl", "enum.jl", "key.jl"]

def read(f):
    return open(f"{SRC}/{f}", encoding="utf-8").read()

def patch(text, pairs, label):
    for old, new in pairs:
        n = text.count(old)
        if n != 1:
            sys.exit(f"{label}: patch matches {n} times: {old[:60]!r}")
        text = text.replace(old, new)
    return text

LAW_U = """    try
        law(Int, m, sol, mask, inactive)
    catch err
        err isa OverflowError || rethrow()
        law(BigInt, m, sol, mask, inactive)
    end
end
"""
# A: one concrete law type, Law{BigInt}: the Int result widened at the boundary.
LAW_A = """    try
        widen_law(law(Int, m, sol, mask, inactive))
    catch err
        err isa OverflowError || rethrow()
        law(BigInt, m, sol, mask, inactive)
    end
end

\"\"\"`L` with every coefficient a `BigInt`, the one law type callers see.\"\"\"
widen_law(L::Law) = Law(L.vars, convert(Poly{BigInt}, L.H), convert(Poly{BigInt}, L.Q),
                        convert(Poly{BigInt}, L.D), [convert(Poly{BigInt}, p) for p in L.sigma],
                        [convert(Poly{BigInt}, p) for p in L.trees], L.shift)
"""
# C: one concrete law type, Law{Int}: checked Int, no retry (overflow is an error).
LAW_C = """    law(Int, m, sol, mask, inactive)
end
"""
ALLO_RETRY = """    if LI !== nothing && typeof(LA) != typeof(LI)
        LA, LI = law(BigInt, m, sol, mask, false), law(BigInt, m, sol, imask, true)
    end
"""
variants = {
    "U": {},
    "A": {"engine.jl": [(LAW_U, LAW_A)]},
    "C": {"engine.jl": [(LAW_U, LAW_C)], "allosteric.jl": [(ALLO_RETRY, "")]},
    "B": {
        "key.jl": [("values(comps[i].terms)", "values(comps[i])"),
                   ("for (mo, c) in comps[i].terms))", "for (mo, c) in comps[i]))")],
        "codegen.jl": [("collect(p.terms)", "collect(p)")],
        "print.jl": [("for (mo, c) in p.terms]", "for (mo, c) in p]"),
                     ("foldl(gcd, values(p.terms))", "foldl(gcd, values(p))"),
                     ("minimum(mo -> mexp(mo, v), keys(p.terms))", "minimum(mo -> mexp(mo, v), keys(p))"),
                     ("for mo in keys(p.terms) for t in mo", "for mo in keys(p) for t in mo"),
                     ("for (mo, x) in p.terms))", "for (mo, x) in p))")],
    },
}
for v, patches in variants.items():
    d = f"{OUT}/{v}"
    os.makedirs(d, exist_ok=True)
    for f in os.listdir(d):
        if f.endswith(".jl"):
            os.remove(f"{d}/{f}")
    for f, pairs in patches.items():
        open(f"{d}/{f}", "w", encoding="utf-8").write(patch(read(f), pairs, f"{v}/{f}"))
    if v == "B":
        b = open(f"{SP}/phase45/research5_interrupted/law/variants/B/poly.jl", encoding="utf-8").read()
        open(f"{d}/poly.jl", "w", encoding="utf-8").write(b)
    incl = ", ".join(f'"{f}"' for f in FILES)
    open(f"{d}/Slice{v}.jl", "w", encoding="utf-8").write(f'''# ABOUTME: Law-type variant {v} of the slice as module Slice{v}: the slice's files, with this
# ABOUTME: directory's patched copies (made by ../mk.py) in place of the originals.
module Slice{v}
for file in ({incl})
    local p = joinpath(@__DIR__, file)
    include(isfile(p) ? p : joinpath("{SRC}", file))
end
end # module
''')
    print(v, sorted(f for f in os.listdir(d)))
