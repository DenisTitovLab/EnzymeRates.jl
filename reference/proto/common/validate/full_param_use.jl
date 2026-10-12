# ABOUTME: Which Full-mode parameters never reach the Full-mode polynomials, per population.
# ABOUTME: Stdlib only.
include(joinpath(@__DIR__, "..", "Population.jl"))
syms(p) = Set(s for (mono, _) in p for (s, _) in mono)
for name in Pop.populations()
    plain = [m for m in Pop.load_population(name) if !m.allosteric && m.ok]
    isempty(plain) && continue
    unused = 0; unused_open = 0; open = 0; examples = String[]
    for m in plain
        used = union(syms(m.pkg_num), syms(m.pkg_den))
        miss = [p for p in m.full_params if p !== :E_total && p ∉ used]
        o = Pop.re_cycle_closure(m)[2] > 0
        open += o
        isempty(miss) && continue
        unused += 1; unused_open += o
        length(examples) < 3 && push!(examples, "$(m.id) $(miss) open_re=$o pass2=$(m.pass2_rename)")
    end
    println("$name: plain $(length(plain)); open RE cycle $open; a Full parameter absent from ",
            "the Full polynomials $unused (of them with an open RE cycle $unused_open)")
    foreach(e -> println("    ", e), examples)
end
