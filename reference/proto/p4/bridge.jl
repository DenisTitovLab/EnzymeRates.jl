# ABOUTME: P4 part 4: Cleland constants as coefficient ratios checked against textbook rate equations by
# ABOUTME: exact identity, Haldane relations, and the Segel transcriptions of the fixtures file.
include(joinpath(@__DIR__, "P4.jl"))
include(joinpath(@__DIR__, "Build.jl"))
using .Poly: Mono, mono, P, term, nterms
const PI = P{Int}
const RF = P4.RF

byid(id) = (pi = split(id, "/"); Pop.load_population(pi[1])[parse(Int, pi[2])])
rfvar(L, v) = RF(term(PI, mono(v => 1)))
concs(L) = (; (x => rfvar(L, v) for (x, v) in L.conc)...)
lawrf(L) = RF(L.num) / RF(L.den)

include(joinpath(@__DIR__, "bridge_text.jl"))

re_random_m = Build.mech("re_random_bibi", [:A, :B], [:P, :Q],
    ["E + A ⇌ E(A)", "E + B ⇌ E(B)", "E(A) + B ⇌ E(A,B)", "E(B) + A ⇌ E(A,B)",
     "E(A,B) <--> E(P,Q)", "E(P,Q) ⇌ E(Q) + P", "E(P,Q) ⇌ E(P) + Q", "E(P) ⇌ E + P",
     "E(Q) ⇌ E + Q"])
cases = [("uni-uni (Uni-Uni spec)", byid("specs/1"), :uni, [uniuni, uniuni_v1v2]),
         ("uni-uni (Segel IX-8, with EA ⇌ EP)", byid("specs/2"), :uni, [uniuni, uniuni_v1v2]),
         ("RE uni-uni", byid("specs/18"), :uni, [uniuni]),
         ("ordered bi-bi (Segel IX-87)", byid("specs/5"), :bi, [ordered]),
         ("Theorell-Chance (Segel IX-122)", byid("specs/6"), :bi, [theorell_chance]),
         ("ping-pong bi-bi (Segel IX-140)", byid("specs/7"), :bi, [pingpong]),
         ("RE random bi-bi (built)", re_random_m, :bi, [re_random]),
         ("ordered bi-bi, RE bindings (spec)", byid("specs/19"), :bi, [ordered])]

function fmt(L, r::RF)
    a, b = P4.reduce_ratio(r.n, r.d)
    b == one(PI) ? P4.fmt_factored(a, L.vars) : P4.fmt_factored(a, L.vars) * " / (" * P4.fmt_factored(b, L.vars) * ")"
end

println("== Part 4a: Cleland constants as coefficient ratios vs textbook equations (exact identity) ==")
for (label, m, kind, texts) in cases, (mode, cyc) in (("AllParameters", :re), ("IndependentParameters", :all))
    L = P4.law(m; cycles = cyc)
    k = kind === :uni ? P4.cleland_uniuni(L; A = m.substrates[1], P = m.products[1]) :
        P4.cleland_bibi(L; A = m.substrates[1], B = m.substrates[2], P = m.products[1], Q = m.products[2])
    undefined = [n for (n, v) in pairs(k) if v === nothing]
    x = concs(L)
    x = kind === :uni ? (A = x[m.substrates[1]], P = x[m.products[1]]) :
        (A = x[m.substrates[1]], B = x[m.substrates[2]], P = x[m.products[1]], Q = x[m.products[2]])
    D = P4.cmap(L, L.den); N = P4.cmap(L, L.num)
    supp = join(sort([P4.fmt_mono(mo, L.vars) for mo in keys(D)]; by = s -> (count(==('*'), s), s)), " ")
    print(rpad(label, 40), rpad(mode, 23), "D support {", supp, "}; undefined ", undefined, "; ")
    for t in texts
        ok = try P4.same(t(k, x), lawrf(L)) catch e; "error: " * sprint(showerror, e)[1:min(end, 60)] end
        print(nameof(t), " identity: ", ok, "; ")
    end
    hs = [n for (n, f) in haldanes if (try P4.same(f(k), k.Keq) catch; false end)]
    println("Haldane identities: ", hs)
    if mode == "AllParameters" && !(label in ("RE uni-uni", "ordered bi-bi, RE bindings (spec)"))
        for (n, v) in pairs(k)
            v === nothing || n === :Keq || println("      ", n, " = ", fmt(L, v))
        end
    end
end

# ─── Part 4b: Segel transcriptions of the fixtures file ────────────────
const FIX = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/er_oracle/test/mechanism_definitions_for_test_enzyme_derivation.jl"
ast = Meta.parseall(read(FIX, String))

walk(f, ex) = (f(ex); ex isa Expr && foreach(a -> walk(f, a), ex.args))
"""Each `let` block's (spec name, macro block, function definitions, analytical_rate_fn expr)."""
function specs_from_ast(ast)
    out = []
    for blk in ast.args
        (blk isa Expr && blk.head === :let) || continue
        mac, funs, name, fn = nothing, Expr[], nothing, nothing
        walk(blk) do ex
            ex isa Expr || return
            ex.head === :macrocall && ex.args[1] === Symbol("@_testhelper_enzyme_mechanism_src") &&
                (mac = ex.args[end])
            ex.head === :function && push!(funs, ex)
            if ex.head === :call && ex.args[1] === :_testhelper_MechanismTestSpec
                for a in ex.args[2:end]
                    a isa Expr && a.head === :kw || continue
                    a.args[1] === :name && (name = a.args[2])
                    a.args[1] === :analytical_rate_fn && (fn = a.args[2])
                end
            end
        end
        name === nothing || push!(out, (name, mac, funs, fn))
    end
    out
end

"""Source steps (flat, in written order) of a macro block: (lhs, rhs, re) with sides as
(conformation, bound (name, inh) pairs, metabolites (name, inh) pairs)."""
function source_steps(mac)
    regs = Symbol[]
    lines = Any[]
    for a in mac.args
        a isa Expr || continue
        items = a.head === :tuple ? a.args : [a]
        for it in items
            it isa Expr && it.head === :call && it.args[1] === :(:) || continue
            it.args[2] === :regulators && push!(regs, it.args[3])
            it.args[2] === :steps && (lines = it.args[3].args)
        end
        a.head === :tuple && a.args[1] isa Expr && a.args[1].args[2] === :regulators &&
            append!(regs, a.args[2:end])
    end
    met(x) = x isa Symbol ? (x, x in regs) : (x.args[1], true)          # x::Inh
    form(f) = f isa Symbol ? (f, Tuple{Symbol, Bool}[]) : (f.args[1], sort(met.(f.args[2:end])))
    side(s) = s isa Expr && s.head === :call && s.args[1] === :+ ?
              (form(s.args[2])..., sort(met.(s.args[3:end]))) : (form(s)..., Tuple{Symbol, Bool}[])
    steps = []
    for l in lines
        l isa Expr || continue
        for st in (l.head === :tuple ? l.args : [l])
            push!(steps, (side(st.args[2]), side(st.args[3]), st.args[1] === :⇌))
        end
    end
    steps
end

pside(m, f, ms) = (m.species[f].conformation, sort([(b.name, b.role === :Inh) for b in m.species[f].bound]),
                   sort([(x.name, x.role === :Inh) for x in ms]))

"""Positional constants (`K{i}`, `k{i}f`, `k{i}r`) of the source steps as RF values from the law's
step monomials, the bridge's swap rule disagreements, and unmatched steps."""
function positional(m, src, L)
    names, vals, mism, unmatched = Symbol[], RF[], Int[], Int[]
    used = falses(length(m.steps))
    for (i, (l, r, re)) in enumerate(src)
        j = findfirst(jj -> !used[jj] && m.steps[jj].re == re &&
                      (pside(m, m.steps[jj].from, m.steps[jj].consumed) in (l, r)) &&
                      (pside(m, m.steps[jj].to, m.steps[jj].released) in (l, r)), eachindex(m.steps))
        j === nothing && (push!(unmatched, i); continue)
        used[j] = true
        s = m.steps[j]
        same_dir = pside(m, s.from, s.consumed) == l
        mo = L.stepmono[j]
        if re
            push!(names, Symbol("K", i)); push!(vals, RF(term(PI, mo[1])))
        else
            f, b = same_dir ? (mo[1], mo[2]) : (mo[2], mo[1])
            push!(names, Symbol("k", i, "f"), Symbol("k", i, "r"))
            push!(vals, RF(term(PI, f)), RF(term(PI, b)))
            bridge_swap = length(s.consumed) == 1 && isempty(s.released) && s.consumed[1].role === :P
            bridge_swap == !same_dir || push!(mism, i)
        end
    end
    NamedTuple{Tuple(names)}(Tuple(vals)), mism, unmatched
end

println("\n== Part 4b: Segel and textbook transcriptions of the fixtures file vs derived laws (exact identity) ==")
specs = Pop.load_population("specs")
sandbox = Module(:SegelSandbox)
nall = 0; npass = Dict("AllParameters" => 0, "IndependentParameters" => 0)
for (name, mac, funs, fn) in specs_from_ast(ast)
    (mac === nothing || fn === nothing) && continue
    k = findfirst(mm -> mm.source == "spec:" * name, specs)
    (k === nothing || specs[k].allosteric) && continue
    m = specs[k]
    global nall += 1
    for f in funs
        Core.eval(sandbox, f)
    end
    rate = Core.eval(sandbox, fn)
    src = source_steps(mac)
    print(rpad(name, 40))
    for (mode, cyc) in (("AllParameters", :re), ("IndependentParameters", :all))
        L = P4.law(m; cycles = cyc)
        pos, mism, unm = positional(m, src, L)
        p = merge(pos, (Et = RF(1),))
        c = (; (x => rfvar(L, v) for (x, v) in L.conc)...)
        ok = try
            P4.same(Base.invokelatest(rate, p, c), lawrf(L))
        catch e
            "error: " * first(sprint(showerror, e), 80)
        end
        ok === true && (npass[mode] += 1)
        # discrimination: doubling any one positional constant the transcription uses must break it
        used = Set{Symbol}()
        walk(ex -> ex isa Symbol ? push!(used, ex) : ex isa QuoteNode && ex.value isa Symbol && push!(used, ex.value), Expr(:block, fn, funs...))
        keysused = [kk for kk in keys(pos) if kk in used]
        broke = count(kk -> !P4.same(Base.invokelatest(rate, merge(p, NamedTuple{(kk,)}((p[kk] * 2,))), c),
                                     lawrf(L)), keysused)
        print(mode, " ", ok, " (doubling each of $(length(keysused)) used constants breaks it: $(broke))",
              isempty(mism) ? "" : " bridge-swap mismatch $(mism)", isempty(unm) ? "" : " unmatched $(unm)", "; ")
    end
    println()
end
println("transcriptions checked: ", nall, "; identical: ", npass)
