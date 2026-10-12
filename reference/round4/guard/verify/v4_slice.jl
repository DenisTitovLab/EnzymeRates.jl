# ABOUTME: Slice-side checks: the slice's verdict and law on RE-turnover mechanisms, and a census
# ABOUTME: of v ≡ 0 and D = 0-on-a-face laws over enumerated parents, their children and mutants.
const SP = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad"
include(joinpath(SP, "proto/common/Population.jl"))
include(joinpath(SP, "slice/src/Slice.jl"))
const S = Slice
using Serialization, Random, Printf

errline(e) = replace(first(sprint(showerror, e), 200), '\n' => ' ')
mets(xs) = S.Met[S.Met(x, false) for x in xs]
Fm(conf, bound...) = S.RawForm(conf, mets(collect(bound)), Symbol[], Symbol[])
st(f, t, c, r, re) = S.RawStep(f, t, mets(c), mets(r), re)

const HAND = [
    ("A RE turnover beside an SS route", [:S], [:P], [
        st(Fm(:E), Fm(:E, :S), [:S], [], true), st(Fm(:E, :S), Fm(:E, :P), [], [], true),
        st(Fm(:E), Fm(:E, :P), [:P], [], true), st(Fm(:E, :S), Fm(:F, :P), [], [], false),
        st(Fm(:F, :P), Fm(:E), [], [:P], false)], [1, 2, 3, 4, 5]),
    ("B mut/expand_bibi/453/flip2", [:A, :B], [:P, :Q], [
        st(Fm(:E), Fm(:E, :A), [:A], [], true), st(Fm(:E), Fm(:E, :P), [:P], [], true),
        st(Fm(:E), Fm(:E, :Q), [:Q], [], true), st(Fm(:E, :P), Fm(:E, :P, :Q), [:Q], [], true),
        st(Fm(:E, :A), Fm(:E, :P, :Q), [:B], [], true), st(Fm(:E, :Q), Fm(:E, :P, :Q), [:P], [], false)],
     [1, 2, 3, 3, 4, 5]),
    ("C all-RE uni-uni (turnover only through RE)", [:S], [:P], [
        st(Fm(:E), Fm(:E, :S), [:S], [], true), st(Fm(:E, :S), Fm(:E, :P), [], [], true),
        st(Fm(:E, :P), Fm(:E), [], [:P], true)], [1, 2, 3]),
]

function part1(rng)
    println("## the slice on RE-turnover mechanisms")
    for (name, subs, prods, steps, group) in HAND
        m = S.Mechanism(subs, prods, steps, group)
        cb = try S.check_boundary(m); "passes" catch e; errline(e) end
        d = try S.derive(m) catch e; println("  $name: derive error ", errline(e)); continue end
        L = d.law
        names = S.fitted_params(m)
        vals = String[]
        for _ in 1:3
            params = (; (x => 10.0^(2rand(rng) - 1) for x in names)..., Keq = 10.0^(2rand(rng) - 1),
                      E_total = 1.0)
            concs = (; (x => 10.0^(2rand(rng) - 1) for x in S.metabolites(m))...)
            push!(vals, @sprintf("%.6g", S.rate_equation(m, concs, params)))
        end
        println("  $name: check_boundary $cb; admit $(something(S.admit(m), "nothing")); Q ≡ 0: $(iszero(L.Q)); v at 3 random points: ",
                join(vals, ", "))
        println("    printed: ", strip(last(split(d.printed, '\n'; keepempty = false))))
    end
end

"""`(nzero, dfwd, drev)` of a derived law: N ≡ 0, D vanishing on the forward / reverse face
(every other variable at a random positive value; coefficients are positive)."""
function law_flags(m, rng)
    L = S.rate_law(m, S.solve(m))
    V = L.vars
    x0 = [10.0^(2rand(rng) - 1) for _ in V.names]
    face(zs) = (x = copy(x0); foreach(z -> x[V.index[z]] = 0.0, zs); x)
    if L isa S.Law
        Dz(x) = S.evaluate(L.D, x) == 0
        return (iszero(L.Q), Dz(face(m.products)), Dz(face(m.substrates)))
    end
    QI, DI = S.inactive_terms(L)
    n = L.n
    function Dc(x)
        WA = prod((S.evaluate(s[1], x)^s[3] for s in L.sites); init = 1.0)
        WI = prod((S.evaluate(s[2], x)^s[3] for s in L.sites); init = 1.0)
        S.evaluate(L.fI, x)^n * S.evaluate(L.active.D, x)^n * WA +
            x[V.index[:L]] * S.evaluate(L.fA, x)^n * S.evaluate(DI, x)^n * WI
    end
    (iszero(L.active.Q) && iszero(QI), Dc(face(m.products)) == 0, Dc(face(m.substrates)) == 0)
end

function census(label, ms, rng)
    cnt = Dict{String, Int}()
    ex = Dict{String, Any}()
    for m in ms
        k = try
            z, f, r = law_flags(m, rng)
            z ? "N≡0" : f ? "D=0 fwd" : r ? "D=0 rev" : "ok"
        catch e
            "error: " * first(errline(e), 40)
        end
        cnt[k] = get(cnt, k, 0) + 1
        get!(ex, k, m)
    end
    println("  $label ($(length(ms))): ", sort!(collect(cnt)))
    ex
end

const PLAIN_ENUM = ["init_uniuni", "init_unibi", "init_bibi", "init_bibi_pp", "expand_bibi", "p1_enum_inh"]
const ALLO_ENUM = ["allo_uniuni", "allo_bibi", "allo_uniuni_inh", "allo_unibi", "allo_bibi_pp", "allo_terbi"]

"""A slice mechanism from the probe's plain-data spec."""
function from_spec(sp)
    met(x) = S.Met(x[1], x[2])
    form(f) = S.RawForm(f[1], S.Met[met(x) for x in f[2]], collect(Symbol, f[3]), collect(Symbol, f[4]))
    steps = [S.RawStep(form(s[1]), form(s[2]), S.Met[met(x) for x in s[3]], S.Met[met(x) for x in s[4]], s[5])
             for s in sp.steps]
    S.Mechanism(sp.substrates, sp.products, steps, collect(Int, sp.group); tags = collect(Symbol, sp.tags),
                multiplicity = sp.multiplicity, sites = S.Site[S.Site(collect(Symbol, s[1]), collect(Symbol, s[2]), s[3])
                                                                for s in sp.sites])
end

function main()
    rng = Xoshiro(5)
    part1(rng)
    println("\n## census of N ≡ 0 and D = 0 on a face")
    parents = [c.mech for set in [PLAIN_ENUM; ALLO_ENUM]
               for c in deserialize(joinpath(SP, "slice/core/data/$set.jls"))]
    census("enumerated parents", parents, rng)
    kids = S.Mechanism[]
    for m in parents
        for f in (S.flip_children, S.split_children)
            try append!(kids, f(m)) catch e; println("  move error: ", errline(e)) end
        end
    end
    unique!(kids)
    census("their flip and split children", kids, rng)
    specs = deserialize(joinpath(SP, "round4/guard/data/specs_full.jls"))
    for tag in ("/del_cat", "/to_ss", "/to_re", "/flip2", "/all_re")
        sel = [r for r in specs if endswith(r.key, tag)]
        ms = S.Mechanism[]
        nerr = 0
        for r in sel
            try push!(ms, from_spec(r.spec)) catch; nerr += 1 end
        end
        ex = census("mutants *$tag (slice ctor errors $nerr)", ms, rng)
        if tag == "/del_cat" && haskey(ex, "N≡0")
            m = ex["N≡0"]
            println("    example N≡0: ", join([S.step_text(m, e) for e in m.edges], "; "))
        end
    end
    sets = ["specs", "handwritten"]
    census("specs + handwritten", [c.mech for set in sets for c in deserialize(joinpath(SP, "slice/core/data/$set.jls"))], rng)
end

main()
