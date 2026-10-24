# ABOUTME: Research library for phase 5's unit rule and printer (R5.4; engine.md Q2-Q7, Q12, Q13):
# ABOUTME: census text to slice mechanisms, the normalizer rules, a printer, value and unit checks.

const SP = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad"
isdefined(Main, :Slice) || include(joinpath(SP, "slice", "src", "Slice.jl"))
const S = Slice
using Printf, Random, Serialization
setprecision(BigFloat, 512)

const QQ = Rational{BigInt}

# ─── Census neutral text → slice mechanism ─────────────────────────────────────────────

"""A census label as a slice metabolite: `Xmod` is reactant X bound as a modifier, a name
that is no reactant is a modifier, any other name is the reactant itself."""
function parse_label(lab::AbstractString, reactants)
    if endswith(lab, "mod") && Symbol(lab[1:(end - 3)]) in reactants
        return S.Met(Symbol(lab[1:(end - 3)]), true)
    end
    x = Symbol(lab)
    S.Met(x, !(x in reactants))
end

"""A census form `Conf(labels; residual=plusA_minusP)` as a slice raw form."""
function parse_form(t::AbstractString, reactants)
    i = findfirst('(', t)
    conf, inner = Symbol(t[1:(i - 1)]), t[(i + 1):(end - 1)]
    res = ""
    occursin("; residual=", inner) && ((inner, res) = split(inner, "; residual="))
    bound = isempty(inner) ? S.Met[] : [parse_label(l, reactants) for l in split(inner, ",")]
    added, subtracted = Symbol[], Symbol[]
    for r in (isempty(res) ? SubString{String}[] : split(res, "_"))
        startswith(r, "plus") ? push!(added, Symbol(r[5:end])) :
            push!(subtracted, Symbol(r[6:end]))
    end
    S.RawForm(conf, bound, added, subtracted)
end

"""A census side `form + label + …` as `(form, labels)`."""
function parse_side(t::AbstractString, reactants)
    parts = split(t, " + ")
    parse_form(parts[1], reactants), S.Met[parse_label(l, reactants) for l in parts[2:end]]
end

"""
    parse_neutral(text, substrates, products) -> slice Mechanism

The mechanism of a census neutral text (`phase3/census/README.md`, "Neutral text").
"""
function parse_neutral(text::AbstractString, subs::Vector{Symbol}, prods::Vector{Symbol})
    reactants = Set([subs; prods])
    allo = startswith(text, "allosteric")
    mult, body, sitetext = 1, "", "none"
    if allo
        rest = text[length("allosteric n=") + 1:end]
        sp = findfirst(' ', rest)
        mult = parse(Int, rest[1:(sp - 1)])
        body, sitetext = split(rest[(sp + 1):end], " sites "; limit = 2)
    else
        body = text[length("plain ") + 1:end]
    end
    steps, group, tags = S.RawStep[], Int[], Symbol[]
    for (g, gt) in enumerate(split(body, " / "))
        tag = ""
        if !endswith(gt, "}")
            gt, tag = split(gt, "}:")
            gt *= "}"
        end
        allo && push!(tags, Symbol(tag))
        for st in split(gt[2:(end - 1)], " | ")
            flag, core = st[(end - 1):end], st[1:(end - 3)]
            a, b = split(core, " <-> ")
            (fa, ca), (fb, cb) = parse_side(a, reactants), parse_side(b, reactants)
            push!(steps, S.RawStep(fa, fb, ca, cb, flag == "RE"))
            push!(group, g)
        end
    end
    sites = S.Site[]
    if allo && sitetext != "none"
        for s in split(sitetext, " / ")
            mm = match(r"^site\((\d+): (.*)\)$", s)
            ligs = [split(l, "::") for l in split(mm[2], ", ")]
            push!(sites, S.Site(Symbol.(first.(ligs)), Symbol.(last.(ligs)),
                                parse(Int, mm[1])))
        end
    end
    allo ? S.Mechanism(subs, prods, steps, group; tags, multiplicity = mult, sites) :
           S.Mechanism(subs, prods, steps, group)
end

"""Roles field `S=A,B;P=P,Q;M=I` of a census law row as `(substrates, products)`."""
function parse_roles(r::AbstractString)
    f = Dict(split(x, "=")[1] => split(x, "=")[2] for x in split(r, ";"))
    sym(s) = isempty(s) ? Symbol[] : Symbol.(split(s, ","))
    sym(f["S"]), sym(f["P"])
end

# ─── Exact affine algebra ──────────────────────────────────────────────────────────────

"""Row-reduce `M` in place; returns the pivot columns."""
function rref!(M::Matrix{QQ})
    r, piv = 1, Int[]
    for c in axes(M, 2)
        r > size(M, 1) && break
        p = findfirst(i -> !iszero(M[i, c]), r:size(M, 1))
        p === nothing && continue
        p += r - 1
        M[[r, p], :] = M[[p, r], :]
        M[r, :] ./= M[r, c]
        for i in axes(M, 1)
            i != r && !iszero(M[i, c]) && (M[i, :] .-= M[i, c] .* M[r, :])
        end
        push!(piv, c)
        r += 1
    end
    piv
end

exact_rank(cols::Vector{Vector{QQ}}) = isempty(cols) ? 0 : length(rref!(hcat(cols...)))

"""Exponent vector of a concentration monomial over the first `nc` variables."""
expvec(mo::S.Mono, nc) = [Int(S.mexp(mo, i)) for i in 1:nc]

"""
    affine(monos, nc) -> (basis, λ)

The first affinely independent monomials of `monos` (in the given order) and the unique
weights `λ` over them with `Σ λ = 1` and `Σ λ·a = 0`; `λ === nothing` when no weights
exist (zero lies outside the affine hull of the exponent vectors).
"""
function affine(monos::Vector{S.Mono}, nc)
    vec(mo) = QQ[1; expvec(mo, nc)]
    basis, B = S.Mono[], Vector{QQ}[]
    for mo in monos
        v = vec(mo)
        exact_rank([B; [v]]) > length(B) && (push!(basis, mo); push!(B, v))
    end
    k = length(B)
    M = hcat(B..., QQ[1; zeros(QQ, nc)])
    piv = rref!(M)
    (k + 1) in piv && return basis, nothing
    λ = zeros(QQ, k)
    for (i, c) in enumerate(piv)
        λ[c] = M[i, end]
    end
    basis, λ
end

# ─── Normalizers ───────────────────────────────────────────────────────────────────────

"""
A normalizer `N = Π_m c(m)^λ_m · x^shift` over denominator monomials `m` (`c(m)` the
coefficient polynomial of the concentration monomial `m`): `kind` names the rule that chose
it, `shift` is the concentration part that stays in `N` (the ratio form's first term).
"""
struct Norm
    kind::Symbol
    lambda::Vector{Pair{S.Mono, QQ}}
    shift::S.Mono
end

"""Indices (into `L.vars`) of the substrates of `m`."""
subidx(m, L) = [L.vars.index[x] for x in m.substrates]

"""
    orderkey(mo, L, nc, subs; order) -> sort key

The term order of the unit rule (spec lines 502-503): total degree, then the count of
non-substrate metabolites, then metabolite order (`metabolites(m)`, larger exponent of an
earlier metabolite first). `order` picks a variant: `:spec` as above; `:nsdeg` counts
non-substrates with multiplicity; `:text` breaks the last tie by C4's monomial text (names
sorted as strings, `A*B^2`); `:nons` drops the non-substrate key; `:sqfree` puts fewer
copies of one metabolite first (`A·B` before `A²`); `:subfirst` orders by the non-substrate
count before the degree.
"""
function orderkey(mo::S.Mono, L, nc, subs; order::Symbol = :spec)
    ex = expvec(mo, nc)
    deg = sum(ex)
    ns = order === :nsdeg ? sum(ex[i] for i in 1:nc if !(i in subs); init = 0) :
         count(i -> ex[i] != 0 && !(i in subs), 1:nc)
    order === :nons && (ns = 0)
    order === :sqfree && return (deg, ns, maximum(ex; init = 0), Tuple(-ex), "")
    order === :subfirst && return (ns, deg, 0, Tuple(-ex), "")
    order === :text || return (deg, ns, 0, Tuple(-ex), "")
    names = sort([(string(L.vars.names[i]), ex[i]) for i in 1:nc if ex[i] != 0])
    (deg, ns, 0, (), join((e == 1 ? n : "$n^$e" for (n, e) in names), "*"))
end

const ORDERS = (:text, :nsdeg, :nons, :sqfree, :subfirst)

"""The unit rule of the spec (lines 490-509): the constant term, else the virtual constant
term on the first affinely independent terms in `orderkey` order (fractional weights
allowed), else the ratio form over the first term."""
function spec_norm(m, L; order::Symbol = :spec)
    D, nc = S.split_conc(L.D, L.vars), length(S.metabolites(m))
    haskey(D, S.Mono()) && return Norm(:constant, [S.Mono() => one(QQ)], S.Mono())
    subs = subidx(m, L)
    monos = sort!(collect(keys(D)); by = mo -> orderkey(mo, L, nc, subs; order))
    basis, λ = affine(monos, nc)
    λ === nothing && return Norm(:ratio, [monos[1] => one(QQ)], monos[1])
    kind = all(isinteger, λ) ? :virtual_int : :virtual_frac
    Norm(kind, [b => l for (b, l) in zip(basis, λ) if !iszero(l)], S.Mono())
end

"""Cleland's `N = Π_X c(ΠS/X) / c(ΠS)^(n−1)` when every factor is a term of the
denominator, else `nothing`."""
function cleland_norm(m, L)
    D, idx = S.split_conc(L.D, L.vars), subidx(m, L)
    whole = S.mono((i => 1 for i in idx)...)
    parts = [S.mmul(whole, S.mono(i => -1)) for i in idx]
    all(mo -> haskey(D, mo), [parts; [whole]]) || return nothing
    n = length(idx)
    lam = [[mo => one(QQ) for mo in parts]; [whole => QQ(1 - n)]]
    Norm(:cleland, filter(p -> !iszero(p.second), lam), S.Mono())
end

"""The slice's normalizer (`slice/src/normalize.jl`) as a `Norm`; its `:coefficient`
fallback divides by the first term's coefficient alone."""
function slice_norm(m, L)
    rule, lam = S.normalizer(m, L)
    Norm(rule, [mo => QQ(e) for (mo, e) in lam], S.Mono())
end

"""`log N` without its concentration part, at the variable values `vals` (BigFloat)."""
logN(nz::Norm, D, vals) = sum(BigFloat(l) * log(S.evaluate(D[mo], vals)) for (mo, l) in nz.lambda)

"""Whether two normalizers give one `N` (relative 1e-100 on `log N` at the given points)."""
function same_N(a::Norm, b::Norm, D, points)
    a.shift == b.shift || return false
    all(points) do vals
        x, y = logN(a, D, vals), logN(b, D, vals)
        abs(x - y) <= big(10.0)^-100 * (1 + abs(x))
    end
end

lamtext(nz::Norm, names) = join(("$(isempty(mo) ? "1" : join(("$(names[v])$(e == 1 ? "" : "^$e")" for (v, e) in mo), "*")):$(l)"
                                 for (mo, l) in nz.lambda), " ")

# ─── Law values ────────────────────────────────────────────────────────────────────────

"""BigFloat values over `vars` from `params` and `concs` (NamedTuples)."""
varvals(vars, params, concs) = BigFloat[vars.conc[i] ? concs[n] : params[n]
                                        for (i, n) in enumerate(vars.names)]

"""`v/E_total` of a plain law at `vals`."""
lawval(L::S.Law, vals) = S.evaluate(L.H, vals) * S.evaluate(L.Q, vals) / S.evaluate(L.D, vals)

"""`v/E_total` of an MWC law (`slice/src/allosteric.jl:10-11`) at `vals` over `A.vars`."""
function lawval(A::S.AlloLaw, vals)
    ev(p) = S.evaluate(p, vals)
    QI, DI = S.inactive_terms(A)
    n, Lv = A.n, vals[A.vars.index[:L]]
    qa, da, qi, di, fa, fi = ev(A.active.Q), ev(A.active.D), ev(QI), ev(DI), ev(A.fA), ev(A.fI)
    wa = prod((ev(s[1])^s[3] for s in A.sites); init = big(1.0))
    wi = prod((ev(s[2])^s[3] for s in A.sites); init = big(1.0))
    num = fi^n * qa * da^(n - 1) * wa + Lv * fa^n * qi * di^(n - 1) * wi
    den = fi^n * da^n * wa + Lv * fa^n * di^n * wi
    ev(A.H) * num / den
end

"""Random parameters (log-uniform 1e-2..1e2) and concentrations (1e-1..1e1) of `m`."""
function random_point(m, rng)
    lu(a, b) = big(10.0)^(a + (b - a) * rand(rng))
    params = (; (n => lu(-2, 2) for n in S.parameters(m))...)
    concs = (; (x => lu(-1, 1) for x in S.metabolites(m))...)
    params, concs
end

# ─── Units ─────────────────────────────────────────────────────────────────────────────

const Unit = Tuple{QQ, QQ}          # (concentration exponent, time exponent)

"""Units of every name a printed law of `m` reads: concentrations, `E_total`, `Keq`, the
independents (`slice/engine/lib.jl:35-49`), regulator constants and `L`."""
function name_units(m, sol)
    u = Dict{Symbol, Unit}()
    for x in S.metabolites(m)
        u[x] = (1, 0)
    end
    u[:E_total] = (1, 0)
    u[:Keq] = (length(m.products) - length(m.substrates), 0)
    fwd, names = S.orient(m), S.column_names(m, sol)
    for c in sol.indep
        col, o = sol.cols[c], S.oriented(m, sol.cols[c].edge, fwd)
        u[names[c]] = !o.re ? (-(col.slot == 1 ? S.nibsum(o.cons) : S.nibsum(o.rel)), -1) :
            S.is_binding(o) ? (1, 0) : (S.nibsum(o.rel) - S.nibsum(o.cons), 0)
    end
    if S.is_allosteric(m)
        for (_, ns) in S.regulator_names(m), n in ns
            u[n] = (1, 0)
        end
        u[:L] = (0, 0)
    end
    u
end

# ─── Printer under a given normalizer ──────────────────────────────────────────────────

const SUP = Dict(zip("0123456789-", "⁰¹²³⁴⁵⁶⁷⁸⁹⁻"))
supers(e::Integer) = e == 1 ? "" : String([SUP[c] for c in string(e)])

"""A coefficient's name for the concentration monomial `mo`: `A_B`, `A²_B`, `P_B⁻¹` (positive
exponents first, each part in metabolite order), `0` for none."""
function conc_label(mo::S.Mono, names)
    isempty(mo) && return "0"
    parts = [(e < 0, v, Int(e)) for (v, e) in mo]
    join((string(names[v], supers(e)) for (_, v, e) in sort(parts)), "_")
end

"""Add `p^e` to the factors printed apart, merging an equal polynomial."""
function addapart!(apart, p, e)
    i = findfirst(x -> x[1] == p, apart)
    i === nothing ? push!(apart, (p, e)) : (apart[i] = (p, apart[i][2] + e))
    filter!(x -> !iszero(x[2]), apart)
end

"""
    printed(m, L, nz, polys) -> (; N, groups)

The polynomials `polys` of law `L` as printed over the normalizer `nz`: its integer and
monomial contents divide every coefficient, its polynomial factors print apart as `N` (their
text, else `nothing`). `groups[j]` lists each concentration monomial of `polys[j]` in the
unit rule's order, shifted by `nz.shift`, as `(printed monomial, coefficient)` with the
coefficient `(c, mo, rest, over)`: `c·mo` times the polynomial `rest` (`nothing` for 1),
over `N` when `over`.
"""
function printed(m, L, nz::Norm, polys)
    D, names, nc = S.split_conc(L.D, L.vars), L.vars.names, length(S.metabolites(m))
    k, mu, apart = one(QQ), S.Mono(), Tuple{Any, QQ}[]
    for (mo, e) in nz.lambda
        c, u, rest = S.content(D[mo])
        if isinteger(e)
            k *= QQ(c)^Int(e)
        elseif c != 1
            addapart!(apart, S.term(typeof(c), S.Mono(), c), e)
        end
        mu = S.mmul(mu, S.mpow(u, Rational{Int}(e)))
        rest === nothing || addapart!(apart, rest, e)
    end
    paren(p) = "(" * S.sum_text(p, names) * ")"
    N = isempty(apart) ? nothing :
        length(apart) == 1 && apart[1][2] == 1 ? S.sum_text(apart[1][1], names) :
        S.ratio_text([S.power_text(paren(p), Rational{Int}(e)) for (p, e) in apart if e > 0],
                     [S.power_text(paren(p), Rational{Int}(-e)) for (p, e) in apart if e < 0])
    function coef(q)
        c, u, rest = S.content(q)
        whole = length(apart) == 1 && apart[1][2] == 1 && rest !== nothing && apart[1][1] == rest
        (QQ(c) / k, S.mmul(u, S.mpow(mu, -1)), whole ? nothing : rest, !isempty(apart) && !whole)
    end
    subs = subidx(m, L)
    shift = S.mpow(nz.shift, -1)
    groups = map(polys) do p
        gs = sort!([(mo, coef(q)) for (mo, q) in S.split_conc(p, L.vars)];
                   by = g -> orderkey(g[1], L, nc, subs))
        [(S.mmul(mo, shift), x) for (mo, x) in gs]
    end
    (; N, groups)
end

"""A group's term as text: its coefficient times its concentration monomial."""
group_text((cm, (c, mo, rest, over)), names) =
    S.term_text(c, S.mmul(cm, mo), names; rest, over)

group_sum(gs, names) = join((group_text(g, names) for g in gs), " + ")

"""`label` times the concentration monomial `cm`: `d_A * A`, `d_P_B⁻¹ * P / B`, `n_B⁻¹ / B`."""
function label_times(label, cm, names)
    t = S.term_text(1, cm, names)
    startswith(t, "1 / ") ? label * " / " * t[5:end] : label * " * " * t
end

"""
    coefficient_sum(prefix, gs, names) -> (sum, legend)

The groups `gs` as a sum of named coefficients (`prefix_A_B * A * B`, `prefix_0` alone for
the constant monomial, `1` for a coefficient of exactly 1 on it) and one legend line per
named coefficient in the fitted constants.
"""
function coefficient_sum(prefix, gs, names)
    unit((cm, x)) = isempty(cm) && x == (one(QQ), S.Mono(), nothing, false)
    label(cm) = "$(prefix)_$(conc_label(cm, names))"
    terms = [unit(g) ? "1" : isempty(g[1]) ? label(g[1]) : label_times(label(g[1]), g[1], names)
             for g in gs]
    legend = ["$(label(g[1])) = $(S.term_text(g[2][1], g[2][2], names; rest = g[2][3], over = g[2][4]))"
              for g in gs if !unit(g)]
    join(terms, " + "), legend
end

"""A reference line's right side: `refstyle = :ratio` prints a ratio of products
(`a * b / (c * Keq)`, `sqrt`, `^(3/2)`); `:product` prints today's style
(`(1 / Keq) * a * (1 / c)`, non-integer powers as Float64)."""
function ref_rhs(mo::S.Mono, names, refstyle::Symbol)
    refstyle === :ratio && return S.term_text(1, mo, names)
    fs = String[]
    keqi = findfirst(==(:Keq), names)
    ord = sort(collect(mo); by = p -> (p.first != keqi, string(names[p.first])))
    for (v, e) in ord
        x = string(names[v])
        push!(fs, e == 1 ? x : e == -1 ? "(1 / $x)" : isinteger(e) ? "$x ^ $(Int(e))" :
                  "$x ^ $(Float64(e))")
    end
    isempty(fs) ? "1.0" : join(fs, " * ")
end

"""The header of a printed law: the `params` and `concs` destructures and the Wegscheider
and Haldane reference lines (each dependent in the independents and `Keq`)."""
function header_lines(m, sol, vars; refstyle = :ratio)
    cn, names = S.column_names(m, sol), vars.names
    lines = ["(; $(join(S.parameters(m), ", "))) = params",
             "(; $(join(S.metabolites(m), ", "))) = concs"]
    for (label, title) in ((:wegscheider, "Wegscheider"), (:haldane, "Haldane"))
        ds = sort!(filter(d -> d.label === label, sol.deps); by = d -> cn[d.col])
        isempty(ds) || push!(lines, "# $title constraints:")
        for d in ds
            mo = S.mono((vars.index[cn[f]] => e for (f, e) in d.terms)...,
                        vars.index[:Keq] => d.keq)
            push!(lines, "$(cn[d.col]) = $(ref_rhs(mo, names, refstyle))")
        end
    end
    lines
end

"""The `v` line: `layout = :spec` prints `E_total * (num) * (H) / (den)`; `:slice` prints
`E_total * (H) * num / (den)`, parenthesizing a numerator of more than one term."""
v_line(num, nnum, H, den; layout = :spec) =
    layout === :spec ? "v = E_total * ($num) * ($H) / ($den)" :
    "v = E_total * ($H) * $(nnum == 1 ? num : "($num)") / ($den)"

"""
    plain_texts(m, L, sol, nz; layout, refstyle) -> (full, coefficients)

A plain law printed over `nz`: the full law (header, `N = …` when a factor prints apart,
the `v` line in the independents) and the coefficient form (`v` line in named
coefficients, then the legend).
"""
function plain_texts(m, L, sol, nz::Norm; layout = :spec, refstyle = :ratio)
    P = printed(m, L, nz, (L.Q, L.D))
    names = L.vars.names
    q, d = P.groups
    H = S.sum_text(L.H, names)
    lines = header_lines(m, sol, L.vars; refstyle)
    P.N === nothing || push!(lines, "N = $(P.N)")
    push!(lines, v_line(group_sum(q, names), length(q), H, group_sum(d, names); layout))
    (qs, lq), (ds, ld) = coefficient_sum("n", q, names), coefficient_sum("d", d, names)
    coef = [v_line(qs, length(q), H, ds; layout);
            "# The coefficients are functions of the fitted constants, not parameters:"]
    P.N === nothing || push!(coef, "N = $(P.N)")
    join(lines, "\n"), join([coef; lq; ld], "\n")
end

"""Conformation `tag`'s law `L` printed over its own normalizer `nz`: the groups of `Q`,
`D` and the flip form's weight, each over the normalizer's factor printed apart, which is
the extra last variable of `names`, `N_tag`."""
function conformation(m, L::S.Law, tag, nz::Norm)
    names = [L.vars.names; Symbol("N_", tag)]
    over(gs) = [(cm, (c, o ? S.mmul(mo, S.mono(length(names) => -1)) : mo, r, false))
                for (cm, (c, mo, r, o)) in gs]
    P = printed(m, L, nz, (L.Q, L.D, L.trees[1] * L.shift))
    (; tag, names, N = P.N, q = over(P.groups[1]), d = over(P.groups[2]),
     f = over(P.groups[3]))
end

"""
    allo_texts(m, A, sol, norms; layout, refstyle) -> (full, coefficients)

An MWC law printed with each conformation over its own normalizer (`norms`, active then
inactive): `Q_s`, `D_s`, `N_s` when a factor prints apart, the flip weights `F_s` only
when they differ and are not 1, `W_s` the regulatory-site factors, then the `v` line
(`slice/src/print_allo.jl:29-63`).
"""
function allo_texts(m, A::S.AlloLaw, sol, norms; layout = :spec, refstyle = :ratio)
    cs = Any[conformation(m, A.active, "A", norms[1])]
    A.inactive === nothing || push!(cs, conformation(m, A.inactive, "I", norms[2]))
    H = S.sum_text(A.H, A.vars.names)
    function build(coefficients)
        local defs = String[]
        local legend = String[]
        function define!(name, gs, c)
            coefficients || return push!(defs, "$name = $(group_sum(gs, c.names))")
            prefix = string(name[1] == 'Q' ? 'n' : lowercase(name[1]), c.tag)
            text, lines = coefficient_sum(prefix, gs, c.names)
            push!(defs, "$name = $text")
            append!(legend, lines)
        end
        for c in cs
            c.N === nothing || push!(coefficients ? legend : defs, "N_$(c.tag) = $(c.N)")
            isempty(c.q) || define!("Q_$(c.tag)", c.q, c)
            define!("D_$(c.tag)", c.d, c)
        end
        flips = [group_sum(c.f, c.names) for c in cs]
        flips[1] == get(flips, 2, "1") ||
            foreach(((c, t),) -> t == "1" || define!("F_$(c.tag)", c.f, c), zip(cs, flips))
        for (i, tag) in ((1, "A"), (2, "I"))
            ws = [string("(", S.sum_text(s[i], A.vars.names), ")", s[3] == 1 ? "" : "^$(s[3])")
                  for s in A.sites if S.nterms(s[i]) > 1]
            isempty(ws) || push!(defs, "W_$tag = " * join(ws, " * "))
        end
        has(x) = any(l -> startswith(l, x * " = "), defs)
        pw(x, e) = !has(x) || e == 0 ? nothing : e == 1 ? x : "$x^$e"
        product(fs...) = (f = filter(!isnothing, [fs...]); isempty(f) ? "1" : join(f, " * "))
        branch(s, o, e) = product(pw("F_$o", A.n), pw("Q_$s", e), pw("D_$s", A.n - e),
                                  pw("W_$s", 1))
        num = filter(!isnothing, [has("Q_A") ? branch("A", "I", 1) : nothing,
                                  has("Q_I") ? "L * " * branch("I", "A", 1) : nothing])
        rest = branch("I", "A", 0)
        den = branch("A", "I", 0) * " + L" * (rest == "1" ? "" : " * " * rest)
        numt = isempty(num) ? "0" : join(num, " + ")
        defs, legend, v_line(numt, length(num), H, den; layout)
    end
    head = header_lines(m, sol, A.vars; refstyle)
    pdefs, _, v = build(false)
    cdefs, clegend, cv = build(true)
    join([head; pdefs; v], "\n"),
    join([cv; cdefs; "# The coefficients are functions of the fitted constants, not parameters:";
          clegend], "\n")
end

# ─── Reading printed text back: values and units ───────────────────────────────────────

struct UnitMismatch <: Exception
    msg::String
end

"""An exponent literal (`2`, `3 / 2`, `0.5`) as an exact rational."""
expo(x::Integer) = QQ(x)
expo(x::AbstractFloat) = QQ(rationalize(BigInt, x))
function expo(ex::Expr)
    ex.head === :call || error("exponent $ex")
    op, a = ex.args[1], expo.(ex.args[2:end])
    op === :/ ? a[1] / a[2] : op === :- ? (length(a) == 1 ? -a[1] : a[1] - a[2]) :
        op === :* ? prod(a) : op === :+ ? sum(a) : error("exponent $ex")
end

"""The BigFloat value of a printed expression under the bindings `env`."""
function value(ex, env)
    ex isa Symbol && return env[ex]
    ex isa Number && return BigFloat(ex)
    ex.head === :call || error("unexpected $(ex.head) in $ex")
    op, args = ex.args[1], ex.args[2:end]
    op === :^ && return value(args[1], env)^BigFloat(expo(args[2]))
    a = [value(x, env) for x in args]
    op === :+ ? sum(a) : op === :* ? prod(a) : op === :/ ? a[1] / a[2] :
        op === :- ? (length(a) == 1 ? -a[1] : a[1] - a[2]) :
        op === :sqrt ? sqrt(a[1]) : error("unexpected operator $op")
end

"""The units `(concentration, time)` of a printed expression; a sum of unequal units throws
`UnitMismatch`."""
function units(ex, env)
    ex isa Symbol && return env[ex]::Unit
    ex isa Number && return (zero(QQ), zero(QQ))
    op, args = ex.args[1], ex.args[2:end]
    if op === :^
        e, u = expo(args[2]), units(args[1], env)
        return (u[1] * e, u[2] * e)
    end
    us = [units(x, env) for x in args]
    if op === :+ || op === :-
        allequal(us) || throw(UnitMismatch("$(join(us, " vs ")) in $(first(string(ex), 120))"))
        return us[1]
    end
    op === :* && return (sum(first, us), sum(last, us))
    op === :/ && return (us[1][1] - us[2][1], us[1][2] - us[2][2])
    op === :sqrt && return (us[1][1] / 2, us[1][2] / 2)
    error("unexpected operator $op")
end

"""The statements of a printed text, comments and blank lines dropped."""
statements(text) = [Meta.parse(l) for l in split(text, "\n")
                    if !isempty(strip(l)) && !startswith(lstrip(l), "#")]

"""
    run_text(text, f, env, bind) -> env

Every statement of `text` evaluated with `f` (`value` or `units`) in an order that defines
each name before use (the coefficient form prints `v` before its legend). A destructure
`(; a, b) = x` binds each name through `bind(name)`.
"""
function run_text(text, f, env, bind)
    todo = statements(text)
    while !isempty(todo)
        n = length(todo)
        for ex in copy(todo)
            ex.head === :(=) || error("not an assignment: $ex")
            lhs, rhs = ex.args
            if lhs isa Expr && lhs.head === :tuple
                foreach(s -> env[s] = bind(s), lhs.args[1].args)
            else
                syms = Set{Symbol}()
                collect_syms!(syms, rhs)
                all(s -> haskey(env, s), syms) || continue
                env[lhs] = f(rhs, env)
            end
            deleteat!(todo, findfirst(x -> x === ex, todo))
        end
        length(todo) == n && error("undefined names in: $(join(string.(todo), "; "))")
    end
    env
end

collect_syms!(out, ex::Symbol) = push!(out, ex)
collect_syms!(out, ex::Expr) = foreach(a -> collect_syms!(out, a), ex.head === :call ? ex.args[2:end] : ex.args)
collect_syms!(out, ex) = out

"""The top-level terms of a sum expression."""
sumterms(ex) = ex isa Expr && ex.head === :call && ex.args[1] === :+ ? ex.args[2:end] : [ex]

"""
    check_units(text, u) -> (ok, note)

Units of a printed law read back from its text: every top-level term of each denominator
sum (`D_…` definitions and, when the `v` line divides by a sum, that sum) and of each flip
weight (`F_…`) is dimensionless, every `W_…` is dimensionless, and `v` has the units of
`E_total` per time.
"""
function check_units(text, u::Dict{Symbol, Unit})
    env = Dict{Symbol, Unit}(u)
    try
        run_text(text, units, env, s -> u[s])
        env[:v] == (QQ(1), QQ(-1)) || return false, "v units $(env[:v])"
        for ex in statements(text)
            lhs, rhs = ex.args
            lhs isa Symbol || continue
            name = string(lhs)
            dens = name == "v" ? (rhs.args[1] === :/ ? [rhs.args[3]] : []) :
                   (startswith(name, "D_") || startswith(name, "F_")) ? [rhs] : []
            for den in dens, t in sumterms(den)
                ut = units(t, env)
                iszero(ut[1]) && iszero(ut[2]) || return false, "$name term $(first(string(t), 80)) has units $ut"
            end
            startswith(name, "W_") && env[lhs] != (0, 0) && return false, "$name units"
        end
        true, ""
    catch err
        err isa UnitMismatch || rethrow()
        false, err.msg
    end
end

"""Relative error of `v` read back from `text` against `expected`, with `params` and
`concs` bound (NamedTuples)."""
function text_error(text, params, concs, expected)
    bind(s) = haskey(params, s) ? BigFloat(params[s]) : BigFloat(concs[s])
    env = Dict{Symbol, Any}(s => bind(s) for s in [collect(keys(params)); collect(keys(concs))])
    run_text(text, value, env, bind)
    Float64(abs(env[:v] - expected) / abs(expected)), env
end

const ARTIFACTS = (r"\^1(?![\d/.)])", r"(?<![\w.)])1 \^", r"\* 1(?![\d./])", r"L \* 1(?!\d)",
                   r"\(1\)", r"\^\(1/1\)")

"""Count of printing artifacts (`^1`, `1 ^`, `* 1`, `L * 1`, `(1)`) in `text`."""
artifacts(text) = sum(length(collect(eachmatch(r, text))) for r in ARTIFACTS)

# ─── One law through every rule and check ──────────────────────────────────────────────

"""Random positive values over `vars` (log-uniform 1e-2..1e2), for comparing normalizers."""
randvals(vars, rng) = [big(10.0)^(-2 + 4 * rand(rng)) for _ in vars.names]

"""
    rules_row(m, L, rng) -> NamedTuple

The normalizer rules on one conformation's law `L`: the spec rule and its order variants,
the slice rule and Cleland's N, with whether each gives the spec rule's `N`.
"""
function rules_row(m, L, rng)
    D, names = S.split_conc(L.D, L.vars), L.vars.names
    pts = [randvals(L.vars, rng) for _ in 1:2]
    spec = spec_norm(m, L)
    variants = Dict(o => spec_norm(m, L; order = o) for o in ORDERS)
    sl = slice_norm(m, L)
    cle = cleland_norm(m, L)
    samepair(a, b) = a.kind === :ratio || b.kind === :ratio ? a.kind === b.kind && a.shift == b.shift :
                     a.kind === :coefficient || b.kind === :coefficient ? false : same_N(a, b, D, pts)
    same(a) = samepair(a, spec)
    (; nD = length(D), constant = haskey(D, S.Mono()), spec_kind = spec.kind,
     spec_lambda = lamtext(spec, names), spec_basis = length(spec.lambda),
     same_order = Dict(o => same(v) for (o, v) in variants),
     kind_order = Dict(o => v.kind for (o, v) in variants),
     cle_order = Dict(o => cle === nothing ? missing : samepair(v, cle) for (o, v) in variants),
     slice_kind = sl.kind, same_slice = same(sl),
     cleland = cle !== nothing, same_cleland = cle === nothing ? missing : same(cle),
     cleland_lambda = cle === nothing ? "" : lamtext(cle, names), spec)
end

"""Each printed dependent's value against `all_parameters`; worst relative error."""
function ref_error(m, sol, env, params)
    allp = S.all_parameters(m, params)
    cn = S.column_names(m, sol)
    worst = 0.0
    for d in sol.deps
        d.label === :kinetic_group && continue
        n = cn[d.col]
        haskey(env, n) || return Inf
        worst = max(worst, Float64(abs(env[n] - allp[n]) / abs(allp[n])))
    end
    worst
end

"""
    analyze(m, rng) -> NamedTuple

One mechanism through the derivation, every normalizer rule (per conformation), the
printer under the spec rule (full law and coefficient form; layout and reference-line
variants), the slice's own printer, and the value, unit, reference and artifact checks.
"""
function analyze(m, rng)
    sol = S.solve(m)
    law = S.rate_law(m, sol)
    allo = law isa S.AlloLaw
    confs = allo ? filter(x -> x[2] !== nothing, [("A", law.active), ("I", law.inactive)]) :
            [("plain", law)]
    rules = [(tag, rules_row(m, L, rng)) for (tag, L) in confs]
    norms = [r.spec for (_, r) in rules]
    texts(; kw...) = allo ? allo_texts(m, law, sol, norms; kw...) :
                            plain_texts(m, law, sol, norms[1]; kw...)
    full, coef = texts()
    full2, coef2 = texts(; layout = :slice, refstyle = :product)
    sfull, scoef = S.rate_equation_string(m, law, sol), S.coefficient_form(m, law)
    u = name_units(m, sol)
    params, concs = random_point(m, rng)
    vars = allo ? law.vars : law.vars
    expected = params.E_total * lawval(law, varvals(vars, params, concs))
    e_full, env = text_error(full, params, concs, expected)
    e_coef, _ = text_error(coef, params, concs, expected)
    e_full2, env2 = text_error(full2, params, concs, expected)
    e_coef2, _ = text_error(coef2, params, concs, expected)
    e_sfull, _ = text_error(sfull, params, concs, expected)
    e_scoef, _ = text_error(scoef, params, concs, expected)
    (uf, nf), (uc, nc) = check_units(full, u), check_units(coef, u)
    (usf, nsf), (usc, nsc) = check_units(sfull, u), check_units(scoef, u)
    refden = maximum((denominator(e) for d in sol.deps if d.label !== :kinetic_group
                      for e in [last.(d.terms); d.keq]); init = 1)
    nq = allo ? 0 : length(S.split_conc(law.Q, law.vars))
    (; allo, n = allo ? law.n : 1, refden, nq, rules = [(tag, Base.structdiff(r, NamedTuple{(:spec,)})) for (tag, r) in rules],
     e_full, e_coef, e_full2, e_coef2, e_sfull, e_scoef,
     e_ref = ref_error(m, sol, env, params), e_ref2 = ref_error(m, sol, env2, params),
     units_full = uf, units_coef = uc, units_note = nf * " " * nc,
     units_sfull = usf, units_scoef = usc, units_snote = nsf * " " * nsc,
     len_full = length(full), len_coef = length(coef), len_sfull = length(sfull),
     len_scoef = length(scoef), len_full2 = length(full2),
     art = artifacts(full) + artifacts(coef), art2 = artifacts(full2) + artifacts(coef2),
     sqrt = occursin("sqrt(", full), ratio_names = occursin("⁻", coef),
     Napart = occursin("\nN", full), full = length(full) <= 4000 ? full : "",
     coef = length(coef) <= 4000 ? coef : "")
end
