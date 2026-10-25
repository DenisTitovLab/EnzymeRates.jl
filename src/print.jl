# ABOUTME: Printed rate laws: the law in the independent constants, with reference lines and
# ABOUTME: the normalizer's factor, its coefficient form, and MWC laws per conformation.

# ─── Text of powers, terms and sums ─────────────────────────────────────────

"""`x^e` as text for a positive exponent `e`: `x`, `x^2`, `sqrt(x)`, `x^(3/2)`."""
_power_text(x, e::_Exp) = e == 1 ? string(x) : e == 1 // 2 ? "sqrt($x)" :
    isinteger(e) ? "$x^$(numerator(e))" : "$x^($(numerator(e))/$(denominator(e)))"

"""Factors as a ratio, `a * b / (c * d)`; `1` when no factor is above the bar."""
_ratio_text(up::Vector{String}, dn::Vector{String}) =
    string(isempty(up) ? "1" : join(up, " * "),
           isempty(dn) ? "" :
           " / " * (length(dn) == 1 ? only(dn) : "(" * join(dn, " * ") * ")"))

"""
    _term_text(c, mo, names; rest = one(_Poly), over = false) -> String

The positive term `c·mo` as a ratio of products over the variable names `names`, times the
polynomial `rest` and over the normalizer's printed factor `N` when `over`:
`2 * A * k1 * (k2 + k3) / (k4 * N)`, and `1` for the empty term.
"""
function _term_text(c, mo::_Mono, names::Vector{Symbol}; rest::_Poly = one(_Poly),
                    over::Bool = false)
    up = String[_power_text(names[v], e) for (v, e) in mo if e > 0]
    dn = String[_power_text(names[v], -e) for (v, e) in mo if e < 0]
    numerator(c) == 1 || pushfirst!(up, string(numerator(c)))
    denominator(c) == 1 || pushfirst!(dn, string(denominator(c)))
    rest == one(_Poly) || push!(up, "(" * _sum_text(rest, names) * ")")
    over && push!(dn, "N")
    _ratio_text(up, dn)
end

"""The nonzero polynomial `p` over the variable names `names` as a sum, `a + 2 * b - c`,
its terms in monomial order."""
function _sum_text(p::_Poly, names::Vector{Symbol})
    s = join(string(c < 0 ? " - " : " + ", _term_text(abs(c), mo, names))
             for (mo, c) in sort!(collect(p.terms); by = first))
    (startswith(s, " + ") ? "" : "-") * chop(s; head = 3, tail = 0)
end

const _SUPERSCRIPTS = Dict(zip("0123456789", "⁰¹²³⁴⁵⁶⁷⁸⁹"))

"""
A coefficient's name for the concentration monomial `mo` over the names `names`: `0` for
none, else each concentration in variable order with its exponent as a superscript, joined
by `_`: `A²_B`.
"""
function _conc_label(mo::_Mono, names::Vector{Symbol})
    isempty(mo) && return "0"
    sup(e) = e == 1 ? "" : String([_SUPERSCRIPTS[d] for d in string(numerator(e))])
    join((string(names[v], sup(e)) for (v, e) in mo), "_")
end

# ─── A law's polynomials over its normalizer ─────────────────────────────────

"""
One printed term of a law's polynomial: the concentration monomial `conc` times its
coefficient `c·mo·rest`, over the normalizer's printed factor `N` when `over`. `c` and `mo`
have the normalizer's integer and monomial contents divided out; `rest` is the
coefficient's polynomial factor, `one(_Poly)` when there is none.
"""
struct _Term
    conc::_Mono
    c::Rational{Int}
    mo::_Mono
    rest::_Poly
    over::Bool
end

"""Multiply the factors printed apart, `apart` (polynomial and exponent), by `p^e`: an
equal polynomial adds its exponent, and a factor whose exponent sums to 0 goes."""
function _apart!(apart::Vector{Tuple{_Poly, _Exp}}, p::_Poly, e::_Exp)
    i = findfirst(x -> x[1] == p, apart)
    i === nothing ? push!(apart, (p, e)) : (apart[i] = (p, apart[i][2] + e))
    filter!(x -> !iszero(x[2]), apart)
end

"""
    _printed(L, nz, polys, subs) -> (; N, terms)

The polynomials `polys` of the conformation law `L` divided by its normalizer `nz`
(`_normalizer`, `subs` the substrates' variables). The normalizer's integer contents to
integer powers and monomial contents divide every coefficient; the other factors, an
integer content to a fractional power included, print apart: `N` is their text, empty
when there is none. `terms[j]` holds a `_Term` for each concentration monomial of
`polys[j]`, in the unit rule's order (`_term_key`). A coefficient whose polynomial factor
is `N`'s one factor drops it and is not over `N`.
"""
function _printed(L::_Law, nz::_Norm, polys::Vector{_Poly}, subs::Vector{Int})
    D, names = _split_conc(L.D, L.vars), L.vars.names
    k, mu, apart = 1 // 1, _Mono(), Tuple{_Poly, _Exp}[]
    for (mo, e) in nz.lambda
        c, u, rest = _content(D[mo])
        if isinteger(e)
            k *= (c // 1)^numerator(e)
        elseif c != 1
            _apart!(apart, _term(_Mono(), c), e)
        end
        mu = _mmul(mu, _mpow(u, e))
        rest == one(_Poly) || _apart!(apart, rest, e)
    end
    paren(p) = "(" * _sum_text(p, names) * ")"
    N = isempty(apart) ? "" :
        length(apart) == 1 && apart[1][2] == 1 ? _sum_text(apart[1][1], names) :
        _ratio_text(String[_power_text(paren(p), e) for (p, e) in apart if e > 0],
                    String[_power_text(paren(p), -e) for (p, e) in apart if e < 0])
    whole = length(apart) == 1 && apart[1][2] == 1 ? apart[1][1] : zero(_Poly)
    # `term` captures only names assigned once, so none is boxed (the loop reassigns `k`).
    scale, inv_mu = k, _mpow(mu, -1)
    function term(conc, q)
        c, u, rest = _content(q)
        w = rest == whole
        _Term(conc, (c // 1) / scale, _mmul(u, inv_mu), w ? one(_Poly) : rest,
              !isempty(apart) && !w)
    end
    terms = Vector{_Term}[
        _Term[term(conc, q) for (conc, q) in sort!(collect(_split_conc(p, L.vars));
                                                   by = x -> _term_key(x[1], L.vars, subs))]
        for p in polys]
    (; N, terms)
end

"""The terms `ts` (`_printed`) as a sum, each concentration monomial times its
coefficient."""
_group_sum(ts::Vector{_Term}, names::Vector{Symbol}) =
    join((_term_text(t.c, _mmul(t.conc, t.mo), names; rest = t.rest, over = t.over)
          for t in ts), " + ")

"""
    _coefficient_sum(prefix, ts, names) -> (sum, lines)

The terms `ts` (`_printed`) as a sum of named coefficients, `prefix_A_B * A * B`, with
`prefix_0` alone for the constant monomial and `1` when its coefficient is 1, and one
line per named coefficient giving its value in the fitted constants.
"""
function _coefficient_sum(prefix::String, ts::Vector{_Term}, names::Vector{Symbol})
    unit(t) = isempty(t.conc) && t.c == 1 && isempty(t.mo) && t.rest == one(_Poly) &&
              !t.over
    label(t) = string(prefix, "_", _conc_label(t.conc, names))
    sum = join((unit(t) ? "1" : isempty(t.conc) ? label(t) :
                string(label(t), " * ", _term_text(1, t.conc, names)) for t in ts),
               " + ")
    lines = String[string(label(t), " = ",
                          _term_text(t.c, t.mo, names; rest = t.rest, over = t.over))
                   for t in ts if !unit(t)]
    sum, lines
end

"""The `v` line, `v = E_total * (Q) * (H) / (D)`: `Q` bare when it is one name."""
function _v_line(q::AbstractString, h::AbstractString, d::AbstractString)
    lead = occursin(r"^[^ ()*/+^-]+$", q) ? q * " * " : "(" * q * ") * "
    "v = E_total * $(lead)($h) / ($d)"
end

const _COEFFICIENTS_COMMENT =
    "# The coefficients are functions of the fitted constants, not parameters:"

# ─── Reference lines and header ──────────────────────────────────────────────

"""
The reference lines of the solve `s` over the variables `vars`: a `# Wegscheider
constraints:` section, then a `# Haldane constraints:` section, each dependent constant by
name as a ratio of products of `Keq` and the independents with exact exponents (`sqrt`,
`^(3/2)`), `= 1` for a constant that thermodynamics pins to 1. A section without a line is
left out; a kinetic group's other members equal its kept constant and print no line.
"""
function _reference_lines(s::_Solve, vars::_Vars)
    lines, value = String[], _column_values(vars, s)
    for (label, title) in ((:wegscheider, "Wegscheider"), (:haldane, "Haldane"))
        ds = sort!(filter(d -> d.label === label, s.dependents); by = d -> s.names[d.col])
        isempty(ds) || push!(lines, "# $title constraints:")
        for d in ds
            push!(lines, "$(s.names[d.col]) = $(_term_text(1, value[d.col], vars.names))")
        end
    end
    lines
end

"""The `params` and `concs` destructuring lines of `m`'s printed law over the solve `s`."""
_header_lines(m::Union{Mechanism, AllostericMechanism}, s::_Solve) =
    String["(; $(join((s.fitted..., :Keq, :E_total), ", "))) = params",
           "(; $(join(metabolites(m), ", "))) = concs"]

# ─── Plain and MWC laws ──────────────────────────────────────────────────────

"""
    _law_text(m, s, L; coefficients = false) -> String

The printed law `L` of `m` over the solve `s` (`rate_equation_string`), in the order it
evaluates. In full: the header (`_header_lines`), the reference lines (`_reference_lines`),
`N = …` when the normalizer has a factor that is not a single term, or an integer to a
fractional power, and the `v` line (`_v_line`) with `Q` and `D` over the normalizer. With
`coefficients`: the header, the coefficient comment, `N` when it prints, one line per
coefficient giving its value, and the `v` line in the named coefficients.
"""
function _law_text(m::Mechanism, s::_Solve, L::_Law; coefficients::Bool = false)
    subs = _substrate_vars(m, L.vars)
    P = _printed(L, _normalizer(L.D, L.vars, subs), [L.Q, L.D], subs)
    names, (q, d) = L.vars.names, P.terms
    h = _sum_text(L.H, names)
    norms = isempty(P.N) ? String[] : ["N = $(P.N)"]
    if coefficients
        (qs, lq) = _coefficient_sum("n", q, names)
        (ds, ld) = _coefficient_sum("d", d, names)
        return join([_header_lines(m, s); _COEFFICIENTS_COMMENT; norms; lq; ld;
                     _v_line(qs, h, ds)], "\n")
    end
    join([_header_lines(m, s); _reference_lines(s, L.vars); norms;
          _v_line(_group_sum(q, names), h, _group_sum(d, names))], "\n")
end

"""
The conformation law `L` of an MWC law, tagged `tag` (`"A"` or `"I"`), printed over its own
normalizer: the terms (`_printed`) of its `Q`, its `D` and its free-E weight (the tree sum
of free E's segment times the law's shift), each over the normalizer's printed factor,
which reads as the name `N_<tag>` appended to the variable names.
"""
function _conformation(L::_Law, tag::String, subs::Vector{Int})
    names = [L.vars.names; Symbol("N_", tag)]
    n_var = _monomial(length(names) => -1)
    over(ts) = _Term[_Term(t.conc, t.c, t.over ? _mmul(t.mo, n_var) : t.mo, t.rest, false)
                     for t in ts]
    P = _printed(L, _normalizer(L.D, L.vars, subs), [L.Q, L.D, L.trees[1] * L.shift], subs)
    (; tag, names, N = P.N, q = over(P.terms[1]), d = over(P.terms[2]),
     f = over(P.terms[3]))
end

"""Push `W_<tag> = (1 + …)^k * …`, the product of the regulatory sites' binding
polynomials `ws` (polynomial and multiplicity) in one conformation, onto `defs`; a site
whose polynomial is 1 is left out, and so is the line when every site is."""
function _site_line!(defs::Vector{String}, tag::String, ws::Vector{Tuple{_Poly, Int}},
                     names::Vector{Symbol})
    parts = [string("(", _sum_text(w, names), ")", k == 1 ? "" : "^$k")
             for (w, k) in ws if _nterms(w) > 1]
    isempty(parts) || push!(defs, "W_$tag = " * join(parts, " * "))
    defs
end

"""
The printed MWC law `A` of `am` over the solve `s`. Each conformation `c` (`A`, then `I`
unless it is free E alone) defines `Q_c` (unless it turns nothing over) and `D_c`, over its
own normalizer `N_c`, which prints as a line when it has a factor that is not a single
term, or an integer to a fractional power; `F_A` and `F_I`, the free-E weights over the
normalizers, print only when they differ and each only when it is not 1; `W_A` and `W_I`
multiply the regulatory sites' factors. The `v` line is
`E_total·H·(F_I^n·Q_A·D_A^(n−1)·W_A + L·F_A^n·Q_I·D_I^(n−1)·W_I) /
(F_I^n·D_A^n·W_A + L·F_A^n·D_I^n·W_I)` without its factors of 1. Both forms open with the
header and end with `v`, in the order they evaluate. In full, the reference lines and the
definitions come between; with `coefficients`, the coefficient comment, the normalizers,
one line per coefficient, named by conformation (`nA_0`, `dA_S`, `fI_B`), and the
definitions: `Q_c`, `D_c` and `F_c` read in those coefficients, `W_A` and `W_I` in the
fitted constants.
"""
function _law_text(am::AllostericMechanism, s::_Solve, A::_AlloLaw;
                   coefficients::Bool = false)
    subs = _substrate_vars(am, A.vars)
    confs = [_conformation(A.active, "A", subs)]
    _flip_only(A.inactive) || push!(confs, _conformation(A.inactive, "I", subs))
    defs, norms, coefs = String[], String[], String[]
    function define!(name, ts, c)
        if coefficients
            prefix = string(name[1] == 'Q' ? 'n' : lowercase(name[1]), c.tag)
            text, lines = _coefficient_sum(prefix, ts, c.names)
            push!(defs, "$name = $text")
            append!(coefs, lines)
        else
            push!(defs, "$name = $(_group_sum(ts, c.names))")
        end
    end
    for c in confs
        isempty(c.N) || push!(coefficients ? norms : defs, "N_$(c.tag) = $(c.N)")
        isempty(c.q) || define!("Q_$(c.tag)", c.q, c)
        define!("D_$(c.tag)", c.d, c)
    end
    flips = [_group_sum(c.f, c.names) for c in confs]
    if flips[1] != (length(flips) == 2 ? flips[2] : "1")
        for (c, t) in zip(confs, flips)
            t == "1" || define!("F_$(c.tag)", c.f, c)
        end
    end
    _site_line!(defs, "A", Tuple{_Poly, Int}[(w, k) for (w, _, k) in A.sites], A.vars.names)
    _site_line!(defs, "I", Tuple{_Poly, Int}[(w, k) for (_, w, k) in A.sites], A.vars.names)
    has(x) = any(l -> startswith(l, x * " = "), defs)
    pw(x, e) = !has(x) || e == 0 ? "" : e == 1 ? x : "$x^$e"
    product(fs...) = (f = filter(!isempty, String[fs...]);
                      isempty(f) ? "1" : join(f, " * "))
    branch(c, o, e) = product(pw("F_$o", A.n), pw("Q_$c", e), pw("D_$c", A.n - e),
                              pw("W_$c", 1))
    num = String[]
    has("Q_A") && push!(num, branch("A", "I", 1))
    has("Q_I") && push!(num, "L * " * branch("I", "A", 1))
    rest = branch("I", "A", 0)
    den = branch("A", "I", 0) * " + L" * (rest == "1" ? "" : " * " * rest)
    v = _v_line(join(num, " + "), _sum_text(A.H, A.vars.names), den)
    head = _header_lines(am, s)
    coefficients ? join([head; _COEFFICIENTS_COMMENT; norms; coefs; defs; v], "\n") :
                   join([head; _reference_lines(s, A.vars); defs; v], "\n")
end

# ─── rate_equation_string ────────────────────────────────────────────────────

"""
    rate_equation_string(m; coefficients = false) -> String

Return the rate law of mechanism `m` as a multi-line `String` (it returns the text; it does
not print). Pass a concrete `Mechanism` / `AllostericMechanism` or its compiled
[`EnzymeMechanism`](@ref) singleton.

The text runs as Julia code once `params` and `concs` are bound. It destructures
[`parameters`](@ref)`(m)` from `params` and the metabolites from `concs`, lists each
dependent constant under `# Wegscheider constraints:` or `# Haldane constraints:` as a
reference line in `Keq` and the independent constants, and ends with
`v = E_total * (Q) * (H) / (D)` in the independent constants. `H = Π[S] − Π[P]/Keq`
vanishes at equilibrium; `Q` and `D` have positive terms. Numerator and denominator are
divided by a normalizer that makes every term of `D` dimensionless, when one exists: `D`'s
constant term when it has one, so `D` reads `1 + S/K…`; otherwise Cleland's normalizer when
`D` holds the terms `ΠS/X` for every substrate `X` and `ΠS`, `ΠS` the product of the
substrates, so each substrate term reads as its concentrations over their Michaelis
constants; otherwise a product of powers of `D`'s coefficients with `D`'s unit, each
concentration counting as a dimension of its own. When no such product exists, some
weighting of the concentrations gives every term of `D` the same nonzero degree. In
`E + A <--> E(A)`, `E(A) + B ⇌ E(Q) + P`, `E(Q) <--> E + Q`, for example, the
rapid-equilibrium Theorell–Chance step makes every term carry exactly one of `B` and `P`.
The law then prints as derived, and its terms keep their units.
A factor of the normalizer that is not a single term, or an integer to a fractional power,
prints once as `N = …`.
An allosteric mechanism prints each conformation's numerator `Q_A`, `Q_I`, denominator
`D_A`, `D_I`, normalizer `N_A`, `N_I`, free-enzyme weight `F_A`, `F_I` (when they differ)
and regulatory-site factor `W_A`, `W_I`, then the MWC combination.

With `coefficients = true` the law reads in named coefficients,
`v = E_total * n_0 * (H) / (1 + d_S * S + …)`, each named by the concentrations it
multiplies (`d_A²_B` for `A^2 * B`). That text runs as Julia code too: it destructures
`params` and `concs`, prints the coefficient comment, `N` (`N_A`, `N_I`) when the
normalizer prints, each coefficient's value in the fitted constants, one line each, for an
allosteric law `Q_c`, `D_c` and `F_c` in those coefficients and the site factors `W_c`,
and ends with the `v` line; the reference lines print only in the full law.

```jldoctest
julia> using EnzymeRates

julia> m = @enzyme_mechanism begin
           substrates: S
           products: P
           steps: begin
               E + S ⇌ E(S)
               E(S) <--> E(P)
               E(P) ⇌ E + P
           end
       end;

julia> print(rate_equation_string(m))
(; K_EP_to_E_P, K_ES_to_E_S, k_ES_to_EP, Keq, E_total) = params
(; S, P) = concs
# Haldane constraints:
k_EP_to_ES = K_EP_to_E_P * k_ES_to_EP / (K_ES_to_E_S * Keq)
v = E_total * (k_ES_to_EP / K_ES_to_E_S) * (S - P / Keq) / (1 + S / K_ES_to_E_S + P / K_EP_to_E_P)

julia> print(rate_equation_string(m; coefficients = true))
(; K_EP_to_E_P, K_ES_to_E_S, k_ES_to_EP, Keq, E_total) = params
(; S, P) = concs
# The coefficients are functions of the fitted constants, not parameters:
n_0 = k_ES_to_EP / K_ES_to_E_S
d_S = 1 / K_ES_to_E_S
d_P = 1 / K_EP_to_E_P
v = E_total * n_0 * (S - P / Keq) / (1 + d_S * S + d_P * P)
```
"""
function rate_equation_string(m::Union{Mechanism, AllostericMechanism};
                              coefficients::Bool = false)
    _check_well_posed(m)
    _rate_text(m, coefficients)
end

rate_equation_string(@nospecialize(m::AbstractEnzymeMechanism);
                     coefficients::Bool = false) =
    _rate_text(_derivation_input(typeof(m)), coefficients)

"""The printed law of the well-posed mechanism `m` (`_law_text`), once class 5 passes on
its derived law (`_with_law`)."""
_rate_text(m::Union{Mechanism, AllostericMechanism}, coefficients::Bool) =
    _with_law((m, s, L) -> _law_text(m, s, L; coefficients), m)
