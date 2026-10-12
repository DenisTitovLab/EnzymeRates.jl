# ABOUTME: Joins the slice-side and package-side verdicts of a run and prints the agreement
# ABOUTME: matrices, the disagreements with examples as step lists, and the cost figures.
# Usage: report.jl <mode> [examples per finding]
using Serialization, Printf, Statistics
const R = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/round4/guard"
const MODE = isempty(ARGS) ? "smoke" : ARGS[1]
const NEX = length(ARGS) > 1 ? parse(Int, ARGS[2]) : 3

sl = deserialize(joinpath(R, "data", "slice_$MODE.jls"))
pk = Dict(r.key => r for r in deserialize(joinpath(R, "data", "pkg_$MODE.jls")))
sp = Dict(r.key => r.spec for r in deserialize(joinpath(R, "data", "specs_$MODE.jls")))

# ─── reasons ────────────────────────────────────────────────────────────────────────────
const PKG_REASONS = [
    r"no form free of one side" => "bottomless", r"close a catalytic cycle" => "RE turnover",
    r"not proportional to net reaction" => "cycle stoichiometry",
    r"differ in kind or in RE/SS" => "group not uniform",
    r"different conformational changes" => "group conformations",
    r"holds the reaction .* twice|both hold the reaction" => "reaction twice",
    r"both render the name" => "form names", r"changes a metabolite's role" => "role change",
    r":OnlyA" => "OnlyA guard", r"more than 5000" => "complexity cap",
    r"reserved name" => "reserved name", r"joins a form to itself|both ends are the form" => "self loop",
    r"both consumed and released" => "both sides"]
function reason(msg, table)
    msg === nothing && return nothing
    i = findfirst(p -> occursin(p.first, msg), table)
    i === nothing ? first(msg, 50) : table[i].second
end
const SLICE_REASONS = [
    r"one reaction" => "reaction twice", r"connected graph" => "disconnected",
    r"15 copies|32 metabolite" => "caps", r"different conformational changes" => "group conformations",
    r"differ in kind or in RE/SS" => "group not uniform", r"joins a form to itself" => "self loop",
    r"takes up and gives off one" => "both sides", r"neither a substrate nor a product" => "role",
    r"holds no step" => "empty group"]

function pkg_verdict(p)
    p === nothing && return "missing"
    p.ctor !== nothing && return "ctor: " * reason(p.ctor, PKG_REASONS)
    p.derive !== nothing && return "derive: " * reason(p.derive, PKG_REASONS)
    occursin(r"E_total \* \(0\)", p.vline) && return "accept (pkg N≡0)"
    occursin(r"/ \(0\)$", p.vline) && return "accept (pkg D≡0)"
    "accept"
end

"""The slice's own boundary verdict (what `derive` runs before the solve) and its errors."""
function slice_verdict(s)
    s.ctor_err !== nothing && return "ctor: " * reason(s.ctor_err, SLICE_REASONS)
    haskey(s, :fatal) && return "fatal"
    st = s.struct_
    st.bottomless !== nothing && return "bottomless"
    st.role !== nothing && return "role change"
    st.names !== nothing && return "names"
    st.admit !== nothing && return "admit"
    "pass"
end
"""Errors of the slice's later stages (solve, reduced law)."""
function slice_later(s)
    haskey(s, :struct_) || return "-"
    s.struct_.solve !== nothing && return "solve: " * first(s.struct_.solve, 40)
    haskey(s.law, :reduced_err) && return "law: " * first(s.law[:reduced_err], 40)
    "ok"
end

"""Law-based verdict: the first failing condition of (a)-(c), or ok; `nothing` if no law."""
function law_verdict(t)
    t === nothing && return nothing
    t.a && return "(a) D≡0"
    t.b && return "(b) N≡0"
    t.dfwd && t.drev && return "(c) D=0 both faces"
    t.dfwd && return "(c) D=0 fwd face"
    t.drev && return "(c) D=0 rev face"
    "ok"
end
law_red(s) = haskey(s, :law) ? law_verdict(get(s.law, :reduced, nothing)) : nothing
law_bare(s) = haskey(s, :law) ? law_verdict(get(s.law, :bare, nothing)) : nothing
lawname(v) = v === nothing ? "no law" : v

"""The oracle's worst verdict over the three points, `-` when it did not run."""
function oracle_verdict(s)
    haskey(s, :num) && s.num isa Dict && !haskey(s.num, :err) || return "-"
    vs = [s.num[f].verdict for f in (:pos, :fwd, :rev)]
    for (v, f) in zip(vs, ("pos", "fwd", "rev"))
        v === :diverges && return "diverges@$f"
    end
    for (v, f) in zip(vs, ("pos", "fwd", "rev"))
        v in (:pattern, :mismatch, :singular) && return "$v@$f"
    end
    any(==(:law_nan), vs) && return "law 0/0, oracle unique"
    "agree"
end

"""The structural checks this study adds: RE turnover cycle, unused reactant, face traps."""
function extra_verdict(s)
    haskey(s, :struct_) || return "-"
    st = s.struct_
    st.re_cycle !== nothing && return "RE turnover"
    st.unused !== nothing && return "unused reactant"
    get(st, :no_turnover, nothing) !== nothing && return "no turnover cycle"
    f = something(st.face_fwd, st.face_rev, "")
    isempty(f) && return "pass"
    startswith(f, "bottomless") || startswith(f, "inactive: bottomless") ? "face: bottomless" :
        occursin("trap", f) ? "face: trap" : occursin("split", f) ? "face: split" : "face: " * first(f, 30)
end

# ─── rendering a spec as a step list ────────────────────────────────────────────────────

mtext(x) = string(x[1], x[2] ? "::Mod" : "")
function ftext(f)
    res = isempty(f[3]) && isempty(f[4]) ? "" :
        "; residual = " * join(string.(f[3]), " + ") * join(" - $x" for x in f[4])
    isempty(f[2]) && isempty(res) ? string(f[1]) : "$(f[1])(" * join(mtext.(f[2]), ", ") * res * ")"
end
function steps_text(x)
    lines = String[]
    for (j, s) in enumerate(x.steps)
        a = join([ftext(s[1]); mtext.(s[3])], " + ")
        b = join([ftext(s[2]); mtext.(s[4])], " + ")
        tag = isempty(x.tags) ? "" : " ::$(x.tags[x.group[j]])"
        push!(lines, "      g$(x.group[j]) $a $(s[5] ? "⇌" : "<-->") $b$tag")
    end
    head = "      [$(join(x.substrates, " + ")) → $(join(x.products, " + "))" *
           (isempty(x.tags) ? "" : ", n=$(x.multiplicity), sites=$(x.sites)") * "]"
    join([head; lines], "\n")
end

function example(s; extra = "")
    println("    ", s.key, "  (", first(s.prov, 160), ")", extra)
    println(steps_text(sp[s.key]))
end

# ─── tables ─────────────────────────────────────────────────────────────────────────────

function table(title, rows, cols, f, g; recs = sl)
    println("\n## ", title)
    counts = Dict{Tuple{String, String}, Int}()
    for s in recs
        k = (string(f(s)), string(g(s)))
        counts[k] = get(counts, k, 0) + 1
    end
    rs = sort!(unique(first.(collect(keys(counts)))))
    cs = sort!(unique(last.(collect(keys(counts)))))
    w = max(18, maximum(length, rs) + 1)
    @printf("%s", rpad(rows * " \\ " * cols, w))
    foreach(c -> @printf("%s", lpad(first(c, 22), 24)), cs)
    println(lpad("total", 8))
    for r in rs
        @printf("%s", rpad(r, w))
        foreach(c -> @printf("%s", lpad(get(counts, (r, c), 0) == 0 ? "." : string(counts[(r, c)]), 24)), cs)
        println(lpad(string(sum(get(counts, (r, c), 0) for c in cs)), 8))
    end
end

countmap_(v) = (d = Dict{Any, Int}(); foreach(x -> d[x] = get(d, x, 0) + 1, v); sort!(collect(d); by = last, rev = true))

group_of(s) = startswith(s.set, "mut_") ? "mutants" :
              s.set == "skeptic_adv" ? "skeptic adv" :
              s.set == "terter_sample" ? "ter-ter sample" :
              any(p -> startswith(s.set, p), ("p1_", "p2_", "p5a_")) ? "adversarial core" : "core"

function main()
    println("# Guard study, mode $MODE: $(length(sl)) records")
    for gname in ("core", "adversarial core", "skeptic adv", "ter-ter sample", "mutants")
        n = count(s -> group_of(s) == gname, sl)
        println("  $gname: $n")
    end
    for kind in ("mut_to_re", "mut_to_ss", "mut_flip2", "mut_all_re", "mut_del_cat")
        println("    $kind: ", count(s -> s.set == kind, sl))
    end
    P(s) = pkg_verdict(get(pk, s.key, nothing))
    table("package verdict × law-based (reduced) verdict", "package", "law", P, s -> lawname(law_red(s)))
    table("slice boundary verdict × law-based (reduced) verdict", "slice", "law", slice_verdict, s -> lawname(law_red(s)))
    table("law-based bare × reduced", "bare", "reduced", s -> lawname(law_bare(s)), s -> lawname(law_red(s)))
    table("added structural checks × law-based (reduced)", "added", "law", extra_verdict, s -> lawname(law_red(s)))
    table("oracle × law-based (reduced)", "oracle", "law", oracle_verdict, s -> lawname(law_red(s)))
    table("oracle × added structural checks", "oracle", "added", oracle_verdict, extra_verdict)
    table("package verdict × oracle", "package", "oracle", P, oracle_verdict)
    table("slice later-stage errors × package verdict", "slice later", "package", slice_later, P)
    for gname in ("core", "adversarial core", "skeptic adv", "ter-ter sample", "mutants")
        table("[$gname] package × law-based (reduced)", "package", "law", P, s -> lawname(law_red(s));
              recs = [s for s in sl if group_of(s) == gname])
    end
    for kind in ("mut_to_re", "mut_to_ss", "mut_flip2", "mut_all_re", "mut_del_cat")
        table("[$kind] package × law-based (reduced)", "package", "law", P, s -> lawname(law_red(s));
              recs = [s for s in sl if s.set == kind])
    end

    println("\n## disagreements, up to $NEX examples each")
    buckets = Dict{String, Vector{Any}}()
    for s in sl
        p, l, o, x = P(s), law_red(s), oracle_verdict(s), extra_verdict(s)
        add(k) = push!(get!(buckets, k, Any[]), s)
        l === nothing && continue
        illlaw = l != "ok"
        pkgrej = !startswith(p, "accept")
        illo = o != "agree" && o != "-"
        illlaw && !pkgrej && add("law flags, package accepts: $l")
        pkgrej && !illlaw && add("package rejects, law ok: $p")
        illo && !illlaw && add("oracle ill-posed, law ok: $o | pkg $p | added $x")
        illlaw && o == "agree" && add("law flags, oracle agrees: $l")
        x != "pass" && !illlaw && add("added check flags, law ok: $x | oracle $o | pkg $p")
        x == "pass" && illo && add("oracle ill-posed, added checks pass: $o | law $l | pkg $p")
    end
    for k in sort!(collect(keys(buckets)))
        v = buckets[k]
        println("\n### $k: $(length(v))")
        sets = Dict{String, Int}()
        foreach(s -> sets[s.set] = get(sets, s.set, 0) + 1, v)
        println("    by set: ", join(("$a $b" for (a, b) in sort!(collect(sets); by = last, rev = true)), ", "))
        for s in v[1:min(NEX, end)]
            pr = get(pk, s.key, nothing)
            msg = pr === nothing ? "" : something(pr.ctor, pr.derive, "")
            st = haskey(s, :struct_) ? s.struct_ : nothing
            ex = st === nothing ? "" : "\n      added: " * string(something(st.re_cycle, st.unused, st.face_fwd, st.face_rev, "-"))
            example(s; extra = isempty(msg) ? ex : "\n      package: " * first(msg, 200) * ex)
            if haskey(s, :num) && s.num isa Dict && !haskey(s.num, :err)
                for f in (:pos, :fwd, :rev)
                    n = s.num[f]
                    @printf("      oracle %s: %s law=%.6g oracle(Λ=1e10,1e20,1e20')=%s\n", f, n.verdict, n.vlaw, string(round.(n.o; sigdigits = 7)))
                end
            end
        end
    end

    println("\n## size of the face discrepancies (oracle pattern/mismatch at a face)")
    for f in (:fwd, :rev)
        law_vs, pat = Float64[], Float64[]
        for s in sl
            haskey(s, :num) && s.num isa Dict && !haskey(s.num, :err) || continue
            n = s.num[f]
            n.verdict in (:pattern, :mismatch) || continue
            o2, o3 = n.o[2], n.o[3]
            push!(pat, abs(o2 - o3) / max(abs(o2), abs(o3), 1e-300))
            isnan(n.vlaw) || push!(law_vs, abs(n.vlaw - o2) / max(abs(o2), abs(n.vlaw), 1e-300))
        end
        isempty(pat) && continue
        @printf("  %s face: n=%d; |o(pattern1)-o(pattern2)|/|o|: median %.2e max %.2e; |law-o|/|o|: median %.2e max %.2e\n",
                f, length(pat), median(pat), maximum(pat),
                isempty(law_vs) ? NaN : median(law_vs), isempty(law_vs) ? NaN : maximum(law_vs))
    end

    println("\n## package laws with a zero numerator or denominator (accepted and derived)")
    for (pat, what) in ((r"E_total \* \(0\)", "N≡0"), (r"/ \(0\)$", "D≡0"))
        ks = [s for s in sl if haskey(pk, s.key) && occursin(pat, pk[s.key].vline)]
        println("  package $what: $(length(ks)); law-based (reduced) verdicts: ",
                join(("$(k)=$(v)" for (k, v) in countmap_(lawname.(law_red.(ks)))), ", "))
    end

    println("\n## costs (seconds, medians and 95th percentiles over records with a law)")
    q(v) = isempty(v) ? "n/a" : @sprintf("median %.2e  p95 %.2e  max %.2e  (n=%d)", median(v), quantile(v, 0.95), maximum(v), length(v))
    for (gname, filt) in (("plain", s -> !s.allo), ("allosteric", s -> s.allo))
        rs = [s for s in sl[50:end] if haskey(s, :law) && haskey(s.law, :t_reduced_test) && filt(s)]
        println("  $gname:")
        println("    reduced law (solve + rate_law): ", q([s.law[:t_reduced_law] for s in rs]))
        println("    support test on reduced law:    ", q([s.law[:t_reduced_test] for s in rs]))
        println("    ratio test/law:                 ", q([s.law[:t_reduced_test] / s.law[:t_reduced_law] for s in rs]))
        rb = [s for s in rs if haskey(s.law, :t_bare_test)]
        println("    bare law(s):                    ", q([s.law[:t_bare_law] for s in rb]))
        println("    support test on bare law:       ", q([s.law[:t_bare_test] for s in rb]))
        pd = [pk[s.key].t_derive for s in rs if haskey(pk, s.key) && pk[s.key].derive === nothing]
        println("    package rate_equation_string:   ", q(pd))
    end
    tt = [s for s in sl if s.set == "terter_sample" && haskey(s, :law) && haskey(s.law, :t_reduced_test)]
    println("  ter-ter sample: law ", q([s.law[:t_reduced_law] for s in tt]), "\n                  test ", q([s.law[:t_reduced_test] for s in tt]))
end

main()
