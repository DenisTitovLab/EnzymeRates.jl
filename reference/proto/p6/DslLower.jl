# ABOUTME: Prototype lowering of a mechanism macro body (the DSL's step grammar) to P6's raw GM
# ABOUTME: input: labels, sides, forms, residuals, X::Inh copies, groups, state tags and sites.
module DslLower

import ..GM

const STATES = (:OnlyA, :OnlyI, :EqualAI, :NonequalAI)

"""The lines of a `begin … end` block."""
lines(b::Expr) = [x for x in b.args if !(x isa LineNumberNode)]

"""A `label: v1, v2, …` line as `(label, [v1, v2, …])`."""
function labelled(line)
    if line isa Expr && line.head === :tuple
        label, first = labelled(line.args[1])
        return label, [first; line.args[2:end]]
    end
    line isa Expr && line.head === :call && line.args[1] === :(:) ||
        error("expected `label: values`, got $line")
    line.args[2], Any[line.args[3]]
end

terms(x) = x isa Expr && x.head === :call && x.args[1] === :+ ? reduce(vcat, terms.(x.args[2:end])) : Any[x]

"""Residual `A + B - P` as `(added, subtracted)`."""
function residual(x, add = Symbol[], sub = Symbol[], sign = 1)
    if x isa Symbol
        push!(sign > 0 ? add : sub, x)
    elseif x.args[1] === :+
        foreach(y -> residual(y, add, sub, sign), x.args[2:end])
    else                                         # binary or unary minus
        length(x.args) == 3 && residual(x.args[2], add, sub, sign)
        residual(x.args[end], add, sub, -sign)
    end
    add, sub
end

"""Lower a mechanism macro body to `(substrates, products, steps, group, kw)`."""
function lower(body::Expr)
    subs, prods, inhs = Symbol[], Symbol[], Symbol[]
    states, sitedefs = Dict{Symbol, Symbol}(), Tuple{Vector{Symbol}, Int}[]
    mult, stepblock, allosteric = 1, nothing, false
    for line in lines(body)
        label, vals = labelled(line)
        if label === :substrates; append!(subs, vals)
        elseif label === :products; append!(prods, vals)
        elseif label in (:regulators, :catalytic_inhibitors); append!(inhs, vals)
        elseif label in (:steps, :catalytic_steps)
            stepblock, allosteric = only(vals), allosteric || label === :catalytic_steps
        elseif label === :catalytic_multiplicity; mult, allosteric = only(vals), true
        elseif label === :allosteric_regulators
            for v in vals; states[v.args[1]] = v.args[2]; end
        elseif label isa Expr && label.args[1] === :regulatory_site
            ligs = only(labelled(l)[2] for l in lines(only(vals)) if labelled(l)[1] === :ligands)
            n = label.args[2]
            push!(sitedefs, (Symbol.(ligs), n isa Expr ? n.args[2] : n))
        else
            error("unknown label $label")
        end
    end
    met(x) = x isa Expr ? GM.Met(x.args[1], true) : GM.Met(x, !(x in subs || x in prods))
    ismet(x) = x isa Expr ? x.head === :(::) : (x in subs || x in prods || x in inhs)
    function form(x)
        x isa Symbol && return GM.RawForm(x, GM.Met[], Symbol[], Symbol[])
        args = x.args[2:end]
        par = findfirst(a -> a isa Expr && a.head === :parameters, args)
        add, sub = par === nothing ? (Symbol[], Symbol[]) : residual(only(args[par].args).args[2])
        GM.RawForm(x.args[1], [met(a) for a in args if !(a isa Expr && a.head === :parameters)], add, sub)
    end
    function side(x)
        ts = terms(x)
        f = only(filter(!ismet, ts))
        form(f), [met(t) for t in ts if t !== f]
    end
    function peel(rhs)                        # a state tag binds to the right side's last term
        ts = terms(rhs)
        last = ts[end]
        last isa Expr && last.head === :(::) && last.args[2] in STATES || return rhs, nothing
        ts[end] = last.args[1]
        (length(ts) == 1 ? ts[1] : Expr(:call, :+, ts...)), last.args[2]
    end
    steps, group, tags = GM.RawStep[], Int[], Symbol[]
    for (g, line) in enumerate(lines(stepblock))
        tag = nothing
        if line isa Expr && line.head === :(::) && line.args[2] in STATES
            line, tag = line.args
        end
        for s in (line isa Expr && line.head === :tuple ? line.args : [line])
            op, lhs, rhs = s.args
            rhs, t = peel(rhs)
            t === nothing || (tag = t)
            (from, c), (to, r) = side(lhs), side(rhs)
            push!(steps, GM.RawStep(from, to, c, r, op === :⇌))
            push!(group, g)
        end
        allosteric && (tag === nothing ? error("group $g has no state tag") : push!(tags, tag))
    end
    sited = Set(l for (ligs, _) in sitedefs for l in ligs)
    sites = [GM.Site(ligs, [states[l] for l in ligs], m) for (ligs, m) in sitedefs]
    for (l, st) in states
        l in sited || push!(sites, GM.Site([l], [st], mult))
    end
    (Symbol.(subs), Symbol.(prods), steps, group, (; tags, multiplicity = mult, sites))
end

end # module
