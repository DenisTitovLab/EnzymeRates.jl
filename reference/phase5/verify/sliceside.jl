# ABOUTME: Verifier helpers on the slice side: build a slice mechanism (or its plain catalytic
# ABOUTME: mechanism, tags and sites dropped) from the plain data `rawrec` writes.
mkmet(S, t) = S.Met(t[1], t[2])
mkform(S, f) = S.RawForm(f[1], S.Met[mkmet(S, x) for x in f[2]], f[3], f[4])

function mkmech(S, r; plain::Bool = false)
    subs, prods, steps, group, tags, n, sites = r
    raws = S.RawStep[S.RawStep(mkform(S, s[1]), mkform(S, s[2]),
                               S.Met[mkmet(S, x) for x in s[3]],
                               S.Met[mkmet(S, x) for x in s[4]], s[5]) for s in steps]
    plain || isempty(tags) ? S.canonicalize(subs, prods, raws, group)[1] :
        S.canonicalize(subs, prods, raws, group; tags, multiplicity = n,
                       sites = S.Site[S.Site(a, b, c) for (a, b, c) in sites])[1]
end
