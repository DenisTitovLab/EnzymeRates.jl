# R5.2 and R5.3: L phantoms

Paths are relative to `phase45/research5/phantom/` unless they start with `SPEC`
(`phase3/proto3_wt/docs/superpowers/specs/2026-10-12-reaction-graph-rewrite-design.md`) or
`key.md` (`phase45/understand/key.md`). `logs/X:N` is line N of a log.

## 0. Answer

- R5.2 passes for the spec's rule (pre-filter plus rank) on every population checked.
- The rule flags 926 variants. Each one's L→0 law equals its catalytic mechanism's plain law. Every flagged variant lies inside that law's family to ≤2.7e-33.
- The modular and BigFloat ranks agree on 37,830 of 37,830 ranks.
- An exact rescaling certificate proves the families equal for 923 of the 926. The other 3 equal the L→0 family only as a union over L (section 2).
- The rank test is not a family test. Three of 60 sampled variants have rank equal to their L→0 law yet a strictly larger family. The pre-filter is what keeps them out (logs/contain_summary.log:14,15,18).
- After the key's own deduplication, the spec's rule removes few fits: 76 of 2,398 in bi-bi round 1, 179 of 7,248 in round 2 and 76 of 2,535 in ping-pong round 1, about 3% (logs/fits.log:2,3,5). The spec expected 10–15% (SPEC:961).
- The certificate needs no pre-filter. Used alone, it removes 512, 1,386 and 512 fits (about 20%), each backed by a proof (logs/fits.log:2,3,5).
- R5.3, fit sources of the 926 phantoms:
  - 660 are `to_allosteric` children whose source is their enumeration parent.
  - 248 have their source at an ancestor or at an ancestor's plain child, always fitted in an earlier round.
  - 18 have no fitted source. All 18 are `re_to_ss` chain-flank flips at n = 1, whose L→0 law is itself rank-deficient.
- Recommendation for the 18: b2 (fit the variant itself, unflagged).

## 1. Inputs and method

- Populations come from the phase-2 package (`census/env`, 4f9a34a) through `scripts/gen.jl`. That script mirrors the census's `expand_labeled` and identify's base tier.
  - bi-bi `A[C], B[N] → P[C], Q[N]`:
    - round 1: 2,752 candidates, 1,833 of them one-subunit variants (logs/gen_bibi.log:2);
    - round 2: 1,000 of 2,752 parents sampled with `Xoshiro(60)`, giving 9,225 candidates, 8,075 of them variants (logs/gen_bibi.log:3).
  - ping-pong `A[CX], B[N] → P[C], Q[NX]`, round 1: 3,220 candidates, 2,164 of them variants (logs/gen_pingpong.log:2).
  - Two regulator populations test the pre-filter's site scope:
    - `uniuni_R`, 85 candidates (logs/gen_uniuni_R.log:1-3);
    - `bibi_R`, 4,144 seeds plus 2,614 round-1 candidates from 300 sampled parents (logs/gen_bibi_R.log:1-2).
- Laws come from the slice (`slice/src`). Its flip at `forms[1]` is phase 4's flip at free E. Admission is the slice's `admit`.
- `scripts/analyze.jl` computes, for each variant:
  - the pre-filter;
  - the modular rank: a Jacobian of log v at n + 6 random points, one draw each mod 2^61−1 and 2^62−57, keeping the larger rank;
  - the BigFloat rank: 192-bit, 2n + 6 points, pivot cut 1e-30, as the census's `slice_rank`;
  - absolute and relative modes;
  - the certificate;
  - the slot map;
  - the keys.
- Comparators: the variant against its L→0 law, which is the law of its catalytic mechanism (`plain = true` in `parse_neutral`), and against its enumeration parent.
- Certificate (`scripts/plib.jl` `certificate`), for n = 1 with QI = 0:
  - write v = N/(Y + L·X) with N = H·fI·QA·WA, Y = fI·DA·WA and X = fA·DI·WI;
  - require X to be a sub-sum of Y with equal coefficients;
  - find weights w over the constants with e·w = s on N and Y−X, and e·w = s+1 on X.
  - Then v(θ, L) = v₀(θ ⊙ (1+L)^w) for all positive θ and L, so the two families are equal.
  - The map was checked at 5 points in 256-bit arithmetic on all 3,971 certificates: worst relative error 8.6e-77 (logs/summary.log:24).
- `scripts/contain.jl` checks numeric containment both ways: Float64 Levenberg–Marquardt with 6 restarts, then a 256-bit polish, at concentrations from 10^-2 to 10^2. The pass mark is 1e-12. The negative controls fail as they should, at 4.0e-2 and 7.3e-3 (logs/contain_summary.log:76-78).
- Runs: the slowest took 168 s; the peak was 1.26 GB (logs/gen_bibi.log, logs/analyze_bibi_1.log).

## 2. R5.2: the pre-filter plus the rank

| population | variants | admitted | pass pre-filter | rank = L→0 rank | flagged | flagged with certificate |
|---|---|---|---|---|---|---|
| bi-bi r1 | 1,833 | 1,833 | 441 | 724 | 188 | 188 |
| bi-bi r2 | 8,075 | 8,050 | 1,464 | 3,030 | 550 | 547 |
| ping-pong r1 | 2,164 | 1,833 | 539 | 897 | 188 | 188 |
| regulator populations | 6,843 | 6,843 | 621 (1,333 bind nothing in the catalytic graph) | 3 | 0 | 0 |

Sources: logs/summary.log:3-10. Ping-pong would flag 228 if admission were ignored (logs/summary.log:5).

- Rank methods:
  - modular equals BigFloat on 18,915 variants × 2 laws (logs/summary.log:14);
  - the BigFloat gap is wide: kept pivots ≥5.0e-8, dropped ≤4.3e-57 (logs/summary.log:15);
  - the slice's fitted count equals the package's `_independent_param_count` on all 18,915 (logs/summary.log:16).
- Flagged against the L→0 law, 926 of 926 on each check (logs/summary.log:23):
  - the edge order is the same;
  - the independent slot sets are the same;
  - the laws are identical under the slot map (6 modular points);
  - key(active law) = key(catalytic plain law);
  - np = np₀ + 1;
  - no inactive independent reaches the law.
- Flagged against the enumeration parent: 660 of 660 plain parents share the key (logs/summary.log:23). The 266 allosteric parents are the wrong comparator. Flagged #8940 is certified equal to its L→0 law yet "contains" its allosteric parent (logs/contain_summary.log:57).
- Containment of the flagged:
  - 40 sampled certified variants pass both ways at fixed L, ≤1.6e-31 (logs/contain_summary.log:4).
  - The 3 uncertified ones (#6010, #6018, #6716; slow steps join free E's segment):
    - variant ⊂ L→0 law to ≤2.7e-33 (logs/contain_summary.log:2);
    - L→0 law ⊂ variant at L = 0.01 in 9 of 9 trials, ≤7.3e-33 (logs/deep.log:1-3,13-15,25-27);
    - at L = 10 it fails in 9 of 9 trials, 9.7e-3 to 2.2e-1 (logs/deep.log:10-12,22-24,34-36);
    - so the families are equal as a union over L, not at fixed L. Inheritance stays right: the variant cannot beat its L→0 law, and it reaches that law's loss as L → 0.
- Re-test of the 14 of 99 (key.md:308):
  - The originals cannot be identified. The verifier's log and the census data were lost in the outage (RECOVERY.md).
  - Regenerated class: 6,567 numeric L-phantom candidates, the census's label: L in the span of the other columns (logs/summary.log:11).
  - Against the L→0 law:
    - 3,971 are certified equal;
    - 1,913 gain rank over it, so their families are strictly larger. Sampled 30 of 30 "contain" it (logs/contain_summary.log:20). This is the "11 contain the parent's law" kind.
    - 683 have equal rank and no certificate. Of 60 sampled, 47 are equal at fixed L, 10 only as L → 0, and 3 strictly contain the L→0 law (logs/contain_summary.log:5-18).
    - logs/summary.log:28-35.
  - Against allosteric enumeration parents, as the original verifier compared them (40 sampled): 18 contain the parent, 20 are the same, 2 differ at 1e-12 (logs/contain_summary.log:51).
  - Against plain enumeration parents: 20 of 20 are the same (logs/contain_summary.log:74).
  - The spec's rule flags none of the rank-gaining or containing cases. Every "contains" or "differs" verdict on a flagged variant comes from comparing with the enumeration parent.
  - No "inside" case recurred in 60 parent comparisons. Phase 2's fixes of known issues 1 and 2 may explain this; not checked.

## 3. R5.3: which fit a phantom inherits

| move | phantoms | source | where |
|---|---|---|---|
| to_allosteric | 660 | enumeration parent | earlier round (logs/sources.log:5-6) |
| change_allo_state | 81 | tagging ancestor | earlier round (logs/sources.log:1) |
| re_to_ss | 110 | the parent's plain sibling (the ancestor's `re_to_ss` child) | earlier round (logs/sources.log:2) |
| split | 57 | the parent's plain sibling | earlier round (logs/sources.log:4) |
| re_to_ss | 18 | never fitted; 8 distinct L→0 keys, np 7 | (logs/summary.log:41,45) |

- No source sits only in the same round or at a non-relative (logs/summary.log:38-44).
- Under a beam, an ancestor was selected and expanded, so its children were fitted. The source therefore always exists before the phantom does.
- The 18 never-fitted phantoms:
  - L→0 law: np₀ 6, rank 5. Variant: np 7, rank 5.
  - The plain move never emits the catalytic mechanism. Example lineage: #3268 ← #303 ← #4. Plain `re_to_ss(#4)` gives only #299 and #300.
  - Cause: `_chain_flank_groups` skips allosteric mechanisms (proto3_wt src/mechanism_enumeration.jl:459-468, "An allosteric mechanism has none"), yet at n = 1 the flank flip stays invisible.
- Under the certificate rule: 3,818 phantoms. 3,756 have a relative source; 62 are never fitted, on 16 keys, all the same chain-flank kind (logs/extra.log:12-19).
- Reachability: 93 non-phantom allosteric round-2 candidates, all `change_allo_state`, have only phantom parents among the 1,000 sampled (logs/extra.log:29). Phantoms must stay expandable (SPEC:650).

## 4. Proposed rulings (key.md:240-275)

1. **Rank method (item 1).** Use the certificate in `src/` instead of the modular rank. Use the modular or BigFloat rank and numeric containment only as test oracles.
   - Why:
     - Inheritance needs equal families, and equal rank does not imply them (logs/contain_summary.log:14,15,18). The certificate proves equal families.
     - It is deterministic, takes 25 µs median against 1.3 ms for one modular draw (logs/timing.log:2-3), and takes about 45 lines in `scripts/plib.jl` against about 135 for the modular rank.
     - It confirms 923 of the spec rule's 926 flags.
   - Cost if wrong: the 3 union-over-L phantoms get fitted. If Denis keeps the rank, the code needs a prime field, the branch Jacobian and 2 draws.
   - This amends SPEC:651-653.
2. **Which parent (item 2).** Option (b): the source is whoever holds the key of the variant's L→0 law, the law of `catalytic(am)`. When that key has no fit, use b2: fit the variant as an ordinary unflagged candidate.
   - Why:
     - The enumeration parent is right only for `to_allosteric` (660 of 926).
     - b2 needs no fit outside the candidate set and no `fit_source_id` without a CSV row.
     - b2 costs at most 18 fits per 8,075 variants (8 keys) under the spec's rule, or 62 (16 keys) under the certificate rule.
   - Cost if wrong (Denis picks b1): a PASS2 path that fits `catalytic(am)` as the key's representative, saving at most 10 fits per round-2 sample.
   - Separately for Denis, a ROADMAP candidate: apply the chain-flank exclusion at one catalytic subunit, which removes all 18 (62).
3. **Name map (item 3).** Option (a), canonical labelling: read the stored canonical-position values through the canonical order of the variant's active law, and report L = 0.
   - Why: the key cache may hold another member of the key. Keys match on 3,818 of 3,818. The slot map (b) also holds on 3,818 of 3,818 but maps only onto `catalytic(am)` (logs/extra.log:2-4).
   - Cost if wrong: none observed.
4. **Beam across sweeps (item 4).** A phantom ranks under its own key and its own `length(fitted_params)`, with its source's loss. Selected keys persist across sweeps:
   - a late member of a selected key expands without spending budget;
   - a late member of a rejected key stays rejected.
   - Why: sources always come first (section 3), phantoms must expand (93 dependants), and this matches the monotone cutoff.
   - Cost if wrong: a subtree lost to a late member. Not measured; no identify run here.
5. **Cross-validation for phantoms (item 5).** Fold losses are computed once per source key and shared. A phantom whose source is in the pool loses the 1-SE tie, having one more parameter and the same fold losses. It is chosen only when its source is absent, and is then reported with L = 0 and `fit_source_id`.
   - Why: equal families give equal fold optima (SPEC:657-658).
   - Cost if wrong: a slot at n + 1 holds a phantom.
6. **Pre-filter scope (item 9).**
   - Under the spec, "binds nothing" means the inactive graph is empty and every site ligand is `:OnlyA`:
     - the slice's `inactive_mask` is empty (slice/src/thermo.jl:44-60);
     - at phase 3's tip f8eb0df that is `isempty(_inactive_groups(am))` (proto3_wt src/rate_eq_derivation.jl:607);
     - after S3f it is `!any(_inactive_mask(am))` (S3_design.md §7).
   - Iso free forms count as binding, so they are fitted, until iso mechanisms land (P4-2).
   - Why:
     - The inactive mask is empty exactly when no inactive step binds, on all 18,915 variants (logs/summary.log:20).
     - The site clause changes no verdict: 0 of 6,843 site-bearing variants are flagged, and the only 3 rank-equal ones fail the pre-filter anyway (logs/summary.log:6-10).
   - Proposed instead, which Denis must approve: drop the pre-filter and flag every admitted variant with a certificate. The certificate already requires QI = 0 and handles sites.
     - Why: about 20% of fits removed instead of about 3% (logs/fits.log:2,3,5), each with a proof, and 20 sampled variants contain both ways (logs/contain_summary.log:19).
     - Cost if wrong: if Denis keeps the pre-filter, the certificate still serves as the exact test inside it (ruling 1). This amends SPEC:659-660.
7. **Rank mode (item 10).** Absolute.
   - Why: relative-mode flags equal absolute-mode flags on all 18,915 variants (logs/summary.log:19). Under ruling 1 the mode does not arise.
   - Cost if wrong: none.

## 5. On phase 3's types

- L→0 law: the law of `catalytic(am)` (proto3_wt src/reaction_graph.jl:756). The slice had to strip tags from the census text instead.
- Slot map: tags change neither form nor edge order. Phase 3 stores the catalytic `Mechanism` as a field, so edge k is the same edge in both.
- Certificate inputs, all in phase 5's allosteric law: the active Q and D, the flip form's weight (free E is `forms[1]`; slice `trees[1]·shift`), fA, fI, QI, DI and the site polynomials.
- Return type: return a `Bool` (R18).
- Not covered:
  - ter reactions, n > 1, and iso mechanisms;
  - phase 4's B v3 guard (the slice's admission stood in for it);
  - children as phase 4's tagging (P4-1) will produce them;
  - round 2 of ping-pong;
  - the 1,752 bi-bi round-1 parents outside the sample.
