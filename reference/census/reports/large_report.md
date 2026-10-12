I ran the census on ter-bi, ter-ter and pyruvate carboxylase. Pyruvate dehydrogenase was skipped: its `init_mechanisms` alone took 397 s and 3.6 GB. Julia time was about 45 min.

## The answer in brief
- **Seeds (round 0):** 99.2–99.4% of fitted seeds fall in no class. Classes D, N and Z are empty there; W and P are each under 1%.
- **Expansion rounds:** 37–57% of fitted children are W, D, P or Z. Almost all of them are one-subunit K-type allosteric variants from `_expand_to_allosteric`:
  - these variants are 73–79% of fitted round-1 children;
  - 51–75% of them are W/D/P/Z;
  - for plain flip and split children the rate is 0.6–2.8%.
- **Worth filtering, cheapest first:**
  - W_admit, the island rule: 6.6–35% of expansion fits, a few µs per check.
  - D, the key duplicates: 22–40% of expansion fits, about 0.6 ms per candidate for law plus key.
  - P with a one-subunit L phantom: 14–27% of expansion fits. In all 47 sampled cases a smaller candidate in the run has the same law (46 times the plain parent).
- **Not worth a filter:** N is 0 everywhere. All 18 Z are also D. The remaining P and W kinds are each under 3%.
- **Today's eq_hash dedup does almost nothing:** it removes 8 of 12,069 compiled candidates. The key would remove 2,046 more. No eq_hash merge is split by the key (0 of 8).

## What ran
- **Pipeline.** I mirrored `_beam_search` on branch `onlya-guard-stranded-segment` (11cac53):
  - Seeds: a deterministic sample of 2,000 (`Xoshiro(20261011)`), then `_base_tier`.
  - Batch stage as in `_process_batch`: structural `seen`, `_eq_complexity > 337`, fitted count > 20, `rate_equation_string` → `_rate_eq_dedup_key`. One representative per new eq_hash is "fitted".
  - Rounds 1–2: random parents drawn from the compiled candidates. Each parent went through `expand_mechanisms` alone inside try/catch, and children were deduplicated with the first parent winning.
- **Fidelity checks.**
  - `seen` also held every base-tier structure of the full seed set (5,158; 271,561; 444,920 structures). No child matched one.
  - For ter-bi only, I rendered the full base tier (4,952 eq_hashes). No fitted child shared an eq_hash with it.
- **Classes**, computed on the slice law (`round4/guard` checks plus slice `admit`, `canonical`, `edge_flux`):
  - **W:** RE turnover; bottomless segment, split or trap (`face_check`) on the faces "all products 0", "all substrates 0" and each single reactant 0; D ≡ 0 on a face; N ≡ 0; `admit`.
  - **D:** same key as a candidate fitted earlier in the run.
  - **P:** rank of dv/dlogθ from exact derivatives at 2n+6 random points, SVD tolerance 1e-9. I report it twice:
    - "P" uses the slice law and slice count (the rewrite).
    - "Ptoday" uses the package's own law, evaluated from its Exprs with complex-step derivatives, against the package count. I computed it on all 1,812 candidates where the slice law is deficient or the counts differ, plus 1,161 random others; all random ones agreed with the slice.
  - **N:** no Vmax, forward or reverse.
  - **Z:** same key as the parent.

## Counts per reaction and round
Percentages are shares of fitted candidates. The last two rows weight the W/D/P/Z union by summed fitted parameters and by summed denominator terms.

| | ter-bi R0 | ter-bi R1 | ter-bi R2 | ter-ter R0 | ter-ter R1 | PC R0 | PC R1 |
|---|---|---|---|---|---|---|---|
| Input | 2,000 of 4,885 seeds | 60 parents → 1,298 | 30 parents → 584 (581 distinct) | 2,000 of 250,855 | 60 → 2,390 | 2,000 of 393,996 | 60 → 2,021 |
| Seed screen (degenerate seeds removed / flip children added) | −65 / +179 | – | – | −95 / +243 | – | −174 / +430 | – |
| Structural `seen` (run / full base) | 0 / – | 0 / 0 | 0 / 0 | 0 / – | 0 / 0 | 0 / – | 0 / 0 |
| Complexity filter (>337) | −91 | −3 | −1 | −294 | −30 | −288 | −32 |
| Parameter cap and compile errors | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| eq_hash dedup | 0 | 0 | −8 | 0 | 0 | 0 | 0 |
| **Fitted** (allosteric) | **2,023** (0) | **1,295** (942) | **572** (498) | **1,854** (0) | **2,360** (1,854) | **1,968** (0) | **1,989** (1,497) |
| W | 0 | 105, 8.1% (all admit) | 39, 6.8% (38 admit, 1 split) | 8, 0.4% (split + trap) | 434, 18.4% (433 admit) | 7, 0.4% (split + trap) | 697, 35.0% (admit) |
| D | 0 | 288, 22.2% | 172, 30.1% | 0 | 786, 33.3% | 0 | 800, 40.2% |
| P (slice) / Ptoday | 17 / 17, 0.8% | 309 / 319, 23.9% / 24.6% | 156 / 169, 27.3% / 29.5% | 6 / 6, 0.3% | 382 / 449, 16.2% / 19.0% | 4 / 4, 0.2% | 275 / 314, 13.8% / 15.8% |
| N | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| Z | 0 | 0 | 18, 3.1% (all also D) | 0 | 0 | 0 | 0 |
| **OK** | **99.2%** | **62.7%** | **50.3%** | **99.2%** | **52.8%** | **99.4%** | **43.5%** |
| W/D/P/Z share by fitted parameters | 1.0% | 38.6% | 50.3% | 0.8% | 49.5% | 0.6% | 58.6% |
| W/D/P/Z share by denominator terms | 0.4% | 25.5% | 32.8% | 1.1% | 33.1% | 0.7% | 41.6% |

Other details:
- **Fitted sizes:** median parameter count is 6–9 and median denominator terms 22–42. The largest law has 4,896 denominator terms (PC round 1).
- **Expansion errors:** `expand_mechanisms` threw on 0 of 210 parents, so I did not reproduce the 5,000-term crash. The parents were within the complexity filter, as in identify.
- **Exclusive combinations, ter-ter round 1:** OK 1,246; D 334; D+W 198; D+P 183; P 163; W 131; D+P+W 71; P+W 34. P alone, which would survive W and D filters, is 163 (6.9%). It is 12.9% in ter-bi R1, 18.0% in ter-bi R2 and 4.6% in PC R1.

**Plain versus allosteric children** (W/D/P/Z share):

| | ter-bi R1 | ter-bi R2 | ter-ter R1 | PC R1 |
|---|---|---|---|---|
| Plain children | 4 of 353 | 2 of 74 | 14 of 506 | 3 of 492 |
| Allosteric children | 479 of 942 (50.8%) | 282 of 498 (56.6%) | 1,100 of 1,854 (59.3%) | 1,121 of 1,497 (74.9%) |

In ter-bi round 2, the 59 fitted allosteric-state-change children were 68% P.

## Class details
**W.**
- The island rule makes up 1,273 of the 1,290 W cases. All come from tagging ping-pong parents, with the message "forms E(; residual …) reach the free enzyme only through :OnlyA steps that are not bindings".
- Face split/trap: 17 of 12,061 fitted. 13 are original seeds (7 ter-ter, 6 PC), 2 are seed-screen flip children and 2 are expansion children.
- Not found anywhere: RE turnover, bottomless segments, D ≡ 0 on a face, v ≡ 0, findings only on single-reactant faces.

**D.**
- 2,025 of 2,046 are same-round sibling K-type variants of one parent, allosteric with allosteric.
- They differ only in the tag of a group the inactive conformation cannot reach anyway.

**P.**
- **Seeds:** 27 cases (16 RE-isomerization gauge in ter-bi, 6 in ter-ter, 4 in PC; 1 other). All 27 are flip children that the seed screen puts in place of degenerate seeds.
- **Expansion:** 1,047 L phantoms, 60 RE-isomerization gauge, 14 flux-free steady-state, 45 other, plus 132 candidates deficient only in the package law.
  - The 14 flux-free cases are :NonequalAI steady-state groups in a catalytically dead inactive conformation. Only the kf/kr ratio of their inactive copy is identifiable.
- **Same-law sample** (Levenberg–Marquardt fit of coefficient ratios):
  - L phantom: 47 of 47 have the same law as a smaller candidate in the run (46 the parent).
  - RE-isomerization gauge and other: 0 of 23 have one.

**N.**
- Zero among fitted candidates.
- The seed screen's removed seeds are all condition C: 65 of 65 in ter-bi and 95 of 95 in ter-ter. None lacks a Vmax.
- I did not finish the PC split because the lock was busy. The law test still finds a Vmax in all 174 PC degenerate seeds, and none of them fails any W check.

**Package law versus slice law (side finding).**
- In 902 fitted children (8–17% of expansion fits), today's package law has one more identifiable dimension than the declared law. An example of a package-only name is `K_EBQ_to_EQ_B`; this is ROADMAP known issue 1.
- The package count exceeds the slice count in 250 fitted children.

## Examples (step lists)
- **W_admit** (ter-bi (1,1157), 7 parameters): `E + P <--> E(; residual = A + D - Q) + B [OnlyA]; (E + A <--> E(A), E(Q; residual = A + D - Q) + A <--> E(A, Q; residual = A + D - Q)) [OnlyA]; (E(; residual = A + D - Q) + Q <--> E(Q; residual = A + D - Q), E(A) + Q <--> E(A, Q), E(Q; residual = A + D - Q) + Q <--> E(A, D, Q)) [OnlyA]; (E(A) + D ⇌ E(Q; residual = A + D - Q), E(Q; residual = A + D - Q) + D ⇌ E(D, Q; residual = A + D - Q), E(A, Q) + D ⇌ E(A, D, Q)) [OnlyA]`
- **W split** (ter-bi (2,460), plain; verdict "S=0 and A=0: surviving forms E and E(P) split"): `(E + A ⇌ E(A), E(D) + A ⇌ E(A, D)); (E + D ⇌ E(D), E(A) + D ⇌ E(A, D), E(P) + D ⇌ E(D, P), E(A, P) + D ⇌ E(A, D, P)); (E + P <--> E(P), E(D) + P <--> E(D, P), E(Q) + P <--> E(P, Q)); (E + Q <--> E(Q), E(P) + Q <--> E(P, Q)); (E(A) + P ⇌ E(A, P), E(A, D) + P ⇌ E(A, D, P)); (E(P) + A ⇌ E(A, P), E(D, P) + A ⇌ E(A, D, P)); E(A, D) + B <--> E(P, Q)`
- **W trap/split in a seed** (ter-ter (0,473)): "P=0: trap {E(A, P)} entered; R=0: E(A, B) and E(A, P) split". The full step list is in `examples.log`.
- **D** (ter-bi (1,12) duplicates (1,8)): `(E + A ⇌ E(A), E(Q) + A ⇌ E(A, Q)) [OnlyA]; (E + Q ⇌ E(Q), E(A) + Q ⇌ E(A, Q)) [EqualAI]; E(A) + D ⇌ E(A, D) [OnlyA]; (E(Q) + B ⇌ E(B, Q), E(A, D) + B ⇌ E(A, B, D)) [EqualAI]; E(Q) + P ⇌ E(P, Q) [EqualAI]; E(P, Q) <--> E(A, B, D) [OnlyA]`. The two differ only in the tag of `E(A) + D ⇌ E(A, D)`, and E(A) is unreachable in the inactive state.
- **P, L phantom** (ter-bi (1,9), 7 parameters, rank 6 in both laws): `(E + A ⇌ E(A), E(Q) + A ⇌ E(A, Q)) [EqualAI]; (E + Q ⇌ E(Q), E(A) + Q ⇌ E(A, Q)) [OnlyA]; E(A) + D ⇌ E(A, D) [EqualAI]; (E(Q) + B ⇌ E(B, Q), E(A, D) + B ⇌ E(A, B, D)) [EqualAI]; E(Q) + P ⇌ E(P, Q) [EqualAI]; E(P, Q) <--> E(A, B, D) [OnlyA]`
- **P, RE-isomerization gauge** (ter-bi (0,1963), a seed-screen flip child, 8 parameters, rank 7): `(E + A ⇌ E(A), E(Q; residual = A + D - Q) + A ⇌ E(A, Q; residual = A + D - Q)); E + P <--> E(P); E(; residual = A + D - Q) + B ⇌ E(B; residual = A + D - Q); (E(; residual = A + D - Q) + Q ⇌ E(Q; residual = A + D - Q), E(A) + Q ⇌ E(A, Q)); E(A) + D ⇌ E(A, D); E(P) ⇌ E(B; residual = A + D - Q); E(Q; residual = A + D - Q) <--> E(A, D)`
- **P, flux-free steady-state step** (ter-bi (2,214)): `(E + A <--> E(A), E(Q) + A <--> E(A, Q)) [NonequalAI]; (E + Q ⇌ E(Q), E(A) + Q ⇌ E(A, Q)) [OnlyA]; E(A) + D ⇌ E(A, D) [EqualAI]; E(Q) + P <--> E(P, Q) [EqualAI]; E(A, D) + B ⇌ E(P, Q) [OnlyA]`
- **Z** (ter-bi (2,95)): it splits `E(Q) + P ⇌ E(P, Q)` out of its parent (1,68)'s P-binding group. The key is unchanged, but the package's parameter count and law dimension both rise from 7 to 8.

## Cost of each check
Warm timings in µs on 120 random fitted ter-ter candidates (median / p90). ter-bi is about half; PC is similar.

| Check | Median | p90 |
|---|---|---|
| `admit` (0 for plain mechanisms) | 0–2 | 1.9 |
| RE turnover | 2 | 2.5 |
| Edge flux | 2 | 2.5 |
| Face check, 8 faces | 24 | 47 |
| `check_boundary` | 50 | 92 |
| Slice solve and law | 295 | 1,539 |
| Key | 318 | 3,443 (max 10,196) |
| Support tests | 8.5 | 60 |
| Saturation | 9 | 65 |
| Rank | 18 | 101 |
| Package `rate_equation_string` (today's cost before every fit) | 2,025 | 6,590 |

A fit is 20 optimizer restarts over the data, so seconds or more. I did not time fits.

## Not checked
- Parents were chosen at random, not by loss as the beam chooses them. Round 2 was run only for ter-bi.
- The full-run eq_hash memo was checked for ter-bi only.
- The 17 face split/trap verdicts are structural only; I did not run the RE-limit oracle on them. The round-4 validation of that check had 6 false alarms in 292.
- I checked no corner faces (for example one substrate and all products at zero), which is where condition-C seeds would show trouble.
- Rank uses one parameter point per candidate. 21 candidates have a weakest direction under 1e-6 (minimum 4.2e-9) and count as full rank; dropped directions were at most 4.1e-16.
- The same-law result rests on a sample of 70.
- No fits were run.
- Your question 1 about the two specs is not part of this census track.

## Files
All under `/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/census/large`:
- **Scripts** (`scripts/`): `census_lib.jl`, `driver.jl`, `run_one.jl`, `summarize.jl`, `report_saved.jl`, `moves.jl`, `examples.jl`, `samelaw.jl`, `costs.jl`, `degenerate.jl`, `extra.jl`, `analyze.sh`.
- **Logs** (`logs/`):
  - Runs: `ter_bi.log`, `ter_ter.log`, `pc.log`.
  - Analyses: `report_all.log`, `moves.log`, `examples.log`, `samelaw_*.log`, `costs.log`, `degenerate.log`, `extra.log`.
  - Smoke test and bi-bi check: `probe.log`, `bibi.log`.
- **Data** (`data/`): `ter_bi_rxn.jls`, `ter_ter_rxn.jls`, `pyruvate_carboxylase_rxn.jls`.