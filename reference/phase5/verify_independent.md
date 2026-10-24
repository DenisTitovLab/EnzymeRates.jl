# Independent recheck of research5 (law, unit, phantom, keycost)

Paths are relative to `scratchpad/phase45/research5/`. My code and logs are in `verify/`
(`verify/logs/X:N` is line N of a log). I wrote all of it; I reused only the subject code (the
slice, `law/variants/{U,B,C}`), the census adapters and `phase3/gcprobe/hoistscan.py`.
Environment: `verify/env` develops `graph_mechanism_wt` 7205ba4, whose src equals 4f9a34a
(`diff -rq`). Julia run time was about 20 min; the logs' `elapsed` values include lock waits.

## 0. Answer

- **Blocking: law type C (law.md L1) is not GC-safe by construction.** One caller loop crashes
  it on both Julia versions. Making `Poly` a `mutable struct` (E) or a `Dict` alias (D) removes
  the crash at the same speed.
- **Major:** the study's scan workload does not type its loops, so its 0-hoist result does not
  exercise the caller pattern behind R18. L7 as worded would pass the crashing code.
- **Confirmed with fresh inputs:** the phantom numbers, the unit rule, the R5.5 count and
  keycost's F9 mechanism.
- **Minor:**
  - phantom.md counts one population twice;
  - unit.md's allosteric coverage is thin (no disagreement).

## 1. Blocking

**B1. The C `Poly` keeps the bug's layout, and a loop that merges it with another type aborts the GC.**

- Crashing loop: `verify/crash_c.jl:19-28`, on 989 mechanisms from `verify/data/scan_pop.tsv`:
  ```julia
  for m in ms
      L = S.rate_law(m)
      y = L isa S.Law ? L.H : L.sites      # Union{Poly{Int}, Vector{…}}
      t += y isa S.Poly ? S.nterms(y) : length(y)
  end
  ```
- C, Julia 1.12.5:
  - aborts with exit 134 and "GC error (probable corruption)" (`logs/crash_C_112.log:3034, 3089`);
  - the trace shows `gc_mark_outrefs` at gc-stock.c:2364, the probe's signature (`:3058`);
  - hoistscan finds 1 hoist with a bits store into a reference field (`logs/hoistscan_all.log:28`).
- C, Julia 1.13.1: exit 134 (`logs/crash_C_113.log:3300, 3356`); 1 dangerous hoist (`hoistscan_all.log:32`).
- The same loop, 3 rounds clean and 0 hoists on both versions:
  - D (B's Dict-alias `poly.jl` with C's Int-only engine): `crash_D_112.log:2998-3002`,
    `crash_D_113.log:3263-3267`;
  - E (C with `mutable struct Poly`, `const terms`): `crash_E_112.log:3272-3276`,
    `crash_E_113.log:3551-3555`.
- Cost: law time E 0.992× and D 0.999× of C on 989 laws, best of 5 interleaved (`logs/time_ce.log:1-3`).
- In the full workload, typed C has 1 dangerous hoist of 496 on 1.12 and 1 of 569 on 1.13
  (`hoistscan_all.log:4, 14`).
  - It sits only in the merge loop: compiling the loops one at a time puts it between the
    `merge_vec` markers (`logs/scan_loc_C_112.log:9065, 9101, 9102`).
  - Derivation, `merge_nothing`, `merge_same` and the identify pass hoist nothing under C.
- The tested C is not L1's proposed type:
  - `rate_law` returns `Union{AlloLaw{Int64}, Law{Int64}}` (`logs/scan_C_112.log:474`), so L1's
    "deletes the retry and the only Union" is false;
  - `Law`, `AlloLaw` and `Poly` are still parametric;
  - `inactive::Union{Nothing, Law{T}}` remains (`law/variants/C/allosteric.jl:19`).
- R18's second clause already forbids this loop. C's safety therefore rests on every future
  caller and on a scan that never compiles such a loop (M1). D and E remove the box the bug
  needs.
- Reproduce:
  ```
  JULIA_LLVM_ARGS="-pass-remarks=julia-licm -pass-remarks-analysis=julia-licm" \
    phase45/jl_research.sh <log> 300 [+1.13] --startup-file=no verify/crash_c.jl C|D|E
  ```
  C exits 134; D and E exit 0.
- **Proposed ruling L1'.**
  - Ruling: keep C's arithmetic (checked `Int64`, no `BigInt` retry, L2 unchanged), but make
    `Poly` a `mutable struct` with a `const` field (E).
  - Why:
    - same speed;
    - no immutable reference-holding wrapper is left to box;
    - E passes the crashing loop on 1.12 and 1.13.
    - D also passes, but a `Dict` alias puts the arithmetic methods on Base's `Dict`, which
      Aqua's piracy check would flag (inferred, not run).
  - Cost if wrong: one heap wrapper per `Poly`, with no time cost measured.

## 2. Major

**M1. The scan workload is untyped, so L7 cannot catch B1.**

- `law/scripts/scan.jl:38` iterates `deserialize(IOBuffer(bytes))`, whose inferred type is
  `Any`. `L` is therefore untyped, every field read is dynamic, and no union-split phi of
  `Poly` values is compiled.
- law.md 1.2 says the loop reads its fields "as the probe's loop did". That holds for the text
  but not for the types.
- My untyped run gives C 0 hoists of 493 (`hoistscan_all.log:2`). The typed run gives 1 of 496
  (`:4`). The typed U control gives 6 of 634 on 1.12 and 6 of 750 on 1.13 (`:8, 18`).
- **Proposed ruling L7'.**
  - Ruling: the PR's scan types its loops: a function barrier over a
    `Vector{Union{Mechanism, AllostericMechanism}}`.
  - It adds one adversarial loop that merges a law field with a field of another type
    (`verify/scan_v.jl:66-90`) as an in-workload positive control.
  - It runs `crash_c.jl`'s loop on the shipped type.
  - Why: an untyped workload reported 0 for code that crashes.
  - Cost if wrong: about 1 minute per Julia version.

## 3. Minor

- **m1. phantom.md counts bi-bi round 1 twice.**
  - The 188 flagged ping-pong round-1 variants have the same neutral texts as the 188 bi-bi
    ones: "identical texts 188", by a Python comparison of `phantom/data/{recs,facts}_{bibi,pingpong}.tsv`.
  - So 738 flagged variants are distinct, not 926.
  - The ping-pong fit-removal row (76/512, `phantom/logs/fits.log:5`) repeats bi-bi's (`:2`).
  - No allosteric ping-pong variant is admitted under the pre-P4-1 tagging, so "every
    population checked" carries no ping-pong evidence.
  - Ruling: the PR body reports 738 distinct.
  - Cost if wrong: none.
- **m2. unit.md: populations barely test the per-conformation rule on allosteric laws.**
  - All 92 conformations of 50 bi-bi round-2 allosteric laws have a constant term
    (`logs/ana_unit2.log:3-4`).
  - Only the study's ping-pong n = 2 population (170 graded of 2,593) and the fixtures exercise
    the virtual and ratio rules there.
  - No disagreement.

## 4. Confirmed (my own code; different inputs)

- **R5.5 (keycost §1): reproduced exactly.**
  - bi-bi uniform round 2: 19 candidates newly over by the sum, 0 by the max, of 34,574
    admitted (`logs/gen_a.log:5`; study `keycost/logs/r55_table.log:48`).
  - The slice recount of both conformations' V×τ matches on 22,069 allosteric candidates with
    0 mismatches, and also gives 19 (`logs/gen_b.log:1`; the study compared 1,808).
  - A fresh seed gives 16 of 33,804, 0.047% (`gen_a.log:9`).
- **R5.2 and R5.3 (phantom), fresh sample.**
  - Inputs: bi-bi round 2 from 300 parents disjoint from the study's 1,000 (overlap 0,
    `gen_a.log:10`). Rank: 256-bit, two θ draws, complete pivoting. The certificate is
    re-derived.
  - Counts (`logs/summ_phantom.log:2-7`):
    - 2,612 variants, 168 flagged, 167 with a certificate;
    - the uncertified one is #4168 (`re_to_ss`);
    - the study had 3 of 550 uncertified.
  - Flagged checks:
    - key(active) = key(L→0) for 168 of 168;
    - np = np₀ + 1 for 168 of 168.
  - Certificate maps hold to 5.2e-77 on 657 certificates.
  - The two θ draws agree on every variant. Pivots kept are at least 1.6e-7; pivots dropped
    are at most 1.1e-76.
  - Fits removed (`summ_phantom.log:10`):
    - spec rule: 71 of 2,525 new keys (2.8%; study 2.5%);
    - certificate rule: 435 (17.2%; study 19.1%).
  - Sources:
    - 80 of 80 `to_allosteric` phantoms at the parent;
    - 2 never fitted (#3222, #3223, `re_to_ss`, np 7), the study's chain-flank class;
    - the rest at relatives (`:11-17`).
    - #4772 and #4773 look like non-relatives under first-parent lineage. With every parent
      edge counted, their holders #1499 and #1509 are also children of their ancestor seed 72
      (`logs/gen_d.log:6-7`), so no disagreement.
- **R5.4 (unit).**
  - 200 ter-ter init laws drawn uniformly outside the study's 2,000 (`gen_c.log:1`;
    `gen_b.log:3` confirms that pick).
  - Rules: 124 constant, 76 graded with integer λ, 0 fractional, 0 ratio. The study had 61%
    and 39%.
  - 0 unit failures.
  - 28 laws are Cleland-eligible; N differs from Cleland's in 23, 82% (`ana_unit2.log:2`).
    The study had 257 of 292, 88%.
  - 150 more ter-ter laws and 92 allosteric conformations: 0 unit failures (`ana_unit.log:2-4`).
  - That first ter-ter draw was biased toward early init order, an error in my sampling, so I
    use it only for the unit check.
- **keycost F9 mechanism.**
  - The helper's loop shape reads 0.625 ns per call in all 5 batches for an inlined law with a
    division, and 1.25 ns behind `@noinline` (`logs/gate_loop.log:1-2`).
  - So the call is hoisted out of the loop, not removed batch by batch, as F9 says.
- **Not rechecked:**
  - law.md: speed and live bytes, Q9-Q11, Q14;
  - keycost: R5.6, R5.7, F4-F6;
  - unit.md: fixtures and printing.
