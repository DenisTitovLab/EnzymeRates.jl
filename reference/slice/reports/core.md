I've built the core stage of the slice. Every check passes on all 4,917 converted mechanisms. One deliverable is missing: the harness blocked writing `slice/core/REPORT.md`, because subagents may not write report files, so the report is this message.

**Files** (code lines from `proto/p5b/codelines.py`; I hand-counted `convert.jl` and got the same 43):

| File | Code lines | Contents |
|---|---|---|
| `slice/src/Slice.jl` | 10 | Loads the core files, then each later stage file that exists, in pipeline order |
| `slice/src/types.jl` | 434 | Mechanism type and structure-only constructor 207, accessors 28, orientation 57, RE forest 27, boundary checks 63, Sig 52 |
| `slice/src/names.jl` | 88 | Naming chokepoint |
| `slice/src/poly.jl` | 151 | Polynomials: rational exponents, checked Int with a BigInt fallback, exact division by the Haldane factor |
| `slice/src/convert.jl` | 43 | Pop record → mechanism, plus the constant map to the oracles' convention |
| **Total** | **726** | |

Also `slice/INTERFACES.md`, the converted sets in `slice/core/data/` (with `MANIFEST.tsv`), and the check scripts and logs in `slice/core/`.

**Results**
- **Conversion:** 4,917 of 4,919 records in 27 sets. That covers all 9 shared populations, the 4 skeptic P5a populations and its adversarial cases, P5a's grid cases, P1's sets plus D1–D5 rebuilt from its script, and P2's sets plus its hand-built cases. The two that fail are hand-built cases the package also rejects, and my constructor gives the same reasons.
- **Not converted:** two P1 files are result tables, not mechanisms. P2's 1,200 fuzz mechanisms were never stored (only 15 were dumped).
- **Orientation:** matches the package's stored direction on 38,263 of 38,263 package-built steps (tier 2 decides 16, tier 3 by form index 55).
- **Equality and hash classes:** 4,862 package-built records give 3,947 classes on both sides, with 0 false splits, 0 false merges and 3,947 distinct hashes.
- **Variants and mutations:**
  - 29,502 rewritten variants (reordered, reversed, relabelled) are all equal with equal hashes, and equal variants compile to the identical type.
  - 39,133 single mutations of ten kinds: 7,746 rejected, 31,387 differ, none equal or sharing a hash.
- **Hash stability:** 4,917 of 4,917 hashes are identical across two processes. A control struct without an explicit hash does change between them, so the check can detect instability.
- **Sig:**
  - The round trip `decode(sig(m)...) == m` holds for all 4,917. The Sig is a median 125 integers, at most 818.
  - Encoding the 239 init_bibi Sigs compiles 556 method instances and 143 KB (budget 4,000 and 1,000 KB; the package today 1,888 and 508 KB).
  - Encoding 2,369 more mechanisms with 49 new Sig lengths adds 1 instance.
  - A generated function that decodes its Sig compiles exactly one instance per fresh mechanism, its own; the decoder is not recompiled.
- **Names:**
  - 77,289 names compared with the package's per-step names: 0 differences unexplained.
  - Each difference comes from a decision:

| Change | Names |
|---|---|
| Residual `_res_+A_-P` → `_res_plusA_minusP` | 1,213 |
| Inhibitor-only names lose their suffix (`K_EI_to_E_I`) | 1,635 |
| A reactant's modifier copy prints `Amod` instead of `Ainh` | 38 |
| Inactive `:OnlyA` members name themselves (untied, decision 4) | 4,636 |
| `:EqualAI` regulator names become bare (`K_Rreg`) | 3 of 149 fitted |

- **Constant convention against the oracle:** on its thermodynamically consistent constants, every edge's forward equilibrium ratio matches (38,449 edges, worst relative error 3.1e-154), and the edge ↔ Pop step map round-trips exactly.
- **Boundary checks:** no false rejections on package-built records. All 16 hand-built bad inputs are rejected as intended, and atom conservation holds on all 4,612 enumerated records.
- **Polynomials:** 2,000 random trials with 0 failures. BigFloat error at most 6.9e-151; Int overflow throws and the BigInt retry works.
- **Costs (median):** constructor 1.78 µs (p90 2.95, max 8.8), orientation 83 ns, hash 174 ns, equality 21 ns, encoding 2.9 µs, decoding 3.7 µs. P6 measured 1.07 µs for its constructor; I did not re-measure the package.

**Deviations and choices for you to confirm**
1. **Form order:** forms sort by bound count, then residual count, then content. Free forms come first, so each RE segment's root has the fewest bound metabolites, matching today's root rule, and the first free form is the natural flip-form candidate. Tier 3 by index still matches the package everywhere.
2. **Name formats I propose:** residuals print `_res_plusA_minusP`, and a reactant bound as a modifier prints `Amod` (`A::Mod` in the DSL).
3. **Naming representative:** I kept today's tie-break by rendered form names. It changed no representative on this data.
4. **Group check:** the constructor computes the orientation so it can check group uniformity in each member's physical-forward direction, as the package does (decision 2). P6 only checked "up to direction".
5. **Additions beyond the decision list:**
   - The constructor rejects disconnected steps; the package does not.
   - `check_atoms(m, atoms)` takes an atom table, since the type has no reaction field.
   - `pop_constants`/`edge_constants` give later stages one tested bridge to the oracles.

**My mistakes, all fixed:**
- My first decoder captured the Sig tuple in closures, which compiled 72 instances per new Sig length.
- One measurement counted the compilation of its own print arguments as Sig cost (2,288 instances).
- Two check-script bugs:
  - The scripts compared retagged mechanisms with un-retagged Pop records, giving 192 false variant failures.
  - A merge mutation only relabelled groups, giving 400 false "equal" results.
- The first constructor took 4.8 µs because of boxed captured variables.

No Julia process of mine is left running.

**What later stages must know**
- Include `Slice.jl` into `Main` before deserializing. Key records by (set, index), because ids repeat across sets; `MANIFEST.tsv` gives each set's Pop file.
- Edge `k`'s constants are `constant_names(m)[k]`: `[K]` or `[kf, kr]` in the physical-forward direction of `oriented(m, k, orient(m))`.
  - `kr` is the "against" constant.
  - A binding's `K` is its dissociation constant.
  - Groups tie slot by slot; in the inactive state `:EqualAI` has no copy, `:NonequalAI` one per group, and `:OnlyA` one per member.
- Compare against the oracles through `pop_constants`/`edge_constants`.
  - The hand-built sets store steps in arbitrary directions and have no package outputs.
  - In rejected P5a cases (`ok = false`), the Pop names and outputs come from a `:NonequalAI` twin.
- `re_forest` roots each segment at its lowest form and scans edges in content order; the engine should reuse it.
- In generated functions, call `decode(S...)` only, and never pass the Sig to another function: anything that sees the concrete tuple type compiles once per length.

Files are in `/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/slice/`:
- `src/Slice.jl`
- `src/types.jl`
- `src/names.jl`
- `src/poly.jl`
- `src/convert.jl`
- `INTERFACES.md`
- `core/data/MANIFEST.tsv`
- `core/logs/`