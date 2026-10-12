I ported the three moves to `slice/src/enum.jl`, which has 357 code lines by `proto/p5b/codelines.py` (I hand-checked `flip_children` at 17). On every plain parent all three moves produce exactly the package's children. On allosteric parents the flip and the dead-end move match too; the split differs on 143 of 760, and every difference traces to the package's own known count issues. I did not write `slice/enum/REPORT.md`, because subagents may not write report files, so this message is the report. The numbers are in `slice/enum/logs/compare_final.log`, `coverage_final.log` and `diag_split_b.log`.

**What I built**
- `flip_children(m)` ports `_expand_re_to_ss`, the package's only flip move (rapid equilibrium to steady state; the reverse appears only inside the split).
- `split_children(m)` and `dead_end_children(m, declared)` port the split and the dead-end move.
- Shared helpers cover flux, rapid-equilibrium segment offsets, the minimal-set search, chain flanks, the hyperbolic check, the copy rule and the split's divisions.
- The split's gain test calls thermo.jl's `solve`; the copy rule uses its `inactive_mask`.
- The harness is in `slice/enum/`: `pkg_moves.jl` runs the package's own moves in er_oracle, and `compare.jl` diffs the results. The `diag_*`, `pkg_counts.jl`, `pkg_issues.jl` and `coverage.jl` scripts do the diagnosis and coverage.
- I documented the interfaces in INTERFACES.md under "Enum".

**Agreement: parents with identical child multisets**

| Parents | Flip | Split | Dead end |
|---|---|---|---|
| init_bibi_pp, 264 plain (includes the 60 expand_bibi parents) | 264/264 | 264/264 | 264/264 (4,330 children) |
| dead1, 300 plain level-1 dead ends | 300/300 | 300/300 | 300/300 |
| p1_enum_inh, 300 plain | 300/300 | 300/300 | 300/300 |
| expand_bibi, 700 (460 allosteric, one subunit) | 700/700 | 632/700 | 700/700 |
| allo_bibi, 300 allosteric (two subunits) | 300/300 | 225/300 | 300/300 |

- The dead-end move declared I plus copies of A and Q.
- All 1,564 parents that have a core conversion convert identically to it (dead1 has none).
- Every child passed `check_boundary` and `check_atoms`, and rebuilding it through `raw`/`canonicalize` reproduced it.
- Every filter was exercised and agrees with the package:
  - the hyperbolic check removed 692 flips;
  - the chain-flank rule removed 74 flips;
  - the copy rule removed 13 split children and 501 dead-end children.
- The solve's independent-parameter count equals the package's count on all plain parents and plain split children.

**Diagnosis of the split differences**
- The 143 parents give 213 package-only and 235 slice-only children.
- I computed the gain verdict of every set of up to three units under both counts. Each count's verdicts rebuild its own observed children exactly.
- 906 verdicts differ. In every one the package's count exceeds the solve's by 1–2, and the P5b-style count agrees with the solve.
- Per-mechanism flags tie every excess to a ROADMAP known issue:
  - **Spurious package gains (773):** the excess is on the child, from issue 2 (a name folded in one state but still used by the other) in 761 and from issues 1 and 2 together in 12. Example: expand_bibi/163, where the fold is `K_EP_to_E_P => K_EBP_to_EB_P` in the active state only.
  - **Missed gains (133):** the excess is on the parent, from issue 1 (a second name from the pruned inactive graph) in 87 and issue 2 in 46. Example: expand_bibi/24, where the parent already fits `K_EBP_to_EP_B`, the very constant the split would add.

**Time per move (median per parent, same timing protocol on both sides)**
- **Flip:** 10–27 µs against the package's 45–316 µs, 4–12× faster.
- **Dead end:** up to 178 µs against up to 1,802 µs, 3.6–12× faster.
- **Split:** 1.2–4.8 ms against 0.09–1.0 ms, so 0.1–0.6× the package's speed.
  - The cause is in thermo.jl, not enum.jl. `cycles` calls `isassigned(path, e.a)` with an Int32 index, which takes Base's exception-based fallback: 40 µs per check against 6 ns with an Int. That makes up about 95% of `solve`.
  - Patched in my test process only (`Int(e.a)`), the split takes 65–340 µs. That is 0.8–1.0× the package's speed on plain parents and 5.7–12× faster on allosteric ones, with the same children on all 1,864 parents.

**Code lines, package (er_oracle) → slice**

| Piece | Package | Slice |
|---|---|---|
| Flip, with the hyperbolic check and chain flanks | 113 | 89 |
| Flux | 71 | 74 |
| Minimal-set search | 20 | 16 |
| Split, with the package's partition counter | 113 | 56 |
| Copy rule | 79 | 73 |
| Dead end | 92 | 49 |
| **Total** | **488** | **357 (27% fewer)** |

The package also leans on `_re_segment_extras` (43 lines in types.jl); I did not count it.

**Deviations from the package code**
- The copy rule checks every copy group of a dead-end child; the package checks only the new group of a plain child. This made no difference on the data.
- When a complex has several twins and none matches its composition, I take the first by content, not by rendered name. No difference on the data.
- The moves do not call `admit`, because admission.jl did not exist yet.
- The allosteric count is the decided joint solve (decision 4), which is why the split differs.

**Effect of the new type on the moves**
- **Easier:**
  - A flag or a group label never moves forms or edges, so flip and split children come from the trusted constructor in O(E).
  - Analyses index forms by integer instead of the package's `Dict{Species}`.
  - The derived physical-forward orientation reproduces the package's stored direction for context forms and the gauge.
- **Harder:**
  - Dead-end children add forms and a table entry, which shifts indices, so they must go through `raw`/`canonicalize` with a form-index-to-RawForm map.
  - With no reaction field, the declared inhibitor names and the atoms must come from the caller.
  - Int32 form indices are a trap, as the thermo.jl issue shows.

**My mistakes, all fixed**
- My first `copy_gauge` draft contained a garbled expression; I caught it on review before it ran.
- `@printf` rejected concatenated format strings.
- My function line counter double-counted adjacent one-line methods.
- A shared helper read `ARGS`.
- Julia's sampling profiler hung on this machine, so I killed it and benchmarked directly.
- In `coverage_final.log`, the allo_bibi "chain-flank 692" line comes from the hyperbolic override still being active; the true value is 0.

No Julia process of mine is left running.

**What later stages must know**
1. Changing `isassigned(path, e.a)` to `isassigned(path, Int(e.a))` (and the same for `e.b`) in thermo.jl's `cycles` would make every `solve` about 20× faster.
2. The driver must filter allosteric dead-end children with `admit`. That flips and splits cannot change admission is my inference; I did not test it.
3. The driver must also run `check_atoms` and `check_names` on children.
4. `bottomless` in types.jl can be rewritten on `segment_offsets`, saving about 10 lines.
5. INTERFACES.md lists the top-level names enum.jl takes, so later files avoid clashes.

Files are in `/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/slice/`:
- src/enum.jl
- INTERFACES.md
- enum/compare.jl
- enum/pkg_moves.jl
- enum/logs/compare_final.log
- enum/logs/coverage_final.log
- enum/logs/diag_split_b.log