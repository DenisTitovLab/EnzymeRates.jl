**Allosteric stage: decision 6 is built, checks (a)–(f) pass, and admission rejects exactly the 81 `:OnlyA` ping-pong mechanisms.** I did not write `REPORT.md` because subagents may not write report files; this message is the report, and the numbers are in the logs.

**What I built** (code lines from `codelines.py`; my hand count of `admission.jl` also gives 61)

| File | Code lines | Contents |
|---|---|---|
| `admission.jl` | 61 | `admit`: one flip form 12, components and arcs 8, G1 26, islands 15 |
| `allosteric.jl` | 100 | `AlloLaw` 10, `allo_law` 26, `cross_weights` 8, `kcat_mask` 2, MWC `kcat_terms` 54 |
| `print_allo.jl` | 61 | each conformation printed in full and in coefficient form |
| **Total** | **222** | |

I also wrote the "Allosteric" section of `INTERFACES.md`; the codegen stage already builds on `AlloLaw`. `run_all.sh` reruns every check.

**Results**
- **(a) Oracle:** all 2,995 allosteric records (2,865 admitted) match P5a's two-conformation oracle. Worst error is 3.8e-32 in BigFloat (the oracle's own floor) and 1.7e-15 in Float64. On 400 small mechanisms, an explicit 1–2-protomer network agrees to 1.6e-33.
- **Limit check:** for each admitted `:OnlyA` mechanism (2,811), I built a direction from component potentials and certified it exactly on the joint system's rows (2,979 of 2,979). Along it the full oracle converges to the deletion law linearly in ε (worst 3.4e-22 at ε = 1e-30). All 85 island rejections have a certified direction whose limit differs from deletion by 100%.
- **(b) Verdicts:** mine equal B1's structural class on all 2,841 population mechanisms. Rejected: expand_bibi 16, allo_bibi_pp 27, allo_terbi 38; everything else is admitted.
  - B1's full-direction runs give "limit = deletion" for all 2,738 admitted mechanisms they covered.
  - The one disagreement, allo_terbi/131, is a sampling artifact. B1 found "limit = deletion" in 13 sampled directions, but its covalent island hangs only off chemistry steps; with island potential −1 the limit differs from deletion by 100%.
- **(c) G1:** the graph test equals an exact rational LP on the joint rows in 2,995 of 2,995. It equals today's guard on P5a's grid (108 of 108). On the AD cases:
  - AD1 ×3 and AD11 are rejected (today's guard admits them); AD15 is admitted (today's guard rejects it).
  - **AD14 and AD14b are admitted, not rejected as you expected.** Their cones are empty only if group members are tied; with decision 6's per-step classes the LP finds a direction, and the skeptic's per-step limit equals deletion.
  - The one-flip-form rule rejects AD2a–e and AD5.
- **(d) Known issues 1 and 2:** none of the 2,995 laws reads a constant outside the parameters, and no two independents share a name. Three fitted constants (specs/40, specs/42) appear in no polynomial: they are speeds of steady-state steps that carry no flux, which decision 12 keeps.
- **(e) Fitted lists:** lengths match the solve stage (1,844 equal, 123 one fewer, 26 two fewer; the 848 set 791/47/10). Compared by kinetic class, lists move in 411 of 1,844: 373 exactly as in P5b, and 38 more because the slice's solve priority picks other dependents. They also move in 101 of 791.
- **(f) Time:** `admit` takes a median 1.0 µs (max 3.5). Solve plus `allo_law` takes a median 104 µs, p90 158 µs, max 0.78 ms.
- **Printing:** both forms evaluate to the law on all 2,865 admitted laws (worst 9.8e-154). Every printed denominator term and flip weight is dimensionless, and there are no `^1` or `L * 1` artifacts. The MWC dimer prints `(Q_A * D_A + L * Q_I * D_I) / (D_A^2 + L * D_I^2)` with `D_A = 1 + S / K_A_ES_to_E_S + P / K_A_EP_to_E_P`.
- **kcat:** all 3,044 corners match the law evaluated near saturation. It equals the package's sampled kcat on 303 of 305 same-model mechanisms.
- **Enumeration:** the three ported moves produced 24,503 children from 1,345 admitted parents, and none was rejected. The 81 come from the package's tagging move, which the slice does not port.

**Deviations**
1. G1 is a graph test (no `:OnlyA` binding inside a component, no directed cycle of them), not Fourier–Motzkin. It is provably equivalent and LP-checked on every record.
2. A second free form must sit in the flip form's rapid-equilibrium segment using only steps outside `:OnlyA` groups. That segment exists in both conformations; no record's verdict changes.
3. **kcat:**
   - A corner saturates at most one ligand per regulatory site. The maximum is unchanged, and this removes the ambiguous corners of specs/34 and specs/38.
   - `kcat_mask` keeps rapid-equilibrium steps between free forms, because in AD13 the flip form sits on a dead-end branch.
4. **Printing:** each conformation is divided by its own normalizer, and the printer uses the raw flip weights (the reduced ones carry units). They print as `F_A`/`F_I` only when they differ and are not 1.

**What later stages must know**
- `kcat(m, params)`, `rate_equation_string(m)` and `coefficient_form(m)` must route an allosteric mechanism to the `AlloLaw` methods.
- The enumeration driver must filter every allosteric child with `admit`.
- Integration can merge two duplicated pieces: the kcat leading-face test (about 12 lines) and the coefficient labels (about 6). The printed locals (`Q_A`, `D_A`, `nA_0`) are not checked against metabolite names.

**Decisions for Denis**
- PFK (specs/35) has no kcat: ATP is both a substrate and an `:OnlyI` regulator, so the saturating rate depends on the ATP:F6P ratio. The package treats the two roles separately.
- AD3b and AD3d get kcat 0, against the package's 0.00004–37. Their substrate-inhibition copy closes a square with the catalytic forms, so the dead-end filter keeps it.
- AD8 errors correctly: it never saturates.

**My mistakes, all fixed**
- The printer first normalized the reduced cross weights, which gave wrong units in 950 laws, and it printed `F_I = 1` raised to a power.
- Several check-script bugs: E_total left out of a comparison, a sign error that left 87 directions uncertified, relative errors on zero limits, and a malformed hand-built case.
- I first miscounted the enumeration children as 25,701 instead of 24,503.

No Julia process of mine is left running.

Everything is in `/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/slice/`:
- `src/admission.jl`
- `src/allosteric.jl`
- `src/print_allo.jl`
- `INTERFACES.md`
- `allo/run_all.sh`
- `allo/logs/`