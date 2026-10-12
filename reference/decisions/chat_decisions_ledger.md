# Decisions made in chat (2026-10-09 to 2026-10-12), not or not fully on the decision cards

Source of truth for reviewers alongside decisions/*.json (Denis's card choices and notes) and
answers/*.json (replies on the cards). Later entries supersede earlier ones.

1. Brief confirmed; ONE umbrella spec, then one writing-plans plan per phase PR, executed with
   subagent-driven development, test-first.
2. Build approach: test-first in the package; the vertical slice (scratchpad/slice) is a reference
   for algorithms and interfaces only; oracles the tests need move into test/.
3. Phases after PR #85 (Horner): fitting scale; joint solve on today's types; graph Mechanism +
   DSL + show/CSV; admission + guard + flip at E; engine/printing/key/codegen. The tests-only
   phase (E11/critic) was DROPPED: each phase updates the pins of the behavior it changes in the
   same PR, and adds an oracle value check before changing any test that checks only a string.
   PRs #85, #86, #87 and a fix for the ter-bi expansion abort merge before phase 1.
4. DSL (2026-10-09, revised 2026-10-11): TWO macros, @enzyme_mechanism and @allosteric_mechanism,
   sharing labels substrates:/products:/modifiers:/steps: and step syntax; only the allosteric one
   takes catalytic_multiplicity:, allosteric_regulators:, regulatory_site(n): lines (one per site,
   multiplicity + per-ligand tags; every declared regulator in exactly one site; no default sites)
   and a tag on every step or group. modifiers: replaces regulators:/catalytic_inhibitors:;
   ::Mod marks only a reactant bound as a non-reacting ligand (dual role); modifier-only names get
   no suffix (K_EI_to_E_I). Removed: catalytic_steps:, catalytic_inhibitors:, dead_end_inhibitors:
   (-> competitive_inhibitors: in @enzyme_reaction), allowed_catalytic_multiplicities:
   (-> oligomeric_state:), multi-argument atom brackets B[N, P], the parenthesized single-step tag
   form, default regulatory sites. Kept: steps either direction, groups as tuples, bare tags,
   residuals, conformation labels, A(1, 2). show(m) prints the DSL; a parser reads it back without
   eval; identify CSV stores DSL text + short stable mechanism_id + parent_id.
5. Types: COMPOSITION. Mechanism = reaction graph; AllostericMechanism wraps a catalytic Mechanism
   + tags + multiplicity + sites; future KNF wraps the same graph. Compiled types
   EnzymeMechanism{Sig}, AllostericEnzymeMechanism{Sig}.
6. E10: conformations flip at free E (in the rewrite, with B1's island rule); EVERY mechanism
   (plain and allosteric) must contain a free form named E (nothing bound, no residual), stored as
   form 1; the E_c test fixture is renamed.
7. Section 3a approved: joint solve (columns, rows, exact RREF, rational exponents, no rename
   pass), fixed priority (kind, occupied, SS before RE, against = reaction direction, ties by
   content names), names after the kept member, parameter API (parameters / fitted_params /
   all_parameters, alphabetical), Full/Reduced removed. "Inactive :OnlyA members name themselves"
   is an internal detail (not user-facing).
8. Section 3b approved with: :OnlyI allowed ONLY on groups whose steps bind a modifier
   (modifiers: names and reactants' ::Mod copies), never on a step that takes up or releases a
   substrate or product as a reactant; admission mirrored in the active conformation, which must
   still turn over. Per-step classes reworded: each member's constant vanishes at its own rate;
   cycle rows still bind them; finite constants outside the law are exactly consistent via the
   solve.
9. E2 refined ("well_defined_only"): a mechanism is valid iff its RE/SS law equals the unique
   fast-RE limit of mass action wherever data sit (faces: every product zero, every substrate
   zero, each single reactant zero). Reject: malformed (incl. no catalytic cycle through an SS
   step, v = 0), RE turnover, a law that is not the unique fast-RE limit on a face (bottomless,
   split, trap; the trap error names the step to make SS), :OnlyA/:OnlyI with no limit. Fit:
   substrate inhibition, abortive complexes, flux-free SS steps, product activation, hand-written
   laws that never saturate. Enumeration keeps not emitting never-saturating mechanisms
   (_degenerate). Specs "Numerator: redundant SS-bind" and "LDH i-state NonequalAI 6-group" are
   invalid -> rejection tests + valid variants. Face check made exact first.
10. E9: no-anchor scale = peak FORWARD rate over ALL concentrations, products included (Denis
    chose this over products-at-zero); deterministic multi-start search in log concentrations,
    24 decades; a law still rising at the box edge stays unscaled and flagged (never a failure).
    kcat applied AFTER the fit: FittingProblem(...; Keq, absolute = false);
    rescale_parameter_values(m, params, kcats), kcats = [(A = 10, B = 5) => 42.0, ...], unlisted
    metabolites zero, geometric mean of measured/predicted, spread reported; identify has no kcat
    keyword. scale_k_to_kcat and _kcat_forward removed.
11. Normalization: "complete unit rule": constant term; else virtual constant term from the units
    equation (sum lambda = 1, sum lambda * a_m = 0; Cleland's N is one solution; fractional
    exponents allowed, printed as roots; canonical basis by lowest-degree terms and fixed order);
    else ratio form (divide by the lowest-degree term). The plain-coefficient fallback goes.
12. Numerator: today's sum over SS steps of net flux times copies of the first substrate moved
    (RE-weight exponents included), then exact division by the Haldane factor H.
13. Enumeration: islands rejected by admission; duplicate tag sets and one-subunit L phantoms are
    KEPT in the frontier but not fitted separately (duplicates share one fit via the rename-proof
    key; L phantoms inherit the parent's fit and are flagged); L phantom = inactive conformation
    binds nothing and does not catalyze AND the identifiable rank shows L is a phantom; the beam
    counts a key once but expands all its members; the "dedup tag sets in today's package now"
    PR is dropped; the tagging move emits no :OnlyI.
14. identify's maxtime docstring is wrong (maxtime applies per restart): fix in the fitting phase.

Spec author's choices NOT explicitly confirmed by Denis (reviewers: flag if risky):
- parameters and all_parameters exported; fitted_params stays internal.
- constructors enforce representation invariants only; guard and admission run at the DSL
  boundary, on enumerated children and at the derivation entry.
- planned source layout; peak search box 1e-12..1e12 in the data's units; the size-scaled speed
  bound's constant calibrated on CI in the last phase.
