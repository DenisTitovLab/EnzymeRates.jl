# Literature and algorithms for a graph-centric derivation engine

**Summary**
- Hill's cycle-flux formula gives the numerator already in Haldane form, (ΠS − ΠP/Keq) times a polynomial with only positive coefficients. Nothing cancels, and Keq appears in one place.
- Haldane relations, Wegscheider relations and kinetic groups are all linear equations in log-parameters. One exact RREF with a fixed column order picks the dependent parameters deterministically. There are always μ−1 Wegscheider relations plus exactly one Haldane relation, where μ is the number of independent cycles.
- Splitting the graph at its articulation points (Volkenstein–Goldstein rule d, generalized) factors both numerator and denominator. A dead-end block always reduces to an equilibrium partition function, which removes one ROADMAP known issue.
- The :OnlyA/:OnlyI limits, the rapid-equilibrium (Cha) limit and the 0/0 case can all be handled by one routine that keeps the lowest-order face of the Newton polytope. The scaling is applied to independent parameters only.
- Normalizing means dividing by one reference monomial. The constant term exists only when the concentration-independent steps form a graph with exactly one terminal strongly connected component. Ping-pong has two, so for it the all-substrate monomial (Cleland/Dalziel) is the natural reference.

---

## 1. Hill's diagram method: cycle fluxes and a Haldane-ready numerator

**Key result**
- **Graph.** Nodes are enzyme forms. Edges are steps, labelled with pseudo-first-order constants (rate constant × concentration).
  - Hill's *directional diagrams* are the King–Altman patterns. D_r is the sum over arborescences directed to node r, and Σ = Σ_r D_r.
- **Net flux of a cycle κ:** J_κ = (Π_κ⁺ − Π_κ⁻)·Σ_κ / Σ.
  - Π_κ^± are the label products around κ in each direction.
  - Σ_κ is the sum of *flux diagrams*: κ plus a forest that directs every off-cycle node into κ. For implementation, Σ_κ = D_c(G/κ), the tree sum of the graph with κ contracted to a single node c.
- **Edge flux and rate.** The net flux of an edge is the signed sum of the cycle fluxes through it. The reaction rate is v = Σ_κ n_κ·J_κ, where n_κ is the number of net reactions per turn of κ (n_κ = 0 for internal loops).
- **Thermodynamics.** Under the cycle constraints, Π_κ⁻ = Π_κ⁺·(ΠP/(Keq·ΠS))^{n_κ}. Ligands that are bound and released within κ cancel. For n_κ = 1:
  - Π_κ⁺ − Π_κ⁻ = (Πk⁺)_κ · g_κ · (ΠS − ΠP/Keq), where g_κ is the shared bystander monomial.
  - Hence **v = (ΠS − ΠP/Keq) · Σ_{κ: n=1} (Πk⁺)_κ g_κ Σ_κ / Σ**.
  - When n_κ ≥ 2 the factor still appears, because xⁿ − yⁿ = (x − y)(…).
  - Beard & Qian prove the general one-way-flux identity J⁺/J⁻ = exp(−ΔG/RT).
- **Example.** Reversible Michaelis–Menten has two parallel E–ES edges (S binding and P release). Its only cycle is the 2-cycle on that parallel pair:
  - v = (k₁k₂S − k₋₁k₋₂P)/(k₋₁ + k₂ + k₁S + k₋₂P) = k₁k₂(S − P/Keq)/(…).
- **Relation to King–Altman.** The KA numerator k_ij·D_i − k_ji·D_j expands into unicyclic terms. Terms whose cycle is the back-and-forth i⇄j cancel in pairs; the rest regroup into Hill's cycle sum. Hill never generates the cancelling pairs. Chou (1989) reports that his rules avoid the 96 mutually cancelling terms of the random Bi Bi product-formation numerator.
- **Cost.**
  - Every simple cycle is a distinct nonzero element of the GF(2) cycle space, so there are at most 2^μ − 1 of them (μ = E − V + 1).
  - To enumerate them, XOR-combine the fundamental cycles of a spanning tree and keep the edge sets that are connected and 2-regular. For μ ≤ 12 that is at most 4095 sets.
  - Flux diagrams reuse the same arborescence routine on G/κ.
  - KDA (Python, NetworkX + SymPy) automates Hill's method. It was validated against a numerical solver and Kirchhoff tree counts on 875 random diagrams. It generates diagrams by "test and select" over edge combinations, which is simple but scales badly.

**Citations**
- Hill TL (1966) J Theor Biol 10:442–459, [10.1016/0022-5193(66)90137-8](https://doi.org/10.1016/0022-5193(66)90137-8).
- Hill TL, *Free Energy Transduction in Biology* (Academic Press 1977); *Free Energy Transduction and Biochemical Cycle Kinetics* (Springer 1989, [10.1007/978-1-4612-3558-3](https://doi.org/10.1007/978-1-4612-3558-3); Dover reprint 2005).
- Schnakenberg J (1976) Rev Mod Phys 48:571, [10.1103/RevModPhys.48.571](https://doi.org/10.1103/RevModPhys.48.571).
- Beard DA, Qian H (2007) PLoS ONE 2:e144, [10.1371/journal.pone.0000144](https://doi.org/10.1371/journal.pone.0000144).
- Chou KC (1989) J Biol Chem 264:12074–12079 (PMID 2745429).
- Awtrey NC, Beckstein O (2024) J Chem Theory Comput, [10.1021/acs.jctc.4c00688](https://doi.org/10.1021/acs.jctc.4c00688); preprint [10.1101/2024.05.27.596119](https://doi.org/10.1101/2024.05.27.596119).

**Relevance to the redesign**
- **G5, G6:** The numerator arrives as (ΠS − ΠP/Keq)·P⁺(k), with P⁺ a polynomial with positive coefficients. No substitute-then-simplify pass is needed and nothing cancels.
- **G11:** Near equilibrium, v is a small factor times a positive sum, so there is no catastrophic cancellation. There are also fewer terms than in KA's numerator, so the generated code is smaller.
- **G1, G4, G8:** Cycles and n_κ come directly from steps stored with their consumed and released metabolites. The same cycle routine labels each cycle Wegscheider (n = 0) or Haldane (n ≠ 0), which section 3 needs.
- **Pitfall:** For cycle enumeration, parallel steps must stay as separate edges; the MM catalytic cycle *is* such a pair. Volkenstein–Goldstein's "add parallel branches" is valid only inside tree sums.
- **G12:** The one-way flux J⁺ = v/(1 − Γ/Keq) and the per-cycle fluxes are natural quantities to report.

## 2. Directed spanning trees and King–Altman patterns

**Key result**
- **Matrix-tree theorem (Tutte).** D_r is the determinant of the Laplacian with row and column r removed.
  - Numerically it costs O(n³), which makes it an ideal test oracle.
  - Expanded symbolically it is a poor generator: the diagonal entries are sums, so the expansion produces cancelling terms.
- **Orienting undirected trees.** When every step is reversible, each undirected spanning tree orients in exactly one way toward any root.
  - So enumerate the undirected trees once (Hill's "partial diagrams") and orient them per root ("directional diagrams").
  - An orientation that needs a missing direction (an irreversible step) contributes 0.
  - Enumeration cost is O(V+E+VN) for undirected and O(V+E+EN) for directed graphs (Gabow & Myers), with later refinements by Kapoor & Ramesh and by Uno.
  - For 12 nodes or fewer, recursive include/exclude over the edge list with union-find is enough.
- **Cornish-Bowden (1977)**, the algebraic form of King–Altman:
  - Pick one outgoing constant per non-target form, reject products that contain a cycle, and get the other numerators by transposing along the path to the target.
  - Parallel steps must be summed first.
  - Repeated constants must be told apart by their position in the table, not by symbol.
  - For the general Bi Bi (8 forms, 12 steps) it produced 288 products per numerator, the same list as Lam & Priest.
- **Other classical methods.**
  - Lam & Priest generate only valid patterns and check them against an independent count.
  - Indge & Childs use Wang algebra (x·x = 0, x + x = 0), which is valid only when every edge has its own variable.
  - Fromm's schematic method produces duplicates that must be deleted.
  - Chou & Forsén and Chou give term-count formulas (useful as a check) and avoid cancelled terms.
- **Volkenstein & Goldstein's simplifications:**
  - (a) Parallel branches add.
  - (b) Symmetric branches merge with statistical factors (two identical sites give 2k₁S).
  - (c) Expansion on an auxiliary node m: D_r = Σ_i P_mr⁽ⁱ⁾·D_mr⁽ⁱ⁾, where each path from m to r is compressed to a point.
  - (d) Parts that share one node: every basic determinant is the product of the parts' determinants, each part directed to the shared node.
- **Software.** KAPattern (per-step steady-state, rapid-equilibrium or irreversible; MathML/MATLAB output), KDA, and GRASP, which uses King–Altman internally.

**Citations**
- King EL, Altman C (1956) J Phys Chem 60:1375, [10.1021/j150544a010](https://doi.org/10.1021/j150544a010).
- Volkenstein MV, Goldstein BN (1966) BBA 115:471–477; companion paper on allosteric models, BBA 115:478–485, [10.1016/0304-4165(66)90446-6](https://doi.org/10.1016/0304-4165(66)90446-6).
- Gabow HN, Myers EW (1978) SIAM J Comput 7:280, [10.1137/0207024](https://doi.org/10.1137/0207024).
- Lam CF, Priest DG (1972) Biophys J 12:248, [10.1016/S0006-3495(72)86084-3](https://doi.org/10.1016/S0006-3495(72)86084-3).
- Fromm HJ (1970) BBRC 40:692, [10.1016/0006-291X(70)90959-9](https://doi.org/10.1016/0006-291X(70)90959-9).
- Indge KJ, Childs RE (1976) Biochem J 155:567, [10.1042/bj1550567](https://doi.org/10.1042/bj1550567).
- Cornish-Bowden A (1977) Biochem J 165:55, [10.1042/bj1650055](https://doi.org/10.1042/bj1650055).
- Chou KC, Forsén S (1980) Biochem J 187:829, [10.1042/bj1870829](https://doi.org/10.1042/bj1870829).
- Qi F, Dash RK, Han Y, Beard DA (2009) BMC Bioinformatics 10:238, [10.1186/1471-2105-10-238](https://doi.org/10.1186/1471-2105-10-238).

**Relevance to the redesign**
- **G2, G4: the simplest correct core.** Enumerate undirected spanning trees per biconnected block (section 6), orient them per root, and give every edge its own symbol.
- **G3: a cheap first prototype.** Check the enumerator against matrix-tree determinants at random numeric rates and against Kirchhoff counts.
- **G1: groups after enumeration.** Every working method needs each edge to have its own identity; applying kinetic groups first breaks Wang algebra and symbol-based cycle rejection.
  - Substituting groups after enumeration only merges positive monomials (giving integer coefficients such as 2k₁k₂), so it can never cancel a term in D.
  - KDA works this way: generic edge rates first, then substitution of the assumptions.
- **G7:** Rule (b) is how n identical subunits collapse into statistical factors when MWC is assembled from parts.

## 3. Thermodynamically consistent parameterizations

**Key result**
- **Generalized Wegscheider condition.** Schuster & Schuster proved it necessary and sufficient for detailed balance, for any stoichiometry and with some reactants held at fixed concentration. In linear-algebra form: ln K ∈ im Nᵀ, equivalently zᵀ ln K = 0 for every z ∈ ker N.
- **Applied to a mechanism.** Let B be the form-incidence matrix and M the metabolite stoichiometry.
  - Every cycle z ∈ ker B satisfies Mz = n_z·r, where r is the overall reaction.
  - It imposes Σ_j z_j (ln k⁺_j − ln k⁻_j) = n_z ln Keq: a Wegscheider relation if n_z = 0, a Haldane relation otherwise.
  - Because z ↦ n_z is linear, all constraints are spanned by **μ−1 Wegscheider relations plus exactly one Haldane relation**.
  - The fundamental cycles of a spanning tree give a basis (Schnakenberg); Vlad & Ross give matrix tools.
- **Liebermeister & Klipp; Liebermeister, Uhlendorf & Klipp.**
  - Independent parameters: energy constants k^G_i = exp(G_i/RT), or independent equilibrium constants taken from a maximal set of linearly independent columns of N; velocity constants k^V = (k⁺_cat·k⁻_cat)^{1/2}; and the K values.
  - Dependent parameters: k^±_cat = k^V·Π_i (k^G_i·k^M_i)^{∓n_i/2}.
  - Every dependent log is linear in the independent logs. Lubitz et al.'s parameter balancing uses this dependency matrix.
- **Ederer & Gilles (thermodynamic-kinetic modelling, TKM).**
  - Each compound has a capacity, and its thermokinetic potential is proportional to its concentration. Each reaction has a resistance, and rate = force/resistance.
  - Detailed balance holds for all parameter values, because step equilibrium constants are ratios of capacities and every cycle telescopes.
- **Saa & Nielsen (GRASP).**
  - Parameters are per-step reversibilities plus a reference state, with path-consistency constraints for random-order mechanisms.
  - The same catalytic mechanism is sampled for the R and T states, and Cleland's rules map the results to macroscopic constants.
- **Alberty (1953)** extended Haldane relations to multi-substrate kinetic constants.

**Citations**
- Schuster S, Schuster R (1989) J Math Chem 3:25, [10.1007/BF01171883](https://doi.org/10.1007/BF01171883).
- Vlad MO, Ross J (2009) WIREs Syst Biol Med 1:348, [10.1002/wsbm.50](https://doi.org/10.1002/wsbm.50).
- Liebermeister W, Klipp E (2006) Theor Biol Med Model 3:41, [10.1186/1742-4682-3-41](https://doi.org/10.1186/1742-4682-3-41).
- Liebermeister W, Uhlendorf J, Klipp E (2010) Bioinformatics 26:1528, [10.1093/bioinformatics/btq141](https://doi.org/10.1093/bioinformatics/btq141).
- Lubitz T et al. (2010) J Phys Chem B 114:16298, [10.1021/jp108764b](https://doi.org/10.1021/jp108764b).
- Ederer M, Gilles ED (2007) Biophys J 92:1846, [10.1529/biophysj.106.094094](https://doi.org/10.1529/biophysj.106.094094).
- Saa P, Nielsen LK (2015) PLoS Comput Biol 11:e1004195, [10.1371/journal.pcbi.1004195](https://doi.org/10.1371/journal.pcbi.1004195).
- Alberty RA (1953) J Am Chem Soc 75:1928 (DOI not verified).

**Relevance to the redesign**
- **G5: deterministic dependent parameters with groups.**
  - All constraints have the form C·θ = n·ln Keq, with θ = G·φ, where G maps step constants to group parameters.
  - Do one exact (rational) RREF of [C·G | n]. Fix the column order from step identity alone, never from the grouping. The pivot columns are then the dependent parameters, and the choice does not depend on how the groups happen to be mixed.
  - A zero row is a constraint the grouping already satisfies.
  - A row that reads 0 = ln Keq means the grouping is thermodynamically infeasible; reject it.
  - A row with two ±1 entries is a Wegscheider rename. The same solve produces it, so the separate rename pass goes away (ROADMAP "Build the Wegscheider rename into the constraint columns").
- **G7:** Solving one system over both conformations' parameters gives the I-state its own Haldane relation. It also makes the conformational constants of liganded states Wegscheider-dependent (MWC's L_X = L·c_X). This removes the root cause of "the two conformations fold Wegscheider ties separately".
- **G12: `dependent_parameter_expression()` is the RREF rows.**
  - Exactly one row contains ln Keq.
  - Keq stays out of the Wegscheider rows if the Haldane row is pivoted on a constant that appears in no Wegscheider cycle, typically a catalytic step's reverse constant (when one exists).
  - The IndependentParameters and AllParameters modes correspond to Liebermeister's "system parameters" versus dependent quantities.
- **G9:** Symmetric per-step parameterizations (k^V plus an equilibrium constant, as in Liebermeister and TKM) make the written direction irrelevant: reversing a step maps K to 1/K. The internal representation can be direction-free even if the displayed names stay textbook names.
- **Caveat:** Energy and capacity parameterizations carry gauge redundancy and are far from textbook names. Use them, if at all, only internally.

## 4. Normal forms of rate equations

**Key result**
- **Cleland's coefficient form.** Cleland writes v = (num₁·AB − num₂·PQ)/(const + coef_A·A + coef_B·B + …) and defines kinetic constants as ratios of coefficients.
  - For ordered Bi Bi (I checked this against a King–Altman derivation): V₁ = num₁/coef_AB, K_a = coef_B/coef_AB, K_b = coef_A/coef_AB, K_ia = const/coef_A, and Keq = num₁/num₂.
  - Dividing by coef_AB gives V₁·AB/(K_ia·K_b + K_b·A + K_a·B + AB).
  - Ping-pong has the same form without the constant term.
  - Cleland lists several Haldane relations per mechanism. They are consistent because the kinetic constants are an over-complete set of coordinates for fewer rate constants.
- **Different authors, different constants.** Alberty (1953), Dalziel (1957), Bloomfield–Peller–Alberty (1962) and Cleland (1963) divide the same ordered Bi Bi polynomial by different rate-constant groups, so they name different constants.
- **Dalziel's form.** E/v = φ₀ + φ_A/A + φ_B/B + φ_AB/(AB) is the Laurent normalization by the all-substrate monomial. Ping-pong is simply φ_AB = 0.
- **Normalization is choosing a reference monomial m\* in supp(D).**
  - Divide N and D by c_{m\*}: a polynomial form, used for evaluation.
  - Or divide by c_{m\*}·x^{m\*}: a dimensionless Laurent form, used for display. Every term becomes a ratio of concentrations to constants, and the reference term becomes 1.
- **When a constant term exists** (my derivation from the matrix-tree theorem plus the terminal-SCC result):
  - The constant term is the tree sum of the subgraph of concentration-independent steps.
  - That sum is nonzero iff the subgraph has a spanning arborescence, i.e. exactly one terminal strongly connected component. The Laplacian kernel dimension equals the number of terminal SCCs (Gunawardena; Mirzaev & Gunawardena).
  - In ping-pong, E and F both leave only by binding, so there are two terminal components and no constant term.
  - In general the lowest degree of D is at least (number of terminal SCCs − 1).
- **Newton-polytope view.**
  - In any concentration limit, the dominant terms of D form the face of Newton(D) that the limit selects.
  - The constant term is the vertex at the origin; ΠS is the vertex reached at forward saturation.
  - An apparent K_M for X is a ratio of coefficients along a polytope edge in direction e_X. For example, K_a = coef_B/coef_AB is the edge from B to AB.
  - Wong & Hanes (and Chou & Forsén) predict which terms appear in D.

**Citations**
- Cleland WW (1963) BBA 67:104, [10.1016/0006-3002(63)91800-6](https://doi.org/10.1016/0006-3002(63)91800-6); 67:173, [10.1016/0006-3002(63)91815-8](https://doi.org/10.1016/0006-3002(63)91815-8); 67:188, [10.1016/0006-3002(63)91816-X](https://doi.org/10.1016/0006-3002(63)91816-X).
- Dalziel K (1957) Acta Chem Scand 11:1706, [10.3891/acta.chem.scand.11-1706](https://doi.org/10.3891/acta.chem.scand.11-1706).
- Wong JT-F, Hanes CS (1962) Can J Biochem Physiol 40:763.
- Nomenclature comparison across these authors: <https://astro.prao.ru/erm/comparison.html>.
- Gunawardena J (2012) PLoS ONE 7:e36321, [10.1371/journal.pone.0036321](https://doi.org/10.1371/journal.pone.0036321).
- Mirzaev I, Gunawardena J (2013) Bull Math Biol 75:2118, [10.1007/s11538-013-9884-8](https://doi.org/10.1007/s11538-013-9884-8).
- Segel IH, *Enzyme Kinetics* (Wiley 1975).

**Relevance to the redesign**
- **G6: one global rule.**
  - Option A: use the origin if it is in supp(D) ("1 + …"). Otherwise use ΠS (the Cleland/Dalziel choice; displayed as "1 + K_a/A + K_b/B + …"). As a last resort, use the lowest monomial in a fixed graded order.
  - Option B: normalize every mechanism by ΠS, so the numerator constant is always V₁.
  - Both are principled. What matters is using one rule everywhere, because mixing rules changes what the fitted constants mean.
  - Evaluate the polynomial form; the Laurent form is singular at zero concentration.
- **Bridge to textbook equations.**
  - Compute the named constants as coefficient ratios.
  - Check a textbook equation by exact polynomial identity.
  - List its dependent named constants (its extra Haldane relations) explicitly.
- **G14:** Normalizing a map from monomial to coefficient, rather than printed expressions, removes artifacts like "1 ^ 2".

## 5. Limits in rational rate laws

**Key result**
- **Scaling and leading terms.** Write each parameter as p_i = ε^{a_i}·p̂_i. K → ∞ means a_K < 0; removing a steady-state step by sending both k → 0 means a = 1 on both constants.
  - A monomial with parameter-exponent vector e has order ⟨a, e⟩.
  - As ε → 0: N/D → N_lead/D_lead if ord N = ord D; N/D → 0 if ord N > ord D; otherwise it diverges.
  - N_lead and D_lead are the faces of the Newton polytopes that minimize ⟨a, ·⟩ (min-plus, i.e. tropical).
  - D has positive coefficients, and Hill's N is (ΠS − ΠP/Keq) times a positive polynomial, so nothing can cancel on the leading face.
  - No Puiseux series are needed; those arise only when solving for unknown roots.
- **Tree sums.** The order of D_r is the cost of a minimum-cost arborescence (Edmonds). This is the W-graph asymptotics of Freidlin–Wentzell and Gorban–Radulescu's "dominant systems".
- **Rapid equilibrium is such a limit.** Scale both constants of the fast steps by λ → ∞ at fixed K, and the steady-state law becomes Cha's law.
  - Topham & Brocklehurst showed Cha's method is exact in this sense.
  - They traced Segel & Martin's objection to parallel pathways through the thermodynamic box that were left out.
- **MWC exclusive binding** (Rubin & Changeux) is the limit c = K_R/K_T → 0.
  - Popova & Sel'kov generalized MWC to a reversible reaction; each conformation's constants must then satisfy the reaction's Haldane relation.
- **Thermodynamic consistency of limits** (Gorban & Yablonsky):
  - Irreversible limits of detailed-balanced systems obey "extended detailed balance": detailed balance on the reversible part, plus a structural condition that no irreversible step lies on an oriented cycle.
  - The limiting mechanism depends on the relative rates at which constants vanish (Gorban, Mirkes & Yablonsky 2013).

**Citations**
- Samal SS et al. (2015) Bull Math Biol 77:2180, [10.1007/s11538-015-0118-0](https://doi.org/10.1007/s11538-015-0118-0).
- Gorban AN, Radulescu O, Zinovyev AY (2010) Chem Eng Sci 65:2310, [10.1016/j.ces.2009.09.005](https://doi.org/10.1016/j.ces.2009.09.005).
- Edmonds J (1967) J Res NBS 71B:233.
- Freidlin MI, Wentzell AD, *Random Perturbations of Dynamical Systems* (Springer), chapter 6.
- Cha S (1968) J Biol Chem 243:820 (PMID 5638598).
- Topham CM, Brocklehurst K (1992) Biochem J 282:261, [10.1042/bj2820261](https://doi.org/10.1042/bj2820261).
- Segel IH, Martin RL (1988) J Theor Biol 135:445, [10.1016/S0022-5193(88)80269-8](https://doi.org/10.1016/S0022-5193(88)80269-8).
- Rubin MM, Changeux JP (1966) J Mol Biol 21:265, [10.1016/0022-2836(66)90097-0](https://doi.org/10.1016/0022-2836(66)90097-0).
- Popova SV, Sel'kov EE (1975) FEBS Lett 53:269, [10.1016/0014-5793(75)80034-2](https://doi.org/10.1016/0014-5793(75)80034-2).
- Gorban AN, Yablonsky GS (2011) Chem Eng Sci 66:5388, [10.1016/j.ces.2011.07.054](https://doi.org/10.1016/j.ces.2011.07.054) (arXiv:1101.5280).
- Gorban AN, Mirkes EM, Yablonsky GS (2013) Physica A 392:1318.

**Relevance to the redesign**
- **G7:** Implement :OnlyA/:OnlyI as a scaling vector on independent parameters only.
  - Each dependent parameter's order then follows linearly from its RREF row: a_dep = Σ c_i·a_i. These are the ROADMAP's "partner limits".
  - The limit of a family of consistent laws is consistent by construction.
- **ROADMAP 0/0 case:** Removing a bridge step zeroes every tree, but N and D both become O(ε), and the leading terms give a finite limit.
- **Check tag combinations.** If a dependent constant on a surviving edge diverges, or a step becomes one-way on a cycle, the tags describe a different mechanism (Gorban–Yablonsky). Reject or report such combinations.
- **One routine, several uses.** The same leading-face routine handles the allosteric tags and the rapid-equilibrium limit. "Steady-state law → Cha law" is a good G3 oracle. :EqualAI and :NonequalAI need only substitution or relabelling.

## 6. Other pieces for a compact, correct engine

- **Factoring at articulation points.** Rule (d) generalizes through Hopcroft–Tarjan biconnected components (linear time).
  - D_r = Π_B D^B_{a(B,r)}, where a(B,r) is the articulation point of block B nearest to r. Σ is a dynamic program over the block-cut tree.
  - Cycles live inside blocks, so the numerator factors too.
  - A dead-end block (no reaction cycle) carries no net flux. With consistent constants it is therefore at equilibrium relative to its attachment point, by Gunawardena's cycle condition.
  - Its factor is then a partition function Z_B, and each dead-end block adds D^core_a·(Z_B − 1) to the denominator.
  - This gives the classic (1 + I/K_I) factors structurally and shrinks generated code (G11).
  - It also fixes the ROADMAP issue "a steady-state step that only reaches dead-end forms fits both rate constants": only the ratio k_off/k_on can enter the law.
- **Cha contraction combined with Hill.** This is my consolidation of Cha's f_i after clearing denominators; I checked it on RE ping-pong and RE MM.
  - Inside each rapid-equilibrium segment, N_i is the path product of K's from the segment's reference form, and Z_s = Σ_{i∈s} N_i.
  - Reduced edge labels are k_e·N_source(e).
  - Then v = Σ_κ n_κ·(Π̃⁺ − Π̃⁻)·Σ̃_κ / Σ_s Z_s·D̃_s. Slow steps inside one segment are self-loop cycles.
  - Once dependent constants are substituted, every path gives the same Laurent monomial for N_i. Cha's "any one path" rule becomes a syntactic identity, and the output is canonical (G5).
- **Node elimination.** Rule (c) is the combinatorial side of Schur-complement (Kron) reduction: k′_ik = k_ik + k_ij·k_jk/Σ_l k_jl. GTH does this without subtractions.
  - It gives O(n³) nested evaluation with no tree expansion.
  - At n = 12 that is about 600 operations, which probably exceeds 120 ns. Treat it as a fallback or an oracle.
- **Oracles for the prototypes (G3).**
  - Matrix-tree determinant versus enumerated tree sums.
  - Numerical null-space steady state versus the derived v.
  - v = 0 at Γ = Keq.
  - Convergence to Cha's law in the λ-limit.
  - Kirchhoff tree counts on random mechanisms, as KDA does.
- **Canonical rate laws (G10).**
  - After substituting dependent parameters, cancelling the structural common factors (from blocks and Cha segments) and normalizing by the fixed reference monomial, equal rate laws give identical coefficient maps. These can be hashed and used to deduplicate in identify.
  - Enzyme forms already carry intrinsic labels (bound metabolites plus residual), so a sorted, orientation-free edge list is already a canonical mechanism key. Canonical graph labelling (McKay & Piperno, nauty) is needed only if those labels were dropped.
- **Coefficient form as the code IR (G11, G2).**
  - Cleland's coefficient form splits evaluation into coefficients(params) and a sparse polynomial in concentrations.
  - With block and Cha factorization, generated code then scales with the number of coefficients plus monomials, not the number of trees.
  - The whole engine reduces to about six generic routines: tree enumeration, cycle-space enumeration, block decomposition, rational RREF, leading-face extraction, and polynomial maps.

**Citations**
- Hopcroft J, Tarjan R (1973) Commun ACM 16:372, [10.1145/362248.362272](https://doi.org/10.1145/362248.362272).
- Grassmann WK, Taksar MI, Heyman DP (1985) Oper Res 33:1107, [10.1287/opre.33.5.1107](https://doi.org/10.1287/opre.33.5.1107).
- Dörfler F, Bullo F (2013) IEEE TCAS-I 60:150 (Kron reduction).
- Meyer CD (1989) SIAM Rev 31:240 (stochastic complementation).
- McKay BD, Piperno A (2014) J Symb Comput 60:94, [10.1016/j.jsc.2013.09.003](https://doi.org/10.1016/j.jsc.2013.09.003).

**Sources:** Bibliographic details were checked against PubMed records (DOIs linked above), publisher, arXiv and KDA documentation pages, and the three PDFs Denis supplied (Volkenstein–Goldstein rules a–d, Cornish-Bowden's notes on repeated constants and parallel steps, and Cha's "any one pathway" rule). Items marked "not verified" are the Alberty 1953 DOI and the Gorban–Mirkes–Yablonsky 2013 DOI. The constant-term criterion and the Cha+Hill consolidation are my own derivations from the cited results, not published statements.