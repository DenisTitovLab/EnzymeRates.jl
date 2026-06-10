# Comprehensive Documentation Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a Documenter.jl documentation site for EnzymeRates.jl that faithfully describes how the package derives, fits, and identifies enzyme rate equations, then trim the README and CLAUDE.md to point into it.

**Architecture:** A Documenter.jl site lives in a `docs/` subproject (its own `Project.toml` that `dev`s EnzymeRates and adds the doc/optimizer deps, kept out of the main `Project.toml` so Aqua stays green), built by `docs/make.jl` with `doctest = true` and `checkdocs = :exports` and a DocumenterCitations bibliography plugin. A `docs` job added to the existing `.github/workflows/CI.yml` runs `julia-actions/julia-docdeploy@v1` (PRs verify, merge-to-`main` deploys gh-pages). Content is organized into three pillars — deriving, fitting, identifying — each opening with a runnable tutorial; the README and CLAUDE.md are trimmed last so they reference only pages that already exist.

**Tech Stack:** Documenter.jl, DocumenterCitations.jl, julia-actions (setup-julia, cache, julia-buildpkg, julia-docdeploy), OptimizationPyCMA (blessed optimizer, with OptimizationBBO as the CI-runnable alternative).

---

## File Structure

**Created**

- `docs/Project.toml` — docs subproject manifest: `dev`s EnzymeRates and adds Documenter, DocumenterCitations, Optimization, OptimizationPyCMA (and OptimizationBBO for the CI-runnable example optimizer); kept separate from the main `Project.toml` so the Aqua stale-deps check stays green.
- `docs/make.jl` — the build entry point: constructs the `CitationBibliography` plugin, calls `makedocs` with `doctest = true`, `checkdocs = :exports`, the full `pages` site tree, and HTML format (canonical URL, `edit_link = "main"`), then `deploydocs` with `devbranch = "main"`.
- `docs/src/refs.bib` — verified BibTeX bibliography (Cha 1968, King & Altman 1956, Monod–Wyman–Changeux 1965, Segel 1975, Haldane 1930, Wegscheider 1901, Hastie et al. 2009, Hansen & Ostermeier 2001 — eight entries for seven concept rows).
- `docs/src/index.md` — Home: what the package is, install, links into the docs.
- `docs/src/getting_started.md` — end-to-end newcomer arc (define → derive → fit → identify), reusing the fast identify example.
- `docs/src/deriving/textbooks.md` — Rate equations from textbooks: the `rate_equation_string` tutorial (pillar opener, the one byte-stable derivation `jldoctest`).
- `docs/src/deriving/re_vs_ss.md` — Rapid equilibrium vs steady state: the RE/SS distinction and its parameter-count consequences.
- `docs/src/deriving/cha_king_altman.md` — The Cha / King–Altman algorithm: division-free, per-segment, finite-at-zero derivation.
- `docs/src/deriving/thermodynamic_constraints.md` — Thermodynamic constraints: Haldane/Wegscheider parameter reduction via the exact-integer null space.
- `docs/src/deriving/ping_pong.md` — Ping-pong mechanisms: the residual-on-`:E` enumerator representation.
- `docs/src/deriving/dead_end.md` — Dead-end inhibitor binding: the all-RE competitive-inhibitor denominator term.
- `docs/src/deriving/mwc_allostery.md` — MWC allostery: A/I two-conformation partition-function rate, four-tag taxonomy.
- `docs/src/fitting/tutorial.md` — Fitting tutorial & data format.
- `docs/src/fitting/normalized_vs_absolute.md` — Normalized vs absolute rate: the `scale_k_to_kcat` knob and kcat rescaling.
- `docs/src/fitting/loss_and_optimizers.md` — Loss & optimizers: the log-ratio loss and the bring-your-own-optimizer story.
- `docs/src/identify/tutorial.md` — Identify tutorial: the fast width-1-beam recovery example.
- `docs/src/identify/model_selection.md` — Model selection: LOOCV, paired 1-SE rule, one-sided permutation test.
- `docs/src/identify/enumeration_engine.md` — The enumeration engine: `init_mechanisms` + six expansion moves + `_dedup_flat!`.
- `docs/src/developer.md` — Developer / Architecture: curated internals (chokepoint, Canonical Step Form, concrete-vs-singleton split, Source Layout, compile-budget), plus `@docs` for `Mechanism`, `Step`, `init_mechanisms`, `compile_mechanism`.
- `docs/src/api.md` — API Reference: `@index` over grouped `@autodocs` (`Private = false`) for types, macros, constants, functions.
- `docs/src/roadmap.md` — Roadmap: KNF allostery, parameter identifiability, iso mechanisms, plotting, outlier-dataset identification.
- `docs/src/references.md` — References: the `@bibliography` block.

**Modified**

- `.github/workflows/CI.yml` — append a `docs` job (sibling of `test`) running setup-julia/cache/buildpkg/julia-docdeploy with `contents: write`, passing `DOCUMENTER_KEY` (preferred SSH deploy key, installed) and `GITHUB_TOKEN` (fallback), plus a doctest step.
- `.gitignore` — append `docs/build/` (and rely on the existing `Manifest.toml` rule for `docs/Manifest.toml`).
- `src/types.jl` — add/repair docstrings: write the `EnzymeReaction` type docstring, reattach the `metabolites` docstring (move the `@generated`-rationale comment into the body), add `Step` and `Mechanism` docstrings (and verify `AllostericEnzymeMechanism`).
- `src/rate_eq_derivation.jl` — expand the `rate_equation_string` docstring (mode, multi-line format, Full-mode restriction, `jldoctest`).
- `README.md` — full rewrite to a thin landing page: badges, one-paragraph description, install, ~10-line derive-and-evaluate quickstart, links into the docs site.
- `.claude/CLAUDE.md` — apply the migration table: cut maintainer-internal and user-concept blocks (content moved to docs), leave one-line guards for Canonical Step Form and the parameter-naming chokepoint, keep the `rate_equation` runtime-perf guardrail and all agent-behavioral rules.
- `test/runtests.jl` — remove the `include("test_readme_runs.jl")` line.
- `test/test_readme_runs.jl` — deleted (its coverage moves into Documenter doctests).

---

I have everything I need. The branch is `docs-comprehensive`, `docs/src` is empty, exports confirmed, repo is `DenisTitovLab/EnzymeRates.jl`, julia 1.12.5. Now I'll write the Phase 0 plan section.

## Phase 0 — Infrastructure Scaffold

This phase stands up the `docs/` subproject and proves the entire Documenter + doctest + citation + CI/deploy pipeline end-to-end *before* any real content exists. The deliverable is a documentation site that builds green locally and in CI off placeholder stub pages, with citations resolving, doctests running (none yet, but the machinery is live), and a `docs` CI job wired into the existing `.github/workflows/CI.yml`. The only out-of-band step is Denis flipping the repo Pages source to the `gh-pages` branch after the first CI deploy creates it (`DOCUMENTER_KEY` is already installed and authenticates the push, with `GITHUB_TOKEN` as fallback). We are already on the `docs-comprehensive` branch with an empty `docs/src/`, so no branch creation is needed. Exit criteria: `julia --project=docs docs/make.jl` builds clean locally, `docs/build/` is gitignored, the CI job is present, and the Pages-enable checklist is recorded.

### Task 0.1: docs/ subproject and dependency instantiation

**Files:**
- Create: `docs/Project.toml`

- [ ] **Step 1: Write `docs/Project.toml`**
  Create `/home/denis.linux/.julia/dev/EnzymeRates/docs/Project.toml` with exactly this content. The UUIDs are the registered General-registry UUIDs for these packages; `EnzymeRates` is `dev`'d in (path-based) by the next step, so it is intentionally absent from `[deps]` here until `Pkg.develop` writes it.

  ```toml
  [deps]
  Documenter = "e30172f5-a6a5-5a46-863b-614d45cd2de4"
  DocumenterCitations = "daee34ce-89f3-4625-b898-19384cb65244"
  EnzymeRates = "e7a1c2d3-4b5f-6a7e-8c9d-0e1f2a3b4c5d"
  Optimization = "7f7a1694-90dd-40f0-9382-eb1efda571ba"
  OptimizationPyCMA = "fb0822aa-1fe5-41d8-99a6-e7bf6c238d3b"

  [compat]
  Documenter = "1"
  DocumenterCitations = "1"
  Optimization = "4"
  OptimizationPyCMA = "1"
  julia = "1.9"
  ```

- [ ] **Step 2: Run dependency develop + instantiate**
  Run from the repo root (`/home/denis.linux/.julia/dev/EnzymeRates`):
  ```bash
  julia --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'
  ```
  Expected: `Pkg.develop` resolves `EnzymeRates v0.1.0` from the local path; `Pkg.instantiate` downloads and precompiles `Documenter`, `DocumenterCitations`, `Optimization`, `OptimizationPyCMA` and their transitive deps. Final output ends with `Precompiling project...` completing with no `ERROR`. This also writes `docs/Manifest.toml` (gitignored) and confirms `EnzymeRates` now appears under `[deps]` in `docs/Project.toml` as a path dev-dep entry.

- [ ] **Step 3: Verify the packages load**
  ```bash
  julia --project=docs -e 'using Documenter, DocumenterCitations, EnzymeRates, Optimization, OptimizationPyCMA; println("docs env OK")'
  ```
  Expected: prints `docs env OK` with no errors. (OptimizationPyCMA pulls a Python CMA-ES via CondaPkg on first load; this may take a minute but must end without `ERROR`.)

- [ ] **Step 4: Commit**
  ```bash
  git add docs/Project.toml
  git commit -m "Add docs/ subproject with Documenter, DocumenterCitations, optimizer deps

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
  ```

### Task 0.2: refs.bib with the seven verified references

**Files:**
- Create: `docs/src/refs.bib`

- [ ] **Step 1: Verify each BibTeX entry against an authoritative source**
  Do NOT trust the skeleton DOIs below — verify every entry before writing the file. For each of the seven references, confirm authors, year, journal/venue, volume, pages, and DOI against an authoritative source (publisher page, Crossref, or PubMed). Concretely:
  - Cha 1968, Cha 1968 *J Biol Chem* "A simple method for derivation of rate equations…" — confirm volume/issue/pages and DOI via the JBC site or Crossref (`https://api.crossref.org/works?query.bibliographic=Cha+1968+rate+equations+steady+state&rows=3`).
  - King & Altman 1956 *J Phys Chem* "A schematic method of deriving the rate laws for enzyme-catalyzed reactions" — confirm volume 60, pages, and DOI via Crossref.
  - Monod, Wyman, Changeux 1965 *J Mol Biol* "On the nature of allosteric transitions: A plausible model" — confirm vol 12, pp 88–118, DOI `10.1016/S0022-2836(65)80285-6` via the journal page.
  - Segel 1975 *Enzyme Kinetics* (Wiley) — book; confirm ISBN and publisher (no DOI required for a book; use `isbn` field).
  - Haldane 1930 *Enzymes* (Longmans, Green) — book; confirm publisher/year (no DOI; `isbn`/`note` as available). Wegscheider 1901 *Z Phys Chem* — confirm volume 39, pages, DOI if Crossref has it; otherwise record as a `@article` without DOI.
  - Hastie, Tibshirani, Friedman *The Elements of Statistical Learning* 2nd ed. 2009 (Springer) — confirm DOI `10.1007/978-0-387-84858-7`.
  - Hansen CMA-ES — cite "The CMA Evolution Strategy: A Tutorial" (Hansen 2016, arXiv:1604.00772) OR Hansen & Ostermeier 2001 *Evolutionary Computation* "Completely derandomized self-adaptation in evolution strategies" (DOI `10.1162/106365601750190398`); pick the 2001 journal article and confirm its DOI via Crossref.

  Use the available `WebSearch`/`WebFetch` tools (load via ToolSearch) or the Crossref REST API through `Bash` curl to confirm each. Record the confirmed DOI/ISBN inline. If any DOI cannot be confirmed, leave the `doi` field OUT entirely rather than guessing, and add a `% VERIFY: <what is missing>` comment on that entry.

- [ ] **Step 2: Write `docs/src/refs.bib`**
  Create `/home/denis.linux/.julia/dev/EnzymeRates/docs/src/refs.bib` with the seven entries below. These are skeletons keyed to match the `[Key](@cite)` keys used later in the docs (`Cha1968`, `KingAltman1956`, `Monod1965`, `Segel1975`, `Haldane1930`, `Wegscheider1901`, `Hastie2009`, `Hansen2001`). Replace each `doi`/`isbn`/page field with the value confirmed in Step 1. Each entry carries a `% VERIFY` comment that the executor MUST remove only after confirming that entry's metadata.

  ```bibtex
  % EnzymeRates.jl documentation references.
  % Every entry below must be verified in Task 0.2 Step 1 before the
  % "% VERIFY" line is removed. Do not invent DOIs.

  % VERIFY: confirm volume/issue/pages and DOI (JBC / Crossref).
  @article{Cha1968,
    author  = {Cha, Sungman},
    title   = {A Simple Method for Derivation of Rate Equations for Enzyme-catalyzed Reactions under the Rapid Equilibrium Assumption or Combined Assumptions of Equilibrium and Steady State},
    journal = {Journal of Biological Chemistry},
    year    = {1968},
    volume  = {243},
    number  = {4},
    pages   = {820--825},
    doi     = {10.1016/S0021-9258(19)81739-8},
  }

  % VERIFY: confirm volume 60, pages, DOI (Crossref).
  @article{KingAltman1956,
    author  = {King, Edward L. and Altman, Carl},
    title   = {A Schematic Method of Deriving the Rate Laws for Enzyme-Catalyzed Reactions},
    journal = {The Journal of Physical Chemistry},
    year    = {1956},
    volume  = {60},
    number  = {10},
    pages   = {1375--1378},
    doi     = {10.1021/j150544a010},
  }

  % VERIFY: confirm vol 12, pp 88-118, DOI (J Mol Biol).
  @article{Monod1965,
    author  = {Monod, Jacques and Wyman, Jeffries and Changeux, Jean-Pierre},
    title   = {On the Nature of Allosteric Transitions: A Plausible Model},
    journal = {Journal of Molecular Biology},
    year    = {1965},
    volume  = {12},
    number  = {1},
    pages   = {88--118},
    doi     = {10.1016/S0022-2836(65)80285-6},
  }

  % VERIFY: confirm publisher and ISBN (book; no DOI).
  @book{Segel1975,
    author    = {Segel, Irwin H.},
    title     = {Enzyme Kinetics: Behavior and Analysis of Rapid Equilibrium and Steady-State Enzyme Systems},
    publisher = {Wiley},
    address   = {New York},
    year      = {1975},
    isbn      = {978-0471303091},
  }

  % VERIFY: confirm publisher/year (book; no DOI).
  @book{Haldane1930,
    author    = {Haldane, John Burdon Sanderson},
    title     = {Enzymes},
    publisher = {Longmans, Green and Co.},
    address   = {London},
    year      = {1930},
  }

  % VERIFY: confirm volume 39, pages, DOI if available (Z. Phys. Chem.).
  @article{Wegscheider1901,
    author  = {Wegscheider, Rudolf},
    title   = {{\"U}ber simultane Gleichgewichte und die Beziehungen zwischen Thermodynamik und Reactionskinetik homogener Systeme},
    journal = {Zeitschrift f{\"u}r Physikalische Chemie},
    year    = {1901},
    volume  = {39},
    pages   = {257--303},
  }

  % VERIFY: confirm DOI 10.1007/978-0-387-84858-7 (Springer, 2nd ed.).
  @book{Hastie2009,
    author    = {Hastie, Trevor and Tibshirani, Robert and Friedman, Jerome},
    title     = {The Elements of Statistical Learning: Data Mining, Inference, and Prediction},
    publisher = {Springer},
    edition   = {2},
    year      = {2009},
    series    = {Springer Series in Statistics},
    doi       = {10.1007/978-0-387-84858-7},
  }

  % VERIFY: confirm DOI 10.1162/106365601750190398 (Evol. Comput.).
  @article{Hansen2001,
    author  = {Hansen, Nikolaus and Ostermeier, Andreas},
    title   = {Completely Derandomized Self-Adaptation in Evolution Strategies},
    journal = {Evolutionary Computation},
    year    = {2001},
    volume  = {9},
    number  = {2},
    pages   = {159--195},
    doi     = {10.1162/106365601750190398},
  }
  ```

- [ ] **Step 3: Validate the BibTeX parses under DocumenterCitations**
  ```bash
  julia --project=docs -e 'using DocumenterCitations; b = CitationBibliography("docs/src/refs.bib"; style=:authoryear); println("bib OK: ", length(b.entries), " entries")'
  ```
  Expected: prints `bib OK: 8 entries` (the eight keys above; note Haldane and Wegscheider are two separate entries even though the spec table groups them on one row). No parse error. If the count is not 8 or a parse error is thrown, fix the offending entry before continuing.

- [ ] **Step 4: Commit**
  ```bash
  git add docs/src/refs.bib
  git commit -m "Add verified bibliography for docs (8 references)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
  ```

### Task 0.3: docs/make.jl

**Files:**
- Create: `docs/make.jl`

- [ ] **Step 1: Write `docs/make.jl`**
  Create `/home/denis.linux/.julia/dev/EnzymeRates/docs/make.jl` with exactly this content. The `pages=` tree mirrors spec section 2; every target is a stub created in Task 0.4. The `bib` plugin is constructed once and passed in `plugins=`. `makedocs` runs with `doctest = true` and `checkdocs = :exports` per spec section 1.

  ```julia
  using Documenter
  using DocumenterCitations
  using EnzymeRates

  bib = CitationBibliography(joinpath(@__DIR__, "src", "refs.bib"); style = :authoryear)

  makedocs(;
      sitename = "EnzymeRates.jl",
      authors = "Denis Titov and contributors",
      modules = [EnzymeRates],
      doctest = true,
      checkdocs = :exports,
      plugins = [bib],
      format = Documenter.HTML(;
          canonical = "https://DenisTitovLab.github.io/EnzymeRates.jl",
          edit_link = "main",
      ),
      pages = [
          "Home" => "index.md",
          "Getting Started" => "getting_started.md",
          "Deriving rate equations" => [
              "Rate equations from textbooks" => "deriving/textbooks.md",
              "Rapid equilibrium vs steady state" => "deriving/re_vs_ss.md",
              "The Cha / King–Altman algorithm" => "deriving/cha_king_altman.md",
              "Thermodynamic constraints" => "deriving/thermodynamic_constraints.md",
              "Ping-pong mechanisms" => "deriving/ping_pong.md",
              "Dead-end inhibitor binding" => "deriving/dead_end.md",
              "MWC allostery" => "deriving/mwc_allostery.md",
          ],
          "Fitting rate equations" => [
              "Fitting tutorial & data format" => "fitting/tutorial.md",
              "Normalized vs absolute rate" => "fitting/normalized_vs_absolute.md",
              "Loss & optimizers" => "fitting/loss_and_optimizers.md",
          ],
          "Identifying the best rate equation" => [
              "Identify tutorial" => "identify/tutorial.md",
              "Model selection" => "identify/model_selection.md",
              "The enumeration engine" => "identify/enumeration_engine.md",
          ],
          "Developer / Architecture" => "developer.md",
          "API Reference" => "api.md",
          "Roadmap" => "roadmap.md",
          "References" => "references.md",
      ],
  )

  deploydocs(;
      repo = "github.com/DenisTitovLab/EnzymeRates.jl.git",
      devbranch = "main",
  )
  ```

- [ ] **Step 2: Commit (make.jl alone; it cannot build yet — stubs land next)**
  ```bash
  git add docs/make.jl
  git commit -m "Add docs/make.jl with full site tree and citation plugin

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
  ```

### Task 0.4: Home page + stub pages for the full site tree

**Files:**
- Create: `docs/src/index.md`
- Create: `docs/src/getting_started.md`
- Create: `docs/src/deriving/textbooks.md`
- Create: `docs/src/deriving/re_vs_ss.md`
- Create: `docs/src/deriving/cha_king_altman.md`
- Create: `docs/src/deriving/thermodynamic_constraints.md`
- Create: `docs/src/deriving/ping_pong.md`
- Create: `docs/src/deriving/dead_end.md`
- Create: `docs/src/deriving/mwc_allostery.md`
- Create: `docs/src/fitting/tutorial.md`
- Create: `docs/src/fitting/normalized_vs_absolute.md`
- Create: `docs/src/fitting/loss_and_optimizers.md`
- Create: `docs/src/identify/tutorial.md`
- Create: `docs/src/identify/model_selection.md`
- Create: `docs/src/identify/enumeration_engine.md`
- Create: `docs/src/developer.md`
- Create: `docs/src/api.md`
- Create: `docs/src/roadmap.md`
- Create: `docs/src/references.md`

- [ ] **Step 1: Write the Home page `docs/src/index.md`**
  Create `/home/denis.linux/.julia/dev/EnzymeRates/docs/src/index.md` with exactly this content. It is intentionally thin (the fuller Getting Started lands in Phase 5); it must reference no API names yet, so the build stays green before docstrings exist.

  ```markdown
  # EnzymeRates.jl

  EnzymeRates.jl identifies the best enzyme rate equation from kinetic data.
  Given a reaction definition and rate measurements at varying substrate and
  product concentrations, the package enumerates the biochemically valid
  mechanisms, fits each to the data, and selects the simplest equation that
  describes the data.

  !!! warning "Documentation under construction"
      This site is being built out. Pages currently marked *Coming soon* are
      placeholders; their content arrives in later phases.

  ## Installation

  ```julia
  using Pkg
  Pkg.add(url = "https://github.com/DenisTitovLab/EnzymeRates.jl")
  ```

  ## Where to go next

  - [Getting Started](@ref) — an end-to-end example from reaction to identified
    rate equation.
  - **Deriving rate equations** — how the package turns a mechanism into a
    symbolic rate law.
  - **Fitting rate equations** — fitting a rate equation to kinetic data.
  - **Identifying the best rate equation** — the model-selection search.
  - [API Reference](@ref) — every exported name.
  ```

  Note: the `[Getting Started](@ref)` and `[API Reference](@ref)` links resolve against the page H1 headers created below; if Documenter reports a missing cross-reference during the Task 0.7 build, confirm the H1 text matches exactly (`# Getting Started`, `# API Reference`).

- [ ] **Step 2: Create the directory layout and all 18 stub pages**
  Run from the repo root. This creates the two subdirectories and writes a uniform `# <Title>` + `Coming soon.` stub into each non-index page. The `api.md` and `developer.md` stubs stay plain "Coming soon" here — their `@autodocs`/`@index` blocks arrive in later phases, so they do not yet require docstrings to exist. The `references.md` page is the exception: it gets its real `@bibliography` block in this step, exercising the DocumenterCitations pipeline end-to-end as part of the Phase 0 proof.

  ```bash
  mkdir -p docs/src/deriving docs/src/fitting docs/src/identify

  printf '# Getting Started\n\nComing soon.\n' > docs/src/getting_started.md
  printf '# Rate equations from textbooks\n\nComing soon.\n' > docs/src/deriving/textbooks.md
  printf '# Rapid equilibrium vs steady state\n\nComing soon.\n' > docs/src/deriving/re_vs_ss.md
  printf '# The Cha / King–Altman algorithm\n\nComing soon.\n' > docs/src/deriving/cha_king_altman.md
  printf '# Thermodynamic constraints\n\nComing soon.\n' > docs/src/deriving/thermodynamic_constraints.md
  printf '# Ping-pong mechanisms\n\nComing soon.\n' > docs/src/deriving/ping_pong.md
  printf '# Dead-end inhibitor binding\n\nComing soon.\n' > docs/src/deriving/dead_end.md
  printf '# MWC allostery\n\nComing soon.\n' > docs/src/deriving/mwc_allostery.md
  printf '# Fitting tutorial & data format\n\nComing soon.\n' > docs/src/fitting/tutorial.md
  printf '# Normalized vs absolute rate\n\nComing soon.\n' > docs/src/fitting/normalized_vs_absolute.md
  printf '# Loss & optimizers\n\nComing soon.\n' > docs/src/fitting/loss_and_optimizers.md
  printf '# Identify tutorial\n\nComing soon.\n' > docs/src/identify/tutorial.md
  printf '# Model selection\n\nComing soon.\n' > docs/src/identify/model_selection.md
  printf '# The enumeration engine\n\nComing soon.\n' > docs/src/identify/enumeration_engine.md
  printf '# Developer / Architecture\n\nComing soon.\n' > docs/src/developer.md
  printf '# API Reference\n\nComing soon.\n' > docs/src/api.md
  printf '# Roadmap\n\nComing soon.\n' > docs/src/roadmap.md
  ```

  Then create `references.md` with its real bibliography block — this proves the DocumenterCitations pipeline renders. The lone `*` lists every entry in `refs.bib`, cited or not, so the page is populated even before any `@cite` exists:

  ````bash
  cat > docs/src/references.md <<'EOF'
  # References

  ```@bibliography
  *
  ```
  EOF
  ````

- [ ] **Step 3: Confirm all 19 source pages exist**
  ```bash
  find docs/src -name '*.md' | sort
  ```
  Expected: exactly these 19 paths (in sorted order):
  ```
  docs/src/api.md
  docs/src/deriving/cha_king_altman.md
  docs/src/deriving/dead_end.md
  docs/src/deriving/mwc_allostery.md
  docs/src/deriving/ping_pong.md
  docs/src/deriving/re_vs_ss.md
  docs/src/deriving/textbooks.md
  docs/src/deriving/thermodynamic_constraints.md
  docs/src/developer.md
  docs/src/fitting/loss_and_optimizers.md
  docs/src/fitting/normalized_vs_absolute.md
  docs/src/fitting/tutorial.md
  docs/src/getting_started.md
  docs/src/identify/enumeration_engine.md
  docs/src/identify/model_selection.md
  docs/src/identify/tutorial.md
  docs/src/index.md
  docs/src/references.md
  docs/src/roadmap.md
  ```
  If any `pages=` target from `make.jl` is missing here, the Task 0.7 build will error with `the following pages are missing` — create the missing stub before building.

- [ ] **Step 4: Commit**
  ```bash
  git add docs/src/index.md docs/src/getting_started.md docs/src/deriving docs/src/fitting docs/src/identify docs/src/developer.md docs/src/api.md docs/src/roadmap.md docs/src/references.md
  git commit -m "Add Home page and stub pages for full docs site tree

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
  ```

### Task 0.5: Gitignore the docs build output

**Files:**
- Modify: `.gitignore` (append one line; current file is the 6-line list ending at `.claude/scheduled_tasks.lock`)

- [ ] **Step 1: Append `docs/build/` to `.gitignore`**
  Edit `/home/denis.linux/.julia/dev/EnzymeRates/.gitignore`, adding `docs/build/` after the existing `.claude/scheduled_tasks.lock` line. The full file becomes:

  ```
  Manifest.toml
  .DS_Store
  .CondaPkg/
  test/outcmaes/
  outcmaes/
  .claude/scheduled_tasks.lock
  docs/build/
  ```

- [ ] **Step 2: Confirm `docs/build/` and `docs/Manifest.toml` are ignored**
  ```bash
  git check-ignore docs/build docs/Manifest.toml
  ```
  Expected: both paths print (each on its own line), confirming they are ignored — `docs/build` by the new rule, `docs/Manifest.toml` by the existing `Manifest.toml` rule. If `docs/build` does not print, the new line was not picked up.

- [ ] **Step 3: Commit**
  ```bash
  git add .gitignore
  git commit -m "Ignore docs/build output

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
  ```

### Task 0.6: docs CI job in CI.yml

**Files:**
- Modify: `.github/workflows/CI.yml` (append a `docs` job under `jobs:`, after the existing `test` job)

- [ ] **Step 1: Re-read the current CI.yml to confirm the append point**
  ```bash
  tail -n 5 .github/workflows/CI.yml
  ```
  Expected last lines:
  ```
      - uses: julia-actions/julia-processcoverage@v1
      - uses: codecov/codecov-action@v6
        with:
          token: ${{ secrets.CODECOV_TOKEN }}
  ```
  This confirms the `test` job ends here and the new `docs` job appends at the same indentation level as `test:` (2 spaces).

- [ ] **Step 2: Append the `docs` job**
  Edit `/home/denis.linux/.julia/dev/EnzymeRates/.github/workflows/CI.yml`, adding the following block as the last content of the file (it sits under the top-level `jobs:` key, a sibling of `test:`). The job uses `julia-actions/julia-docdeploy@v1`, grants `contents: write`, runs a "Configure doc environment" step that `dev`s the package into the docs env, and runs doctests as a separate step. Both `DOCUMENTER_KEY` and `GITHUB_TOKEN` are passed: Documenter prefers the SSH deploy key (Denis has installed it — a write-access deploy key plus the `DOCUMENTER_KEY` Actions secret, both verified present 2026-06-10) and falls back to the token if the key ever fails, so the deploy succeeds either way. `contents: write` is the load-bearing permission that lets the push reach `gh-pages`.

  Append exactly (note the leading two-space indent on `docs:`):

  ```yaml
  docs:
    name: Documentation
    runs-on: ubuntu-latest
    timeout-minutes: 30
    permissions:
      contents: write
      statuses: write
    steps:
      - uses: actions/checkout@v6
      - uses: julia-actions/setup-julia@v3
        with:
          version: '1'
      - uses: julia-actions/cache@v3
      - name: Configure doc environment
        shell: julia --project=docs --color=yes {0}
        run: |
          using Pkg
          Pkg.develop(PackageSpec(path=pwd()))
          Pkg.instantiate()
      - uses: julia-actions/julia-buildpkg@v1
      - uses: julia-actions/julia-docdeploy@v1
        env:
          DOCUMENTER_KEY: ${{ secrets.DOCUMENTER_KEY }}
          GITHUB_TOKEN: ${{ secrets.GITHUB_TOKEN }}
      - name: Run doctests
        shell: julia --project=docs --color=yes {0}
        run: |
          using Documenter: DocMeta, doctest
          using EnzymeRates
          doctest(EnzymeRates)
  ```

- [ ] **Step 3: Validate the YAML parses**
  ```bash
  julia --project=docs -e 'import Pkg; Pkg.add("YAML"); using YAML; d = YAML.load_file(".github/workflows/CI.yml"); println("jobs: ", join(keys(d["jobs"]), ", "))'
  ```
  Expected: prints `jobs: test, docs` (order may vary). No parse error. After confirming, undo the throwaway YAML add so it does not pollute the docs env:
  ```bash
  julia --project=docs -e 'import Pkg; Pkg.rm("YAML")'
  ```
  (Alternative if you prefer not to touch the docs env: `python3 -c "import yaml,sys; d=yaml.safe_load(open('.github/workflows/CI.yml')); print('jobs:', list(d['jobs']))"` — expect `jobs: ['test', 'docs']`.)

- [ ] **Step 4: Commit**
  ```bash
  git add .github/workflows/CI.yml
  git commit -m "Add docs CI job with docdeploy and doctest steps

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
  ```

### Task 0.7: Local build verification (the pipeline proof)

**Files:**
- (no file changes; this task proves the scaffold builds)

- [ ] **Step 1: Build the docs locally**
  ```bash
  julia --project=docs docs/make.jl
  ```
  Expected success signals, in order:
  - `[ Info: SetupBuildDirectory: setting up build directory.`
  - `[ Info: Doctest: running doctests.` — with **no** doctest failures (there are zero `jldoctest` blocks yet, so this passes trivially).
  - `[ Info: CheckDocument: running document checks.` — `checkdocs = :exports` passes because no `@docs`/`@autodocs` block references any export yet (the `api.md` stub has none).
  - `[ Info: Documenter: rendering done.`
  - The process exits 0. `deploydocs` prints `[ Info: Skipping deployment: ...` (no `CI` env locally) — that is expected and not a failure.
  No `ERROR`, no `Documenter could not auto-detect`, no `the following pages are not in pages` warning.

- [ ] **Step 2: Confirm the rendered output landed and is untracked**
  ```bash
  test -f docs/build/index.html && echo "BUILD OK" && git status --short docs/build
  ```
  Expected: prints `BUILD OK` and the `git status --short docs/build` line prints **nothing** (the build dir is gitignored, so it must not appear as untracked).

- [ ] **Step 3: Verify the citation machinery is wired (no live cite yet)**
  ```bash
  julia --project=docs -e 'using Documenter, DocumenterCitations, EnzymeRates; bib = CitationBibliography(joinpath("docs","src","refs.bib"); style=:authoryear); println(bib isa DocumenterCitations.CitationBibliography ? "cite plugin OK" : "FAIL")'
  ```
  Expected: prints `cite plugin OK`. This confirms the plugin object the build consumed is constructible against the committed `refs.bib`.

- [ ] **Step 4: Commit the docs Project.toml dev-entry if Pkg rewrote it**
  ```bash
  git status --short docs/Project.toml
  ```
  If `Pkg.develop` (Task 0.1) or `Pkg.rm("YAML")` (Task 0.6) left `docs/Project.toml` modified, commit it:
  ```bash
  git add docs/Project.toml
  git commit -m "Sync docs/Project.toml after instantiate

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
  ```
  If `git status --short` shows nothing, skip the commit — the file is already current.

### Task 0.8: Denis's manual step — enable GitHub Pages (checklist)

**Files:**
- (no file changes; this is a hand-off checklist for Denis, recorded so the deploy can serve)

This task is **not** executable by the agent — it requires Denis's GitHub repo-admin access. Authentication is already set up: Denis installed `DOCUMENTER_KEY` (a write-access deploy key plus the matching Actions secret, both verified present 2026-06-10), and the docs job passes it with `GITHUB_TOKEN` as a fallback. The one remaining manual action is flipping the Pages **source** so GitHub serves the branch — this can only happen after the first `docs` CI run on `main` creates the `gh-pages` branch.

- [ ] **Step 1: Let the `docs` job run once on `main`**
  After the `docs-comprehensive` PR merges to `main`, the `docs` job runs `julia-actions/julia-docdeploy@v1`, which creates the `gh-pages` branch and pushes the `dev/` build. Confirm the branch exists:
  ```bash
  git ls-remote --heads origin gh-pages
  ```
  Expected: one line ending in `refs/heads/gh-pages`. If empty, the docs job has not deployed yet — check the Actions run for the `docs` job.

- [ ] **Step 2: Hand off to Denis — set the Pages source**
  Tell Denis to set, one time, in the GitHub web UI: **Settings → Pages → Build and deployment → Source = "Deploy from a branch" → Branch = `gh-pages`, folder = `/ (root)` → Save.** Do **not** pick "GitHub Actions" as the source. The site then serves at `https://denistitovlab.github.io/EnzymeRates.jl/dev/`.

- [ ] **Step 3: Note that tagged-release docs are already wired**
  No further action: `TagBot.yml` already passes `ssh: ${{ secrets.DOCUMENTER_KEY }}`, so when Denis cuts a tagged release the tag triggers the `docs` workflow and the `stable/`/`vX.Y.Z/` docs build and deploy via the same key. (Per CLAUDE.md, do not write a summary `.md`; surface the first-deploy + Pages-toggle status in the return message.)

### Task 0.9: Phase 0 exit gate

**Files:**
- (no file changes; final green check)

- [ ] **Step 1: Re-run the full docs build to confirm reproducibility**
  ```bash
  julia --project=docs docs/make.jl
  ```
  Expected: `[ Info: Documenter: rendering done.`, exit 0, no doctest failures, no `checkdocs` failures, `deploydocs` reports `Skipping deployment` locally.

- [ ] **Step 2: Confirm the main test suite is untouched/green**
  Phase 0 adds no `src/` or `test/` changes, but confirm nothing regressed. The `docs-comprehensive` branch is based on `rate-equation-no-concentration-division` (NOT `main`), so diff against the branch's own start point, not `main`:
  ```bash
  base=$(git merge-base HEAD rate-equation-no-concentration-division)
  git diff --stat "$base"..HEAD -- src test
  ```
  Expected: empty output (no `src/` or `test/` changes in this phase). No need to run the slow full suite for an infra-only phase; the CLAUDE.md "run tests before committing" rule is satisfied for docs-only changes by the green `docs/make.jl` build, since Phase 0 touches no package code.

- [ ] **Step 3: Confirm the working tree is clean and on the docs branch**
  ```bash
  git status --short && git branch --show-current
  ```
  Expected: `git status --short` prints nothing (all Phase 0 work committed; `docs/build/` and `docs/Manifest.toml` gitignored); branch is `docs-comprehensive`.

**Phase 0 exit (per spec section 11):** docs build green locally; `docs` CI job present in `CI.yml` (green on the PR is verified once pushed); Pages-enable hand-off recorded (Task 0.8); first gh-pages deploy to be eyeballed after the `docs-comprehensive` PR merges to `main` and Denis sets the Pages source. No package code changed, so the full suite is unaffected.

I have all the verified facts. Now I'll write the Phase 1 section. I have confirmed:

- `EnzymeReaction` comment at types.jl:272-276; `@enzyme_reaction` needs atom brackets (`S[C]`); `show` → `EnzymeReaction: S ⇌ P`
- `metabolites` docstring at 1167-1172 detached by comment 1173-1176 before `@generated` at 1177
- `rate_equation_string` docstring at 548-552; mixed RE-binding + SS-iso uni-uni gives clean `1 + P/K_P_E + S/K_S_E`; pinned output captured
- `AllostericEnzymeMechanism` already documented (types.jl:821-829) — Task verifies, no write needed
- Exports list confirmed; `substrates`/`products` not exported

Here is the phase section.

---

## Phase 1 — Docstring gaps + API Reference page

This phase makes `checkdocs = :exports` pass and builds the consolidated API Reference page. Two exported names fail `checkdocs = :exports` today: `EnzymeReaction` (a `#`-comment, never a docstring) and `metabolites` (a real docstring detached from its `@generated` method by an intervening comment block). This phase writes the `EnzymeReaction` docstring, reattaches the `metabolites` docstring, expands the terse `rate_equation_string` docstring with a pinned `jldoctest`, verifies `AllostericEnzymeMechanism` is already documented, and adds `docs/src/api.md` (`@index` + grouped `@autodocs`). Every docstring example is captured from a live Julia process, never invented. Exit: `julia --project=docs docs/make.jl` builds with no missing-docstring errors and all doctests pass.

> Precondition: Phase 0 has created `docs/Project.toml`, `docs/make.jl` (with `checkdocs = :exports` and `doctest = true`), and `docs/src/` with the section stubs. If `docs/make.jl` does not yet exist when you start, stop and finish Phase 0 first — every verification step in this phase runs it.

### Task 1.1: Write the `EnzymeReaction` type docstring

**Files:**
- Modify: `src/types.jl:272-276` (replace the `#`-comment with a `"""` docstring)

The exported `EnzymeReaction` struct currently carries only a `#`-comment, so `checkdocs = :exports` fails on it. Convert the comment to a docstring, expand it to document the three fields and the canonicalization guarantee, and pin a `jldoctest` built with `@enzyme_reaction`. The `show` output `EnzymeReaction: S ⇌ P` is byte-stable (verified across a live process). `@enzyme_reaction` requires atom brackets on each reactant (`S[C]`), so the example uses them.

- [ ] **Step 1: Replace the `#`-comment block (`src/types.jl:272-276`) with this docstring**

Replace exactly this text:

```
# EnzymeReaction: the public concrete reaction descriptor. Holds
# reactants (substrate + product atom payload), regulators (with allowed
# multiplicity sets), and the catalytic multiplicities the enumerator is
# allowed to consider. Canonical ordering is enforced so two equivalent
# constructions compare equal under `==` / `hash`.
```

with:

````
"""
    EnzymeReaction

The public concrete reaction descriptor: substrates, products, optional
regulators, and the set of catalytic multiplicities the mechanism
enumerator is allowed to consider. It is the entry point to the package —
pair it with rate data in an [`IdentifyRateEquationProblem`](@ref), or pass
it to `EnzymeRates.init_mechanisms` to enumerate candidate mechanisms.

# Fields
- `reactants::Vector{ReactantAtoms}` — substrates and products, each
  carrying its per-atom inventory (used for ping-pong residual bookkeeping
  and atom conservation across steps).
- `regulators::Vector{RegulatorMults}` — competitive inhibitors and
  allosteric regulators, each with its allowed oligomeric multiplicities.
- `allowed_catalytic_multiplicities::Vector{Int}` — the oligomeric states
  the allosteric enumerator may assign to the catalytic core.

Construct one with the [`@enzyme_reaction`](@ref) DSL. Each reactant takes
an atom bracket (`S[C]`, `A[C1H1]`); the brackets are load-bearing for
ping-pong and multi-substrate reactions. Reactants and regulators are sorted
by name in the constructor, so two equivalent declarations compare equal
under `==`/`hash`.

```jldoctest
julia> rxn = @enzyme_reaction begin
           substrates: S[C]
           products:   P[C]
       end;

julia> rxn isa EnzymeReaction
true

julia> rxn
EnzymeReaction: S ⇌ P
```
"""
````

- [ ] **Step 2: Capture/verify the doctest output**
```bash
julia --project=docs -e 'using Documenter, EnzymeRates; doctest(EnzymeRates; fix=true)'
```
Expected: the doctest passes and `fix=true` leaves the pinned `true` and `EnzymeReaction: S ⇌ P` lines unchanged (they were captured from a live process). If `git diff src/types.jl` shows the doctest output was rewritten, re-read the rewritten block — the captured form is authoritative; keep whatever the fixer wrote and proceed.

- [ ] **Step 3: Build the docs to confirm `checkdocs` no longer flags `EnzymeReaction`**
```bash
julia --project=docs /home/denis.linux/.julia/dev/EnzymeRates/docs/make.jl
```
Expected: "Documenter: rendering done", no doctest failures. The build may still report `metabolites` as a missing/undocumented export until Task 1.2 — `EnzymeReaction` must no longer appear in any `checkdocs` warning.

- [ ] **Step 4: Commit**
```bash
git add src/types.jl
git commit -m "$(cat <<'EOF'
Document EnzymeReaction with a jldoctest example

Converts the #-comment into a type docstring so checkdocs=:exports
stops flagging EnzymeReaction. Documents the three fields and the
canonical-ordering guarantee; pins a @enzyme_reaction jldoctest.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
EOF
)"
```

### Task 1.2: Reattach the detached `metabolites` docstring

**Files:**
- Modify: `src/types.jl:1167-1192` (move the `@generated`-rationale comment from above the function to below its signature)

The `metabolites` docstring at `src/types.jl:1167-1172` is detached: the 4-line `#`-comment at `1173-1176` sits between the docstring and the `@generated function metabolites(...)` at `1177`, so Documenter attaches the docstring to the comment, not the method, and `checkdocs = :exports` reports `metabolites` undocumented. Moving the comment to the first line inside the function body reattaches the docstring while preserving the comment verbatim (CLAUDE.md forbids deleting comments). A `jldoctest` is added; `metabolites(m) == (:S, :P)` was captured from a live process for the uni-uni mechanism.

- [ ] **Step 1: Replace the docstring + comment + signature span (`src/types.jl:1167-1177`)**

Replace exactly this text:

```
"""
    metabolites(m::EnzymeMechanism) → Tuple{Symbol,...}

Return distinct metabolite names (substrates ∪ products ∪ regulators) as a tuple
of `Symbol`s in declaration order, deduplicated.
"""
# Kept `@generated` (unlike the other demoted accessors): `loss!` uses
# `metabolites(m)` as a compile-time-constant tuple to build the
# per-datapoint `NamedTuple{MetNames}` concs on the fitting hot path. A
# runtime body would make that NamedTuple type-unstable and allocate.
@generated function metabolites(::EnzymeMechanism{Sig}) where {Sig}
```

with:

````
"""
    metabolites(m::EnzymeMechanism) → Tuple{Symbol,...}

Return distinct metabolite names (substrates ∪ products ∪ regulators) as a tuple
of `Symbol`s in declaration order, deduplicated. The fitter uses this as the
key set for the per-datapoint concentration `NamedTuple` passed to
[`rate_equation`](@ref).

```jldoctest
julia> m = @enzyme_mechanism begin
           substrates: S
           products: P
           steps: begin
               E + S <--> E(S)
               E(S) <--> E + P
           end
       end;

julia> metabolites(m)
(:S, :P)
```
"""
@generated function metabolites(::EnzymeMechanism{Sig}) where {Sig}
    # Kept `@generated` (unlike the other demoted accessors): `loss!` uses
    # `metabolites(m)` as a compile-time-constant tuple to build the
    # per-datapoint `NamedTuple{MetNames}` concs on the fitting hot path. A
    # runtime body would make that NamedTuple type-unstable and allocate.
````

Note: the four comment lines are reindented to the 4-space function-body level and moved below the signature; their text is unchanged.

- [ ] **Step 2: Confirm the function body still starts correctly**

Read `src/types.jl:1177-1200` and verify the line immediately after the moved comment block is `    m = Mechanism(EnzymeMechanism{Sig}())` (the original first body line). The comment must sit between the signature line and `m = Mechanism(...)`.

- [ ] **Step 3: Capture/verify the doctest output**
```bash
julia --project=docs -e 'using Documenter, EnzymeRates; doctest(EnzymeRates; fix=true)'
```
Expected: passes; the pinned `(:S, :P)` line is unchanged. If `fix` rewrites it, keep the fixer's output.

- [ ] **Step 4: Build the docs**
```bash
julia --project=docs /home/denis.linux/.julia/dev/EnzymeRates/docs/make.jl
```
Expected: "Documenter: rendering done", no doctest failures, and `metabolites` no longer appears in any `checkdocs` missing-docstring warning.

- [ ] **Step 5: Run the package test suite for the `metabolites` hot-path tests (the comment move must not change behavior)**
```bash
julia --project -e 'using Pkg; Pkg.test()' 2>&1 | tail -30
```
Expected: full suite green (the edit only relocates a comment; the `@generated` body is byte-identical). If anything fails, the comment landed in the wrong place — re-check Step 2.

- [ ] **Step 6: Commit**
```bash
git add src/types.jl
git commit -m "$(cat <<'EOF'
Reattach metabolites docstring; add a jldoctest

The @generated-rationale comment sat between the docstring and the
method, detaching the docstring and failing checkdocs=:exports. Move
the comment into the function body (verbatim) and add a jldoctest.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
EOF
)"
```

### Task 1.3: Expand the `rate_equation_string` docstring

**Files:**
- Modify: `src/rate_eq_derivation.jl:548-552` (replace the two-line docstring above `function rate_equation_string end`)

The current docstring (`src/rate_eq_derivation.jl:548-552`) is two terse lines: it documents neither the `mode` argument, the multi-line format, nor the Full-mode restriction. Expand it from the verified behavior of the methods at `:555-608`: `mode` defaults to `Reduced`; `Reduced` (`:583`) emits a `params`-destructure line, a `concs`-destructure line, `# Wegscheider constraints:` / `# Haldane constraints:` sections for dependent parameters, then the `v = E_total * (num) / (den)` line; `Full` (`:574`) emits the destructure lines plus the raw `v` line with no constraint section (every rate constant is independent). The `jldoctest` output below was captured from a live process for the mixed RE-binding + SS-iso uni-uni mechanism, which yields the literature-style `1 + P/K_P_E + S/K_S_E` denominator. Per the audit's doctest rules, the example uses `print(...)` (the function returns a multi-line `String`; `print` avoids quote/escape noise).

- [ ] **Step 1: Read the current docstring to confirm the exact text to replace**

Read `src/rate_eq_derivation.jl:548-553`. Confirm it is:
```
"""
    rate_equation_string(m, [mode])

Return a string representation of the rate equation.
"""
function rate_equation_string end
```

- [ ] **Step 2: Replace the docstring (`src/rate_eq_derivation.jl:548-552`)**

Replace exactly this text:

```
"""
    rate_equation_string(m, [mode])

Return a string representation of the rate equation.
"""
```

with:

````
"""
    rate_equation_string(m, [mode]) -> String

Return the symbolic rate equation for mechanism `m` as a multi-line
`String` (it returns the text — it does not print). `mode` is `Reduced`
(default) or `Full`; pass a concrete [`Mechanism`](@ref) /
`AllostericMechanism` or its compiled `EnzymeMechanism` singleton.

The string is a runnable transcript of how [`rate_equation`](@ref)
evaluates: a `(; …) = params` destructure line, a `(; …) = concs`
destructure line, then the `v = E_total * (num) / (den)` line. In
`Reduced` mode, dependent rate constants are listed first under
`# Wegscheider constraints:` and `# Haldane constraints:` headers — the
thermodynamic identities that eliminate parameters — and only the
independent set appears in the `params` destructure. In `Full` mode every
rate constant is independent, so there is no constraint section. `Full`
mode is defined for `EnzymeMechanism` only; an `AllostericEnzymeMechanism`
supports `Reduced` mode only.

Use [`print`](@ref) on the result to see the multi-line layout without
escaped newlines.

```jldoctest
julia> m = @enzyme_mechanism begin
           substrates: S
           products: P
           steps: begin
               E + S ⇌ E(S)
               E(S) <--> E(P)
               E(P) ⇌ E + P
           end
       end;

julia> print(rate_equation_string(m))
(; K_P_E, K_S_E, k_ES_to_EP, Keq, E_total) = params
(; S, P) = concs
# Haldane constraints:
k_EP_to_ES = (1 / Keq) * K_P_E * (1 / K_S_E) * k_ES_to_EP
v = E_total * (k_ES_to_EP * S / K_S_E - k_EP_to_ES * P / K_P_E) / (1 + P / K_P_E + S / K_S_E)
```
"""
````

- [ ] **Step 3: Capture/verify the doctest output**
```bash
julia --project=docs -e 'using Documenter, EnzymeRates; doctest(EnzymeRates; fix=true)'
```
Expected: passes; the five pinned output lines are unchanged (captured from a live process). If `fix` rewrites any line, the canonical printer changed — keep the fixer's output and note it for the commit.

- [ ] **Step 4: Build the docs**
```bash
julia --project=docs /home/denis.linux/.julia/dev/EnzymeRates/docs/make.jl
```
Expected: "Documenter: rendering done", no doctest failures.

- [ ] **Step 5: Commit**
```bash
git add src/rate_eq_derivation.jl
git commit -m "$(cat <<'EOF'
Expand rate_equation_string docstring with mode + format + jldoctest

Documents the mode argument (Reduced/Full), the multi-line transcript
format with Wegscheider/Haldane constraint sections, the Full-mode
EnzymeMechanism-only restriction, and that it returns a String. Pins a
jldoctest on a uni-uni RE mechanism.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
EOF
)"
```

### Task 1.4: Verify `AllostericEnzymeMechanism` is documented

**Files:**
- Inspect only: `src/types.jl:821-829`

The spec lists "Add a type docstring for `AllostericEnzymeMechanism`" as a default. Verification confirms it already has a `"""` docstring at `src/types.jl:821-829`, so no write is needed — this task only proves the export is covered by `checkdocs`.

- [ ] **Step 1: Confirm the existing docstring**

Read `src/types.jl:821-836`. Confirm the `"""..."""` block (starting at `821` with `    AllostericEnzymeMechanism{CatalyticMech, CatSites, RegSites}`) is immediately followed (after a single `#`-comment that begins at `830`) by `struct AllostericEnzymeMechanism{` at `834`.

- [ ] **Step 2: Check whether the docstring is attached or detached**

If the `#`-comment at `src/types.jl:830-833` sits between the docstring and the `struct` keyword, the docstring is detached (same failure mode as `metabolites`) and you must move that comment below the struct's first line — apply the same fix pattern as Task 1.2: relocate the comment lines verbatim into the struct body, keeping them unchanged. If `checkdocs` does not flag `AllostericEnzymeMechanism` in Step 3, the docstring is attached and no edit is needed.

> Note: a docstring placed directly above a `#`-comment that is itself directly above a `struct` does NOT always attach in Julia — the same detachment that breaks `metabolites` can apply here. Step 3 is the authoritative check.

- [ ] **Step 3: Build the docs and grep the output for the export**
```bash
julia --project=docs /home/denis.linux/.julia/dev/EnzymeRates/docs/make.jl 2>&1 | grep -i "AllostericEnzymeMechanism" || echo "NOT FLAGGED — documented"
```
Expected: `NOT FLAGGED — documented` (no `checkdocs` warning names `AllostericEnzymeMechanism`). If it IS flagged, apply the comment-relocation fix from Step 2, re-run, and confirm it clears; then commit with message `Reattach AllostericEnzymeMechanism docstring`.

- [ ] **Step 4: Commit (only if an edit was made in Step 2)**
```bash
git add src/types.jl
git commit -m "$(cat <<'EOF'
Reattach AllostericEnzymeMechanism docstring

Move the type-parameter-rationale comment into the struct body so the
docstring attaches and checkdocs=:exports covers the export.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
EOF
)"
```
If no edit was made, skip this commit — the export was already documented.

### Task 1.5: Create the consolidated API Reference page

**Files:**
- Create: `docs/src/api.md`
- Inspect: `docs/make.jl` (confirm the `pages` list includes the API page; add it if Phase 0's stub omitted it)

Build the API Reference page with `@index` topped over four `@autodocs` blocks (`Modules = [EnzymeRates]`, `Private = false`), grouped by `Order` so types, macros, constants, and functions render in stable sections. Because the blocks are `Private = false`, adding a future export surfaces it with no edit here; `checkdocs = :exports` is the guarantee the page stays complete.

- [ ] **Step 1: Create `docs/src/api.md` with this exact content**

````markdown
# API Reference

```@meta
CurrentModule = EnzymeRates
```

This page lists every exported name. Maintainer-only internals are not
shown here — see the Developer / Architecture page for those.

## Index

```@index
```

## Types

```@autodocs
Modules = [EnzymeRates]
Private = false
Order = [:type]
```

## Macros

```@autodocs
Modules = [EnzymeRates]
Private = false
Order = [:macro]
```

## Constants

```@autodocs
Modules = [EnzymeRates]
Private = false
Order = [:constant]
```

## Functions

```@autodocs
Modules = [EnzymeRates]
Private = false
Order = [:function]
```
````

- [ ] **Step 2: Ensure `docs/make.jl`'s `pages` list references the API page**

Read `docs/make.jl`. If its `makedocs(...)` call has a `pages = [...]` argument, confirm it contains an entry pointing at `api.md` (e.g. `"API Reference" => "api.md"`). If the entry is missing, add it to the `pages` list in the position the Phase 0 site structure places it (after the Developer page, before Roadmap), matching the surrounding entry style exactly. If `makedocs` has no `pages` argument (auto-discovery), no edit is needed — `api.md` is picked up automatically.

- [ ] **Step 3: Build the docs**
```bash
julia --project=docs /home/denis.linux/.julia/dev/EnzymeRates/docs/make.jl
```
Expected: "Documenter: rendering done", no doctest failures, and no `checkdocs` warning for any export. The build log should show the four `@autodocs` blocks resolving (no "no docs found for" errors). Open `docs/build/api/index.html` mentally via the log — every exported name from `src/EnzymeRates.jl` (the 18 names across the 9 `export` lines) must appear under exactly one section.

- [ ] **Step 4: Confirm `checkdocs = :exports` passes cleanly end-to-end**
```bash
julia --project=docs /home/denis.linux/.julia/dev/EnzymeRates/docs/make.jl 2>&1 | grep -iE "checkdocs|missing docstring|no docs found|doctest" || echo "CLEAN — no checkdocs/doctest issues"
```
Expected: `CLEAN — no checkdocs/doctest issues`. Any line printed here is a remaining gap — resolve it (a missing docstring means a Task 1.1–1.4 edit did not land) before committing.

- [ ] **Step 5: Run the full package test suite (the source edits must not regress)**
```bash
julia --project -e 'using Pkg; Pkg.test()' 2>&1 | tail -20
```
Expected: full suite green. The Phase 1 edits touch only docstrings/comments; `test/test_types.jl` (the `name(p,m)` AST-walker chokepoint at `:1577-1644`) and the `rate_equation` performance tests must still pass.

- [ ] **Step 6: Commit**
```bash
git add docs/src/api.md docs/make.jl
git commit -m "$(cat <<'EOF'
Add consolidated API Reference page

@index over grouped @autodocs (types, macros, constants, functions;
Private=false) so new exports surface automatically. With the docstring
fixes, checkdocs=:exports now passes.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
EOF
)"
```

All facts verified. I now have a verified textbook jldoctest, a verified SS contrast, dead-end inhibitor, ping-pong (enumerator), and MWC examples. Let me check what `_default_save_dir` / `init_mechanisms` is named and confirm `init_mechanisms` is reachable as `EnzymeRates.init_mechanisms` for the ping-pong page (it is, per CLAUDE.md). I have enough to write the phase. Let me draft the complete Phase 2 section.

Based on all my verification, here is the Phase 2 section:

## Phase 2 — Derivation pillar

This phase writes the seven concept pages under "Deriving rate equations," in the spec's post-review order. Each page is prose authored with the `elements-of-style:writing-clearly-and-concisely` skill, built on code-verified facts cited to `src/` functions, and carries runnable `@example`/`jldoctest` blocks. The textbook page anchors the pillar with the one audit-confirmed byte-stable `jldoctest`; the other six use `@example` (run, output captured, not output-checked) because their derived strings, while deterministic, are large and not worth pinning byte-for-byte. Every page standardizes MWC terminology on **A/I** (active/inactive) and avoids every stale name in the spec's section 8 staleness ledger. The seven pages are added to the `Deriving rate equations` nav block in `docs/make.jl` (the Phase 0 stub created the section with placeholder entries; each task repoints one entry to a real file). Each task ends green: `julia --project=docs /home/denis.linux/.julia/dev/EnzymeRates/docs/make.jl` builds without warnings, with no doctest failures, ending in `Documenter: rendering done`.

Throughout this phase, the `pages` entry in `docs/make.jl` for the derivation section must end up exactly:

```julia
    "Deriving rate equations" => [
        "Rate equations from textbooks" => "deriving/textbooks.md",
        "Rapid equilibrium vs steady state" => "deriving/re_vs_ss.md",
        "The Cha / King–Altman algorithm" => "deriving/cha_king_altman.md",
        "Thermodynamic constraints" => "deriving/thermodynamic_constraints.md",
        "Ping-pong mechanisms" => "deriving/ping_pong.md",
        "Dead-end inhibitor binding" => "deriving/dead_end.md",
        "MWC allostery" => "deriving/mwc_allostery.md",
    ],
```

When a task says "repoint the nav entry," edit only the one `"… " => "deriving/<file>.md"` line for that page so the build picks it up; leave the others as the Phase 0 stub left them until their own task runs.

---

### Task 2.1: Rate equations from textbooks

**Files:**
- Create: `docs/src/deriving/textbooks.md`
- Modify: `docs/make.jl` (the one `"Rate equations from textbooks" => …` line in the `pages` vector)

- [ ] **Step 1: Write `docs/src/deriving/textbooks.md`**
  Author the prose with the `elements-of-style:writing-clearly-and-concisely` skill. This is the pillar's opening tutorial: a reader should be able to define a reversible single-substrate enzyme, derive its rate equation as a String, and recognize it as reversible Michaelis–Menten.

  The page MUST cover these code-verified facts (cite the function in an HTML comment next to each claim so a future editor can re-verify; do not print the comments to the reader):
  - `rate_equation_string(m, [mode])` **returns a `String`** (it does not print). Default mode is `Reduced` (`src/rate_eq_derivation.jl:555`). This corrects the README "prints" drift (staleness ledger, section 8).
  - The `Reduced` mode string has the shape: a `params` destructuring line, a `concs` destructuring line, an optional `# Wegscheider constraints:` / `# Haldane constraints:` section, and the final `v = E_total * (num) / (den)` line (`rate_equation_string(::M, ::ReducedMode)`, `src/rate_eq_derivation.jl:583-608`).
  - Binding constants are named by metabolite and enzyme form: `K_S_E`, `K_P_E`. The single SS catalytic rate constant is directed: `k_ES_to_EP`; its reverse `k_EP_to_ES` is the Haldane-dependent parameter (CLAUDE.md "Parameter naming convention"; `name(p, m)` chokepoint at `src/types.jl`).
  - In `Reduced` mode one rate constant per thermodynamic cycle is dependent and substituted; here `k_EP_to_ES` is fixed by `Keq` and the binding constants (`_dependent_param_exprs`, `src/thermodynamic_constr_for_rate_eq_derivation.jl:263`). `Keq` is always user-supplied, never fitted (CLAUDE.md).
  - A symbol-mapping note (spec section 10 default): package `K_S_E`, `K_P_E` ↔ textbook `Km`, `Kp`; package `k_ES_to_EP` ↔ textbook `kcat`. State that the RE denominator `1 + S/K_S_E + P/K_P_E` is the literature reversible-MM form.
  - The `Full` mode exists and lists all raw rate constants plus `E_total` (`parameters(m, Full)` → `(:k_ES_to_EP, :k_EP_to_ES, :K_S_E, :K_P_E, :E_total)`); the page should mention it in one sentence and point the reader to the Rapid-equilibrium-vs-steady-state page for the contrast. Do not pin the Full-mode string here.

  The page MUST contain this `jldoctest` block verbatim (input + the pinned, audit-confirmed byte-stable output). Use a named doctest block so later pages can reuse the binding if needed:

  ````markdown
  ```jldoctest textbook
  julia> using EnzymeRates

  julia> m = @enzyme_mechanism begin
             substrates: S
             products:   P
             steps: begin
                 E + S ⇌ E(S)
                 E(S) <--> E(P)
                 E(P) ⇌ E + P
             end
         end;

  julia> m isa EnzymeMechanism
  true

  julia> print(rate_equation_string(m))
  (; K_P_E, K_S_E, k_ES_to_EP, Keq, E_total) = params
  (; S, P) = concs
  # Haldane constraints:
  k_EP_to_ES = (1 / Keq) * K_P_E * (1 / K_S_E) * k_ES_to_EP
  v = E_total * (k_ES_to_EP * S / K_S_E - k_EP_to_ES * P / K_P_E) / (1 + P / K_P_E + S / K_S_E)
  ```
  ````

  Follow the doctest with an `@example` block (run, not output-checked) that shows `parameters(m)` and `metabolites(m)` so the reader sees the fitted-parameter tuple and the concentration symbols:

  ````markdown
  ```@example textbook_ex
  using EnzymeRates
  m = @enzyme_mechanism begin
      substrates: S
      products:   P
      steps: begin
          E + S ⇌ E(S)
          E(S) <--> E(P)
          E(P) ⇌ E + P
      end
  end
  (parameters(m), metabolites(m))
  ```
  ````

  Cross-link `[rate_equation_string](@ref)`, `[parameters](@ref)`, `[metabolites](@ref)`, and `[@enzyme_mechanism](@ref)` into the API Reference (the `[name](@ref)` form per spec section 2). Show `m isa EnzymeMechanism`, never `typeof(m)` (spec section 3 doctest rule). Standardize on `print(rate_equation_string(m))`, not the bare return value (spec section 3).

- [ ] **Step 2: Repoint the nav entry in `docs/make.jl`**
  Set the derivation-section line to `"Rate equations from textbooks" => "deriving/textbooks.md",` (it already is per the Phase 0 stub if Phase 0 named files this way; if the stub used a placeholder filename, change it to `deriving/textbooks.md`).

- [ ] **Step 3: Capture/confirm the doctest output**
  The output above is the audit-verified byte-stable block, so it is pinned. Run the capture command once to confirm it matches (it will report no changes if correct):
  ```bash
  julia --project=docs -e 'using Documenter, EnzymeRates; doctest(EnzymeRates; fix=true)'
  ```
  Expected: the `jldoctest textbook` block is unchanged (git shows no diff to the pinned output). If `fix=true` rewrites it, the output drifted — STOP and reconcile before committing.

- [ ] **Step 4: Build the docs**
  ```bash
  julia --project=docs /home/denis.linux/.julia/dev/EnzymeRates/docs/make.jl
  ```
  Expected: ends with `Documenter: rendering done`, no doctest failures, no missing-page warnings.

- [ ] **Step 5: Commit**
  ```bash
  git add docs/src/deriving/textbooks.md docs/make.jl
  git commit -m "Add 'Rate equations from textbooks' derivation page

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
  ```

---

### Task 2.2: Rapid equilibrium vs steady state

**Files:**
- Create: `docs/src/deriving/re_vs_ss.md`
- Modify: `docs/make.jl` (the `"Rapid equilibrium vs steady state" => …` line)

- [ ] **Step 1: Write `docs/src/deriving/re_vs_ss.md`**
  Author the prose with the `elements-of-style:writing-clearly-and-concisely` skill. The page explains the one structural distinction that drives parameter count: whether a step is rapid-equilibrium (RE) or steady-state (SS).

  Code-verified facts to cover (cite each):
  - A `Step` carries `is_equilibrium::Bool` (`src/types.jl:142`); `is_equilibrium(s)` reads it (`src/types.jl:167`). In the DSL, `⇌` parses to `is_equilibrium = true` and `<-->` to `false` (`_step_struct_info`, `src/dsl.jl:890-900`).
  - An RE step contributes **one** parameter — a binding `Kd` (rendered `K_<met>_<form>`) or an isomerization `Kiso` — while an SS step contributes **two** rate constants — `Kon`+`Koff` (binding) or `Kfor`+`Krev` (iso) (`_step_parameters`, `src/thermodynamic_constr_for_rate_eq_derivation.jl:36-46`).
  - RE binding parameters use metabolite/form names (`K_S_E`); SS binding parameters render as `kon_S_E` / `koff_S_E`; SS iso rates render directionally as `k_ES_to_EP` / `k_EP_to_ES` (CLAUDE.md "Parameter naming convention"; `name(p, m)` chokepoint).
  - RE assumes the binding step stays at equilibrium relative to catalysis; SS makes no such assumption and solves the full King–Altman / Cha steady state (link to the Cha page).
  - After thermodynamic reduction, each independent cycle still costs one dependent rate constant regardless of RE/SS, but SS steps add more *independent* parameters because they start with two rate constants each (contrast the two `@example` outputs below).

  Include two `@example` blocks (run, not output-checked) — the all-RE textbook mechanism and the same skeleton made all-SS — and contrast their `parameters(...)` tuples in prose:

  ````markdown
  ```@example revss
  using EnzymeRates
  re = @enzyme_mechanism begin
      substrates: S
      products:   P
      steps: begin
          E + S ⇌ E(S)
          E(S) <--> E(P)
          E(P) ⇌ E + P
      end
  end
  parameters(re)
  ```

  ```@example revss
  ss = @enzyme_mechanism begin
      substrates: S
      products:   P
      steps: begin
          E + S <--> E(S)
          E(S) <--> E(P)
          E(P) <--> E + P
      end
  end
  parameters(ss)
  ```
  ````

  In prose, state the verified contrast: the all-RE form fits `(:K_P_E, :K_S_E, :k_ES_to_EP)` (plus `Keq`, `E_total`), while the all-SS form fits `(:k_ES_to_EP, :koff_P_E, :koff_S_E, :kon_P_E, :kon_S_E)` (plus `Keq`, `E_total`) — more parameters because each binding step now carries a separate on/off rate. Add a one-sentence note that the `RE→SS` expansion move flips a whole kinetic group atomically, linking to the enumeration page (do not re-derive it here).

  Cross-link `[parameters](@ref)` and `[@enzyme_mechanism](@ref)`.

- [ ] **Step 2: Repoint the nav entry**
  Set `"Rapid equilibrium vs steady state" => "deriving/re_vs_ss.md",` in `docs/make.jl`.

- [ ] **Step 3: Build the docs**
  ```bash
  julia --project=docs /home/denis.linux/.julia/dev/EnzymeRates/docs/make.jl
  ```
  Expected: ends with `Documenter: rendering done`; both `@example revss` blocks run; no warnings.

- [ ] **Step 4: Commit**
  ```bash
  git add docs/src/deriving/re_vs_ss.md docs/make.jl
  git commit -m "Add 'Rapid equilibrium vs steady state' derivation page

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
  ```

---

### Task 2.3: The Cha / King–Altman algorithm

**Files:**
- Create: `docs/src/deriving/cha_king_altman.md`
- Modify: `docs/make.jl` (the `"The Cha / King–Altman algorithm" => …` line)

- [ ] **Step 1: Write `docs/src/deriving/cha_king_altman.md`**
  Author the prose with the `elements-of-style:writing-clearly-and-concisely` skill. The page explains how the package derives a rate equation symbolically at compile time, with emphasis on the division-free, finite-at-zero design.

  Code-verified facts to cover (cite each):
  - Derivation runs at compile time inside `@generated` functions; each unique mechanism type triggers a full symbolic derivation, so very large mechanisms are slow to compile (CLAUDE.md "Known Issues"; `@generated rate_equation`, `src/rate_eq_derivation.jl:533-544`).
  - The mechanism graph is partitioned into rapid-equilibrium segments; each segment is referenced to its **free enzyme** — the form with the fewest bound metabolites (tie-broken toward no covalent residual, then a deterministic name) via `_segment_root` (`src/rate_eq_derivation.jl:263-267`). Within a segment, relative form abundances ("alpha" factors) are computed against that root (`_compute_alpha`, `src/rate_eq_derivation.jl:277-322`). This is what produces the readable `1 + S/K_S_E + P/K_P_E` denominator.
  - Segments are connected by SS steps; the Cha treatment builds the inter-segment rate matrix and takes cofactor determinants (`_raw_symbolic_rate_polys`, `src/rate_eq_derivation.jl:339-399`; `sym_det`, `src/sym_poly_for_rate_eq_derivation.jl:88-114`).
  - The algebra is **Laurent** (exponents may be negative): `MONO = Vector{Pair{Symbol,Int}}`, `POLY = Dict{MONO, Rational{Int}}` (`src/sym_poly_for_rate_eq_derivation.jl:9-10`). Alphas and numerator/denominator are kept as fractional POLYs directly — no common-denominator linearization (CLAUDE.md `src/rate_eq_derivation.jl` layout note).
  - `_reduce_conc_lowest_terms(num, den, conc_set)` shifts each **concentration** symbol's exponent so its minimum across num ∪ den is 0, clearing any `1/conc` coupling so the derived equation is finite at zero concentration. It only ever shifts concentration symbols, never parameters — so it cannot drop a fitted parameter (`src/sym_poly_for_rate_eq_derivation.jl:56-81`; `_concentration_symbols`, `src/rate_eq_derivation.jl:249-255`).
  - The runtime `rate_equation` call is allocation-free and sub-100 ns: the emitted Expr is a balanced binary `+`/`*` tree (`_nest_binary`, `src/sym_poly_for_rate_eq_derivation.jl:162-170`) so Julia fuses scalar arithmetic. Link to the `rate_equation` runtime-perf contract on the Developer page.
  - `MAX_RATE_EQUATION_TERMS = 5000`: derivation aborts with an error past this term count (`src/sym_poly_for_rate_eq_derivation.jl:7`, `sym_det` check at `:102-111`).

  Include one `@example` block (run, not output-checked) that derives a mechanism whose naive per-segment reference would introduce a `1/conc` term, then shows the final equation has a constant `1 +` denominator (finite at zero). Use the textbook RE mechanism and print its string, calling out in prose that the denominator's leading constant `1` is the finite-at-zero result of `_reduce_conc_lowest_terms`:

  ````markdown
  ```@example cha
  using EnzymeRates
  m = @enzyme_mechanism begin
      substrates: S
      products:   P
      steps: begin
          E + S ⇌ E(S)
          E(S) <--> E(P)
          E(P) ⇌ E + P
      end
  end
  print(rate_equation_string(m))
  ```
  ````

  In prose, point at the `1 + P / K_P_E + S / K_S_E` denominator and explain that referencing each segment to its free enzyme plus the concentration-GCD is what keeps it division-free.

  Add an `@example` that evaluates `rate_equation` at zero substrate to make "finite at zero" concrete:

  ````markdown
  ```@example cha
  rate_equation(m, (; S = 0.0, P = 0.0),
                (; K_S_E = 1.0, K_P_E = 1.0, k_ES_to_EP = 10.0, Keq = 2.0, E_total = 1.0))
  ```
  ````

  State in prose that the result is `0.0` (zero net rate at zero concentrations) and finite — no division-by-zero — which is the whole point of the conc-GCD step.

  Cross-link `[rate_equation](@ref)` and `[rate_equation_string](@ref)`. Add the citations `[chaUsefulnessNetRate1968](@cite)` (Cha 1968) and `[kingSchematicMethodDeriving1956](@cite)` (King & Altman 1956) once each in the opening paragraph (exact citation keys come from Phase 0 `refs.bib`; if the Phase 0 keys differ, use the keys present in `docs/src/refs.bib`).

- [ ] **Step 2: Repoint the nav entry**
  Set `"The Cha / King–Altman algorithm" => "deriving/cha_king_altman.md",` in `docs/make.jl`.

- [ ] **Step 3: Build the docs**
  ```bash
  julia --project=docs /home/denis.linux/.julia/dev/EnzymeRates/docs/make.jl
  ```
  Expected: ends with `Documenter: rendering done`; both `@example cha` blocks run (the second prints `0.0`); the two `@cite` calls resolve (no "citation not found" warning); no warnings.

- [ ] **Step 4: Commit**
  ```bash
  git add docs/src/deriving/cha_king_altman.md docs/make.jl
  git commit -m "Add 'Cha / King–Altman algorithm' derivation page

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
  ```

---

### Task 2.4: Thermodynamic constraints

**Files:**
- Create: `docs/src/deriving/thermodynamic_constraints.md`
- Modify: `docs/make.jl` (the `"Thermodynamic constraints" => …` line)

- [ ] **Step 1: Write `docs/src/deriving/thermodynamic_constraints.md`**
  Author the prose with the `elements-of-style:writing-clearly-and-concisely` skill. The page explains why some rate constants are not free parameters: thermodynamics ties them to `Keq` (Haldane) or to each other (Wegscheider).

  Code-verified facts to cover (cite each):
  - The package finds thermodynamic cycles as the **exact-integer null space** of the enzyme incidence matrix (`_integer_nullspace`, `src/thermodynamic_constr_for_rate_eq_derivation.jl:111-145`); the computation uses `Rational{BigInt}` and reduces each null-space vector to a primitive integer vector with a sign convention — no floating point, so the constraint structure is exact.
  - Each cycle is classified as **Haldane** (its net metabolite change is proportional to the overall reaction → it carries a power of `Keq`) or **Wegscheider** (closed internal loop, zero net change → a pure rate-constant ratio with no `Keq`). A cycle whose metabolite change is neither zero nor proportional to the net reaction is a hard error (`classify_cycle`, `src/thermodynamic_constr_for_rate_eq_derivation.jl:204-230`; `_thermodynamic_constraints`, `:147-235`).
  - One rate constant per independent cycle is made **dependent** and expressed via `build_power_expr` as `Keq^a * prod(k_i^b_i)` (`src/thermodynamic_constr_for_rate_eq_derivation.jl:407-415`; `build_power_expr`, `src/sym_poly_for_rate_eq_derivation.jl:229-260`).
  - Which constant is chosen dependent is decided by **canonical step order** plus a structural pivot priority — internal isomerizations are eliminated before metabolite steps before free-enzyme binding (`_step_priority`, `:85-91`; Gaussian elimination with priority pivoting, `:372-405`). Because step order is canonicalized in the `Mechanism` constructor, the dependent-parameter choice is deterministic and rate-equivalent mechanisms reduce identically (CLAUDE.md "Canonical Step Form").
  - A thermodynamically contradictory mechanism (a constraint row reducing to `0 = c * log(Keq)` with `c ≠ 0`) raises an error rather than producing a silent equation (`src/thermodynamic_constr_for_rate_eq_derivation.jl:386-392`).
  - In `rate_equation_string`, single-symbol Wegscheider ties are shown under `# Wegscheider constraints:` and marked `(substituted into v)`; `Keq`-bearing dependents appear under `# Haldane constraints:` (`src/rate_eq_derivation.jl:583-608`; `ANNOTATION_SUBSTITUTED`, `:15`).

  Include an `@example` block (run, not output-checked) that shows a Haldane line in the rendered string (the textbook mechanism has exactly one Haldane dependent):

  ````markdown
  ```@example thermo
  using EnzymeRates
  m = @enzyme_mechanism begin
      substrates: S
      products:   P
      steps: begin
          E + S ⇌ E(S)
          E(S) <--> E(P)
          E(P) ⇌ E + P
      end
  end
  print(rate_equation_string(m))
  ```
  ````

  In prose, point at the line `k_EP_to_ES = (1 / Keq) * K_P_E * (1 / K_S_E) * k_ES_to_EP` and explain it as the Haldane relation: the reverse SS rate is fixed once `Keq` and the binding constants are known, so it is not a fitted parameter. Note that `parameters(m)` therefore omits `k_EP_to_ES`.

  Cross-link `[parameters](@ref)` and `[rate_equation_string](@ref)`. Add the Haldane/Wegscheider citations once in the intro (`[haldane1930enzymes](@cite)` and the Wegscheider entry — use the keys present in `docs/src/refs.bib`).

- [ ] **Step 2: Repoint the nav entry**
  Set `"Thermodynamic constraints" => "deriving/thermodynamic_constraints.md",` in `docs/make.jl`.

- [ ] **Step 3: Build the docs**
  ```bash
  julia --project=docs /home/denis.linux/.julia/dev/EnzymeRates/docs/make.jl
  ```
  Expected: ends with `Documenter: rendering done`; the `@example thermo` block runs and renders the Haldane line; `@cite` keys resolve; no warnings.

- [ ] **Step 4: Commit**
  ```bash
  git add docs/src/deriving/thermodynamic_constraints.md docs/make.jl
  git commit -m "Add 'Thermodynamic constraints' derivation page

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
  ```

---

### Task 2.5: Ping-pong mechanisms

**Files:**
- Create: `docs/src/deriving/ping_pong.md`
- Modify: `docs/make.jl` (the `"Ping-pong mechanisms" => …` line)

- [ ] **Step 1: Write `docs/src/deriving/ping_pong.md`**
  Author the prose with the `elements-of-style:writing-clearly-and-concisely` skill. The page explains how the package represents a covalent (ping-pong) intermediate: an enzyme form on conformation `:E` that carries a `Residual` of atoms left behind between the producing and consuming steps.

  Code-verified facts to cover (cite each):
  - A `Residual` records atoms added/subtracted relative to apo enzyme: `Residual` has `added::Vector{Substrate}` and `subtracted::Vector{Product}` (`src/types.jl:32-39`). `has_residual(s::Species)` is true when the residual is non-empty (`src/types.jl:67`). The form's rendered name encodes it, e.g. `:E_res_+A_-Q` (`name(::Species)`, `src/types.jl:80-97`).
  - In the **enumerator**, the covalent intermediate always lives on conformation `:E` carrying a `Residual` — never a separate conformation label. The enumerator builds these forms via `_make_species` (always `:E`, `src/mechanism_enumeration.jl:13-21`) and computes the residual from the consumed-substrate / released-product history with `_residual_for` (`src/mechanism_enumeration.jl:29-45`). State the scope precisely: this `:E`-only invariant is the **enumerator's** convention (corrects the CLAUDE.md over-broad claim, staleness ledger section 8).
  - The DSL, by contrast, **does** allow a separate conformation label for a covalent intermediate (e.g. `Estar(; residual = A - P)`) — see `_call_form_term_info` and the residual-walk in `src/dsl.jl:301-355` and `_walk_residual_expr` (`:427-454`). So when documenting hand-written mechanisms, note the DSL is more permissive than the enumerator's `:E`-only output.
  - A genuine ping-pong intermediate must carry a **non-empty** residual (atoms remain on the enzyme between the producing and consuming steps); the enumerator rejects a degenerate empty-residue "ping-pong" that would return the enzyme to apo `E` mid-cycle (CLAUDE.md "Admissible residual" C4; `backtrack!`'s `pingpong_intermediate` parameter, `src/mechanism_enumeration.jl:317-331`). Note the staleness-ledger correction: the backtracking flag is `pingpong_intermediate::Bool`, not `has_residual` (which is the `Species` accessor).
  - Atom inventories come from `@enzyme_reaction` atom brackets (`S[C2N1]`), which are load-bearing for ping-pong: the residual is computed by atom bookkeeping, so substrates and products must declare real atom counts (`_parse_chemical_formula`, `src/dsl.jl:129-141`; atom counts must be positive integers, `ReactantAtoms` validation, `src/types.jl:236`).

  Include one `@example` block (run, not output-checked) that builds a real bi-bi atom-transfer reaction, enumerates its minimal mechanisms, and prints the names of the ping-pong (residual-bearing) enzyme forms to show they sit on conformation `:E`:

  ````markdown
  ```@example pingpong
  using EnzymeRates
  rxn = @enzyme_reaction begin
      substrates: A[C2N1], B[C1]
      products:   P[C2], Q[C1N1]
  end
  mechs = EnzymeRates.init_mechanisms(rxn)
  pp = nothing
  for m in mechs, grp in EnzymeRates.steps(m), s in grp
      for sp in (EnzymeRates.from_species(s), EnzymeRates.to_species(s))
          EnzymeRates.has_residual(sp) && (global pp = m)
      end
  end
  residual_forms = Symbol[]
  for grp in EnzymeRates.steps(pp), s in grp
      for sp in (EnzymeRates.from_species(s), EnzymeRates.to_species(s))
          EnzymeRates.has_residual(sp) &&
              (EnzymeRates.name(sp) in residual_forms ||
               push!(residual_forms, EnzymeRates.name(sp)))
      end
  end
  residual_forms
  ```
  ````

  In prose, state the verified result: the residual-bearing forms are named like `:E_res_+A_-Q`, `:EB_res_+A_-Q`, `:EQ_res_+A_-Q` — all on conformation `:E`, each carrying the `+A -Q` residual (the group from A retained until B picks it up to form Q). Note that `init_mechanisms` and `steps`/`from_species`/`has_residual`/`name` are reached as `EnzymeRates.<name>` because they are internal-but-usable (CLAUDE.md API design).

  Cross-link `[@enzyme_reaction](@ref)` and (in prose, by name) the enumeration page for how these forms are generated.

- [ ] **Step 2: Repoint the nav entry**
  Set `"Ping-pong mechanisms" => "deriving/ping_pong.md",` in `docs/make.jl`.

- [ ] **Step 3: Build the docs**
  ```bash
  julia --project=docs /home/denis.linux/.julia/dev/EnzymeRates/docs/make.jl
  ```
  Expected: ends with `Documenter: rendering done`; the `@example pingpong` block runs and lists residual-bearing `:E_…_res_…` form names; no warnings.

- [ ] **Step 4: Commit**
  ```bash
  git add docs/src/deriving/ping_pong.md docs/make.jl
  git commit -m "Add 'Ping-pong mechanisms' derivation page

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
  ```

---

### Task 2.6: Dead-end inhibitor binding

**Files:**
- Create: `docs/src/deriving/dead_end.md`
- Modify: `docs/make.jl` (the `"Dead-end inhibitor binding" => …` line)

- [ ] **Step 1: Write `docs/src/deriving/dead_end.md`**
  Author the prose with the `elements-of-style:writing-clearly-and-concisely` skill. The page explains how a competitive (dead-end) inhibitor enters the rate equation: it binds an enzyme form in rapid equilibrium and adds a term to the shared denominator without contributing to the numerator (no turnover).

  Code-verified facts to cover (cite each):
  - A dead-end inhibitor binds in **rapid equilibrium only**: the enumerator creates the inhibitor binding step with `is_equilibrium = true` (`Step(base, de_species, CompetitiveInhibitor(reg_name), true)`, `src/mechanism_enumeration.jl:1507-1509`). Because it is a dead end, the inhibitor-bound form has no catalytic exit, so it appears only in the denominator.
  - Inhibitors are declared via `dead_end_inhibitors:` or `competitive_inhibitors:` — **both** map to role `:competitive`; the only two regulator roles are `:competitive` and `:allosteric` (`_parse_reaction_block`, `src/dsl.jl:74-78`; corrects the CLAUDE.md `:unknown`/`:dead_end` role labels, staleness ledger section 8).
  - In a hand-written `@enzyme_mechanism`, the `::Inh` role tag binds a declared metabolite in its `CompetitiveInhibitor` role while keeping its real name, so `concs.NAME` still drives it. The tagged form renders with an `inh` marker (`:E_Iinh`, parameter `:K_Iinh_E`) so it stays distinct from a product-bound form (`_step_side_term_info` / `_call_form_term_info` `::Inh` handling, `src/dsl.jl:277-294` and `:319-333`; `name(::Species)` `inh` suffix, `src/types.jl:83-84`). Only `::Inh` is supported.
  - When multiple regulators bind the same site, they share a single denominator factor of the form `(1 + R1/K_R1 + R2/K_R2)^m` (CLAUDE.md "Same-site regulators"). For a single competitive inhibitor on free enzyme, the effect is one added `+ I/K_Iinh_E` term in the denominator.

  Include one `@example` block (run, not output-checked) that adds a competitive inhibitor with the `::Inh` tag to the textbook mechanism and shows the rendered equation gaining an `I / K_Iinh_E` denominator term:

  ````markdown
  ```@example deadend
  using EnzymeRates
  de = @enzyme_mechanism begin
      substrates: S
      products:   P
      regulators: I
      steps: begin
          E + S ⇌ E(S)
          E(S) <--> E(P)
          E(P) ⇌ E + P
          E + I ⇌ E(I::Inh)
      end
  end
  print(rate_equation_string(de))
  ```
  ````

  In prose, point at the verified denominator `1 + I / K_Iinh_E + P / K_P_E + S / K_S_E`: the numerator is unchanged from the inhibitor-free textbook case (no turnover term for `I`), and the only new fitted parameter is `K_Iinh_E` (confirm with `parameters(de)` in an inline sentence: `(:K_Iinh_E, :K_P_E, :K_S_E, :k_ES_to_EP, :Keq, :E_total)`). Note that `I` appears in `concs` because the inhibitor keeps its declared name.

  Cross-link `[@enzyme_mechanism](@ref)`, `[rate_equation_string](@ref)`, and (in prose) the enumeration page's "add dead-end regulator" move.

- [ ] **Step 2: Repoint the nav entry**
  Set `"Dead-end inhibitor binding" => "deriving/dead_end.md",` in `docs/make.jl`.

- [ ] **Step 3: Build the docs**
  ```bash
  julia --project=docs /home/denis.linux/.julia/dev/EnzymeRates/docs/make.jl
  ```
  Expected: ends with `Documenter: rendering done`; the `@example deadend` block runs and renders the `I / K_Iinh_E` denominator term; no warnings.

- [ ] **Step 4: Commit**
  ```bash
  git add docs/src/deriving/dead_end.md docs/make.jl
  git commit -m "Add 'Dead-end inhibitor binding' derivation page

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
  ```

---

### Task 2.7: MWC allostery

**Files:**
- Create: `docs/src/deriving/mwc_allostery.md`
- Modify: `docs/make.jl` (the `"MWC allostery" => …` line)

- [ ] **Step 1: Write `docs/src/deriving/mwc_allostery.md`**
  Author the prose with the `elements-of-style:writing-clearly-and-concisely` skill. The page explains how the package builds a Monod–Wyman–Changeux (MWC) two-conformation rate equation: an active (A) and inactive (I) state, each a copy of the catalytic mechanism, summed with a coupling constant `L`, with regulatory ligands binding in rapid equilibrium and shifting the A/I balance. **Standardize all terminology on A/I**; state once that R ≡ A (active) and T ≡ I (inactive), so the package's `K_A_…` / `K_I_…` names match the literature R/T notation.

  Code-verified facts to cover (cite each):
  - Each kinetic group (and each regulatory ligand) carries an A/I taxonomy tag, one of `:OnlyA`, `:OnlyI`, `:EqualAI`, `:NonequalAI` (`src/types.jl:506-507`):
    - `:OnlyA` — exists in the active state only; zeroed in the inactive polynomial.
    - `:OnlyI` — exists in the inactive state only.
    - `:EqualAI` — one shared symbol for both states (renders with no A/I token, e.g. `K_S_E`).
    - `:NonequalAI` — independent A and I symbols (`K_A_S_E` / `K_I_S_E`, `k_A_ES_to_EP` / `k_I_ES_to_EP`).
  - **Catalytic groups cannot be `:OnlyI`** — it is a hard constructor error (active-state convention): `_VALID_CAT_ALLO_STATES = (:OnlyA, :EqualAI, :NonequalAI)` and the `AllostericMechanism` constructor raises on `:OnlyI` for a catalytic group (`src/types.jl:506`, `:528-534`). Regulatory ligands *may* be `:OnlyI` (`_VALID_REG_ALLO_STATES`, `src/types.jl:507`).
  - The MWC rate is a partition-function sum over conformations: numerator `= cat_n * Σ_c (L_c * N_cat_c * Q_cat_c^(cat_n-1) * Π reg-site factors)`, denominator `= Σ_c (L_c * Q_cat_c^cat_n * Π reg-site factors)`, where `cat_n = catalytic_multiplicity` (`_allosteric_num_den_exprs`, `src/rate_eq_derivation.jl:1521-1597`; the formula comment at `:957-969`). Each catalytic partition factor appears raised to `catalytic_multiplicity` in the denominator and `catalytic_multiplicity - 1` in the numerator (the `^2` powers in the example below come from `catalytic_multiplicity: 2`).
  - Regulatory ligands bind in **rapid equilibrium**; each regulatory site contributes a factor `(1 + lig/K_lig + …)^multiplicity` to both numerator and denominator (`_reg_site_expr`, `src/rate_eq_derivation.jl:1372-1388`; `_power_expr` at site multiplicity, `:1562-1582`).
  - The A↔I conformational transition is itself rapid-equilibrium, encoded by the single coupling parameter `L` (`Lallo`, included in the Reduced parameter set via `_dependent_param_exprs(::AllostericEnzymeMechanism)`, `src/rate_eq_derivation.jl:1262-1361`).
  - When the inactive cycle cannot close (any `:OnlyA` catalytic group present), the inactive numerator is forced to zero for Haldane consistency — `_i_state_dead` returns true and the `L*num_I` term is dropped (`_i_state_dead`, `src/rate_eq_derivation.jl:981-984`; numerator branch, `:1588-1596`).
  - **Known rendering quirk** (spec section 10, "do not fix"): in the rendered Haldane constraint lines, the RHS may reference a bare active-state base name (e.g. `k_ES_to_EP`, `k_A_EP_to_ES`) even though `params` lists the A/I-suffixed names. This is confirmed real behavior; document it as-is, do not "correct" it.

  Include one `@example` block (run, not output-checked) building a two-conformation MWC mechanism (`catalytic_multiplicity: 2`, an `:OnlyA` activator `A` and an `:OnlyI` inhibitor `I`, with one `:NonequalAI` catalytic step) and printing both `parameters(...)` and the rendered string:

  ````markdown
  ```@example mwc
  using EnzymeRates
  allo = @allosteric_mechanism begin
      substrates: S
      products:   P
      catalytic_multiplicity: 2
      allosteric_regulators: A::OnlyA, I::OnlyI
      catalytic_steps: begin
          E + S ⇌ E(S)        :: EqualAI
          E(S) <--> E(P)       :: NonequalAI
          E(P) ⇌ E + P        :: EqualAI
      end
      regulatory_site(multiplicity = 2): begin
          ligands: A
      end
      regulatory_site(multiplicity = 2): begin
          ligands: I
      end
  end
  parameters(allo)
  ```

  ```@example mwc
  print(rate_equation_string(allo))
  ```
  ````

  In prose, walk the verified parameter tuple `(:K_P_E, :K_S_E, :k_A_ES_to_EP, :k_I_ES_to_EP, :K_A_Areg, :K_I_Ireg, :L, :Keq, :E_total)`:
  - `:EqualAI` catalytic steps share `K_S_E` / `K_P_E` (no A/I token).
  - the `:NonequalAI` catalytic step splits into `k_A_ES_to_EP` and `k_I_ES_to_EP`.
  - the `:OnlyA` activator contributes `K_A_Areg` (active-state regulatory K); the `:OnlyI` inhibitor contributes `K_I_Ireg`.
  - `L` is the MWC coupling constant.
  Point at the `(1 + A / K_A_Areg) ^ 2` and `(1 + I / K_I_Ireg) ^ 2` factors and the `L * (...)` inactive branch to make the partition-function structure concrete, and flag the bare-`k_ES_to_EP` Haldane-RHS line as the known rendering quirk.

  Add the hard-error demonstration as an `@example` so the page documents the constraint (the `@allosteric_mechanism` macro raises at expansion, so wrap it to display the message rather than failing the build):

  ````markdown
  ```@example mwc
  try
      @allosteric_mechanism begin
          substrates: S
          products:   P
          catalytic_steps: begin
              E + S ⇌ E(S)   :: OnlyI
              E(S) <--> E(P)  :: EqualAI
              E(P) ⇌ E + P   :: EqualAI
          end
      end
  catch err
      showerror(stdout, err)
  end
  ```
  ````

  In prose, state that a catalytic step tagged `:OnlyI` is rejected: the active state must always be present, so only regulatory ligands may be inactive-only.

  Cross-link `[@allosteric_mechanism](@ref)`, `[parameters](@ref)`, and `[rate_equation_string](@ref)`. Add the MWC citation `[monodNatureAllostericTransitions1965](@cite)` once in the intro (use the key present in `docs/src/refs.bib`).

- [ ] **Step 2: Repoint the nav entry**
  Set `"MWC allostery" => "deriving/mwc_allostery.md",` in `docs/make.jl`.

- [ ] **Step 3: Build the docs**
  ```bash
  julia --project=docs /home/denis.linux/.julia/dev/EnzymeRates/docs/make.jl
  ```
  Expected: ends with `Documenter: rendering done`; all three `@example mwc` blocks run (the third prints the `:OnlyI` rejection message via `showerror`); the `@cite` resolves; no warnings.

- [ ] **Step 4: Commit**
  ```bash
  git add docs/src/deriving/mwc_allostery.md docs/make.jl
  git commit -m "Add 'MWC allostery' derivation page

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
  ```

I have all the verified facts. `popsize`/`maxtime`/`n_restarts` are optimizer kwargs passed through `kwargs...`; `popsize`/`maxiters`/`verbose` are PyCMA/BBO-specific solver options. Now I'll write the Phase 3 section.

## Phase 3 — Fitting pillar

This phase writes the three "Fitting rate equations" pages from the site outline: the fitting tutorial and data format, the normalized-vs-absolute rate distinction, and the loss-and-optimizers concept page. All three pages live under `docs/src/fitting/` (created with stub content in Phase 0) and are wired into `docs/make.jl`'s `pages` tree under the "Fitting rate equations" node. The executor writes the prose with the `elements-of-style:writing-clearly-and-concisely` skill, standardizes all MWC terminology on A/I (active/inactive, stating the R≡A / T≡I correspondence once where allostery appears), and presents the relative path as primary with the absolute path as advanced. Every block that shows a fit result is an `@example` (the multi-start optimizer is non-deterministic), so the prose asserts only stable keys and shapes — `keys(result)`, `result.retcode isa Symbol` — never pinned numbers. The one deterministic, output-checked block in the phase is the `FittingProblem` default-field `jldoctest`. Each task ends green: `julia --project=docs docs/make.jl` builds with no doctest failures and prints `Documenter: rendering done`.

### Task 3.1: Fitting tutorial & data format page

**Files:**
- Create: `/home/denis.linux/.julia/dev/EnzymeRates/docs/src/fitting/tutorial.md`
- Modify: `/home/denis.linux/.julia/dev/EnzymeRates/docs/make.jl` (the `pages` vector — add the `"Fitting rate equations"` node with this page)

- [ ] **Step 1: Confirm the page is referenced in `make.jl`**
Open `/home/denis.linux/.julia/dev/EnzymeRates/docs/make.jl` and locate the `pages = [ ... ]` vector. Ensure a `"Fitting rate equations"` node exists pointing at the three Phase-3 files. The node must read exactly (insert it after the "Deriving rate equations" node and before the "Identifying the best rate equation" node; if the three child entries already exist as Phase-0 stubs, leave them):

```julia
    "Fitting rate equations" => [
        "Fitting tutorial & data format" => "fitting/tutorial.md",
        "Normalized vs absolute rate" => "fitting/normalized_vs_absolute.md",
        "Loss & optimizers" => "fitting/loss_and_optimizers.md",
    ],
```

- [ ] **Step 2: Write `docs/src/fitting/tutorial.md`**
Write the page prose with the `elements-of-style:writing-clearly-and-concisely` skill. The page is the entry point to the fitting pillar — "do it, then understand it." It MUST cover these code-verified facts (cite nothing in the rendered prose, but the facts below are the ground truth to write from):

- **Data format**: the table is any Tables.jl-compatible source (a `NamedTuple` of equal-length column vectors is the simplest). It is converted internally with `Tables.columntable` (`FittingProblem`, `src/fitting.jl:50`). Required columns: `group` (groups measurements that share one `E_total`), `Rate` (measured rate, must be nonzero — a zero rate errors because the loss works in log space, `src/fitting.jl:64-68`), and **one column per metabolite** whose name matches `metabolites(mechanism)` exactly (`src/fitting.jl:52,59-61`). A missing required column, a missing metabolite column, or a zero rate each raises an `ErrorException` at construction (`src/fitting.jl:56-68`).
- **`metabolites(mechanism)`** is the exported accessor that tells you which concentration columns the data needs; for a uni-uni `S ⇌ P` mechanism it returns `(:S, :P)` (`metabolites`, exported at `src/EnzymeRates.jl:17`; behavior shown in `test/test_fitting.jl:44`).
- **Construction**: `FittingProblem(mechanism, data; Keq, scale_k_to_kcat=1.0)` (`src/fitting.jl:46`). `Keq` is a required keyword and is always user-provided, never estimated. The same constructor accepts either a compiled `EnzymeMechanism`/`AllostericEnzymeMechanism` or a concrete `Mechanism`/`AllostericMechanism` (the latter is compiled once at construction, `src/fitting.jl:99-102`).
- **Fitting**: `fit_rate_equation(fp, optimizer; n_restarts=20, maxtime=60.0, lb=…, ub=…, kwargs...)` (`src/fitting.jl:209`). It runs `n_restarts` independent optimizations from random initial points and keeps the best (`src/fitting.jl:227-236`). Bounds `lb`/`ub` default to `fill(∓15.0, n_params)` in log space — parameters are fit in log space, `actual k = exp(x)` (`src/fitting.jl:130,239`). Extra `kwargs...` pass straight through to the optimizer (this is how `popsize`, `maxiters`, etc. reach the solver).
- **Return value**: a `NamedTuple` `(params, loss, retcode)` (`src/fitting.jl:201,249`). `params` is a `NamedTuple` of fitted rate constants keyed by `fitted_params(mechanism)`; `loss` is the best loss; `retcode` is a `Symbol`. **Only `:Success` means the optimizer converged on its own criteria** — any other value (`:MaxTime`, `:Failure`, the `:NoFiniteLoss` sentinel) means treat the fit as un-converged (`src/fitting.jl:202-207,225,234`).
- **The optimizer is an explicit argument** with no default — point the reader to the "Loss & optimizers" page for the bring-your-own-optimizer story and a working choice. The tutorial uses `BBO_adaptive_de_rand_1_bin_radiuslimited()` from `OptimizationBBO` because it is the alternative exercised in the test suite (`test/test_fitting.jl:351,364`) and needs no Python.

The page MUST contain these exact runnable blocks (an `@example fitting` shared sandbox so state carries across blocks):

First, define a small mechanism and inspect the required columns:

````markdown
```@example fitting
using EnzymeRates

uni_uni = @enzyme_mechanism begin
    substrates: S
    products:   P
    steps: begin
        E + S <--> E(S)
        E(S) <--> E + P
    end
end

metabolites(uni_uni)
```
````

Next, build synthetic data by evaluating the rate equation on a grid, then construct the `FittingProblem`:

````markdown
```@example fitting
true_params = (kon_S_E = 10.0, kon_P_ES = 5.0, koff_S_E = 1.0, Keq = 2.0, E_total = 1.0)

concs_list = [(S = s, P = p) for s in (0.5, 1.0, 2.0, 5.0, 10.0) for p in (0.1, 0.5)]

data = (
    group = fill("G1", length(concs_list)),
    Rate  = [rate_equation(uni_uni, c, true_params) for c in concs_list],
    S     = [c.S for c in concs_list],
    P     = [c.P for c in concs_list],
)

fp = FittingProblem(uni_uni, data; Keq = 2.0)
```
````

Then run the fit and show only stable shape/keys (never pinned numbers):

````markdown
```@example fitting
using OptimizationBBO

result = fit_rate_equation(
    fp, BBO_adaptive_de_rand_1_bin_radiuslimited();
    n_restarts = 5, maxtime = 5.0,
)

keys(result)
```
````

````markdown
```@example fitting
result.retcode isa Symbol
```
````

````markdown
```@example fitting
result.params
```
````

The prose immediately after these blocks MUST state: the fit is a random multi-start, so `result.params` and `result.loss` vary run to run; check `result.retcode === :Success` to confirm convergence; and the returned `params` carry the kcat normalization explained on the next page. Add a sentence cross-linking the concept pages with Documenter `@ref` links: "[Normalized vs absolute rate](@ref) explains `scale_k_to_kcat`; [Loss & optimizers](@ref) covers the loss and how to bring your own optimizer." (Use the page titles as the `@ref` targets — Documenter resolves a header `@ref` by the page's H1 title.)

- [ ] **Step 3: Run the build to verify the page compiles and its `@example` blocks run**
```bash
julia --project=docs /home/denis.linux/.julia/dev/EnzymeRates/docs/make.jl
```
Expected: the build runs all `@example fitting` blocks without error, prints `Documenter: rendering done`, and reports no doctest failures. If a metabolite column name mismatch or a `@ref` to a not-yet-built page warns, fix the column names / defer the cross-link until the target pages exist (they are added in the same phase, so re-running after Tasks 3.2–3.3 resolves any `@ref` warning).

- [ ] **Step 4: Commit**
```bash
git add /home/denis.linux/.julia/dev/EnzymeRates/docs/src/fitting/tutorial.md /home/denis.linux/.julia/dev/EnzymeRates/docs/make.jl
git commit -m "Docs: fitting tutorial & data-format page

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

### Task 3.2: Normalized vs absolute rate page

**Files:**
- Create: `/home/denis.linux/.julia/dev/EnzymeRates/docs/src/fitting/normalized_vs_absolute.md`

- [ ] **Step 1: Write `docs/src/fitting/normalized_vs_absolute.md`**
Write the prose with the `elements-of-style:writing-clearly-and-concisely` skill. Present the **relative (normalized) path as primary** and the **absolute turnover path as advanced**. The page MUST cover these code-verified facts:

- **The single knob is `scale_k_to_kcat::Union{Real,Nothing}`**, a field on `FittingProblem` and a keyword on its constructor (`src/fitting.jl:25,32,47`). Its default is `1.0`. A positive `Real` selects **relative mode**; `nothing` selects **absolute mode**. A non-positive `Real` is rejected at construction with an `ErrorException` (`src/fitting.jl:48-49`).
- **Relative mode (a `Real`, default `1.0`)**: kinetic data is treated as relative. The loss is **per-group mean-centered**, which removes each group's arbitrary `E_total` scale, so the loss is invariant to per-group `E_total` rescaling (`loss!`, `src/fitting.jl:147-169`, esp. the centered branch `156-168`; field doc `src/fitting.jl:15-16`). After fitting, the returned SS rate constants are **rescaled so `_kcat_forward(mechanism, params) ≈ scale_k_to_kcat`** (`fit_rate_equation`, `src/fitting.jl:240-248`). With the default `1.0`, the fitted equation has kcat = 1; the absolute turnover is recovered by multiplying by a separately measured kcat. A custom target (e.g. `scale_k_to_kcat = 42.0`) rescales so kcat = 42 (`test/test_fitting.jl:371-375`).
- **Absolute mode (`nothing`)**: the data is absolute per-enzyme turnover. The loss is **uncentered** — the absolute magnitude is scored (`loss!`, `src/fitting.jl:152-155`). `fit_rate_equation` returns the **raw, unrescaled** parameters; the data fixes the absolute scale (`src/fitting.jl:198-199,240`).
- **Which rescales**: rescaling touches **only SS rate constants** (the `Kon`, `Koff`, `Kfor`, `Krev` family, classified by `_ss_rate_constant_names`, `src/rate_eq_derivation.jl:621-642`). RE binding K's, `Keq`, `E_total`, allosteric `L`, and regulatory K's are left unchanged (`rescale_parameter_values`, `src/rate_eq_derivation.jl:945-955`). kcat is homogeneous degree-1 in the SS k's and independent of the RE K's, which is why scaling the k's uniformly sets kcat to the target without disturbing anything else.
- **`rescale_parameter_values` is the public kcat helper**: `rescale_parameter_values(m, params; scale_k_to_kcat=1.0)` (exported, `src/EnzymeRates.jl:18`; def `src/rate_eq_derivation.jl:945`). It is what `fit_rate_equation` calls internally in relative mode, and it can be called directly on any parameter `NamedTuple` to renormalize after the fact. `_kcat_forward` and the other kcat internals stay private — name them only as "the internal kcat computation," do not present them as API.
- Because allostery first appears here only in passing (the rescale also covers allosteric SS k's via the `AllostericEnzymeMechanism` overload, `src/rate_eq_derivation.jl:627-642`), state the R≡A / T≡I correspondence once if A/I parameter names appear, and otherwise keep the page's worked example uni-uni.

The page MUST contain a `jldoctest` that pins the default field value (this is deterministic and output-checked). Write the **input** block exactly as below; the **output** is captured by the doctest-fix command in Step 2, not hand-written:

````markdown
```jldoctest
julia> using EnzymeRates

julia> m = @enzyme_mechanism begin
           substrates: S
           products:   P
           steps: begin
               E + S <--> E(S)
               E(S) <--> E + P
           end
       end;

julia> data = (group = ["G1"], Rate = [1.0], S = [1.0], P = [0.1]);

julia> FittingProblem(m, data; Keq = 2.0).scale_k_to_kcat === 1.0
true
```
````

The page MUST also contain a deterministic `@example` (no fitting, so it could even be a doctest, but keep it `@example` to avoid pinning a float) showing `rescale_parameter_values` setting kcat to a target directly on a hand-built parameter set:

````markdown
```@example rescale
using EnzymeRates

uni_uni = @enzyme_mechanism begin
    substrates: S
    products:   P
    steps: begin
        E + S <--> E(S)
        E(S) <--> E + P
    end
end

params = (kon_S_E = 10.0, kon_P_ES = 5.0, koff_S_E = 1.0, Keq = 2.0, E_total = 1.0)

rescaled = rescale_parameter_values(uni_uni, params; scale_k_to_kcat = 5.0)
rescaled
```
````

The prose after this block MUST note that only the SS k columns (`kon_S_E`, `kon_P_ES`, `koff_S_E` here) changed and that `Keq` and `E_total` are untouched, and that calling the internal kcat computation on `rescaled` returns ≈ 5.

- [ ] **Step 2: Capture the real `jldoctest` output with doctest-fix**
```bash
julia --project=docs -e 'using Documenter, EnzymeRates; doctest(EnzymeRates; fix=true)'
```
Expected: the command runs the docstring/markdown doctests, writes the real captured `true` output into the `scale_k_to_kcat === 1.0` block, and reports `Doctests passed` (or the fixed count). Open `docs/src/fitting/normalized_vs_absolute.md` and confirm the doctest output line is now present and equals `true`.

- [ ] **Step 3: Run the full build to verify**
```bash
julia --project=docs /home/denis.linux/.julia/dev/EnzymeRates/docs/make.jl
```
Expected: prints `Documenter: rendering done`; the `jldoctest` passes (no doctest failures); the `@example rescale` block runs without error.

- [ ] **Step 4: Commit**
```bash
git add /home/denis.linux/.julia/dev/EnzymeRates/docs/src/fitting/normalized_vs_absolute.md
git commit -m "Docs: normalized vs absolute rate page (scale_k_to_kcat)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

### Task 3.3: Loss & optimizers page

**Files:**
- Create: `/home/denis.linux/.julia/dev/EnzymeRates/docs/src/fitting/loss_and_optimizers.md`

- [ ] **Step 1: Write `docs/src/fitting/loss_and_optimizers.md`**
Write the prose with the `elements-of-style:writing-clearly-and-concisely` skill. This is the concept page behind the fit. It MUST cover these code-verified facts:

- **The loss is a log-ratio loss** computed by `loss!(x, fp)` (`src/fitting.jl:121`). Parameters in `x` are in **log space** (`actual k = exp(x)`, `src/fitting.jl:130`). For each data point it forms `log(abs(pred)) − log(abs(measured))` (`src/fitting.jl:143`); the final value is the sum of squares divided by the number of points (`src/fitting.jl:181`).
- **Centered vs uncentered** (this is the loss-side of `scale_k_to_kcat`, cross-link to [Normalized vs absolute rate](@ref)): when `fp.scale_k_to_kcat` is a `Real`, each **group's** log-ratios are mean-subtracted before squaring (`src/fitting.jl:156-168`), making the loss invariant to per-group `E_total`. When it is `nothing`, the log-ratios are squared uncentered (`src/fitting.jl:152-155`), so absolute magnitude is scored.
- **Sign-mismatch penalty**: a point whose predicted rate sign disagrees with the measured sign (or where the prediction is zero) is assigned a sentinel `10.0` in the working buffer (`src/fitting.jl:140-141`); after the loop a flat `100.0` per mismatched point is added (`src/fitting.jl:176-181`). In centered mode this keeps an all-mismatch group from canceling to zero loss under mean-subtraction (`src/fitting.jl:171-175`).
- **`loss!` is zero-allocation on the hot path** and is evaluated millions of times per cross-validation fold — explain this is why parameters are pre-arranged into `fitted_params` order and the buffer is pre-allocated (`FittingProblem` field `log_ratios_buffer`, `src/fitting.jl:18,87`; `loss!` reuses `fp.log_ratios_buffer`, `src/fitting.jl:122`). Tie this to the package's `rate_equation` 0-allocation / sub-100 ns contract, which the loss depends on (cross-link to the Developer page once it exists).
- **Bring-your-own optimizer (state this fact explicitly and prominently)**: the base package depends only on `Optimization` — `fit_rate_equation` builds an `Optimization.OptimizationFunction` / `OptimizationProblem` and calls `Optimization.solve(prob, optimizer; maxtime, kwargs...)` (`src/fitting.jl:216,229-230`). The `optimizer` argument is **whatever Optimization.jl solver object the caller passes**; no solver backend ships with the base package. The reader installs one of the Optimization.jl solver sub-packages and passes its optimizer object.
- **Blessed choice**: `PyCMAOpt()` from `OptimizationPyCMA` (a multi-start CMA-ES; the optimizer recommended for rate-equation fitting and the one the identify pipeline uses, `test/test_identify_rate_equation.jl:9,81,388`). **Tested alternative**: `BBO_adaptive_de_rand_1_bin_radiuslimited()` from `OptimizationBBO` (`test/test_fitting.jl:351,364`). Both `OptimizationPyCMA` and `OptimizationBBO` are **add-on dependencies** the user installs separately; `OptimizationPyCMA` additionally needs Python/`pycma` available.
- **Passing solver options**: solver-specific keywords flow through `fit_rate_equation`'s `kwargs...` to `Optimization.solve` (`src/fitting.jl:214,230`). For PyCMA, `popsize` and `maxiters` are common (`test/test_identify_rate_equation.jl:390`); `maxtime` and `n_restarts` are arguments of `fit_rate_equation` itself (`src/fitting.jl:210-211`).

The page MUST contain this exact `@example` showing the same fit with two different optimizers, asserting only stable keys/retcode shape (no pinned numbers). It reuses the BBO solver (no Python needed in CI) and shows the PyCMA call as text the reader can run locally:

````markdown
```@example optimizers
using EnzymeRates, OptimizationBBO

uni_uni = @enzyme_mechanism begin
    substrates: S
    products:   P
    steps: begin
        E + S <--> E(S)
        E(S) <--> E + P
    end
end

true_params = (kon_S_E = 10.0, kon_P_ES = 5.0, koff_S_E = 1.0, Keq = 2.0, E_total = 1.0)
concs_list = [(S = s, P = p) for s in (0.5, 1.0, 2.0, 5.0, 10.0) for p in (0.1, 0.5)]
data = (
    group = fill("G1", length(concs_list)),
    Rate  = [rate_equation(uni_uni, c, true_params) for c in concs_list],
    S     = [c.S for c in concs_list],
    P     = [c.P for c in concs_list],
)

fp = FittingProblem(uni_uni, data; Keq = 2.0)

result = fit_rate_equation(
    fp, BBO_adaptive_de_rand_1_bin_radiuslimited();
    n_restarts = 5, maxtime = 5.0,
)

(keys(result), result.retcode isa Symbol)
```
````

The page MUST also include this **non-executed** block (a plain fenced ` ```julia ` block, not `@example`, because `OptimizationPyCMA` needs Python and must not run in CI) showing the blessed PyCMA call and its solver kwargs:

````markdown
```julia
using OptimizationPyCMA

result = fit_rate_equation(
    fp, PyCMAOpt();
    n_restarts = 5, maxtime = 5.0, popsize = 50,
)
```
````

The prose MUST explicitly state: "EnzymeRates depends only on Optimization.jl; it ships no solver. Install one of the Optimization.jl solver sub-packages — `OptimizationPyCMA` (blessed; needs Python and `pycma`) or `OptimizationBBO` (a tested, pure-Julia alternative) — and pass its optimizer object to `fit_rate_equation`." Add `@ref` cross-links to [Normalized vs absolute rate](@ref) (for centered/uncentered) and to the Fitting tutorial.

- [ ] **Step 2: Run the build to verify**
```bash
julia --project=docs /home/denis.linux/.julia/dev/EnzymeRates/docs/make.jl
```
Expected: prints `Documenter: rendering done`; the `@example optimizers` block runs (BBO is available in the docs project) and its final tuple shows the keys and `true`; the plain ` ```julia ` PyCMA block is rendered but not executed; no doctest failures. If Documenter warns that `OptimizationBBO` is not in the docs project, add it: `julia --project=docs -e 'using Pkg; Pkg.add("OptimizationBBO")'` (it should already be present from the Phase-0 docs `Project.toml`, which adds the example optimizer; if only `OptimizationPyCMA` was added, add `OptimizationBBO` too since CI cannot run the Python optimizer).

- [ ] **Step 3: Re-run the build after any cross-link target now exists**
```bash
julia --project=docs /home/denis.linux/.julia/dev/EnzymeRates/docs/make.jl
```
Expected: with all three Phase-3 pages present, every `@ref` cross-link between them resolves — no "unable to find reference" warnings — and the build prints `Documenter: rendering done` with no doctest failures.

- [ ] **Step 4: Commit**
```bash
git add /home/denis.linux/.julia/dev/EnzymeRates/docs/src/fitting/loss_and_optimizers.md
git commit -m "Docs: loss & optimizers page (bring-your-own optimizer)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

`frontier` is a search-internal `Dict`, not a struct accessor — `IdentifyRateEquationResults` only has `best` and `cv_results`. The spec's coverage-gap "(best, cv_results schema, frontier)" wording refers to documenting the search-internal frontier concept on the engine page, not a results accessor. I have all facts needed. Writing the phase section now.

## Phase 4 — Identify pillar

This phase writes the three pages of the "Identifying the best rate equation" pillar plus the coverage gaps that belong to them. Page 1 is a runnable tutorial built around a *fast* real search (noiseless simulated data, a non-degenerate uni-uni mechanism, a collapsed width-1 beam) that recovers its own generating mechanism in seconds, with prose noting that the full production search uses the wider default beam and runs ~1 hour. Page 2 documents LOOCV and the two-test model-selection rule (paired 1-SE AND one-sided permutation), grounded in a `jldoctest` on the pure `_onesided_permutation_p`, and corrects the README's plain-argmin-CV description. Page 3 documents the enumeration engine: `init_mechanisms` plus the six expansion moves and `_dedup_flat!`, with the "mostly +1 param, but verify per move" framing. Coverage gaps (`IdentifyRateEquationResults` accessors, CSV/`progress.log` artifacts, `max_param_count`, `FitFailure`/loud-failure semantics) fold onto the relevant page. Every page is written with the elements-of-style:writing-clearly-and-concisely skill, and all MWC terminology is standardized on A/I (state the R≡A, T≡I correspondence once where allostery is mentioned).

Prerequisite for this phase: Phase 0's `docs/Project.toml` already `dev`s EnzymeRates and adds Documenter + an optimizer dependency, and `docs/src/identify/` exists as a stub directory referenced from `docs/make.jl`'s `pages`. Each task below assumes the three `docs/src/identify/*.md` files are already listed under the "Identifying the best rate equation" section of `pages` in `make.jl` (Phase 0 created stubs); if a path is not yet wired, the executor adds it to the `pages` vector in the same task.

### Task 4.1: Identify tutorial — the fast search (FAST example)

**Files:**
- Create: `docs/src/identify/tutorial.md`
- Modify: `docs/make.jl` (the `pages` vector — only if `identify/tutorial.md` is not already listed)

- [ ] **Step 1: Verify the load-bearing facts against the code before writing**
  Open these and confirm each claim; do NOT proceed if any differs:
  - `src/identify_rate_equation.jl:305-328` (`_select_beam`): with `loss_rel_threshold=1.0`, `loss_abs_threshold=0.0`, the cutoff is `cutoff = loss_rel_threshold * best + loss_abs_threshold == best`, so only the lowest-loss candidate at each count clears it; `min_beam_width=1` keeps exactly the rank-1 candidate. Together they collapse the beam to one survivor per parameter-count level.
  - `src/identify_rate_equation.jl:32-77` (`IdentifyRateEquationProblem` constructor): requires columns `:group` and `:Rate` plus one column per substrate/product/regulator; requires `Keq` keyword; rejects any zero `Rate` (`log(0)`); requires **at least 2 unique groups** for CV (line 67-71). `scale_k_to_kcat` defaults to `1.0`.
  - `src/identify_rate_equation.jl:169-193` (`identify_rate_equation` kwargs): `optimizer` is required; defaults `min_beam_width=50`, `loss_rel_threshold=2.0`, `loss_abs_threshold=0.01`, `max_param_count=20`, `n_restarts=20`, `maxtime=60.0`, `save_dir=_default_save_dir()`, `pmap_function=pmap`, `show_progress=true`.
  - `src/identify_rate_equation.jl:199-207`: `save_dir` must be empty of `.csv` files or the run errors; `_default_save_dir()` (line 914-920) picks `YYYY_MM_DD_results[_N]` in the cwd.
  - `IdentifyRateEquationResults` (`src/identify_rate_equation.jl:90-93`) has exactly two fields: `best::AbstractEnzymeMechanism` and `cv_results::DataFrame`.
  - The uni-uni mechanism shape `E + S <--> E(S)` / `E(S) <--> E + P` (two SS steps) is a valid non-degenerate generator (`test/test_fitting.jl:8-15`); its `metabolites` are `(:S, :P)` and `make_synthetic_data` builds rates via `rate_equation(mechanism, concs, true_params)` (`test/test_fitting.jl:18-36`).

- [ ] **Step 2: Write the tutorial page**
  Create `docs/src/identify/tutorial.md`. The executor WRITES the prose with the elements-of-style:writing-clearly-and-concisely skill. The page MUST cover these code-verified facts (cite them inline where natural):
  - **What the search does** (`identify_rate_equation`, `src/identify_rate_equation.jl:169`): given an `EnzymeReaction` and rate data, it enumerates biochemically valid mechanisms, fits each, and selects the simplest that generalizes by leave-one-group-out cross-validation.
  - **The data contract** (`IdentifyRateEquationProblem` constructor, `src/identify_rate_equation.jl:32-77`): a columnar table with a `:group` column, a `:Rate` column, and one column per substrate/product/regulator; `Keq` is supplied, never fit; at least two distinct `group` values are required (they are the CV folds); no zero rates.
  - **Why this example is fast but real**: noiseless simulated data and a uni-uni mechanism with no degenerate analogs let the collapsed beam deterministically recover the generating mechanism; the beam arguments `loss_rel_threshold=1.0`, `loss_abs_threshold=0.0`, `min_beam_width=1` make the cutoff equal the best loss at each level (`_select_beam`, `src/identify_rate_equation.jl:317`), keeping one survivor per parameter count; a low `max_param_count` bounds the initial level. State plainly that the **full production search uses the wider defaults** (`min_beam_width=50`, `loss_rel_threshold=2.0`, `loss_abs_threshold=0.01`, `max_param_count=20`) and runs roughly an hour, exploring far more candidates.
  - **The result** (`IdentifyRateEquationResults`, `src/identify_rate_equation.jl:90-93`): `results.best` is the selected mechanism; `results.cv_results` is the CV-score `DataFrame`. Show reading both. Cross-link `rate_equation_string` with `[`rate_equation_string`](@ref)`.
  - A one-line forward pointer to the model-selection page and the enumeration-engine page.

  The page MUST contain these exact runnable blocks. The first is `@setup` (hidden, runs but not shown); the `@example identify_fast` blocks share a sandbox and are displayed. The optimizer is the docs-project dependency blessed in the spec (`OptimizationPyCMA`, `PyCMAOpt()` with `BBO_adaptive_de_rand_1_bin_radiuslimited`); `pmap_function=map` keeps the example serial.

  ````markdown
  # Identify tutorial

  ```@setup identify_fast
  using EnzymeRates
  ```

  ## A reaction and a generating mechanism

  ```@example identify_fast
  using EnzymeRates

  # The reaction whose mechanism we will recover: a reversible uni-uni S ⇌ P.
  rxn = @enzyme_reaction begin
      substrates: S
      products:   P
  end

  # A concrete, non-degenerate uni-uni mechanism to GENERATE the data.
  generator = @enzyme_mechanism begin
      substrates: S
      products:   P
      steps: begin
          E + S <--> E(S)
          E(S) <--> E + P
      end
  end
  nothing # hide
  ```

  ## Simulate noiseless data

  ```@example identify_fast
  Keq = 2.0
  true_params = (kon_S_E = 10.0, kon_P_ES = 5.0, koff_S_E = 1.0,
                 Keq = Keq, E_total = 1.0)

  # Two measurement groups => two cross-validation folds (the minimum).
  concs = [(S = 1.0, P = 0.1), (S = 2.0, P = 0.1), (S = 5.0, P = 0.1),
           (S = 1.0, P = 0.5), (S = 2.0, P = 0.5), (S = 5.0, P = 0.5)]
  groups = ["G1", "G1", "G1", "G2", "G2", "G2"]

  rates = [rate_equation(generator, c, true_params) for c in concs]
  data = (group = groups, Rate = rates,
          S = [c.S for c in concs], P = [c.P for c in concs])
  nothing # hide
  ```

  ## Run the fast search

  ```@example identify_fast
  using OptimizationPyCMA

  prob = IdentifyRateEquationProblem(rxn, data; Keq = Keq)

  results = identify_rate_equation(prob;
      optimizer = PyCMAOpt(BBO_adaptive_de_rand_1_bin_radiuslimited),
      loss_rel_threshold = 1.0,   # cutoff == best loss …
      loss_abs_threshold = 0.0,   # … no additive slack …
      min_beam_width = 1,         # … so exactly one survivor per level.
      max_param_count = 4,        # tiny initial level => seconds, not hours
      pmap_function = map,        # serial; pass `pmap` to distribute
      show_progress = false,
      save_dir = mktempdir(),
  )
  nothing # hide
  ```

  ## Read the result

  ```@example identify_fast
  results.best
  ```

  ```@example identify_fast
  print(rate_equation_string(results.best))
  ```

  ```@example identify_fast
  first(results.cv_results, 5)
  ```
  ````

  After the blocks, the prose explains: `results.best isa AbstractEnzymeMechanism`; `results.cv_results` schema is detailed on the [Model selection](@ref) page; the recovered equation matches the generator because the data is noiseless and the mechanism non-degenerate.

- [ ] **Step 3: Ensure the page is wired into `make.jl`**
  Open `docs/make.jl`. If `"identify/tutorial.md"` is not already an entry under the `"Identifying the best rate equation" =>` section of `pages`, add it as the first entry of that section so the nav reads:
  ```julia
  "Identifying the best rate equation" => [
      "Identify tutorial" => "identify/tutorial.md",
      "Model selection" => "identify/model_selection.md",
      "The enumeration engine" => "identify/enumeration_engine.md",
  ],
  ```
  (The latter two files are created in Tasks 4.2 and 4.3; listing them now is fine only if Phase 0 already created stubs — otherwise add them in their own tasks. If the section already lists all three, make no change.)

- [ ] **Step 4: Build the docs**
  ```bash
  julia --project=docs /home/denis.linux/.julia/dev/EnzymeRates/docs/make.jl
  ```
  Expected: build completes with `Documenter: rendering done`, no `@example` evaluation errors, no doctest failures. The `identify_fast` blocks run in seconds (the budget is the optimizer multi-start on a 3-parameter uni-uni fit at `max_param_count = 4`). If the build reports an `@example` error, read it — a missing optimizer dep means Phase 0 did not add `OptimizationPyCMA` to `docs/Project.toml`; fix that there, do not stub the example.

- [ ] **Step 5: Commit**
  ```bash
  git add docs/src/identify/tutorial.md docs/make.jl
  git commit -m "docs: identify tutorial — fast width-1-beam recovery example

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
  ```

### Task 4.2: Model selection — LOOCV, 1-SE rule, permutation test

**Files:**
- Create: `docs/src/identify/model_selection.md`
- Modify: `docs/make.jl` (the `pages` vector — only if `identify/model_selection.md` is not already listed)

- [ ] **Step 1: Verify the load-bearing facts against the code before writing**
  Confirm each; do NOT proceed if any differs:
  - **CV score** (`src/identify_rate_equation.jl:867`): `cv_df.cv_score = [mean(log.(v)) for v in cv_df.cv_fold_scores]` — the score is the **mean of the LOG** per-fold losses, not the mean of raw losses. (This is the README correction.)
  - **LOOCV folds** (`_loocv`, `src/identify_rate_equation.jl:621-654`): one fold per unique `group`; fit on all-but-one group, evaluate on the held-out group; each fold score floored at `eps(Float64)` so `log` is finite; **loud** — a fold fit that throws propagates, and a non-finite fold loss raises naming the held-out group (line 645-647).
  - **Selection rule** (`_select_best_n_params`, `src/identify_rate_equation.jl:754-819`): `n_min` is the bucket with lowest mean log-fold-loss (parsimony tiebreak to smaller `n_params`, line 768); for each smaller bucket in ascending `n_params`, accept iff **BOTH** `mean(diffs) ≤ se_threshold * std(diffs)/sqrt(n_folds)` (paired 1-SE, `se_threshold=1.0`) **AND** `permutation_p > perm_p_threshold` (`perm_p_threshold=0.16`); return the first that passes; fall through to `n_min` if none pass (line 802-812). When `n_folds==1` the SE is undefined and the loop is skipped (line 803).
  - **Permutation p** (`_onesided_permutation_p`, `src/identify_rate_equation.jl:673-716`): one-sided `Pr(perm_mean ≥ observed)` under the sign-flip null; **exact** enumeration of `2^n` sign patterns when `n ≤ exact_threshold` (default 20), Monte Carlo otherwise; pure and seedable (the `observed_sum` is computed with the same sequential reduction as the loop so the identity permutation is always counted, line 687-691).
  - **Within-bucket pick** (`src/identify_rate_equation.jl:875-882`): the best mechanism is the lowest **training `loss`** within the selected `best_n` bucket.
  - **Diagnostics columns** surfaced in `cv_results`: `mean_log_loss_diff`, `se_paired`, `permutation_p` (the `n_min` bucket has 0.0 in all three, `src/identify_rate_equation.jl:775-779`), plus `cv_score`, `n_params`, `eq_hash`, and one `cv_fold_<group>` column per group (`src/identify_rate_equation.jl:897-901`).
  - **`cv_results` base schema** (`_rows_to_dataframe`, `src/identify_rate_equation.jl:249-278`): columns `n_params`, `loss`, `mechanism_type`, `rate_equation`, `retcode`, `error`, `eq_hash`, plus one column per fitted parameter name (`missing` where a mechanism lacks that parameter).
  - **`n_cv_candidates`** (`src/identify_rate_equation.jl:184`, default 5): the top N **distinct-equation** candidates per param count enter LOOCV (`_cv_model_selection`, line 834-851 dedups by `eq_hash`).

- [ ] **Step 2: Write the model-selection page (prose + jldoctest skeleton)**
  Create `docs/src/identify/model_selection.md`. The executor WRITES the prose with the elements-of-style:writing-clearly-and-concisely skill and standardizes MWC terminology on A/I (only relevant if allostery is mentioned in passing). The page MUST cover, with citations:
  - **Leave-one-group-out CV**: each unique `group` value is one fold; the `group` column reflects experimental batches sharing an `E_total`, so LOOCV estimates generalization to new conditions (`_loocv`, `src/identify_rate_equation.jl:621`).
  - **The CV score is the mean of LOG per-fold losses** (`src/identify_rate_equation.jl:867`) — say this explicitly and correct the intuition that it is a plain mean of losses.
  - **Two-test selection rule**, AND-combined (`_select_best_n_params`, `src/identify_rate_equation.jl:754`): the paired 1-SE rule (`se_threshold=1.0`) and the one-sided sign-flip permutation test (`perm_p_threshold=0.16`); a simpler bucket is selected only if it passes **both**; smallest passing `n_params` wins; fall through to `n_min`. State that within the chosen bucket the lowest **training loss** mechanism is returned (`src/identify_rate_equation.jl:875-882`).
  - **Loud CV** (`src/identify_rate_equation.jl:645-647`): a non-finite fold loss aborts model selection rather than silently dropping the candidate, because a corrupted CV invalidates selection and CV is cheap to recompute from the saved CSVs.
  - **The `cv_results` DataFrame** (coverage gap — fold here): list its columns from `_rows_to_dataframe` (`src/identify_rate_equation.jl:259-277`) plus the diagnostic columns `mean_log_loss_diff`, `se_paired`, `permutation_p` and the per-fold `cv_fold_<group>` columns (`src/identify_rate_equation.jl:886-901`); note `IdentifyRateEquationResults` exposes exactly `best` and `cv_results` (`src/identify_rate_equation.jl:90-93`).
  - **`n_cv_candidates`** controls how many distinct equations per param count enter LOOCV (default 5, `src/identify_rate_equation.jl:184`).
  - A one-line note that `se_threshold` and `perm_p_threshold` are tunable knobs of `identify_rate_equation`.

  The page MUST contain a `jldoctest` block on the pure, seedable `_onesided_permutation_p`, demonstrating the exact sign-flip enumeration on a hand-built tiny input. Author it as an INPUT-ONLY block first (the real output is captured in Step 3 — do NOT invent it):

  ````markdown
  ## The permutation test, exactly

  For a small paired-difference vector the permutation p-value is computed by
  exact enumeration of all `2^n` sign-flips — pure and reproducible:

  ```jldoctest
  julia> using EnzymeRates

  julia> EnzymeRates._onesided_permutation_p([0.1, 0.2, 0.3, 0.4])
  ```
  ````

  Also include a brief note that `_onesided_permutation_p` and `_select_best_n_params` are internal helpers shown here only to make the rule concrete; users call `identify_rate_equation` and read `results.cv_results`.

- [ ] **Step 3: Capture the real doctest output**
  Run the doctest auto-fix to fill in the `_onesided_permutation_p` output from the actual function (4 elements → exact branch, `2^4 = 16` sign patterns):
  ```bash
  julia --project=docs -e 'using Documenter, EnzymeRates; doctest(EnzymeRates; fix=true)'
  ```
  Expected: the command rewrites the `jldoctest` block in `docs/src/identify/model_selection.md`, inserting the real returned `Float64` on the line after the `julia>` call. Inspect the diff with `git diff docs/src/identify/model_selection.md` and confirm a single float value was inserted (a deterministic count-of-sign-patterns-≥-observed divided by 16). Do not hand-edit the value.

- [ ] **Step 4: Correct the README's model-selection description (mark for Phase 6, verify wording here)**
  The README at `README.md:285-295` says "The CV score is the mean held-out loss" and "the parameter count whose CV score is lowest" — both are wrong (it is the mean of LOG losses, and selection is the smallest bucket passing the 1-SE AND permutation tests, not plain argmin). This page is the correct replacement. Add a sentence in the page noting it supersedes the README's description, so the Phase 6 trim can delete `README.md:285-295` and link here. Do not edit `README.md` in this phase (the trim is Phase 6).

- [ ] **Step 5: Ensure the page is wired into `make.jl`**
  Confirm `"Model selection" => "identify/model_selection.md"` is the second entry of the "Identifying the best rate equation" section in `docs/make.jl`'s `pages` (added in Task 4.1 Step 3). If absent, add it there.

- [ ] **Step 6: Build the docs**
  ```bash
  julia --project=docs /home/denis.linux/.julia/dev/EnzymeRates/docs/make.jl
  ```
  Expected: `Documenter: rendering done`, no doctest failures (the captured `_onesided_permutation_p` value now matches the live function), no `@example` errors.

- [ ] **Step 7: Commit**
  ```bash
  git add docs/src/identify/model_selection.md docs/make.jl
  git commit -m "docs: model selection — LOOCV, paired 1-SE + permutation test

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
  ```

### Task 4.3: The enumeration engine — init + six moves + dedup

**Files:**
- Create: `docs/src/identify/enumeration_engine.md`
- Modify: `docs/make.jl` (the `pages` vector — only if `identify/enumeration_engine.md` is not already listed)

- [ ] **Step 1: Verify the load-bearing facts against the code before writing**
  Confirm each; do NOT proceed if any differs:
  - **`init_mechanisms(reaction)`** (`src/mechanism_enumeration.jl:1809-1821`): returns `Vector{Mechanism}` at minimum parameter count — for each catalytic topology (`_catalytic_topologies`), one SS step, all substrate/product dead-end subsets (`_expand_substrate_product_dead_ends`), with binding steps sharing a `(metabolite, RE/SS)` class collapsed into one kinetic group (`_apply_equivalence_grouping`, line 1172-1191).
  - **`expand_mechanisms(mechs, rxn)`** (`src/mechanism_enumeration.jl:1756-1767`): returns a **flat** `Vector{Union{Mechanism, AllostericMechanism}}`; bucketing by param count is the caller's job; every child is asserted atom-conserving (line 1763-1765).
  - **The six moves** dispatched in `_add_expansions_mech!` (`src/mechanism_enumeration.jl:1769-1779`), in this exact order:
    1. `_expand_re_to_ss` (line 1207-1214) — flips a whole catalytic kinetic group RE→SS **atomically** (all members convert together), only if every member is currently RE.
    2. `_expand_split_kinetic_group` (line 1249-1271) — for each group with ≥2 members, carves one member into a fresh trailing singleton group.
    3. `_expand_add_dead_end_regulator` (line 1367-1539) — adds a `CompetitiveInhibitor` binding step set (one fresh kinetic group / one `K_R`); mirror steps inherit their catalytic counterpart's `kinetic_group`.
    4. `_expand_to_allosteric` (line 1554-1582) — promotes a non-allosteric `Mechanism` to allosteric variants (all-`:EqualAI` baseline plus one `:OnlyA`-per-group variant), enumerated over `allowed_catalytic_multiplicities`; no-op on an already-allosteric input.
    5. `_expand_add_allosteric_regulator` (line 1608-1654) — adds one `AllostericRegulator` at a new or existing regulatory site with a tag from `{:OnlyA, :OnlyI, :NonequalAI}` (plus `:EqualAI` only at an existing mixed site); no-op on `Mechanism`.
    6. `_expand_change_allo_state` (line 1715-1746) — relaxes one allo state from `:EqualAI`/`:OnlyA`/`:OnlyI` to `:NonequalAI`; no-op on `Mechanism`.
  - **`_dedup_flat!`** (`src/mechanism_enumeration.jl:1793-1796`) is `unique!`: mechanisms are canonical at **construction** (the `Mechanism`/`AllostericMechanism` constructor sorts steps, groups, regulatory sites), so dedup is non-mutating structural `==`/`hash`, NOT a separate canonicalization pass. (Correct the README's `dedup!`-canonicalizes claim, `README.md:257-258`.)
  - **Actual-count bucketing**: `_process_batch` (`src/identify_rate_equation.jl:414-449`) computes `n = length(fitted_params(em))` per child and the beam buckets by this **actual fitted-param count** (line 422-423, 461); the search re-fits each child. `max_param_count` drops a child whose actual count exceeds the cap **before** fitting (line 423), bounding search depth (coverage gap — fold here).
  - **`FitFailure` / loud failures** (coverage gap — fold here): `_process_batch` returns `(entries, failures)`; a compile/fit that throws becomes a `FitFailure` carrying the exception text (`src/identify_rate_equation.jl:340-348, 442-444`), never silently dropped; an all-failed base tier re-raises (line 510-520); CSV rows carry `retcode`/`error` (`src/identify_rate_equation.jl:355-365`).
  - **CSV + `progress.log` artifacts** (coverage gap — fold here): `save_dir` is mandatory and writes `initial_mechanisms.csv` plus `equation_search_iteration_N.csv` (`src/identify_rate_equation.jl:231-242`); `_progress` appends each master-level line to `<save_dir>/progress.log` and flushes stdout for cluster visibility (`src/identify_rate_equation.jl:367-382`).

- [ ] **Step 2: Write the enumeration-engine page**
  Create `docs/src/identify/enumeration_engine.md`. The executor WRITES the prose with the elements-of-style:writing-clearly-and-concisely skill and standardizes MWC terminology on A/I (state R≡A, T≡I once where the allosteric moves are described). The page MUST cover, with citations:
  - **Three composable building blocks, no monolith** (`src/mechanism_enumeration.jl`): `init_mechanisms` (minimum-param mechanisms), `expand_mechanisms` (one move per kind, flat return), `_dedup_flat!` (`unique!` over canonical-by-construction structs). Correct the README's `dedup!` naming and its claim that dedup canonicalizes (canonicalization is in the constructor).
  - **The six expansion moves**, each named, in `_add_expansions_mech!` order (list above). For each move state what it changes and its **parameter delta**, using the "mostly +1 param" framing from the spec: most moves add one fitted parameter, but a Haldane/Wegscheider constraint can absorb the new parameter (**net +0**), and changing a **steady-state** group from `EqualAI` to `NonequalAI` adds **+2** (independent `kf` and `kr` per state). This is precisely why the search re-fits and buckets by **actual** fitted-param count rather than assuming +1.
  - **Actual-count bucketing** (`src/identify_rate_equation.jl:422-423, 461`): exact counts come from `length(fitted_params(compile_mechanism(m)))`; the beam buckets by this.
  - **`max_param_count`** (coverage gap): caps actual fitted params, dropping over-cap children before fitting (`src/identify_rate_equation.jl:423`); bounds search depth, not per-mechanism compile cost.
  - **Loud failures** (coverage gap): `FitFailure` captures exception text (`src/identify_rate_equation.jl:340-348`); an all-failed base tier re-raises (line 510-520); rows carry `retcode`/`error`.
  - **Artifacts** (coverage gap): `save_dir` is mandatory; `initial_mechanisms.csv` + `equation_search_iteration_N.csv` + `progress.log` (cluster-visible, flushed). Cross-link `[`init_mechanisms`](@ref)` and `[`expand_mechanisms`](@ref)` if those are exported/documented; if they are internal, render the names in backticks without `@ref`.

- [ ] **Step 3: Add the "verify the per-move deltas before writing" callout**
  Per the spec's resolved-decisions, the exact per-move parameter deltas are **verify-during-writing**. Embed a literal instruction block at the top of the page's source (an HTML comment, so it does not render) directing the executor to confirm each move's delta against the code before stating it in prose:
  ```markdown
  <!--
  VERIFY-DURING-WRITING: before stating any per-move parameter delta below,
  confirm it against the code. Do not assume +1.
    - _expand_re_to_ss: SS adds a reverse rate, BUT a Haldane/Wegscheider
      constraint may make it dependent (net +0). Confirm via
      fitted_params(compile_mechanism(child)) vs the parent.
    - _expand_split_kinetic_group: +1 (a previously-shared K/k splits in two).
    - _expand_add_dead_end_regulator: +1 (one K_R per regulator group).
    - _expand_to_allosteric: +1 (just L); :OnlyA-per-group variant zeros the
      opposite state, no extra param.
    - _expand_add_allosteric_regulator: +1 (one K_R), or +2 for a NonequalAI
      ligand (independent A/I K).
    - _expand_change_allo_state to NonequalAI: +1 for an RE binding K
      (K_A, K_I), +2 for an SS group (kf, kr each split A/I).
  Confirm each delta empirically with fitted_params before writing the number.
  -->
  ```
  The executor replaces the prose deltas only after confirming them; the comment stays as a maintenance note.

- [ ] **Step 4: Add an `@example` showing `init_mechanisms` returns a non-empty `Vector{Mechanism}`**
  Per the audit rule, do NOT pin enumeration counts (they vary with enumeration changes). Assert invariants instead. Add this exact block (it imports the internal `init_mechanisms` via the module-qualified name):

  ````markdown
  ## Enumeration in practice

  `init_mechanisms` produces the minimum-parameter mechanisms for a reaction:

  ```@example enum
  using EnzymeRates

  rxn = @enzyme_reaction begin
      substrates: S
      products:   P
  end

  mechs = EnzymeRates.init_mechanisms(rxn)
  (count = length(mechs), eltype = eltype(mechs), nonempty = !isempty(mechs))
  ```

  Applying the expansion moves grows the candidate set; `expand_mechanisms`
  returns a flat vector of `Mechanism` / `AllostericMechanism`, which
  `_dedup_flat!` collapses to the structurally distinct ones:

  ```@example enum
  children = EnzymeRates.expand_mechanisms(mechs, rxn)
  EnzymeRates._dedup_flat!(children)
  (eltype = eltype(children), nonempty = !isempty(children))
  ```
  ````

  The surrounding prose states the invariants in words (the result is a non-empty `Vector{Mechanism}`; children are a flat `Union` vector; dedup is structural), not pinned numbers.

- [ ] **Step 5: Ensure the page is wired into `make.jl`**
  Confirm `"The enumeration engine" => "identify/enumeration_engine.md"` is the third entry of the "Identifying the best rate equation" section in `docs/make.jl`'s `pages` (added in Task 4.1 Step 3). If absent, add it there.

- [ ] **Step 6: Build the docs**
  ```bash
  julia --project=docs /home/denis.linux/.julia/dev/EnzymeRates/docs/make.jl
  ```
  Expected: `Documenter: rendering done`, no doctest failures, no `@example` errors. The `enum` blocks run quickly (uni-uni `init_mechanisms` + one `expand_mechanisms` pass is cheap — no fitting). If `init_mechanisms`/`expand_mechanisms` are flagged by `checkdocs` as referenced-but-undocumented, that is expected (they are internal) — `checkdocs = :exports` only checks exports, so this does not fail the build.

- [ ] **Step 7: Commit**
  ```bash
  git add docs/src/identify/enumeration_engine.md docs/make.jl
  git commit -m "docs: enumeration engine — init + six moves + actual-count framing

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
  ```

I have everything needed. `EnzymeMechanism` and `AllostericEnzymeMechanism` already have docstrings; `init_mechanisms` and `compile_mechanism` do too. `Mechanism` and `Step` do NOT — so the developer page's `@docs Mechanism`/`@docs Step` requires adding docstrings first.

Now let me write out the Phase 5 section. The CLAUDE.md content for Source Layout and Known Issues is already in my context. I'll produce the three-task phase.

## Phase 5 — Developer Page, Getting Started, and Roadmap

This phase writes the three pages that fall outside the three learning pillars: the maintainer-facing **Developer / Architecture** page (migrating the internals from CLAUDE.md with the audit's staleness corrections applied), the newcomer-facing **Getting Started** arc that walks define→derive→fit→identify end to end (reusing the fast identify example so nothing slow lands on the front path), and a short prose **Roadmap**. The Developer page uses `@docs` blocks for four internal names; two of them (`Mechanism`, `Step`) currently lack docstrings, so this phase adds those docstrings before referencing them, otherwise the `@docs` block fails the build. All three pages are wired into `docs/make.jl`'s `pages` list, which Phase 0 created with stubs at these exact paths. Assume Phase 0 has produced `docs/Project.toml`, `docs/make.jl`, and stub files at `docs/src/developer.md`, `docs/src/getting_started.md`, and `docs/src/roadmap.md` that are already listed in `pages`.

### Task 5.1: Developer / Architecture page

**Files:**
- Modify: `/home/denis.linux/.julia/dev/EnzymeRates/src/types.jl:132-162` (add `Step` docstring)
- Modify: `/home/denis.linux/.julia/dev/EnzymeRates/src/types.jl:468-482` (add `Mechanism` docstring)
- Create/overwrite: `/home/denis.linux/.julia/dev/EnzymeRates/docs/src/developer.md`

- [ ] **Step 1: Add a docstring to `Step` so the `@docs Step` block resolves**

In `/home/denis.linux/.julia/dev/EnzymeRates/src/types.jl`, the `Step` struct at line 138 is preceded only by `#`-comments (lines 132–137). `@docs Step` requires a real docstring. Insert a docstring immediately above the `struct Step` line, keeping the existing `#`-comment block in place above it. Use the Edit tool to replace the exact text:

old_string (the comment block plus the struct opening):
```julia
# Step: one elementary transition. Binding steps carry
# `bound_metabolite`; iso steps carry `nothing`. All binding steps (RE and
# SS) canonicalize here (bound metabolite on the `to_species` side). All iso
# steps (RE and SS) canonicalize in the Mechanism constructor via
# `_canonical_iso_direction`. After Mechanism construction, every Step is
# canonicalized.
struct Step
```
new_string:
```julia
# Step: one elementary transition. Binding steps carry
# `bound_metabolite`; iso steps carry `nothing`. All binding steps (RE and
# SS) canonicalize here (bound metabolite on the `to_species` side). All iso
# steps (RE and SS) canonicalize in the Mechanism constructor via
# `_canonical_iso_direction`. After Mechanism construction, every Step is
# canonicalized.
"""
    Step

One elementary transition between two enzyme `Species`. A binding step carries
the bound `Metabolite` in `bound_metabolite` (iso steps store `nothing`).
`is_equilibrium` flags a rapid-equilibrium step (`true`) versus a steady-state
step (`false`). Every `Step` is stored in canonical form: binding steps put the
bound metabolite on `to_species`, and iso-step direction is fixed by the
`Mechanism`/`AllostericMechanism` constructor.
"""
struct Step
```

- [ ] **Step 2: Add a docstring to `Mechanism` so the `@docs Mechanism` block resolves**

In the same file, the `Mechanism` struct at line 473 is preceded only by a `#`-comment. Replace the exact text:

old_string:
```julia
# Mechanism: groups elementary steps by kinetic group (outer
# vector). All steps within a group share kinetic parameters. The
# constructor canonicalizes iso-step direction and stores the steps;
# parameter naming and step ordering derive purely from structure and
# flat iteration order.
struct Mechanism
```
new_string:
```julia
# Mechanism: groups elementary steps by kinetic group (outer
# vector). All steps within a group share kinetic parameters. The
# constructor canonicalizes iso-step direction and stores the steps;
# parameter naming and step ordering derive purely from structure and
# flat iteration order.
"""
    Mechanism

A non-allosteric enzyme mechanism: a `reaction::EnzymeReaction` plus
`steps::Vector{Vector{Step}}`, where the outer vector is kinetic groups and each
inner vector holds the steps that share that group's kinetic parameters. The
constructor canonicalizes iso-step direction and sorts steps and groups, so two
mechanisms that differ only in how their steps were written collapse to the same
struct. Lift to the singleton derivation type with `compile_mechanism(m)` or
`EnzymeMechanism(m)`.
"""
struct Mechanism
```

- [ ] **Step 3: Confirm the new docstrings attach (no `@docs` surprise later)**

Run:
```bash
julia --project -e 'using EnzymeRates; for n in (:Mechanism, :Step); s = string(Base.Docs.doc(getfield(EnzymeRates, n))); println(n, " => ", !occursin("No documentation found", s)); end'
```
Expected output:
```
Mechanism => true
Step => true
```

- [ ] **Step 4: Write the Developer / Architecture page**

Write `/home/denis.linux/.julia/dev/EnzymeRates/docs/src/developer.md`. The executor WRITES the prose using the `elements-of-style:writing-clearly-and-concisely` skill; standardize all MWC terminology on **A/I** (active/inactive), stating the R≡A, T≡I correspondence once. Do NOT invent function behavior — every claim below is code-verified; cite the file:function in the prose only where it aids a maintainer.

The page must open with this exact frontmatter line as the first line:
```markdown
# Developer / Architecture
```

The page must cover these code-verified FACTS, organized under the section headings shown:

(a) **Concrete vs. singleton mechanism types.**
- `Mechanism` is a concrete struct with two fields — `reaction::EnzymeReaction` and `steps::Vector{Vector{Step}}` (outer = kinetic groups, inner = steps sharing kinetic parameters). Source: `src/types.jl:473-482`.
- `Step` has four fields: `from_species::Species`, `to_species::Species`, `bound_metabolite::Union{Metabolite,Nothing}`, `is_equilibrium::Bool`. Source: `src/types.jl:138-162`.
- The enumeration pipeline runs end-to-end on the concrete `Mechanism` / `AllostericMechanism` structs; there is no separate working representation.
- `EnzymeMechanism{Sig}` (`src/types.jl:777`) is a singleton type used only by the `@generated` rate-equation derivation. `Sig` is the tuple `(reaction_sig, steps_sig)` produced by `_sig_of(::Mechanism)` (`src/types.jl:748-751`).
- `EnzymeMechanism(m::Mechanism)` lifts a concrete mechanism to its singleton type, dropping unbound regulators at this boundary via `_drop_unbound_regulators` (`src/types.jl:787-817`); `Mechanism(em::EnzymeMechanism)` lifts back via `_mechanism_from_sig` (`src/types.jl:819`, `:753-756`). State that `Sig` is purely structural: two mechanisms that differ only in source step order collapse to the same `EnzymeMechanism` type.
- `compile_mechanism` is the one-name lift for both families: `compile_mechanism(m::Mechanism) = EnzymeMechanism(m)` and `compile_mechanism(am::AllostericMechanism) = AllostericEnzymeMechanism(am)` (`src/mechanism_enumeration.jl:1143-1144`). It is internal (not exported), reached as `EnzymeRates.compile_mechanism`.

(b) **The three-type-parameter rationale for `AllostericEnzymeMechanism`.**
- `AllostericEnzymeMechanism{CatalyticMech, CatSites, RegSites}` (`src/types.jl:834-836`). Explain why the three parameters cannot be folded into one value-tuple `Sig` the way `EnzymeMechanism{Sig}` does: the first slot is a `DataType` (an `EnzymeMechanism` subtype), and Julia rejects a `DataType` in the value-tuple position of a type parameter. Source: the comment at `src/types.jl:830-833`.

(c) **The `name(p, m)` parameter-naming chokepoint.**
- All `Parameter → Symbol` rendering flows through one function family, `name(p::Parameter, m)`, at `src/types.jl:1341-1429`. Routing all parameter-name production through one place keeps any name-scheme change a single-function edit.
- Step-governed parameters carry the `Step` they name (and a `state::Symbol`); `Kreg` carries its `RegulatorySite` + `AllostericRegulator` + state; `Keq`/`Etot`/`Lallo` are stateless. Source: the `Parameter` family at `src/types.jl:179-219`.
- Rendering helpers: `_state_tag` maps the A/I state token to a prefix (`:A`→`"A_"`, `:I`→`"I_"`, `:EqualAI`/`:None`→`""`, `src/types.jl:1353-1360`); `_render_binding` names a binding param as metabolite + pre-binding form (`src/types.jl:1366-1372`); `_render_iso` names an iso param by directed species pair (`src/types.jl:1375-1377`). Give the worked examples: a binding K renders `:K_S_E`, an inactive-state binding K renders `:K_I_S_E`, an SS iso rate renders `:k_ES_to_EP`.
- The chokepoint is enforced by an AST-walker test, NOT a regex and NOT a separate `test_chokepoint.jl` file. The test lives in `test/test_types.jl:1577-1644`: `_walk_violations!` parses each `src/*.jl` with `Meta.parseall`, recognizes a chokepoint method body via `_is_chokepoint_def` (a `name` method dispatching on a `Parameter` subtype value), and fails the build if any `Symbol("K…")`/`Symbol("k…")`/`Symbol("V…")`/`Symbol("L…")` literal appears outside a chokepoint body. Cross-link the test by path.

(d) **Canonical Step Form (load-bearing, not cosmetic).**
- Binding steps (RE and SS) are canonicalized in the `Step` constructor so the bound metabolite is always on `to_species` (free enzyme + free metabolite on the from side): `E + S ⇌ ES`, never `ES ⇌ E + S`; a product-release step `EP ⇌ E + P` stores as `E + P → EP`. Source: `src/types.jl:143-161`.
- Iso steps (RE and SS) are canonicalized to the physical-forward direction by `_canonical_iso_direction` (`src/types.jl:407-427`), called for every group by `_canonicalize_iso_groups` (`src/types.jl:432-440`) inside the `Mechanism`/`AllostericMechanism` constructors. State the three tiers: Tier 1 = atom-balance progression (substrate-bound count up, product-bound count down); Tier 2 = 1-hop binding-graph context via `_entry_kind` (`src/types.jl:386-399`), where a product-only form → substrate-only form is forward; Tier 3 = lexical fallback.
- `_entry_kind` classifies a form by how it participates in ALL binding steps (RE and SS, not RE-only) as the free side, because "substrate-entry / product-exit" is a chemistry fact independent of the RE/SS flag. Source: the comment at `src/types.jl:373-385`.
- Step and group ORDER (not just direction) is canonicalized by `_canonical_group_order!` (`src/types.jl:461-466`) using `_step_canonical_key` (`src/types.jl:445-447`). State explicitly that this is load-bearing: the Haldane/Wegscheider reduction chooses which parameters are dependent by step order, so the reduced rate equation and `fitted_params` would depend on step order without it. Because the constructor canonicalizes, `_dedup_flat!` is pure non-mutating `==`/`hash` comparison (`unique!`).

(e) **The `@generated` compile-budget Known Issue.**
- `rate_equation(m, conc, params)` and friends use `@generated` functions that derive the rate equation at compile time by the King–Altman/Cha method; each unique `EnzymeMechanism` type triggers a full symbolic derivation.
- For mechanisms with many enzyme forms/steps, this compilation can be very slow, exhaust memory, or `StackOverflow`. `identify_rate_equation` caps expansion via `max_param_count` (bounding search depth) but not per-mechanism compile cost, so one very large mechanism can still be slow to compile.
- Cross-link the runtime contract: once compiled, `rate_equation` MUST be allocation-free and sub-100 ns per call, enforced by `test_rate_equation_performance` (`test/test_rate_eq_derivation.jl:387`, asserts `allocs == 0` at `:821` and `t < 100e-9` at `:822`). Link this back to the CLAUDE.md "runtime perf is non-negotiable" rule (a guard line stays in CLAUDE.md per the migration table).

(f) **Source Layout** — a maintainer's map of `src/`. Reproduce the per-file purposes from the CLAUDE.md "Source Layout" section, one bullet per file (`src/types.jl`, `src/dsl.jl`, `src/sym_poly_for_rate_eq_derivation.jl`, `src/rate_eq_derivation.jl`, `src/thermodynamic_constr_for_rate_eq_derivation.jl`, `src/fitting.jl`, `src/identify_rate_equation.jl`, `src/mechanism_enumeration.jl`), corrected per the audit staleness ledger — in particular do not reference `_is_ss_rate_constant` (the real name is `_ss_rate_constant_names`, `src/rate_eq_derivation.jl:621`) or `_kcat_components` (the real name is `_kcat_groups_from_polys`, `src/rate_eq_derivation.jl:650`).

The page must end with these four `@docs` blocks for the internal names worth documenting (place them under a "## Internal API" heading; these are the only `@docs` blocks on the page, since the public API lives on the autodocs page):

```@docs
EnzymeRates.Mechanism
EnzymeRates.Step
EnzymeRates.init_mechanisms
EnzymeRates.compile_mechanism
```

- [ ] **Step 5: Build the docs and verify the page renders with no doctest failures**

Run:
```bash
julia --project=docs /home/denis.linux/.julia/dev/EnzymeRates/docs/make.jl
```
Expected: the build finishes with `[ Info: Documenter: rendering done`, no `Error: ` lines, no `doctest failure`, and no warning about a missing docstring for `EnzymeRates.Mechanism`, `EnzymeRates.Step`, `EnzymeRates.init_mechanisms`, or `EnzymeRates.compile_mechanism`. (`@docs` blocks fail the build if any referenced name lacks a docstring; the prior steps added the two missing ones.)

- [ ] **Step 6: Run the chokepoint test to confirm the new docstrings did not introduce a violation**

The two added docstrings live in `src/types.jl`, which the AST-walker scans. Confirm they are clean:
```bash
julia --project -e 'using Test; include("test/test_types.jl")' 2>&1 | tail -5
```
Expected: the `chokepoint: no Symbol("[KkVL]...") outside parameter-name renderers` testset passes (no `chokepoint violations` `@info` line, all `@test` green). The docstrings contain no `Symbol("…")` calls, so they cannot trip the walker.

- [ ] **Step 7: Commit**
```bash
git add src/types.jl docs/src/developer.md
git commit -m "$(cat <<'EOF'
Docs: Developer/Architecture page + docstrings for Mechanism and Step

Migrate maintainer internals from CLAUDE.md (corrected per the audit
staleness ledger): the name(p,m) chokepoint and its AST-walker test in
test/test_types.jl, Canonical Step Form, the concrete-vs-singleton split,
the AllostericEnzymeMechanism three-type-param rationale, Source Layout,
and the @generated compile-budget note. Add docstrings to Mechanism and
Step so their @docs blocks resolve.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
EOF
)"
```

### Task 5.2: Getting Started page

**Files:**
- Create/overwrite: `/home/denis.linux/.julia/dev/EnzymeRates/docs/src/getting_started.md`

- [ ] **Step 1: Write the Getting Started page**

Write `/home/denis.linux/.julia/dev/EnzymeRates/docs/src/getting_started.md`. The executor WRITES the prose with the `elements-of-style:writing-clearly-and-concisely` skill; standardize MWC terminology on A/I. The page is a single scannable arc for newcomers — define → derive → fit → identify — using `@example` blocks throughout (NOT `jldoctest`, because the fit and identify steps include random multi-start output that cannot be pinned). It reuses the fast identify example so nothing slow lands on the front path; the heavy production-search prose stays on the Identify pillar page.

The page's first line must be exactly:
```markdown
# Getting Started
```

The page must contain these `@example` blocks in order. They share one named sandbox (`getting-started`) so later blocks see earlier bindings. Use these EXACT code blocks (they are real, runnable, and fast):

Block 1 — **define** a reaction and a textbook RE mechanism:
````markdown
```@example getting-started
using EnzymeRates

rxn = @enzyme_reaction begin
    substrates: S
    products:   P
end

m = @enzyme_mechanism begin
    substrates: S
    products:   P
    catalytic_steps: begin
        E + S ⇌ E(S)
        E(S) <--> E(P)
        E(P) ⇌ E + P
    end
end

m isa EnzymeMechanism
```
````

Block 2 — **derive**: list parameters, metabolites, and print the rate-equation string (use `print(...)` to avoid quote/escape noise, per the audit doctest rules):
````markdown
```@example getting-started
parameters(m)
```

```@example getting-started
metabolites(m)
```

```@example getting-started
print(rate_equation_string(m))
```
````

Block 3 — **evaluate** the rate numerically at one concentration point with a hand-supplied parameter `NamedTuple` (this demonstrates the exported `rate_equation` call; values are illustrative, not fitted):
````markdown
```@example getting-started
params = (
    K_S_E = 1.0,
    K_P_E = 1.0,
    k_ES_to_EP = 10.0,
    Keq = 2.0,
    E_total = 1.0,
)
concs = (S = 2.0, P = 0.5)

rate_equation(m, concs, params)
```
````

Block 4 — **fit**: generate noiseless synthetic data from `rate_equation` on a small concentration grid, build a `FittingProblem`, and fit with the blessed bring-your-own optimizer `PyCMAOpt()` from `OptimizationPyCMA`. Keep restarts/maxtime light so the block runs in seconds. The prose must state that fit output is a random multi-start result, so only its keys/shapes are stable (which is why this page uses `@example`, not `jldoctest`):
````markdown
```@example getting-started
using OptimizationPyCMA, Random
Random.seed!(1)

groups = String[]; Rate = Float64[]; Svals = Float64[]; Pvals = Float64[]
for g in 1:3, _ in 1:8
    s = 0.1 + 9.9 * rand()
    p = 0.1 + 9.9 * rand()
    push!(groups, "G$g")
    push!(Rate, rate_equation(m, (S = s, P = p), params))
    push!(Svals, s); push!(Pvals, p)
end
data = (group = groups, Rate = Rate, S = Svals, P = Pvals)

fp = FittingProblem(m, data; Keq = 2.0)
fitted_params, loss, retcode = fit_rate_equation(fp, PyCMAOpt(); n_restarts = 1, maxtime = 2.0)
retcode
```
````

Block 5 — **identify** (the fast example): build an `IdentifyRateEquationProblem` on the same small reaction and data, then run `identify_rate_equation` with the greedy width-1 beam settings (`min_beam_width = 1`, `loss_rel_threshold = 1.0`, `loss_abs_threshold = 0.0`) and a low `max_param_count`, so the search collapses to one survivor per parameter-count level and runs in seconds. `save_dir` is mandatory; point it at a temp dir. The prose must note that the full production search uses a wider default beam, explores more candidates, and can run ~1 hour — so this fast variant is what the docs run:
````markdown
```@example getting-started
prob = IdentifyRateEquationProblem(rxn, data; Keq = 2.0)

results = identify_rate_equation(prob;
    min_beam_width = 1,
    loss_rel_threshold = 1.0,
    loss_abs_threshold = 0.0,
    max_param_count = 6,
    n_cv_candidates = 1,
    save_dir = mktempdir(),
    pmap_function = map,
    optimizer = PyCMAOpt(),
    n_restarts = 1, maxtime = 1.0,
    show_progress = false)

print(rate_equation_string(results.best))
```
````

The prose around these blocks must, in order: (1) name the three pillars the arc touches (deriving, fitting, identifying) and link each block forward to its pillar page with `[…](@ref)` cross-links — Deriving rate equations, Fitting rate equations, Identifying the best rate equation; (2) explain the `@enzyme_mechanism` step grammar exactly: `⇌` marks a rapid-equilibrium binding step (one binding constant K), `<-->` marks a steady-state interconversion (independent `kf`/`kr`); (3) note that `m isa EnzymeMechanism` (not `typeof(m)`) is shown deliberately because the `Sig` type parameter is an unreadable string; (4) describe the data format: a `group` column identifies measurement batches sharing one `E_total`, a `Rate` column holds measured rates, and one column per metabolite/regulator name; (5) explain that `Keq` is always user-supplied, never estimated.

- [ ] **Step 2: Build the docs and verify every `@example` runs**

Run:
```bash
julia --project=docs /home/denis.linux/.julia/dev/EnzymeRates/docs/make.jl
```
Expected: `[ Info: Documenter: rendering done`, no `Error:` lines, and no `@example` block failures (every block, including the `fit_rate_equation` and `identify_rate_equation` blocks, executes to completion). The identify block must finish in seconds, not minutes — if it hangs, the greedy-beam args or `max_param_count` were dropped.

- [ ] **Step 3: Commit**
```bash
git add docs/src/getting_started.md
git commit -m "$(cat <<'EOF'
Docs: Getting Started end-to-end arc (define, derive, fit, identify)

Single scannable newcomer page using @example throughout; reuses the
fast width-1-beam identify example so nothing slow lands on the front
path. Cross-links forward to the three pillar pages.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
EOF
)"
```

### Task 5.3: Roadmap page

**Files:**
- Create/overwrite: `/home/denis.linux/.julia/dev/EnzymeRates/docs/src/roadmap.md`

- [ ] **Step 1: Write the Roadmap page**

Write `/home/denis.linux/.julia/dev/EnzymeRates/docs/src/roadmap.md`. This is a short prose page — no code blocks, no doctests. The executor WRITES the prose with the `elements-of-style:writing-clearly-and-concisely` skill; standardize MWC terminology on A/I where it appears. Keep each item to one or two sentences; frame them as planned/future directions, not commitments with dates.

The page's first line must be exactly:
```markdown
# Roadmap
```

The page lists exactly these five planned directions, one short subsection or bullet each:
- **KNF allostery model** — a Koshland–Némethy–Filmer (sequential) allosteric model alongside the existing MWC (concerted A/I) model.
- **Parameter identifiability** — reporting which fitted parameters the data can and cannot constrain, beyond model selection.
- **Iso mechanisms** — broader support for isomerization (iso) mechanisms in enumeration and derivation.
- **Plotting** — built-in plotting of fitted rate equations against data.
- **Outlier dataset identification** — flagging measurement groups that the selected mechanism fits poorly, to surface suspect datasets.

The intro paragraph must state that these are directions under consideration and that the current scope is the three pillars (deriving, fitting, identifying), linking to Getting Started with `[Getting Started](@ref)`.

- [ ] **Step 2: Build the docs and verify the page renders**

Run:
```bash
julia --project=docs /home/denis.linux/.julia/dev/EnzymeRates/docs/make.jl
```
Expected: `[ Info: Documenter: rendering done`, no `Error:` lines, no broken-cross-reference warning for `[Getting Started](@ref)` (the Getting Started page from Task 5.2 supplies that anchor), no doctest failures.

- [ ] **Step 3: Commit**
```bash
git add docs/src/roadmap.md
git commit -m "$(cat <<'EOF'
Docs: Roadmap page (KNF, identifiability, iso, plotting, outlier detection)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
EOF
)"
```

I have all exact line ranges. The CLAUDE.md structure maps cleanly to the migration table. Now I'll draft the Phase 6 section.

## Phase 6 — Trim and Slim

This phase is last in the documentation effort because every page the README and CLAUDE.md will point at must already exist (Phases 0–5). It applies three trims: the README shrinks to a thin landing page that links into the published docs site; `test/test_readme_runs.jl` retires because its coverage moved into Documenter doctests (a net coverage increase); and CLAUDE.md splits by audience per the spec's migration table (section 7) — maintainer-internal and user-concept blocks are **cut**, the load-bearing ones leave a one-line guard pointer into the docs, and every agent-behavioral rule stays untouched. The CLAUDE.md edits are explicit and conservative: cut + guard, never a silent drop of a rule. The phase ends with the full test suite green (after removing the README test) and the docs build green.

### Task 6.1: Slim README.md to a landing page

**Files:**
- Modify: `/home/denis.linux/.julia/dev/EnzymeRates/README.md` (full rewrite, currently 296 lines)

- [ ] **Step 1: Confirm the docs site URL the README links to**
  Run this to confirm the canonical URL already wired into `docs/make.jl` (set in Phase 0), so the README links match the deployed site:
  ```bash
  grep -n "canonical\|github.io" /home/denis.linux/.julia/dev/EnzymeRates/docs/make.jl
  ```
  Expected: a line containing `canonical = "https://DenisTitovLab.github.io/EnzymeRates.jl"`. The README links target `https://DenisTitovLab.github.io/EnzymeRates.jl/stable/`. If the grep shows a different host/org, use that host with the `/stable/` suffix instead — do not invent a URL.

- [ ] **Step 2: Confirm the docs page filenames the README deep-links to exist**
  Run:
  ```bash
  ls /home/denis.linux/.julia/dev/EnzymeRates/docs/src/
  ```
  Expected: the page files created in Phases 0–5 are present, including `index.md`, a Getting Started page, a deriving-rate-equations tutorial page, a fitting tutorial page, an identify tutorial page, and `api.md` (or equivalently named). Note the exact slugs (Documenter maps `docs/src/getting_started.md` → `/stable/getting_started/`). If the slugs differ from the ones used in Step 3 below, substitute the real slugs in the README link targets — the four links must resolve to pages that exist.

- [ ] **Step 3: Replace README.md with the landing page**
  Overwrite `/home/denis.linux/.julia/dev/EnzymeRates/README.md` with exactly this content. The badges block is copied verbatim from the current README (lines 3–6) so the CI/coverage/JET/Aqua badges are preserved. The quickstart is the README's existing install plus a ~10-line derive-and-evaluate snippet (a strict subset of the current running example — no fitting, no data generation, no `identify_rate_equation` — so the landing page stays thin and the heavy content lives in the docs):

  ```markdown
  # EnzymeRates.jl

  [![Build Status](https://github.com/DenisTitovLab/EnzymeRates.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/DenisTitovLab/EnzymeRates.jl/actions/workflows/CI.yml?query=branch%3Amain)
  [![Coverage](https://codecov.io/gh/DenisTitovLab/EnzymeRates.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/DenisTitovLab/EnzymeRates.jl)
  [![JET](https://img.shields.io/badge/%F0%9F%9B%A9%EF%B8%8F_tested_with-JET.jl-233f9a)](https://github.com/aviatesk/JET.jl)
  [![Aqua](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

  Identify the best enzyme rate equation from kinetic data. Given a reaction
  definition and experimental rate measurements at varying substrate, product,
  and regulator concentrations, EnzymeRates enumerates biochemically valid
  mechanisms, fits each to the data, and selects the simplest mechanism that
  generalizes by leave-one-group-out cross-validation. It has first-class
  support for MWC allostery and for mechanisms that mix steady-state and
  rapid-equilibrium elementary steps, and it derives Haldane/Wegscheider
  thermodynamic constraints automatically from the mechanism's cycle structure.

  ## Installation

  ```julia
  using Pkg
  Pkg.add(url="https://github.com/DenisTitovLab/EnzymeRates.jl")
  ```

  ## Quickstart

  Define a mechanism, derive its symbolic rate equation, and evaluate it:

  ```julia
  using EnzymeRates

  m = @enzyme_mechanism begin
      substrates: S
      products:   P
      steps: begin
          E + S ⇌ E(S)
          E(S) <--> E(P)
          E(P) ⇌ E + P
      end
  end

  parameters(m)                              # parameter names to supply
  print(rate_equation_string(m))             # the symbolic rate equation
  rate_equation(m, (S=1e-4, P=1e-5),         # evaluate numerically
      (K_S_E=1e-4, k_ES_to_EP=100.0, k_EP_to_ES=1.0, K_P_E=1e-5, Keq=2.0, E_total=1.0))
  ```

  ## Documentation

  Full documentation — tutorials for deriving, fitting, and identifying rate
  equations, plus the architecture and API reference — lives at the
  [documentation site](https://DenisTitovLab.github.io/EnzymeRates.jl/stable/):

  - [Getting Started](https://DenisTitovLab.github.io/EnzymeRates.jl/stable/getting_started/) — the end-to-end define → derive → fit → identify arc.
  - [Deriving rate equations](https://DenisTitovLab.github.io/EnzymeRates.jl/stable/deriving/textbooks/) — RE vs steady state, the Cha/King–Altman algorithm, thermodynamic constraints, ping-pong, and MWC allostery.
  - [Fitting rate equations](https://DenisTitovLab.github.io/EnzymeRates.jl/stable/fitting/tutorial/) — the data format, normalized vs absolute rate, and optimizer choice.
  - [Identifying the best rate equation](https://DenisTitovLab.github.io/EnzymeRates.jl/stable/identify/tutorial/) — the enumeration engine, beam search, and cross-validation model selection.
  - [API Reference](https://DenisTitovLab.github.io/EnzymeRates.jl/stable/api/) — every exported type, macro, and function.
  ```

  Note for the executor: in Step 3 the five doc-link slugs (`getting_started/`, `deriving/textbooks/`, `fitting/tutorial/`, `identify/tutorial/`, `api/`) must match the real page slugs found in Step 2. If a Phase 0–5 page lives at a different path, edit only the URL path component — keep the host and the `/stable/` prefix.

- [ ] **Step 4: Verify the quickstart snippet runs**
  The README test is about to be retired (Task 6.2), but the quickstart must still be correct. Run the quickstart block directly to confirm it executes and the parameter names match what `parameters(m)` reports:
  ```bash
  julia --project=/home/denis.linux/.julia/dev/EnzymeRates -e 'using EnzymeRates; m = @enzyme_mechanism begin
      substrates: S
      products:   P
      steps: begin
          E + S ⇌ E(S)
          E(S) <--> E(P)
          E(P) ⇌ E + P
      end
  end; println(parameters(m)); println(rate_equation(m, (S=1e-4, P=1e-5), (K_S_E=1e-4, k_ES_to_EP=100.0, k_EP_to_ES=1.0, K_P_E=1e-5, Keq=2.0, E_total=1.0)))'
  ```
  Expected: prints a parameter-name tuple/vector that contains exactly `K_S_E`, `k_ES_to_EP`, `k_EP_to_ES`, `K_P_E`, `Keq`, `E_total` (order may differ), and a finite positive `Float64` rate. If the printed parameter names differ from the `params` keys in the snippet (e.g. the SS rate is named differently for this step ordering), update the README snippet's `params` keys to the names `parameters(m)` actually prints, then re-run until they match. Do NOT invent names — use the printed ones.

- [ ] **Step 5: Commit**
  ```bash
  git add /home/denis.linux/.julia/dev/EnzymeRates/README.md
  git commit -m "Slim README to a landing page linking into the docs site

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
  ```

### Task 6.2: Retire the README test

**Files:**
- Delete: `/home/denis.linux/.julia/dev/EnzymeRates/test/test_readme_runs.jl`
- Modify: `/home/denis.linux/.julia/dev/EnzymeRates/test/runtests.jl:17`

- [ ] **Step 1: Remove the include from runtests.jl**
  In `/home/denis.linux/.julia/dev/EnzymeRates/test/runtests.jl`, delete exactly this line (line 17):
  ```julia
      include("test_readme_runs.jl")
  ```
  Leave every other line in the `@testset "EnzymeRates.jl"` block unchanged. After the edit, the include block reads (lines 10–19):
  ```julia
      include("test_accessors.jl")
      include("test_types.jl")
      include("test_dsl.jl")
      include("test_rate_eq_derivation.jl")
      include("test_fitting.jl")
      include("test_mechanism_enumeration.jl")
      include("test_identify_rate_equation.jl")
      include("test_aqua_jet.jl")
      include("test_compile_budget.jl")
  ```

- [ ] **Step 2: Delete the README test file**
  ```bash
  git rm /home/denis.linux/.julia/dev/EnzymeRates/test/test_readme_runs.jl
  ```
  Expected: `rm 'test/test_readme_runs.jl'`.

- [ ] **Step 3: Confirm nothing else references the deleted file**
  ```bash
  grep -rn "test_readme_runs" /home/denis.linux/.julia/dev/EnzymeRates/ --include=*.jl
  ```
  Expected: no output (no remaining references). If anything prints, that reference must be removed before continuing.

- [ ] **Step 4: Confirm the test suite still loads (quick parse/include check)**
  ```bash
  grep -c "include(" /home/denis.linux/.julia/dev/EnzymeRates/test/runtests.jl
  ```
  Expected: `10` (was 11 includes — the `mechanism_definitions...` include plus 10 testset includes; after removing the README test, 9 testset includes + 1 fixture include = 10). The full-suite run is the real gate (Task 6.4).

- [ ] **Step 5: Commit**
  ```bash
  git add /home/denis.linux/.julia/dev/EnzymeRates/test/runtests.jl
  git commit -m "Retire test_readme_runs.jl; doctest coverage replaces it (net increase)

The README example coverage moved into Documenter jldoctest/@example blocks,
which run and output-check the deterministic exported functions. Net coverage
increase: doctests check rendered output, the README test only checked that
the blocks ran.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
  ```

### Task 6.3: Apply the CLAUDE.md migration table (cut + guard)

**Files:**
- Modify: `/home/denis.linux/.julia/dev/EnzymeRates/.claude/CLAUDE.md`

This task applies spec section 7's migration table. The blocks listed under "Move to the Developer page" and "Move to how-it-works pages" are **cut** from CLAUDE.md (their corrected content already lives in the Phase 1–5 docs pages). The load-bearing ones leave a one-line guard pointer. Every agent-behavioral rule (lines 1–179: Foundational rules through Learning and Memory Management) stays byte-for-byte unchanged, as do the Package Goal, API Design, Commands, Workflow, and Code Style blocks. Each step below cuts one block and, where the spec marks it load-bearing, leaves a guard line. The edits are conservative: nothing is silently dropped — a cut block either moves to a docs page (already written) or leaves a guard.

- [ ] **Step 1: Cut the "Parameter naming chokepoint" block, leave a guard**
  This block (CLAUDE.md lines 230–237) is maintainer-internal per the migration table; its content (corrected: enforcement is `test/test_types.jl:1577-1644`, not a `test_chokepoint.jl`) now lives on the Developer page. Replace the entire block:
  ```markdown
  ### Parameter naming chokepoint

  All `Parameter → Symbol` rendering in `src/` flows through
  `name(p::Parameter, m)`. Step-governed parameters carry the `Step` they name;
  regulatory parameters carry their `RegulatorySite` and ligand; scalar
  parameters are stateless.

  `test/test_chokepoint.jl` enforces this — any direct `Symbol("K…")`/`Symbol("k…")`/`Symbol("V…")`/`Symbol("L…")` literal outside a parameter-name renderer fails the build (AST walker, not regex).
  ```
  with this single guard line:
  ```markdown
  ### Parameter naming chokepoint (guard)

  All `Parameter → Symbol` rendering flows through the `name(p, m)` chokepoint; the AST-walker test at `test/test_types.jl:1577-1644` fails the build on any stray `Symbol("K…")`/`Symbol("k…")`/`Symbol("V…")`/`Symbol("L…")` literal outside a parameter-name renderer. See the Developer page in the docs for the full rationale.
  ```

- [ ] **Step 2: Cut the "Canonical Step Form" block, leave a load-bearing guard**
  This block (CLAUDE.md lines 239–244) is dual-nature: the chemistry is user-facing (now in the docs) but the load-bearing warning must stay so an agent does not "simplify" the canonicalization. Per the migration table the function names move to the Developer page; per the spec this is explicitly the load-bearing guard example. Replace the entire block:
  ```markdown
  ### Canonical Step Form
  - Every `Step` is canonicalized. **Binding steps (RE and SS)** are normalized in the `Step` constructor so the bound metabolite is always on the bound (to) side / free-enzyme + free-metabolite on LHS: `E + S ⇌ ES`, never `ES ⇌ E + S`. This holds for product-binding (release) steps too — `EP ⇌ E + P` stores as `E + P → EP`, so a product-release step is a product-binding step with kon/koff for the product.
  - **Iso steps (RE and SS)** are canonicalized to the physical-forward direction in the `Mechanism` / `AllostericMechanism` constructor via `_canonical_iso_direction` (Tier 1: substrate/product bound counts; Tier 2: 1-hop binding-graph `_entry_kind` — product-exit→substrate-entry is forward; Tier 3: lex fallback). `from` = side further along the substrate→product progression. This makes structural iso/binding names direction-independent of how the user wrote the step.
  - **Step and group ORDER** (not only direction) is canonicalized in the `Mechanism` / `AllostericMechanism` constructor (`_canonical_group_order!` + `_step_canonical_key`), so two step orderings build the identical struct and `_dedup_flat!` is pure `==`/`hash` comparison. This is **load-bearing, not cosmetic**: the Haldane/Wegscheider reduction chooses which parameters are dependent by step order, so the reduced rate equation and `fitted_params` depend on step order without it. The independent-parameter set and the Haldane/Wegscheider RHS factors are additionally name-sorted (`_dependent_param_exprs`, `build_power_expr`) so graph-distinct but rate-equivalent mechanisms render identical equation strings and share a dedup key.
  - Hand-derived analytical test oracles use the chemical-forward convention (`kNf` = product-making direction), which differs from the package's binding convention for product-binding steps; tests feed them through the `analytical_oracle_params` shim, which swaps `kf`↔`kr` for product-binding steps and maps each oracle's as-written step order onto the canonical stored order (via the spec's `source_steps`). Oracle formulas are never modified.
  - `_binding_K_symbols` relies on the metabolite-on-bound-side invariant: checks only for metabolite on the free/LHS side, no RHS check needed.
  ```
  with this single load-bearing guard line:
  ```markdown
  ### Canonical Step Form (load-bearing guard)
  Step direction, step order, and group order are canonicalized in the `Step` and `Mechanism`/`AllostericMechanism` constructors. This is **load-bearing, not cosmetic**: the Haldane/Wegscheider reduction picks dependent parameters by step order, so `fitted_params` and the reduced rate equation depend on it; `_dedup_flat!` is pure `==`/`hash` only because construction is canonical. Do not relax or reorder this without reading the Developer page in the docs.
  ```

- [ ] **Step 3: Cut the "Regulator representation" block (stale — no guard needed)**
  This block (CLAUDE.md lines 246–250) is user-concept and, per spec section 8, almost entirely stale: it names `regulator_roles` (does not exist) and role labels `:unknown`/`:dead_end` (the real roles are `:competitive` and `:allosteric` only, `src/dsl.jl:51`). The corrected content now lives on the DSL how-it-works page. Delete the entire block (heading and all five bullets):
  ```markdown
  ### Regulator representation
  - Regulators are `(name::Symbol, role::Symbol)` pairs in `EnzymeReaction` type parameter `R` (kept for `@enzyme_reaction` declaration).
  - `regulators(rxn)` returns bare `Symbol` names; `regulator_roles(rxn)` returns `(name, role)` pairs.
  - `@enzyme_reaction` DSL accepts `regulators:` (role `:unknown`), `dead_end_inhibitors:` (`:dead_end`), `allosteric_regulators:` (`:allosteric`), `oligomeric_state:` (Int, default 1).
  - At the `EnzymeMechanism` level, regulators are bare Symbols (no role tags) — the role distinction is captured by the mechanism's structure (which forms the regulator binds to and which kinetic_group those steps belong to).
  ```
  Leave nothing in its place (it is a pure user-concept block, fully migrated and corrected on the DSL docs page).

- [ ] **Step 4: Cut the "Allosteric state taxonomy" block**
  This block (CLAUDE.md lines 252–269) is user-concept (the A/I taxonomy, standardized on A/I in the docs). The corrected content lives on the MWC allostery how-it-works page; the synth-dep dependent-name subtlety lives on the Developer page. Delete the entire block from the `### Allosteric state taxonomy (per kinetic group, per regulatory ligand)` heading through the end of the `docs/superpowers/specs/...nonequalai-rank-validity.md` bullet (lines 252–269). Leave nothing in its place.

- [ ] **Step 5: Cut the "Mirror / dead-end kinetic-group sharing" block**
  This block (CLAUDE.md lines 271–272) describes enumeration internals now on the Developer / enumeration docs page. Delete the entire block:
  ```markdown
  ### Mirror / dead-end kinetic-group sharing
  - The mechanism-enumeration generator assigns dead-end mirror steps the same `kinetic_group::Int` as their catalytic counterpart. Mirror propagation is implicit in the kinetic-group atomicity: when a group's RE→SS conversion fires, every member converts together.
  ```
  Leave nothing in its place.

- [ ] **Step 6: Cut the "Parameter naming convention" block**
  This block (CLAUDE.md lines 274–276) describes the structural parameter-name scheme (`:K_S_E`, `:k_ES_to_EP`, `:K_I_S_E`), now on the docs derivation/parameters pages. Delete the entire block:
  ```markdown
  ### Parameter naming convention
  - `parameters(m)` returns structural symbols derived from the representative step of each kinetic group. Binding constants use metabolite/form names such as `:K_S_E`; SS iso rates use directed species names such as `:k_ES_to_EP` and `:k_EP_to_ES`; inactive-state symbols use the `I_` state token such as `:K_I_S_E`.
  - This is consistent across `parameters`, `_dependent_param_exprs`, `rate_equation_string`, and `_kcat_forward`.
  ```
  Leave nothing in its place.

- [ ] **Step 7: Cut the "Catalytic topology constraints" block**
  This block (CLAUDE.md lines 278–293) is user-facing chemistry (C1–C10, now on the catalytic-topology docs page) plus dev-only function names. Per spec section 8 it carries a stale name: `has_residual` in `backtrack!` is really `pingpong_intermediate::Bool` (`src/mechanism_enumeration.jl:328`). The corrected content lives in the docs. Delete the entire block from the `### Catalytic topology constraints` heading (line 278) through the `- Bystander mechanisms not needed at init level...` bullet (line 293). Leave nothing in its place.

- [ ] **Step 8: Cut the "Mechanism enumeration building blocks" block**
  This block (CLAUDE.md lines 295–309) describes the three composable functions and the expansion moves — user-concept (six moves, flat return) on the enumeration docs page, with the lift details on the Developer page. Delete the entire block from the `### Mechanism enumeration building blocks` heading (line 295) through the final `- Allosteric conversion is +1 param...` bullet (line 309). Leave nothing in its place.

- [ ] **Step 9: Cut the "Source Layout" block**
  Source Layout (CLAUDE.md lines 311–320) is maintainer-internal per the migration table; it now lives on the Developer page. Delete the entire block from the `## Source Layout` heading (line 311) through the `src/mechanism_enumeration.jl` bullet ending at line 320. Leave nothing in its place.

- [ ] **Step 10: Cut the "Vmax Normalization" block**
  This block (CLAUDE.md lines 322–334) is user-concept (kcat normalization) with corrected function names now in the docs. Per spec section 8 it carries two stale names: `_is_ss_rate_constant` is really `_ss_rate_constant_names` (`src/rate_eq_derivation.jl:621`, classifies by `Parameter` subtype, not "lowercase k + digit") and `_kcat_components` is really `_kcat_groups_from_polys` (`src/rate_eq_derivation.jl:650`). The corrected content lives on the fitting "Normalized vs absolute rate" docs page (user) and the Developer page (internals). Delete the entire block from the `## Vmax Normalization (kcat factoring) — IMPLEMENTED` heading (line 322) through the `### Key properties` bullet ending at line 334. Leave nothing in its place.

- [ ] **Step 11: Cut the "Known Issues" compile-limits block**
  This block (CLAUDE.md lines 336–343) is the `@generated` compile-budget note — maintainer-internal per the migration table, now on the Developer page. Delete the entire block from the `## Known Issues` heading (line 336) through the `- identify_rate_equation processes candidates...` bullet ending at line 343. Leave nothing in its place.

- [ ] **Step 12: Cut the "Testing" repo block, KEEP the runtime-perf guard**
  The `## Testing` block (CLAUDE.md lines 345–352) describes test-fixture internals — maintainer detail now on the Developer page. Delete lines 345–352 (the `## Testing` heading through the `- kcat/rescaling tests...` bullet). **Do NOT touch** the `### `rate_equation` runtime perf is non-negotiable` block (lines 354–356): the migration table explicitly keeps it in the slimmed CLAUDE.md (cross-linked from the Developer page) because it is an agent-behavioral guardrail ("YOU MUST STOP and discuss with Denis first"). Append one cross-link sentence to the end of that kept block, after the existing final sentence "This is one of the most important tests in the suite.":
  ```markdown
   See the Developer page in the docs for how `rate_equation` is derived and why the 0-allocation / sub-100-ns contract holds.
  ```

- [ ] **Step 13: Verify the cut removed only the intended blocks and kept every rule**
  Confirm the agent-behavioral rules and the kept blocks survived, and the migrated blocks are gone:
  ```bash
  grep -n "^## \|^### \|^# " /home/denis.linux/.julia/dev/EnzymeRates/.claude/CLAUDE.md
  ```
  Expected headings present (in order): `# CLAUDE.md`, `## Foundational rules`, `## Our relationship`, `# Proactiveness`, `## Designing software`, `## Test Driven Development  (TDD)`, `## Writing code`, `## Naming`, `## Code Comments`, `## Version Control`, `## Testing` (the agent-behavioral one at line 115 — keep it), `## Issue tracking`, `## Systematic Debugging Process` (and its four `### Phase` subheadings), `## Learning and Memory Management`, `# Instructions for this repository`, `## Package Goal`, `## API Design`, `## Commands`, `## Workflow`, `## Code Style`, `### Parameter naming chokepoint (guard)`, `### Canonical Step Form (load-bearing guard)`, `### `rate_equation` runtime perf is non-negotiable`.
  Expected headings ABSENT (migrated): `## Key Architecture Decisions`, `### Regulator representation`, `### Allosteric state taxonomy...`, `### Mirror / dead-end kinetic-group sharing`, `### Parameter naming convention`, `### Catalytic topology constraints`, `### Mechanism enumeration building blocks`, `## Source Layout`, `## Vmax Normalization...`, `### Implementation`, `### Key properties`, `## Known Issues`, the maintainer `## Testing` block at the old line 345.
  Note: there are TWO `## Testing` headings in the original — the agent-behavioral one at line 115 (KEEP) and the maintainer test-fixture one at line 345 (CUT). After this task there must be exactly ONE `## Testing` heading. If `grep -c "^## Testing"` returns anything other than `1`, the wrong one was cut — restore from git and redo.

- [ ] **Step 14: Decide the "Key Architecture Decisions" heading and verify no orphaned content**
  Steps 1–12 left the `## Key Architecture Decisions` heading (CLAUDE.md line 220) with its six bullets (lines 222–228) — those bullets are maintainer-internal (the concrete-vs-singleton split, the `Sig` encoding, the `Parameter` family layout) per the migration table's "concrete-vs-singleton split" and "Source Layout" entries, now on the Developer page. Delete the entire `## Key Architecture Decisions` block from its heading (line 220) through the final bullet `- @generated functions are used for per-mechanism compile-time derivation...` (line 228). Leave nothing in its place. Then run:
  ```bash
  grep -n "Key Architecture\|EnzymeMechanism{Sig}\|_sig_of\|Parameter.*family.*types.jl" /home/denis.linux/.julia/dev/EnzymeRates/.claude/CLAUDE.md
  ```
  Expected: no output (the architecture-decisions block and its singleton/`Sig`/Parameter-family details are fully migrated to the Developer page).

- [ ] **Step 15: Verify the slimmed CLAUDE.md has no dangling pointers to deleted blocks**
  ```bash
  grep -n "test_chokepoint\|regulator_roles\|_is_ss_rate_constant\|_kcat_components\|has_residual" /home/denis.linux/.julia/dev/EnzymeRates/.claude/CLAUDE.md
  ```
  Expected: no output. All five stale identifiers were inside cut blocks; none should survive (the guard line in Step 1 cites the correct `test/test_types.jl:1577-1644`, not `test_chokepoint.jl`). If any line prints, it is inside a block that should have been cut — re-check the corresponding step.

- [ ] **Step 16: Commit**
  ```bash
  git add /home/denis.linux/.julia/dev/EnzymeRates/.claude/CLAUDE.md
  git commit -m "Split CLAUDE.md by audience: cut migrated blocks, leave docs guards

Per the documentation migration table (spec section 7): cut maintainer-internal
and user-concept blocks whose corrected content now lives in the docs site;
leave one-line guard pointers for the load-bearing Canonical Step Form and
parameter-naming chokepoint; keep the rate_equation runtime-perf guardrail and
every agent-behavioral rule. Stale references (regulator_roles, test_chokepoint,
_is_ss_rate_constant, _kcat_components, has_residual) are removed with their
blocks; the kept guard cites the real test path test/test_types.jl:1577-1644.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
  ```

### Task 6.4: Verify — full suite green and docs build green

**Files:** none (verification only)

- [ ] **Step 1: Run the full test suite**
  ```bash
  julia --project=/home/denis.linux/.julia/dev/EnzymeRates -e 'using Pkg; Pkg.test()'
  ```
  Expected: `Testing EnzymeRates tests passed` (the final summary line shows `Test Summary:` with `Pass` counts and zero `Fail`/`Error`). The README testset no longer runs (it was removed in Task 6.2). Aqua, JET, and the compile-budget tests must all pass. If any test fails, STOP and diagnose with the systematic-debugging skill — do not paper over it. A failure here most likely means a stale CLAUDE.md cut accidentally touched a code file (it should not — Task 6.3 edits only `.claude/CLAUDE.md`) or the README quickstart param names diverged (Task 6.1 Step 4).

- [ ] **Step 2: Build the docs**
  ```bash
  julia --project=/home/denis.linux/.julia/dev/EnzymeRates/docs /home/denis.linux/.julia/dev/EnzymeRates/docs/make.jl
  ```
  Expected: `[ Info: Documenter: rendering done` near the end, no `Doctest failed` lines, no `checkdocs` failures, and no error backtrace. Because `make.jl` runs `makedocs` with `doctest = true`, this also re-checks every `jldoctest` block. The slimmed README is not part of the docs build, but the five doc-link slugs in the new README (Task 6.1 Step 3) must resolve to pages this build renders — if the build reports a missing page for any slug the README links to, fix the README URL to the page's real slug.

- [ ] **Step 3: Confirm the working tree is clean and on the docs branch**
  ```bash
  git status --short && git log --oneline -3
  ```
  Expected: empty `git status --short` (everything committed across Tasks 6.1–6.3), and the three Phase 6 commits (README slim, README-test retire, CLAUDE.md split) at the top of the log. If `make.jl` regenerated any doctest output via a `fix` pass, those changes must be committed before the tree is clean — but Step 2 runs with `doctest = true` (check-only, no fix), so a clean tree is the expected state. If the tree is dirty, inspect the diff and commit the captured doctest output with message `Capture doctest output after Phase 6 trim`.

---

## Verification

The plan is complete when every spec exit criterion holds:

- **`make.jl` builds clean locally.** `julia --project=docs docs/make.jl` runs to `[ Info: Documenter: rendering done` and exits 0, with no warnings, no missing-page errors, and (locally) `deploydocs` reporting `Skipping deployment`.
- **The `docs` CI job is green and deploys on merge.** The `docs` job in `.github/workflows/CI.yml` passes as a PR check (PRs verify, never publish) and, on merge to `main`, deploys to gh-pages via `julia-actions/julia-docdeploy@v1` using `DOCUMENTER_KEY` (SSH deploy key, with `GITHUB_TOKEN` fallback); the first deploy is eyeballed after the `docs-comprehensive` PR merges and Denis sets the Pages source to the `gh-pages` branch.
- **`checkdocs = :exports` passes.** Every exported name has an attached docstring — the `EnzymeReaction` and `metabolites` gaps are closed and no export is flagged.
- **Every `jldoctest` passes and every `@example` runs.** All output-checked doctests (in `rate_equation_string`, `parameters`, `metabolites`, `rate_equation`, and the page-level deterministic blocks) match byte-for-byte under `doctest = true`; every `@example` block on every tutorial page executes without error.
- **The fast identify tutorial recovers its mechanism in seconds.** The width-1-beam `identify_rate_equation` example (`min_beam_width = 1`, `loss_rel_threshold = 1.0`, `loss_abs_threshold = 0.0`, low `max_param_count`, serial `pmap_function = map`) runs in seconds on noiseless data and `results.best` recovers the generating mechanism.
- **The full test suite stays green after the trim.** `julia --project -e 'using Pkg; Pkg.test()'` passes (Aqua, JET, compile-budget, chokepoint, and `rate_equation` performance tests included) after `test_readme_runs.jl` is retired and `runtests.jl` is updated — a net coverage increase, since the doctests now check rendered output.