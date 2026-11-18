# CLAUDE.md Rewrite Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Rewrite `.claude/CLAUDE.md` around lowering complexity through deep modules, with a module table generated from `src/` and checked by the boundary test.

**Architecture:** `test/module_map.jl` gains a Markdown table generator beside the Mermaid one; the two generated blocks share one "problems" function and one "replace" function. `test/update_module_map.jl` rewrites both blocks. The prose of `.claude/CLAUDE.md` is then rewritten around the generated table, section by section as the spec lists.

**Tech Stack:** Julia 1.11+, `Test`, the parsing helpers already in `test/module_map.jl`.

**Spec:** `docs/superpowers/specs/2026-11-18-claude-md-rewrite-design.md`

## Global Constraints

- `.claude/CLAUDE.md` is one file under 200 lines; target about 150.
- Calm wording: no capitals for emphasis, no "Rule #1", no "YOU MUST", no "NEVER" in capitals.
- The maintainer is "the maintainer", never named.
- No change under `src/`. No new public name. No new test machinery beyond the table check.
- Every helper under `test/` carries the `_testhelper_` prefix.
- 92-character lines and 4-space indentation in `.jl` files. Markdown table rows may exceed 92.
- The headings `### Derivation tests` and `### Enumeration-engine tests` and their rule numbers stay: comments in `test/Enumeration/test_mechanism_enumeration.jl` and `test/Derivation/test_oracles.jl` cite them.
- These strings stay verbatim in the new file: `50 ns + 0.5 ns × denominator terms`, `1.25 ×`, `8 \`log\` calls`, `20 batches`, and both command blocks under "Commands".
- One Julia process at a time. Full suite in the background under `timeout 1200`.
- Commit messages end with the two attribution lines the session specifies.

## Review Focus

1. A Windows checkout gives `.claude/CLAUDE.md` CRLF line endings; the table check must compare as LF. Test in Task 1.
2. A `CLAUDE.md` with no table markers must yield one problem line, not an exception. Test in Task 1.
3. A module docstring whose first sentence spans two lines must become one table cell on one line. Test in Task 1 (the fixture's `Fitting` docstring).
4. A rule of the current file dropped by accident rather than by the spec's "Removed" list. Task 2, Step 2 walks the current file rule by rule.
5. A perf bound or command altered while compressing. Task 2, Step 5 greps for the verbatim strings.

---

### Task 1: Generated module table, checked against `src/`

**Files:**
- Modify: `test/module_map.jl` (after `_testhelper_module_map`, about lines 223-300)
- Modify: `test/update_module_map.jl`
- Modify: `test/test_boundaries.jl` (testsets "module map diagram" and "the developer docs' module map matches src/")
- Modify: `.claude/CLAUDE.md` (insert the marked block under `### Module boundaries (guard)`)

**Interfaces:**
- Consumes: `_testhelper_api(root, mod)`, `_testhelper_exports(file)`, `_testhelper_type_names(root)`, `_testhelper_map_order(names, types, exports)`, `_testhelper_box_sentence(file, mod)`, `_testhelper_LAYERS`, `_testhelper_MAP_BOXES`, `_testhelper_module_sources(root)`, `_testhelper_module_map(root)`.
- Produces:
  - `_testhelper_module_table(root) -> String`: the Markdown table, ending in `"\n"`.
  - `_testhelper_MAP_BLOCK` and `_testhelper_TABLE_BLOCK`: `NamedTuple`s `(; opening, closing, generate, described)`. `opening` and `closing` are the literal lines around the generated text, `generate` is `_testhelper_module_map` or `_testhelper_module_table`, `described` is the phrase the problem line starts with.
    - Map: `opening = "```mermaid\n"`, `closing = "```\n"`, `described = "the Mermaid module map in docs/src/developer.md"`.
    - Table: `opening = "<!-- module table, generated from src/ by test/update_module_map.jl -->\n"`, `closing = "<!-- end of module table -->\n"`, `described = "the module table in .claude/CLAUDE.md"`.
  - `_testhelper_generated_block_problems(doc, block, root) -> Vector{String}`: empty when the text between `block.opening` and `block.closing` in `doc` (CRLF read as LF) equals `block.generate(root)`; otherwise one line, `"$(block.described) differs from the one generated from src/; rewrite it with \`julia test/update_module_map.jl\`"`. A `doc` without the block gives that same line.
  - `_testhelper_replace_generated_block(doc, block, text) -> String`.
  - These two replace `_testhelper_module_map_block_problems`, `_testhelper_replace_module_map` and `_testhelper_MAP_FENCE`, which are deleted.

Table format. Header `| Module | What it does | Public names (★ exported) | May call |`, separator `|:--|:--|:--|:--|`, then one row per box of `_testhelper_MAP_BOXES` that is the top level or has a folder, in that order. Cells:
- Module: the module name; the top level is `EnzymeRates (top level)`.
- What it does: `_testhelper_box_sentence`, with each `<br/>` replaced by a space.
- Public names: the module's public names in `_testhelper_map_order`, each in backticks, an exported one preceded by `★ `, joined by `, `. The top level lists the exports no module owns.
- May call: `_testhelper_LAYERS[m]` joined by `, `, `none` when empty, `all` for the top level.

- [ ] **Step 1: Write the failing tests**

In `test/test_boundaries.jl`, inside the testset "module map diagram", after `expected` is defined and using the same `tree`:

```julia
    expected_table = join([
        "| Module | What it does | Public names (★ exported) | May call |",
        "|:--|:--|:--|:--|",
        "| EnzymeRates (top level) | Does it all. | ★ `identify`, ★ `Result` | all |",
        "| Fitting | Fits a law to data. | ★ `fit` | Mechanisms, Derivation, Data |",
        "| Mechanisms | Stores things. | ★ `metabolites`, `Mechanism`, " *
        "`compile_mechanism` | none |", ""], "\n")
    _testhelper_with_tree(tree...) do root
        @test _testhelper_module_table(root) == expected_table
        table, map = _testhelper_TABLE_BLOCK, _testhelper_MAP_BLOCK
        doc = "text\n\n$(table.opening)$(expected_table)$(table.closing)\nmore\n"
        @test _testhelper_generated_block_problems(doc, table, root) == String[]
        @test _testhelper_generated_block_problems(
            replace(doc, "\n" => "\r\n"), table, root) == String[]
        stale = replace(doc, "Stores things." => "Stores.")
        @test _testhelper_generated_block_problems(stale, table, root) == [
            "the module table in .claude/CLAUDE.md differs from the one generated " *
            "from src/; rewrite it with `julia test/update_module_map.jl`"]
        @test _testhelper_replace_generated_block(stale, table, expected_table) == doc
        @test !isempty(_testhelper_generated_block_problems("no table", table, root))
        @test !isempty(_testhelper_generated_block_problems("no diagram", map, root))
    end
```

In the same testset, change the existing Mermaid assertions to the shared functions: `_testhelper_module_map_block_problems(x, root)` becomes `_testhelper_generated_block_problems(x, _testhelper_MAP_BLOCK, root)`, and `_testhelper_replace_module_map(stale, expected)` becomes `_testhelper_replace_generated_block(stale, _testhelper_MAP_BLOCK, expected)`.

Replace the testset "the developer docs' module map matches src/" with:

```julia
@testset "the generated module map and table match src/" begin
    root = joinpath(@__DIR__, "..")
    doc = read(joinpath(root, "docs", "src", "developer.md"), String)
    @test _testhelper_module_map_problems(doc, _testhelper_SRC) == String[]
    @test _testhelper_generated_block_problems(
        doc, _testhelper_MAP_BLOCK, _testhelper_SRC) == String[]
    @test _testhelper_generated_block_problems(
        read(joinpath(root, ".claude", "CLAUDE.md"), String), _testhelper_TABLE_BLOCK,
        _testhelper_SRC) == String[]
end
```

- [ ] **Step 2: Run the tests and confirm they fail**

Run: `julia --project -e 'using TestEnv; TestEnv.activate(); using Test, EnzymeRates; include("test/test_boundaries.jl")'`
Expected: an `UndefVarError` for `_testhelper_module_table` in "module map diagram".

- [ ] **Step 3: Implement in `test/module_map.jl`**

Add `_testhelper_module_table(root)`, the two block constants, `_testhelper_generated_block_problems` and `_testhelper_replace_generated_block` with the signatures above. Delete `_testhelper_MAP_FENCE`, `_testhelper_module_map_block_problems` and `_testhelper_replace_module_map`. Find a block by its literal `opening` and `closing` (the first `closing` after the `opening`), not by a regex built from them. Update the file's `ABOUTME` lines to name both generated blocks.

- [ ] **Step 4: Rewrite `test/update_module_map.jl`**

It rewrites the Mermaid block of `docs/src/developer.md` and the table block of `.claude/CLAUDE.md`, each with `_testhelper_replace_generated_block`. Update its `ABOUTME` lines.

- [ ] **Step 5: Insert the block into `.claude/CLAUDE.md` and generate it**

Under `### Module boundaries (guard)`, after its first paragraph, add the `opening` line, an empty table position and the `closing` line, each followed by a blank line as Markdown needs. Then run `julia test/update_module_map.jl`.
Expected: the block holds seven rows (top level, DSL, Enumeration, Fitting, Derivation, Data, Mechanisms) and `git diff docs/src/developer.md` is empty.

- [ ] **Step 6: Run the tests and confirm they pass**

Run the command of Step 2.
Expected: every testset passes, with no warnings in the output.

- [ ] **Step 7: Check that the guard fails when it should**

Add `, zzz_probe` to the `public` statement in `src/Data/Data.jl`, run `julia --project -e 'include("test/module_map.jl"); foreach(println, _testhelper_generated_block_problems(read(".claude/CLAUDE.md", String), _testhelper_TABLE_BLOCK, _testhelper_SRC))'`.
Expected: the "module table in .claude/CLAUDE.md differs" line. Then `git checkout src/Data/Data.jl` and confirm `git status --short src/` is empty.

- [ ] **Step 8: Update the developer docs' sentence about the generator**

In `docs/src/developer.md`, sections "Modules" and "The boundary test": where they say `julia test/update_module_map.jl` rewrites the diagram, say it also rewrites the module table of `.claude/CLAUDE.md`, and that the boundary test checks both.

- [ ] **Step 9: Commit**

```bash
git add test/module_map.jl test/update_module_map.jl test/test_boundaries.jl \
    .claude/CLAUDE.md docs/src/developer.md
git commit -m "Generate a module table in CLAUDE.md and check it against src/"
```

---

### Task 2: Rewrite the prose of `.claude/CLAUDE.md`

**Files:**
- Modify: `.claude/CLAUDE.md` (whole file except the generated block)
- Modify: `docs/src/developer.md` only where Step 6 finds a quoted rule that changed

**Interfaces:**
- Consumes: the generated block from Task 1, kept byte for byte between its two marker lines.
- Produces: the file the spec's "Content of the new file, by section" describes: eleven sections in the spec's order, with the spec's content for each.

- [ ] **Step 1: Save the current file for comparison**

Run: `git show HEAD:.claude/CLAUDE.md > "$SCRATCH/claude-md-before.md"`, with `$SCRATCH` the session scratchpad directory.

- [ ] **Step 2: Account for every rule of the current file**

Walk the saved file top to bottom. For each rule, write one line in `$SCRATCH/rule-ledger.md`: the rule in a few words, then `kept in section N`, `merged into section N` or `removed: <the item of the spec's "Removed from the current file" list>`. A rule that fits none of the three is a gap: keep it in the nearest section and note it for the report.

- [ ] **Step 3: Write the new file**

Sections, in order, with these headings:

1. `# EnzymeRates: instructions for agents`, then `## Goal: low complexity`
2. `## Module map` (the generated block, then the facts the spec lists for it)
3. `## Design rules`
4. `## Changing an interface`
5. `## Playbooks` with `### Feature`, `### Refactor`, `### Bug fix`, `### Optimization`
6. `## Julia style`
7. `## Names, comments and docstrings`
8. `## Guards`, closing with the paragraph on when to propose a guard
9. `## Tests` with `### Derivation tests` and `### Enumeration-engine tests`, rules numbered as now
10. `## Commands and workflow`
11. `## Working with the maintainer`

Write each section from the spec's entry for it. Prose rules: one idea per sentence, the reason beside the rule where the reason is not obvious, no rule stated twice. The Guards section gives each guard at most three lines and names its test file.

- [ ] **Step 4: Check length and banned wording**

Run: `wc -l .claude/CLAUDE.md` (expected: under 200) and
`grep -nE "Denis|Rule #1|YOU MUST|NEVER|ALWAYS|Circle K|replaced\.|sycophant|glaze" .claude/CLAUDE.md` (expected: no output).

- [ ] **Step 5: Check the strings that must survive**

Run: `grep -c -F -e '50 ns + 0.5 ns × denominator terms' -e '1.25 ×' -e '8 `log` calls' -e '20 batches' -e '### Derivation tests' -e '### Enumeration-engine tests' -e 'timeout 1200 julia --project' -e 'TestEnv.activate()' .claude/CLAUDE.md`
Expected: at least 8 matching lines; then confirm each of the eight patterns matches by grepping it alone.

- [ ] **Step 6: Find quotes of changed rules elsewhere**

Run: `grep -rn -i "CLAUDE.md\|inline single-use\|YAGNI" docs/src ROADMAP.md test --include=*.md --include=*.jl`
Expected: the references to the full-suite command in `ROADMAP.md` and to the enumeration and derivation rule headings in `test/`, all still true. Fix any line that quotes a rule this rewrite changed.

- [ ] **Step 7: Run the boundary test**

Run the command of Task 1, Step 2.
Expected: every testset passes; the table block is intact.

- [ ] **Step 8: Read the new file as a newcomer**

Read the whole file once, top to bottom, asking of each line whether an agent would act differently without it, and of each pair of rules whether they can conflict. Delete a line that changes nothing; reword a pair that conflicts. Rerun Steps 4, 5 and 7 after any edit.

- [ ] **Step 9: Commit**

```bash
git add .claude/CLAUDE.md docs/src/developer.md
git commit -m "Rewrite CLAUDE.md around low complexity and deep modules"
```

---

### Task 3: Verify the branch and hand it over

**Files:**
- Delete: `docs/superpowers/specs/2026-11-18-claude-md-rewrite-design.md`
- Delete: `docs/superpowers/plans/2026-11-18-claude-md-rewrite.md`

- [ ] **Step 1: Run the full suite**

Run in the background: `timeout 1200 julia --project -e 'using Pkg; Pkg.test(julia_args=["--heap-size-hint=2500M"])'`
Expected: all tests pass and the output has no warnings. Wait for the exit notice and read the result.

- [ ] **Step 2: Report to the maintainer**

Give the line count before and after, the ledger's removed rules, any gap Step 2 of Task 2 found, and the full-suite result. Ask whether to open a pull request.

- [ ] **Step 3: Delete the spec and this plan**

When the maintainer asks for the pull request, in the same branch:

```bash
git rm docs/superpowers/specs/2026-11-18-claude-md-rewrite-design.md \
    docs/superpowers/plans/2026-11-18-claude-md-rewrite.md
git commit -m "Remove the CLAUDE.md rewrite's spec and plan"
```
