# Concrete-Types Refactor Audit — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use
> `superpowers:subagent-driven-development` (recommended) or
> `superpowers:executing-plans` to implement this plan task-by-task.
> Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Execute the 3-pass audit defined in
`docs/superpowers/specs/2026-05-30-refactor-audit-workflow-design.md`
and produce the findings document at
`docs/superpowers/2026-05-30-refactor-audit-findings.md`. The audit
reports every simplification opportunity the concrete-types and
structural-parameter-names refactors enabled but did not complete,
across all 9 src files plus blocking tests, against a dual goal:
(1) reduce non-comment non-doc src LOC, (2) simplify the code.

**Architecture:** Serial **Pass 1** catalog (the executor reads every
src file top-to-bottom and captures suspects in a gitignored scratch
file). Parallel **Pass 2** verification (Agent-tool dispatches in
parallel, one query per agent, results merged back to the scratch
file). Serial **Pass 3** synthesis (cluster findings, build a
dependency graph, identify blocking tests, sequence, write the
committed findings doc).

**Tech Stack:** Julia 1.x package; bash for git ops; `grep` / `find`
for verification scaffolding; Agent tool (subagent_type=`Explore` for
read-only verifies, `general-purpose` for cross-file comparisons) for
parallel Pass-2 dispatches; Markdown for the scratch + findings docs.

---

## Reference: Suspect-Criteria Table (applies to every Pass-1 task)

Every Pass-1 task in this plan applies the same suspect criteria.
Quoted from the workflow spec — do not modify in-place; if you need to
extend, do it via the spec, not here.

| Pattern | Suspect criterion |
|---|---|
| **Dead code** | Non-exported symbol; visible call sites only in other suspect code, in removed-style branches, or in test-only "internal helper" assertions |
| **Duplication** | Same algorithm, regex, dispatch, or formula appears in 2+ functions or files |
| **Symbol-tuple plumbing** | Code converts Sig → opaque Symbol tuple → back, OR regenerates structure (substrate list, product list, kinetic-group rep, bound metabolite) that `Step` / `Species` / `Mechanism` already carries directly |
| **Compile-time accessor** | `@generated` walking `Sig` where plain `Mechanism` field access would work, AND not on the `rate_equation` hot path (verified by reading the `rate_equation` `@generated` body; the accessor perf test `test/test_accessors.jl` is explicitly negotiable) |
| **Test-private helper** | A `_`-prefixed helper that has a direct test assertion against its return value (the test constrains the helper's signature, blocking refactor) |
| **String-keyed projection** | `name_map` and related — projections built from rendered Symbol strings that structural parameter keys obsoleted |
| **Permissive parser + post-hoc guard** | A parser branch that accepts a wider input language than intended, paired with a later validator rejecting the slack. Example: `_assert_no_opaque_terms` over `dsl.jl:259-262` |
| **Stale spec/stage comment** *(doc)* | Comment mentions "legacy", "old Sig", "previous path", "deprecated", "OLD:", "moved", "renamed from", "from Stage Nx", "see YYYY-MM-DD-…", "Phase N", "per <past spec doc>". Flag for removal; does **not** reduce non-comment non-doc LOC |
| **Comment used as docstring** *(doc)* | Function or struct whose explanation lives in a block of `#`-comments above it rather than a `"""docstring"""` on the definition. Flag for conversion; does **not** reduce non-comment non-doc LOC |

**On-encounter confidence:** `H` (unambiguous from local reading), `M`
(likely, verification needed), `L` (hunch).

**Scratch row format** (one bullet per suspect):

```
- L<start>-L<end> | <category> | <one-line summary> | <H/M/L>
```

Use `intentionally-not-flagged` as the category and `-` as the
confidence for regions you intentionally skip (with a reason in the
summary). Every line of every file must end up either flagged as a
suspect or marked intentionally-not-flagged.

---

## Phase 0 — Setup

### Task 0.1: Gitignore the scratch file and create its scaffold

**Files:**
- Modify: `.gitignore`
- Create: `docs/superpowers/scratch-refactor-audit-notes.md` (NOT committed)

- [ ] **Step 1: Append the scratch path to .gitignore**

Read `.gitignore` and append a new line block at the end. The current
file ends with `2026-05-02_identify_ldh_results/`. Edit it to add:

```
# Refactor audit scratch (Pass 1 catalog; deleted at end of audit)
docs/superpowers/scratch-refactor-audit-notes.md
```

- [ ] **Step 2: Verify the gitignore entry is honored**

Run:
```bash
git check-ignore docs/superpowers/scratch-refactor-audit-notes.md
```
Expected output:
```
docs/superpowers/scratch-refactor-audit-notes.md
```
(The path printed = it is ignored. If empty, the gitignore entry is
wrong.)

- [ ] **Step 3: Commit the gitignore change**

```bash
git add .gitignore
git commit -m "$(cat <<'EOF'
chore: gitignore refactor-audit scratch notes file

Adds docs/superpowers/scratch-refactor-audit-notes.md (the Pass-1
catalog working file for the refactor audit) to .gitignore. The
scratch file is deleted at the end of the audit; only the findings
doc is committed.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

- [ ] **Step 4: Create the scratch file with section scaffolding**

Use Write to create `docs/superpowers/scratch-refactor-audit-notes.md`
with this exact initial content (one section per src file in the
EnzymeRates.jl inclusion order, plus a Pass-2 / Pass-3 staging
section):

```markdown
# Refactor Audit — Pass 1 Scratch (uncommitted)

This file is gitignored. It is the running Pass-1 catalog and is
deleted at the end of the audit; the committed deliverable is
`docs/superpowers/2026-05-30-refactor-audit-findings.md`.

Format for each suspect entry: `- L<start>-L<end> | <category> | <one-line summary> | <H/M/L>`
See `docs/superpowers/specs/2026-05-30-refactor-audit-workflow-design.md` for criteria.

## src/types.jl  (1,550 LOC)
<!-- Pass 1 suspects below. Every line in [1, 1550] must end up either
     flagged or in an `intentionally-not-flagged` entry with reason. -->

## src/dsl.jl  (1,138 LOC)

## src/sym_poly_for_rate_eq_derivation.jl  (315 LOC)

## src/rate_eq_derivation.jl  (1,706 LOC)

## src/thermodynamic_constr_for_rate_eq_derivation.jl  (417 LOC)

## src/fitting.jl  (211 LOC)

## src/mechanism_enumeration.jl  (2,196 LOC)

## src/identify_rate_equation.jl  (886 LOC)

## src/EnzymeRates.jl  (37 LOC)

## Pass 2 — Verification queries
<!-- Filled in Phase 2. Each entry: `Q-NNN | <query template> | <suspects covered>`. -->

## Pass 2 — Verification results
<!-- Filled in Phase 2. Each entry: `Q-NNN result: <summary>`. -->

## Pass 3 — Findings (promoted suspects)
<!-- Filled in Phase 3. Each entry: `F-NNN | <suspect-id> | <category> | <H/M/L>`.
     Full finding bodies go into the committed findings doc, not here. -->

## Pass 3 — Drops (verification failed)
<!-- Filled in Phase 3. Each entry: `<suspect-id> dropped: <one-line reason>`. -->
```

- [ ] **Step 5: Verify the scratch file is not staged**

Run:
```bash
git status --short docs/superpowers/scratch-refactor-audit-notes.md
```
Expected output: empty (no `??` line). If `??` appears, the gitignore
isn't catching the path — fix the .gitignore entry before continuing.

### Task 0.2: Measure the non-comment non-doc src LOC baseline

**Files:**
- Modify: `docs/superpowers/scratch-refactor-audit-notes.md` (record
  baseline near the top)

- [ ] **Step 1: Inspect a representative Julia file to confirm the
      comment / docstring patterns**

Read the first 50 lines of `src/types.jl` to confirm:
  - Single-line comments start with `#` (with or without leading
    whitespace)
  - Docstrings are triple-quoted blocks `"""…"""` attached to
    definitions (multi-line possible)
  - The `ABOUTME:` line is a `#` comment

- [ ] **Step 2: Run the baseline counter**

Use this command (counts lines that are NOT blank, NOT a `#` comment,
and NOT inside a `"""…"""` docstring block):

```bash
for f in src/*.jl; do
    total=$(wc -l < "$f")
    nccnd=$(awk '
        BEGIN { in_ds=0 }
        {
            line=$0
            # detect docstring boundaries
            n = gsub(/"""/, "\"\"\"", line)
            stripped = line
            sub(/^[[:space:]]+/, "", stripped)
            sub(/[[:space:]]+$/, "", stripped)
            if (in_ds) {
                if (n % 2 == 1) in_ds = 0
                next
            }
            if (n % 2 == 1) {
                # docstring opens on this line; only count if there is
                # also code on the line BEFORE the opening """ — for
                # simplicity, treat the whole line as docstring
                in_ds = 1
                next
            }
            if (stripped == "") next
            if (substr(stripped, 1, 1) == "#") next
            count++
        }
        END { print count+0 }
    ' "$f")
    printf "%6d  %6d  %s\n" "$total" "$nccnd" "$f"
done | tee /tmp/refactor_audit_baseline.txt
echo "---"
awk '{ tot+=$1; ncc+=$2 } END { printf "TOTAL: %d total LOC, %d non-comment non-doc LOC\n", tot, ncc }' /tmp/refactor_audit_baseline.txt
```

Expected output: a 9-line table (one row per src file) plus a TOTAL
line. The non-comment non-doc total is the baseline number.

- [ ] **Step 3: Record the baseline in the scratch file**

Use Edit to insert a line near the top of
`docs/superpowers/scratch-refactor-audit-notes.md` (right after the
opening paragraph), inserting:

```
**Baseline (non-comment non-doc src LOC):** <number from Step 2>
**Computed on:** 2026-05-30
```

Replace `<number from Step 2>` with the actual TOTAL number from the
command output.

---

## Phase 1 — Pass 1: Catalog suspects (serial, by Claude)

**For every task in Phase 1**, the procedure is identical:

1. Read the entire src file in one Read call (no `offset`/`limit`)
2. Walk the file top-to-bottom; for each function / method / struct /
   macro, evaluate against the suspect-criteria table at the top of
   this plan
3. For each suspect, append a bullet to the matching `## src/<file>`
   section of the scratch file using the format
   `- L<start>-L<end> | <category> | <one-line summary> | <H/M/L>`
4. For every line range NOT flagged as a suspect, append an entry of
   the form `- L<start>-L<end> | intentionally-not-flagged | <one-line reason> | -`
5. Confirm the union of `L<start>-L<end>` ranges (suspect +
   intentionally-not-flagged) covers `[1, file_LOC]` exactly. If a
   gap exists, re-read the gap region and add an entry.

Do NOT verify any suspect in Pass 1 — verification is Pass 2's job.
The goal here is breadth, not precision.

### Task 1.1: Catalog src/types.jl

**Files:**
- Read: `src/types.jl` (1,550 LOC)
- Modify: `docs/superpowers/scratch-refactor-audit-notes.md` (the
  `## src/types.jl` section)

- [ ] **Step 1: Read src/types.jl in full**

Use the Read tool on `/home/denis.linux/.julia/dev/EnzymeRates/src/types.jl`
with no offset / limit. The tool returns lines 1-2000 by default; this
file is 1,550 LOC, so one Read call is enough.

- [ ] **Step 2: Walk top-to-bottom and flag suspects**

For each function / struct / `@generated` body / macro you encounter,
apply the suspect-criteria table at the top of this plan. Areas to
look at especially carefully given the spec's "Current state" notes:

- **Mechanism ↔ Sig conversion machinery** (lines ~500-816 per the
  prior-agent audit): is the singleton-type bridge collapsible?
- **Generated accessors over Sig** (lines ~974-1321 per prior-agent
  audit): which are on `rate_equation`'s hot path and which could be
  plain Mechanism field access?
- **Parameter family** (lines ~195-219 per CLAUDE.md): are there
  Parameter types that no longer have callers, or constructors that
  duplicate?
- **`_species_name_from_sig`** (around line 1408 per prior memory):
  symbol-tuple plumbing — confirm and capture.
- **`_step_tuple_from_sig`** if present in this file (or
  cross-reference if it lives in rate_eq_derivation.jl).

- [ ] **Step 3: Append suspects to the scratch file**

Append every suspect bullet under the `## src/types.jl  (1,550 LOC)`
header. Example bullets:

```
- L500-L520 | symbol-tuple-plumbing | `_species_name_from_sig` rebuilds opaque form-name via string-join; the Mechanism value carries the structure directly | H
- L974-L1010 | compile-time-accessor | `@generated substrates(::Type{EnzymeMechanism{Sig}})` walks Sig where `m.reaction.reactants` field access works; verify hot-path | M
- L1408-L1430 | stale-spec-comment | comment `# Phase 2 — see 2026-05-26-phase2-enumerator-decomposed-species.md` | H
```

- [ ] **Step 4: Add intentionally-not-flagged entries for skipped
       ranges**

For any contiguous block of lines you read but did not flag as a
suspect (e.g. simple field accessors, exports, public docstrings,
load-bearing constructors), append an entry of the form:

```
- L1-L50 | intentionally-not-flagged | module preamble + imports | -
- L100-L194 | intentionally-not-flagged | public type aliases + RateEquationMode hierarchy — production surface | -
```

These ranges plus the suspect ranges must cover `[1, 1550]` exactly.

- [ ] **Step 5: Run the meta-completeness check for types.jl**

Use bash + awk to verify the ranges cover the file exactly. Save this
check as a reusable script — you will use it for every Pass-1 task.

```bash
python3 - <<'PY'
import re, sys
with open("docs/superpowers/scratch-refactor-audit-notes.md") as fh:
    text = fh.read()
section = re.search(r"## src/types\.jl[^\n]*\n(.*?)(?=\n## |\Z)", text, re.DOTALL)
if not section:
    sys.exit("types.jl section not found")
body = section.group(1)
ranges = []
for line in body.splitlines():
    m = re.match(r"\s*-\s*L(\d+)-L(\d+)\s*\|", line)
    if m:
        ranges.append((int(m.group(1)), int(m.group(2))))
ranges.sort()
total_lines = 1550
covered = [False] * (total_lines + 1)
for lo, hi in ranges:
    for i in range(lo, hi+1):
        covered[i] = True
gaps = []
in_gap = False
gap_start = None
for i in range(1, total_lines+1):
    if not covered[i] and not in_gap:
        in_gap = True; gap_start = i
    elif covered[i] and in_gap:
        in_gap = False
        gaps.append((gap_start, i-1))
if in_gap:
    gaps.append((gap_start, total_lines))
print(f"types.jl: {len(ranges)} ranges, {len(gaps)} gaps")
if gaps:
    for g in gaps:
        print(f"  GAP: L{g[0]}-L{g[1]}")
    sys.exit(1)
print("OK: coverage complete")
PY
```

Expected output: `types.jl: <N> ranges, 0 gaps\nOK: coverage complete`.
If gaps reported, re-read the gap region and add a suspect or
`intentionally-not-flagged` entry to close it, then re-run the check.

### Task 1.2: Catalog src/dsl.jl

**Files:**
- Read: `src/dsl.jl` (1,138 LOC)
- Modify: `docs/superpowers/scratch-refactor-audit-notes.md` (the
  `## src/dsl.jl` section)

Key areas to scrutinize for this file:

- **`_term_bare_enzyme` / `_step_side_term_info`** (lines ~258-280) and
  **`_assert_no_opaque_terms`** (lines ~382-399): the permissive
  parser + post-hoc guard pattern the spec flags. Capture the
  redundancy.
- **`_call_form_term_info`** and the `.sym` synthesis branch (prior
  audit referenced lines 369-378): is the dead `:call`-branch
  Symbol-synthesis tail still alive?
- **Two callers of `_parse_steps_block_with_groups`** (lines ~618,
  ~1204 per prior audit): are both callers symmetric? Is the
  `side_terms_per_step` post-processing duplicated?
- **`_synthesize_species_name`** (around line 389 per prior audit): if
  dead, capture.
- **Stale spec/stage comments**: this file is large and has been edited
  across many phases; expect many doc-category flags.
- The macro emission tail (`EnzymeMechanism(Mechanism(...))` etc.):
  any branches for the deleted legacy 2-arg ctor?

- [ ] **Step 1: Read src/dsl.jl in full**

Use the Read tool on `/home/denis.linux/.julia/dev/EnzymeRates/src/dsl.jl`
with no offset / limit (1,138 LOC fits one Read call).

- [ ] **Step 2: Walk top-to-bottom and flag suspects per the criteria
       table at the top of this plan**

For each function / struct / macro you encounter, evaluate against
every row of the suspect-criteria table. Treat the "Key areas"
bullets above as priorities, but do not skip other parts of the file.

- [ ] **Step 3: Append suspect bullets to the
       `## src/dsl.jl  (1,138 LOC)` section of the scratch file**

Use the format
`- L<start>-L<end> | <category> | <one-line summary> | <H/M/L>`.

- [ ] **Step 4: Append `intentionally-not-flagged` entries for
       contiguous skipped ranges, covering `[1, 1138]` exactly**

Use the format
`- L<start>-L<end> | intentionally-not-flagged | <one-line reason> | -`.

- [ ] **Step 5: Run the meta-completeness check for dsl.jl**

```bash
python3 - <<'PY'
import re, sys
FILE = "src/dsl.jl"; NLINES = 1138
with open("docs/superpowers/scratch-refactor-audit-notes.md") as fh:
    t = fh.read()
sec = re.search(rf"## {re.escape(FILE)}[^\n]*\n(.*?)(?=\n## |\Z)",
                t, re.DOTALL)
assert sec, f"section not found: {FILE}"
rs = [(int(a), int(b)) for a, b in
      re.findall(r"\s*-\s*L(\d+)-L(\d+)\s*\|", sec.group(1))]
cov = [False] * (NLINES + 1)
for lo, hi in rs:
    for i in range(max(1, lo), min(NLINES, hi) + 1):
        cov[i] = True
gaps = []; ig = False; gs = None
for i in range(1, NLINES + 1):
    if not cov[i] and not ig:
        ig = True; gs = i
    elif cov[i] and ig:
        ig = False; gaps.append((gs, i - 1))
if ig:
    gaps.append((gs, NLINES))
print(f"{FILE}: {len(rs)} ranges, {'OK' if not gaps else 'FAIL'}")
for g in gaps:
    print(f"  GAP: L{g[0]}-L{g[1]}")
sys.exit(1 if gaps else 0)
PY
```

Expected output: `src/dsl.jl: <N> ranges, OK`. If `FAIL`, re-read the
gap region and add suspect or intentionally-not-flagged entries to
close it, then re-run.

### Task 1.3: Catalog src/sym_poly_for_rate_eq_derivation.jl

**Files:**
- Read: `src/sym_poly_for_rate_eq_derivation.jl` (315 LOC)
- Modify: scratch `## src/sym_poly_for_rate_eq_derivation.jl` section

Key areas to scrutinize for this file:

- `POLY` / `MONO` algebra: any redundant constructors / accessors?
- `_rename_symbols` and `_zero_symbols_in_poly`: the spec mentions
  these as targets for the A/I rename collapse. Is one a special case
  of the other? Could both fold into a single mapping API?
- Any helpers used only by deleted callers?

- [ ] **Step 1: Read in full**

Use Read on
`/home/denis.linux/.julia/dev/EnzymeRates/src/sym_poly_for_rate_eq_derivation.jl`
with no offset / limit.

- [ ] **Step 2: Walk top-to-bottom and flag suspects per the
       criteria table**

- [ ] **Step 3: Append suspect bullets to the
       `## src/sym_poly_for_rate_eq_derivation.jl  (315 LOC)` section**

- [ ] **Step 4: Append intentionally-not-flagged entries covering
       `[1, 315]` exactly**

- [ ] **Step 5: Run the meta-completeness check for
       sym_poly_for_rate_eq_derivation.jl**

```bash
python3 - <<'PY'
import re, sys
FILE = "src/sym_poly_for_rate_eq_derivation.jl"; NLINES = 315
with open("docs/superpowers/scratch-refactor-audit-notes.md") as fh:
    t = fh.read()
sec = re.search(rf"## {re.escape(FILE)}[^\n]*\n(.*?)(?=\n## |\Z)",
                t, re.DOTALL)
assert sec, f"section not found: {FILE}"
rs = [(int(a), int(b)) for a, b in
      re.findall(r"\s*-\s*L(\d+)-L(\d+)\s*\|", sec.group(1))]
cov = [False] * (NLINES + 1)
for lo, hi in rs:
    for i in range(max(1, lo), min(NLINES, hi) + 1):
        cov[i] = True
gaps = []; ig = False; gs = None
for i in range(1, NLINES + 1):
    if not cov[i] and not ig:
        ig = True; gs = i
    elif cov[i] and ig:
        ig = False; gaps.append((gs, i - 1))
if ig:
    gaps.append((gs, NLINES))
print(f"{FILE}: {len(rs)} ranges, {'OK' if not gaps else 'FAIL'}")
for g in gaps:
    print(f"  GAP: L{g[0]}-L{g[1]}")
sys.exit(1 if gaps else 0)
PY
```

Expected: `src/sym_poly_for_rate_eq_derivation.jl: <N> ranges, OK`.

### Task 1.4: Catalog src/rate_eq_derivation.jl

**Files:**
- Read: `src/rate_eq_derivation.jl` (1,706 LOC) — large; may need
  multiple Read calls if the default 2000-line limit splits it (it
  shouldn't, but verify)
- Modify: scratch `## src/rate_eq_derivation.jl` section

This is one of the highest-leverage files. Per the spec, the
derivation back-end opaque-tuple round-trip is IN SCOPE. Areas:

- **`_raw_symbolic_rate_polys`** and the cleanup note around line 324
  (prior-agent audit): symbol-tuple regeneration.
- **King-Altman / Cha rate derivation** consumer of `reactions(m)` /
  `enzyme_forms(m)` Symbol identity matching (prior audit cited lines
  142-143, 445).
- **AllostericEnzymeMechanism MWC assembly** (`_build_allosteric_rate_body`,
  `rate_equation_string`): duplication with `_kcat_forward`,
  `_kcat_components`?
- **A/I rename helpers** `_onlyA_parameters`, `_I_rename_parameters`,
  `_all_i_state_parameters`, `_T_rename`, `_build_kinetic_rename_map`
  (lines ~992-1202, ~1204-1373, ~1423-1613 per prior audit): can these
  collapse into one structural transform?
- **Per-state policy logic**: catalog every place the A/I taxonomy is
  branched on. These should be one helper.
- **`@generated rate_equation`**: this is the SACRED hot path. Capture
  its line range as `intentionally-not-flagged` with reason "hot path,
  perf-bound" — do not flag for refactor unless the refactor preserves
  0-alloc / <100 ns.
- Any helper that is referenced only by tests in
  `test/test_rate_eq_derivation.jl` (test-private helper).

- [ ] **Step 1: Read in full**

Use Read on
`/home/denis.linux/.julia/dev/EnzymeRates/src/rate_eq_derivation.jl`
with no offset / limit. The file is 1,706 lines; the Read tool's
default 2000-line cap accommodates this.

- [ ] **Step 2: Walk top-to-bottom and flag suspects per the
       criteria table at the top of this plan**

Pay special attention to the "Key areas" list above. The `@generated
rate_equation` body is the SACRED hot path — mark its line range as
`intentionally-not-flagged` with reason "hot path, perf-bound".

- [ ] **Step 3: Append suspect bullets to the
       `## src/rate_eq_derivation.jl  (1,706 LOC)` section**

- [ ] **Step 4: Append intentionally-not-flagged entries covering
       `[1, 1706]` exactly**

- [ ] **Step 5: Run the meta-completeness check for
       rate_eq_derivation.jl**

```bash
python3 - <<'PY'
import re, sys
FILE = "src/rate_eq_derivation.jl"; NLINES = 1706
with open("docs/superpowers/scratch-refactor-audit-notes.md") as fh:
    t = fh.read()
sec = re.search(rf"## {re.escape(FILE)}[^\n]*\n(.*?)(?=\n## |\Z)",
                t, re.DOTALL)
assert sec, f"section not found: {FILE}"
rs = [(int(a), int(b)) for a, b in
      re.findall(r"\s*-\s*L(\d+)-L(\d+)\s*\|", sec.group(1))]
cov = [False] * (NLINES + 1)
for lo, hi in rs:
    for i in range(max(1, lo), min(NLINES, hi) + 1):
        cov[i] = True
gaps = []; ig = False; gs = None
for i in range(1, NLINES + 1):
    if not cov[i] and not ig:
        ig = True; gs = i
    elif cov[i] and ig:
        ig = False; gaps.append((gs, i - 1))
if ig:
    gaps.append((gs, NLINES))
print(f"{FILE}: {len(rs)} ranges, {'OK' if not gaps else 'FAIL'}")
for g in gaps:
    print(f"  GAP: L{g[0]}-L{g[1]}")
sys.exit(1 if gaps else 0)
PY
```

Expected: `src/rate_eq_derivation.jl: <N> ranges, OK`.

### Task 1.5: Catalog src/thermodynamic_constr_for_rate_eq_derivation.jl

**Files:**
- Read: `src/thermodynamic_constr_for_rate_eq_derivation.jl` (417 LOC)
- Modify: scratch
  `## src/thermodynamic_constr_for_rate_eq_derivation.jl` section

Areas:
- **`_thermodynamic_constraints`** and **`_dependent_param_exprs_kernel`**
  (prior audit cited lines 88-90, 207-209): consume `reactions(m)` /
  `enzyme_forms(m)` via Symbol identity — the derivation-back-end
  opaque-tuple cluster.
- Synth-dep machinery (`_dep_inactive_name`, `_add_case_b_renames!`):
  any duplicate of A/I rename in rate_eq_derivation.jl?
- Kinetic-group merge map logic: structural representation that
  Mechanism already carries?

- [ ] **Step 1: Read in full**

Use Read on
`/home/denis.linux/.julia/dev/EnzymeRates/src/thermodynamic_constr_for_rate_eq_derivation.jl`
with no offset / limit.

- [ ] **Step 2: Walk top-to-bottom and flag suspects per the
       criteria table**

- [ ] **Step 3: Append suspect bullets to the
       `## src/thermodynamic_constr_for_rate_eq_derivation.jl  (417 LOC)`
       section**

- [ ] **Step 4: Append intentionally-not-flagged entries covering
       `[1, 417]` exactly**

- [ ] **Step 5: Run the meta-completeness check for
       thermodynamic_constr_for_rate_eq_derivation.jl**

```bash
python3 - <<'PY'
import re, sys
FILE = "src/thermodynamic_constr_for_rate_eq_derivation.jl"
NLINES = 417
with open("docs/superpowers/scratch-refactor-audit-notes.md") as fh:
    t = fh.read()
sec = re.search(rf"## {re.escape(FILE)}[^\n]*\n(.*?)(?=\n## |\Z)",
                t, re.DOTALL)
assert sec, f"section not found: {FILE}"
rs = [(int(a), int(b)) for a, b in
      re.findall(r"\s*-\s*L(\d+)-L(\d+)\s*\|", sec.group(1))]
cov = [False] * (NLINES + 1)
for lo, hi in rs:
    for i in range(max(1, lo), min(NLINES, hi) + 1):
        cov[i] = True
gaps = []; ig = False; gs = None
for i in range(1, NLINES + 1):
    if not cov[i] and not ig:
        ig = True; gs = i
    elif cov[i] and ig:
        ig = False; gaps.append((gs, i - 1))
if ig:
    gaps.append((gs, NLINES))
print(f"{FILE}: {len(rs)} ranges, {'OK' if not gaps else 'FAIL'}")
for g in gaps:
    print(f"  GAP: L{g[0]}-L{g[1]}")
sys.exit(1 if gaps else 0)
PY
```

Expected: `src/thermodynamic_constr_for_rate_eq_derivation.jl: <N> ranges, OK`.

### Task 1.6: Catalog src/fitting.jl

**Files:**
- Read: `src/fitting.jl` (211 LOC)
- Modify: scratch `## src/fitting.jl` section

Small file, mostly Optimization.jl wrapping. Areas:
- `FittingProblem` constructor: any branches for removed shapes?
- `loss!` and `fit_rate_equation`: any dead config knobs?
- Parameter-name handling at the fitter boundary: does it touch
  `name_map` or structural keys?

- [ ] **Step 1: Read in full**

Use Read on `/home/denis.linux/.julia/dev/EnzymeRates/src/fitting.jl`
with no offset / limit.

- [ ] **Step 2: Walk top-to-bottom and flag suspects per the
       criteria table**

- [ ] **Step 3: Append suspect bullets to the
       `## src/fitting.jl  (211 LOC)` section**

- [ ] **Step 4: Append intentionally-not-flagged entries covering
       `[1, 211]` exactly**

- [ ] **Step 5: Run the meta-completeness check for fitting.jl**

```bash
python3 - <<'PY'
import re, sys
FILE = "src/fitting.jl"; NLINES = 211
with open("docs/superpowers/scratch-refactor-audit-notes.md") as fh:
    t = fh.read()
sec = re.search(rf"## {re.escape(FILE)}[^\n]*\n(.*?)(?=\n## |\Z)",
                t, re.DOTALL)
assert sec, f"section not found: {FILE}"
rs = [(int(a), int(b)) for a, b in
      re.findall(r"\s*-\s*L(\d+)-L(\d+)\s*\|", sec.group(1))]
cov = [False] * (NLINES + 1)
for lo, hi in rs:
    for i in range(max(1, lo), min(NLINES, hi) + 1):
        cov[i] = True
gaps = []; ig = False; gs = None
for i in range(1, NLINES + 1):
    if not cov[i] and not ig:
        ig = True; gs = i
    elif cov[i] and ig:
        ig = False; gaps.append((gs, i - 1))
if ig:
    gaps.append((gs, NLINES))
print(f"{FILE}: {len(rs)} ranges, {'OK' if not gaps else 'FAIL'}")
for g in gaps:
    print(f"  GAP: L{g[0]}-L{g[1]}")
sys.exit(1 if gaps else 0)
PY
```

Expected: `src/fitting.jl: <N> ranges, OK`.

### Task 1.7: Catalog src/mechanism_enumeration.jl

**Files:**
- Read: `src/mechanism_enumeration.jl` (2,196 LOC) — largest file
- Modify: scratch `## src/mechanism_enumeration.jl` section

The biggest file. Per the prior-agent audit and CLAUDE.md, key areas:

- **Native expansion moves** (`_expand_re_to_ss`,
  `_expand_split_kinetic_group`, `_expand_add_dead_end_regulator`,
  `_expand_to_allosteric`, `_expand_add_allosteric_regulator`,
  `_expand_change_allo_state`) at lines ~1090-1675 per prior audit:
  do they share copy-and-rebuild patterns that could collapse into
  small internal update constructors (`_with_steps`,
  `_with_cat_states`, `_with_reg_sites`)?
- **`_catalytic_topologies` backtracking**: dead branches? Can the
  multi-tier rules merge?
- **Name-map / canonical-hash projection** at lines ~1840-2030 per
  prior audit: structural-parameter-names refactor obsoleted some of
  this.
- **Mirror-step propagation** logic: implicit per kinetic-group
  atomicity — is there leftover explicit mirror code?
- **`_to_group_list`** and `_make_species`: any branches for
  bound-form representations that were retired?
- Helpers prefixed with `_` that look like leftover scaffolding from
  the spec/scratch removal in earlier phases (e.g. `_RawSpec` was
  deleted per commit `d17192c`, but leftover callers might exist).
- **`_n_fit_params_estimate`**: known under-counter per the
  `project_n_fit_params_estimate_undercounts` memory; do not delete,
  but flag for the deferred fix.

- [ ] **Step 1: Read in full**

Use Read on
`/home/denis.linux/.julia/dev/EnzymeRates/src/mechanism_enumeration.jl`
with no offset / limit. This file is 2,196 lines, larger than the
default 2000-line Read cap; the tool reads the first 2,000 lines.
Then issue a second Read call with `offset=2000, limit=200` to get
lines 2001-2196.

- [ ] **Step 2: Walk top-to-bottom and flag suspects per the
       criteria table**

Expect a large suspect count here. Pay special attention to the "Key
areas" list above.

- [ ] **Step 3: Append suspect bullets to the
       `## src/mechanism_enumeration.jl  (2,196 LOC)` section**

- [ ] **Step 4: Append intentionally-not-flagged entries covering
       `[1, 2196]` exactly**

- [ ] **Step 5: Run the meta-completeness check for
       mechanism_enumeration.jl**

```bash
python3 - <<'PY'
import re, sys
FILE = "src/mechanism_enumeration.jl"; NLINES = 2196
with open("docs/superpowers/scratch-refactor-audit-notes.md") as fh:
    t = fh.read()
sec = re.search(rf"## {re.escape(FILE)}[^\n]*\n(.*?)(?=\n## |\Z)",
                t, re.DOTALL)
assert sec, f"section not found: {FILE}"
rs = [(int(a), int(b)) for a, b in
      re.findall(r"\s*-\s*L(\d+)-L(\d+)\s*\|", sec.group(1))]
cov = [False] * (NLINES + 1)
for lo, hi in rs:
    for i in range(max(1, lo), min(NLINES, hi) + 1):
        cov[i] = True
gaps = []; ig = False; gs = None
for i in range(1, NLINES + 1):
    if not cov[i] and not ig:
        ig = True; gs = i
    elif cov[i] and ig:
        ig = False; gaps.append((gs, i - 1))
if ig:
    gaps.append((gs, NLINES))
print(f"{FILE}: {len(rs)} ranges, {'OK' if not gaps else 'FAIL'}")
for g in gaps:
    print(f"  GAP: L{g[0]}-L{g[1]}")
sys.exit(1 if gaps else 0)
PY
```

Expected: `src/mechanism_enumeration.jl: <N> ranges, OK`.

### Task 1.8: Catalog src/identify_rate_equation.jl

**Files:**
- Read: `src/identify_rate_equation.jl` (886 LOC)
- Modify: scratch `## src/identify_rate_equation.jl` section

Areas:
- **`name_map` and string-keyed projections** (lines ~300-390, 421-496
  per prior audit): structural-parameter-names refactor obsoleted these.
- **`_canonical_rate_eq_hash_data_impl_struct`** (from Stage 6.1): is
  there a residual non-struct version?
- **`_project_cached_params`**: rendered-Symbol projection that
  structural keys could obsolete?
- **Beam search + LOOCV logic**: dead config knobs from earlier
  iterations?
- **`IdentifyRateEquationProblem` single constructor on
  `EnzymeReaction`**: any leftover branches for removed shapes?
- **Stale phase/stage comments** referencing Stage 6, 7, 7d, etc.

- [ ] **Step 1: Read in full**

Use Read on
`/home/denis.linux/.julia/dev/EnzymeRates/src/identify_rate_equation.jl`
with no offset / limit.

- [ ] **Step 2: Walk top-to-bottom and flag suspects per the
       criteria table**

- [ ] **Step 3: Append suspect bullets to the
       `## src/identify_rate_equation.jl  (886 LOC)` section**

- [ ] **Step 4: Append intentionally-not-flagged entries covering
       `[1, 886]` exactly**

- [ ] **Step 5: Run the meta-completeness check for
       identify_rate_equation.jl**

```bash
python3 - <<'PY'
import re, sys
FILE = "src/identify_rate_equation.jl"; NLINES = 886
with open("docs/superpowers/scratch-refactor-audit-notes.md") as fh:
    t = fh.read()
sec = re.search(rf"## {re.escape(FILE)}[^\n]*\n(.*?)(?=\n## |\Z)",
                t, re.DOTALL)
assert sec, f"section not found: {FILE}"
rs = [(int(a), int(b)) for a, b in
      re.findall(r"\s*-\s*L(\d+)-L(\d+)\s*\|", sec.group(1))]
cov = [False] * (NLINES + 1)
for lo, hi in rs:
    for i in range(max(1, lo), min(NLINES, hi) + 1):
        cov[i] = True
gaps = []; ig = False; gs = None
for i in range(1, NLINES + 1):
    if not cov[i] and not ig:
        ig = True; gs = i
    elif cov[i] and ig:
        ig = False; gaps.append((gs, i - 1))
if ig:
    gaps.append((gs, NLINES))
print(f"{FILE}: {len(rs)} ranges, {'OK' if not gaps else 'FAIL'}")
for g in gaps:
    print(f"  GAP: L{g[0]}-L{g[1]}")
sys.exit(1 if gaps else 0)
PY
```

Expected: `src/identify_rate_equation.jl: <N> ranges, OK`.

### Task 1.9: Catalog src/EnzymeRates.jl

**Files:**
- Read: `src/EnzymeRates.jl` (37 LOC)
- Modify: scratch `## src/EnzymeRates.jl` section

Tiny file (37 LOC) — just the module declaration, exports, and
`include` lines. Areas:
- Any exported symbol whose definition is gone or whose only callers
  are tests of removed shapes?
- Any export-list comment that references stages or specs?
- Any `using` that has become dead since a referenced submodule was
  removed?

- [ ] **Step 1: Read in full**

Use Read on `/home/denis.linux/.julia/dev/EnzymeRates/src/EnzymeRates.jl`
with no offset / limit (37 lines).

- [ ] **Step 2: Walk top-to-bottom and flag suspects per the
       criteria table**

- [ ] **Step 3: Append suspect bullets to the
       `## src/EnzymeRates.jl  (37 LOC)` section**

- [ ] **Step 4: Append intentionally-not-flagged entries covering
       `[1, 37]` exactly**

- [ ] **Step 5: Run the meta-completeness check for EnzymeRates.jl**

```bash
python3 - <<'PY'
import re, sys
FILE = "src/EnzymeRates.jl"; NLINES = 37
with open("docs/superpowers/scratch-refactor-audit-notes.md") as fh:
    t = fh.read()
sec = re.search(rf"## {re.escape(FILE)}[^\n]*\n(.*?)(?=\n## |\Z)",
                t, re.DOTALL)
assert sec, f"section not found: {FILE}"
rs = [(int(a), int(b)) for a, b in
      re.findall(r"\s*-\s*L(\d+)-L(\d+)\s*\|", sec.group(1))]
cov = [False] * (NLINES + 1)
for lo, hi in rs:
    for i in range(max(1, lo), min(NLINES, hi) + 1):
        cov[i] = True
gaps = []; ig = False; gs = None
for i in range(1, NLINES + 1):
    if not cov[i] and not ig:
        ig = True; gs = i
    elif cov[i] and ig:
        ig = False; gaps.append((gs, i - 1))
if ig:
    gaps.append((gs, NLINES))
print(f"{FILE}: {len(rs)} ranges, {'OK' if not gaps else 'FAIL'}")
for g in gaps:
    print(f"  GAP: L{g[0]}-L{g[1]}")
sys.exit(1 if gaps else 0)
PY
```

Expected: `src/EnzymeRates.jl: <N> ranges, OK`.

### Task 1.10: Pass-1 meta-completeness check across ALL files

**Files:**
- Read: `docs/superpowers/scratch-refactor-audit-notes.md`

- [ ] **Step 1: Run the cross-file completeness check**

```bash
python3 - <<'PY'
import re, sys
files = {
    "src/types.jl": 1550,
    "src/dsl.jl": 1138,
    "src/sym_poly_for_rate_eq_derivation.jl": 315,
    "src/rate_eq_derivation.jl": 1706,
    "src/thermodynamic_constr_for_rate_eq_derivation.jl": 417,
    "src/fitting.jl": 211,
    "src/mechanism_enumeration.jl": 2196,
    "src/identify_rate_equation.jl": 886,
    "src/EnzymeRates.jl": 37,
}
with open("docs/superpowers/scratch-refactor-audit-notes.md") as fh:
    text = fh.read()
fail = False
for fname, nlines in files.items():
    esc = re.escape(fname)
    sec = re.search(rf"## {esc}[^\n]*\n(.*?)(?=\n## |\Z)", text, re.DOTALL)
    if not sec:
        print(f"  MISSING SECTION: {fname}"); fail = True; continue
    ranges = []
    for line in sec.group(1).splitlines():
        m = re.match(r"\s*-\s*L(\d+)-L(\d+)\s*\|", line)
        if m:
            ranges.append((int(m.group(1)), int(m.group(2))))
    covered = [False] * (nlines + 1)
    for lo, hi in ranges:
        for i in range(max(1,lo), min(nlines, hi)+1):
            covered[i] = True
    gaps = []
    in_gap = False; gs = None
    for i in range(1, nlines+1):
        if not covered[i] and not in_gap:
            in_gap = True; gs = i
        elif covered[i] and in_gap:
            in_gap = False; gaps.append((gs, i-1))
    if in_gap: gaps.append((gs, nlines))
    status = "OK" if not gaps else f"FAIL ({len(gaps)} gap(s))"
    print(f"  {fname}: {len(ranges)} entries, {status}")
    if gaps:
        fail = True
        for g in gaps: print(f"    GAP: L{g[0]}-L{g[1]}")
sys.exit(1 if fail else 0)
PY
```

Expected output: all 9 files report `OK`. If any report `FAIL`, go
back to the corresponding Task 1.x, re-read the gap region, and add
missing entries; then re-run this check.

- [ ] **Step 2: Tally suspect counts per category**

```bash
python3 - <<'PY'
import re, collections
with open("docs/superpowers/scratch-refactor-audit-notes.md") as fh:
    text = fh.read()
cat = collections.Counter()
for line in text.splitlines():
    m = re.match(r"\s*-\s*L\d+-L\d+\s*\|\s*([\w\-]+)", line)
    if m:
        cat[m.group(1)] += 1
for c, n in cat.most_common():
    print(f"  {c}: {n}")
PY
```

Expected output: counts per category. Healthy ranges would be:
- `intentionally-not-flagged`: ~50-200 (file structure entries)
- `dead`: dozens
- `duplication`: dozens
- `symbol-tuple-plumbing`: 5-20
- `compile-time-accessor`: 5-30
- `test-private-helper`: 5-20
- `string-keyed-projection`: 3-10
- `permissive-parser-guard`: 1-5
- `stale-spec-comment`: dozens (Phase-stage / past-spec references)
- `comment-as-docstring`: 5-30

If a category that the spec called out has ZERO entries, recheck the
relevant files — a missing category likely means a suspect was missed.

---

## Phase 2 — Pass 2: Verify (parallel subagents)

### Task 2.1: Batch suspects into verification queries

**Files:**
- Modify: `docs/superpowers/scratch-refactor-audit-notes.md` (the
  `## Pass 2 — Verification queries` section)

- [ ] **Step 1: Read the scratch file**

Use Read on `docs/superpowers/scratch-refactor-audit-notes.md`.

- [ ] **Step 2: Group all (non-intentionally-skipped) suspects by
       verification query type**

For each suspect bullet, decide which verification query template
matches:

| Template | Use for | Query shape |
|---|---|---|
| **A. Caller search** | Dead-code suspects, test-private-helper suspects | "Find all references to the Julia symbol `<NAME>` in `src/` and `test/`. For each match, report `file:line`. Categorize each as: (a) definition site, (b) function call, (c) Symbol literal in code, (d) string mention in comment / docstring / test string. Report under 200 words." |
| **B. Body comparison** | Duplication suspects | "Compare the bodies of two Julia functions: `<NAME_A>` at `<FILE>:<L1>-<L2>` and `<NAME_B>` at `<FILE>:<L3>-<L4>`. List shared lines, differing lines, and any hidden side effects in either. Report whether they could collapse to a single function. Under 250 words." |
| **C. Hot-path call-graph check** | Compile-time-accessor suspects | "Trace which functions are called inside `rate_equation`'s `@generated` body in `src/rate_eq_derivation.jl`, two levels deep. Of these candidate symbols `<LIST>`, report which appear in that call graph. Under 200 words." |
| **D. Struct-native rewrite check** | Symbol-tuple-plumbing suspects | "Read `<FILE>:<L1>-<L2>` (the function `<NAME>`). It currently regenerates `<X>` from `Sig`. Confirm: is `<X>` already on the `Mechanism` value as a field or trivially derivable from one? Report yes/no with the field path." |
| **E. String-key check** | String-keyed-projection suspects | "Read `<FILE>:<L1>-<L2>` (the `name_map` / projection at this location). Identify what rendered-Symbol string is the key. Is the same information available as a structural key on `Parameter` / `Step` / `RegulatorySite`? Report yes/no with the structural-key path." |
| **F. Parser-tighten check** | Permissive-parser-guard suspects | "Read `src/dsl.jl:<L1>-<L2>` (the parser branch at this location). Confirm: if the branch refused non-conformation-shaped bare Symbols (regex `^[A-Z][a-z0-9]*(_[a-z0-9]+)*$`), would the post-hoc guard `_assert_no_opaque_terms` become redundant? Report yes/no with the redundancy claim." |
| **G. Stale-comment confirm** | Stale-spec-comment / comment-as-docstring suspects | (No verification needed; these are direct-read.) |

Doc-category suspects (G) skip Pass 2 — they go straight to promoted.

- [ ] **Step 3: Write the query list to the scratch file**

Under the `## Pass 2 — Verification queries` section, append one
bullet per query:

```
- Q-001 | template A | suspects covered: types.jl L500-L520, types.jl L600-L620 | NAME: `_species_name_from_sig`
- Q-002 | template B | suspects covered: rate_eq_derivation.jl L992-L1010, rate_eq_derivation.jl L1204-L1230 | A=`_onlyA_parameters`, B=`_I_rename_parameters`
- ...
```

Aim to keep each query covering 1-3 closely related suspects (so the
verification result maps cleanly to which suspect it bears on).

- [ ] **Step 4: Sanity check**

Every suspect from the meta-completeness check (Task 1.10) that is
NOT category G (doc) should appear in at least one Q-NNN entry. If a
suspect was skipped, decide whether it should be merged into an
existing query or get its own.

### Task 2.2: Dispatch parallel verifications

**Files:**
- Modify: `docs/superpowers/scratch-refactor-audit-notes.md` (the
  `## Pass 2 — Verification results` section)

- [ ] **Step 1: Pick a batch size**

Dispatch in batches of **at most 10 parallel Agent calls per message**
(Agent-tool runtime budget). If you have N queries, you will dispatch
ceil(N/10) batches.

- [ ] **Step 2: Dispatch a batch**

For each batch of up to 10 queries, send **one message** containing N
parallel Agent tool calls. Use `subagent_type="Explore"` for all
read-only queries (templates A, C, D, E, F); use
`subagent_type="general-purpose"` for body-comparison queries
(template B) since the comparison may need richer reasoning. Each
Agent call's `prompt:` field gets the verbatim filled-in template
text (resolve all `<NAME>` / `<FILE>` / `<LIST>` placeholders before
dispatching).

Example dispatch shape (this is illustrative — your actual queries
will vary by what Pass 1 found):

```
Agent(description="caller search _species_name_from_sig",
      subagent_type="Explore",
      prompt="Find all references to the Julia symbol `_species_name_from_sig` ...")
Agent(description="body compare _onlyA vs _I_rename",
      subagent_type="general-purpose",
      prompt="Compare the bodies of `_onlyA_parameters` at src/rate_eq_derivation.jl:992-1010 and ...")
... up to 10 calls in one message ...
```

- [ ] **Step 3: Record results**

For each batch, append to the `## Pass 2 — Verification results`
section of the scratch file one bullet per query:

```
- Q-001 result: `_species_name_from_sig` is called only by `_step_tuple_from_sig` (types.jl:1430) and one `@generated` body in rate_eq_derivation.jl:142. Both regenerate opaque tuples from Sig. → suspects PROMOTED.
- Q-002 result: bodies of `_onlyA_parameters` and `_I_rename_parameters` share the rep-step → param-name pattern but diverge on the state-rename mapping. Collapsible into one structural transform parameterized by allo-state. → suspects PROMOTED to a Cluster move.
- ...
```

- [ ] **Step 4: Repeat for remaining batches**

Until every Q-NNN has a recorded result.

### Task 2.3: Promote / demote / drop suspects based on results

**Files:**
- Modify: `docs/superpowers/scratch-refactor-audit-notes.md` (the
  `## Pass 3 — Findings (promoted suspects)` and
  `## Pass 3 — Drops (verification failed)` sections)

- [ ] **Step 1: Walk every suspect entry**

For each suspect bullet under any `## src/<file>` section, find its
verification result (via the Q-NNN cross-reference) and decide:

  - **PROMOTE**: evidence supports the suspect → assign a new finding
    ID `F-NNN` (sequential, starting at F-001 in
    Pass-1-section order: types.jl entries first, then dsl.jl, etc.)
  - **DEMOTE**: partial evidence → promote but with confidence `L`
    (Low) and a note like "needs design".
  - **DROP**: verification refuted the suspect → log under
    `## Pass 3 — Drops` with a one-line reason.

- [ ] **Step 2: Update each suspect bullet in place**

Append `| → F-NNN` (or `| → DROPPED: <reason>`) to each suspect
bullet:

```
- L500-L520 | symbol-tuple-plumbing | `_species_name_from_sig` … | H | → F-008
- L600-L620 | dead | `_orphan_helper` looks unused | M | → DROPPED: used by test/test_X.jl:300 as fixture
```

- [ ] **Step 3: Populate `## Pass 3 — Findings (promoted suspects)`**

For each PROMOTED suspect, append a one-line index entry:

```
- F-001 | src/types.jl:L100-L120 | dead | H
- F-002 | src/types.jl:L500-L520 | symbol-tuple-plumbing | H
- ...
```

This is just an index — full finding bodies are written in Phase 3.

- [ ] **Step 4: Populate `## Pass 3 — Drops`**

For each DROPPED suspect, append one bullet with the reason. This
becomes the "drops with one-line reason" record per the spec.

---

## Phase 3 — Pass 3: Synthesize and write the findings doc

### Task 3.1: Cluster findings and build the dependency graph

**Files:**
- Modify: `docs/superpowers/scratch-refactor-audit-notes.md` (add a
  new section `## Pass 3 — Clusters + dep graph`)

- [ ] **Step 1: Read the promoted-findings index**

From the `## Pass 3 — Findings (promoted suspects)` section.

- [ ] **Step 2: Group by shared architectural move**

Walk the index and assign each finding to a cluster (lettered A, B,
C, …). Heuristics:

- All findings that recommend touching the same logical concept go in
  one cluster (e.g. every finding under "symbol-tuple-plumbing" that
  recommends a struct-native context is in **Cluster: derivation
  back-end struct-native rewrite**).
- All findings whose recommendation is "collapse this A/I rename" go
  in **Cluster: A/I rename collapse**.
- Findings whose recommendation is "demote `EnzymeMechanism{Sig}`" go
  in **Cluster: singleton-type demotion**.
- Singleton findings that don't fit a cluster stay solo (referenced as
  themselves in the dep graph).

- [ ] **Step 3: Build the dependency graph**

For each finding, decide: does its recommendation assume another
finding has already landed? Common dependencies:

- Collapsing a `@generated` accessor depends on demoting the singleton
  type that exposes it (Cluster: singleton-type demotion).
- Deleting a `name_map` projection may depend on rewriting consumers
  to use structural keys.
- Replacing symbol-tuple plumbing in the derivation back-end depends on
  introducing a struct-native analysis context first.

Record dependencies in the scratch under `## Pass 3 — Clusters + dep
graph`:

```
### Clusters
- **Cluster A — derivation back-end struct-native rewrite**: F-008, F-015, F-042, F-051
- **Cluster B — A/I rename collapse**: F-023, F-024, F-025, F-026
- **Cluster C — singleton-type demotion**: F-031, F-032, F-033, F-034
- ...

### Dependency edges
- F-042 → depends on → F-031 (collapse accessor requires demoted type)
- Cluster A → depends on → Cluster C (struct-native rewrite assumes singleton demoted first)
- F-076 → depends on → F-023 (string-keyed projection rewrite assumes A/I rename collapse landed)
- ...
```

### Task 3.2: Identify blocking tests per finding

**Files:**
- Read: `test/*.jl`
- Modify: `docs/superpowers/scratch-refactor-audit-notes.md` (extend
  the index entries with a blocking-tests field)

- [ ] **Step 1: For each promoted finding, scan tests for blocking
       coverage**

A test is **blocking** if it would fail or no longer compile after the
recommendation is applied (per the spec). Use grep to scan:

```bash
# For a finding that recommends deleting symbol `_foo_bar`:
grep -n "_foo_bar" test/*.jl
# Inspect each match: does the test assert on the helper's return
# value (blocking) or just call the public API that wraps it (not
# blocking)?
```

Dispatch this in parallel for batches of findings — one Agent call
per finding (Explore subagent) with the prompt:

```
"For the finding that recommends <RECOMMENDATION> at <FILE>:<L1>-<L2>:
identify every test in test/*.jl that would fail or no longer compile
after the change is applied. Report each as test/<file>:<line> with a
one-line reason. Tests that merely call public API are not blocking
and should be excluded. Under 150 words."
```

- [ ] **Step 2: Append blocking-test lists to each finding's index
       entry**

```
- F-001 | src/types.jl:L100-L120 | dead | H | blocking: none
- F-008 | src/types.jl:L500-L520 | symbol-tuple-plumbing | H | blocking: test/test_rate_eq_derivation.jl:L420-L425 (asserts on form-name shape)
- ...
```

### Task 3.3: Sequence findings into waves

**Files:**
- Modify: `docs/superpowers/scratch-refactor-audit-notes.md` (add
  `## Pass 3 — Sequencing` section)

- [ ] **Step 1: Topologically sort findings by dependency graph**

Within each topological layer, prefer high-confidence dead-code first
(easy wins), then duplication collapse, then architectural.

- [ ] **Step 2: Write the proposed sequencing as wave bullets**

```
### Wave 1 — no-deps, high-confidence dead code
- F-001, F-003, F-007, F-012, F-019, ... (N findings, ~M LOC)

### Wave 2 — duplication collapse (depends on Wave 1 verified)
- Cluster B (A/I rename collapse), F-018, F-022, ... (N findings, ~M LOC)

### Wave 3 — architectural (depends on Wave 2 verified)
- Cluster C (singleton-type demotion), then Cluster A (derivation
  back-end struct-native rewrite), then ... (N findings, ~M LOC)

### Wave 4 — doc-hygiene sweep (no behavior change; safe to land last
or alongside other commits)
- All doc-category findings batched (~N findings)
```

### Task 3.4: Write the findings doc

**Files:**
- Create: `docs/superpowers/2026-05-30-refactor-audit-findings.md`

- [ ] **Step 1: Compose the doc header and executive summary**

Use Write to create the file with this initial content (fill in
real numbers from the scratch):

```markdown
# Concrete-Types Refactor Audit — Findings

**Date:** 2026-05-30
**Branch:** refactor-to-concrete-types-instead-of-symbols
**Baseline (non-comment non-doc src LOC):** <number from Phase 0 Task 0.2>
**Workflow:** docs/superpowers/specs/2026-05-30-refactor-audit-workflow-design.md

## §1 Executive summary

- <N total findings (D dead-code, Du duplication, A architectural,
  T test-surface, H doc-hygiene)>
- Estimated savings: ~<Z> non-comment non-doc src LOC (<P>% of baseline)
- Doc-hygiene count (does not reduce LOC): <H>
- Simplification themes:
  1. <one-line top architectural move>
  2. <one-line second move>
  3. <...>
- Hypothesis test (Denis: ~half of src removable): **supported / partially supported / not supported**
  - Verdict reasoning: <one short paragraph>
```

- [ ] **Step 2: Write the §2 Findings section, src-file by src-file**

For each src file in inclusion order, write one `### src/<file>`
subsection. Under each, list every promoted finding with the full
per-finding format from the spec:

```markdown
### src/types.jl

#### F-001  Remove orphan helper `_foo_bar`

**Location:** src/types.jl:100-120
**Category:** Dead code
**Confidence:** High
**LOC saving (non-comment non-doc):** ~21
**Simplification gain:** removes one unused symbol from the namespace; no behavior change
**Depends on:** none
**Blocking tests:** none
**Recommendation:** `_foo_bar` was introduced in Stage 6 to support the
legacy Sig path. Pass-2 caller search confirmed no production callers
remain. Delete the function and the unused `import Foo` it required.

#### F-002  ...
```

Pull the body of each finding from the scratch index plus the
verification result. Each `Recommendation:` paragraph should cite the
verification evidence (the Q-NNN result) in plain English — do NOT
include Q-NNN IDs in the committed doc, since the scratch file is
not committed.

- [ ] **Step 3: Write §3 Dependency clusters**

```markdown
## §3 Dependency clusters

- **Cluster A — derivation back-end struct-native rewrite**: F-008, F-015, F-042, F-051. Total ~<N> LOC + simplification: the King-Altman / Wegscheider consumers stop regenerating opaque Symbol tuples and instead match on `Mechanism` field values directly.
- **Cluster B — A/I rename collapse**: F-023..F-026. Total ~<N> LOC + simplification: one structural transform parameterized by allo-state replaces three near-duplicate helpers.
- ...
```

- [ ] **Step 4: Write §4 Suggested sequencing**

Copy the wave bullets from the scratch `## Pass 3 — Sequencing`
section. Format as a numbered list of waves.

- [ ] **Step 5: Write §5 Hard constraints tracked**

```markdown
## §5 Hard constraints tracked

- **`rate_equation` perf budget**: findings touching the derivation
  back-end (Cluster A) must include perf evidence (0 alloc, <100 ns)
  in their impl plans. Findings F-NNN, F-NNN are the relevant ones.
- **Test coverage**: Findings F-NNN, F-NNN delete tests — each names
  its behavior-test replacement or justifies "truly internal, no
  behavior loss".
- **Front-end struct-family unification**: must not be reintroduced as
  parallel reps. No finding in this audit proposes that, but the
  constraint is recorded for impl-plan reviewers.
```

- [ ] **Step 6: Self-review the findings doc**

Re-read the committed doc top to bottom. Check:
- Every finding has all 8 fields (Location, Category, Confidence,
  LOC saving, Simplification gain, Depends on, Blocking tests,
  Recommendation)
- Every "Depends on" reference points to a real F-NNN in the doc
- Every "Blocking tests" entry cites a real test file:line
- No `<…>` placeholder slipped through
- Executive-summary numbers match the body counts

### Task 3.5: Commit the findings doc and clean up the scratch

**Files:**
- Delete: `docs/superpowers/scratch-refactor-audit-notes.md`
- Commit: `docs/superpowers/2026-05-30-refactor-audit-findings.md`

- [ ] **Step 1: Add and commit the findings doc**

```bash
git add docs/superpowers/2026-05-30-refactor-audit-findings.md
git commit -m "$(cat <<'EOF'
docs: refactor-audit findings

Output of the 3-pass audit defined in
docs/superpowers/specs/2026-05-30-refactor-audit-workflow-design.md.
Lists every simplification opportunity the concrete-types and
structural-parameter-names refactors enabled but did not complete,
organized by src file with dependency clusters and a suggested
implementation sequencing. Tracks both the LOC-reduction target
(non-comment non-doc src LOC) and the simplification target.

Follow-on impl plans (one per dependency cluster) drive the actual
simplifications.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

- [ ] **Step 2: Delete the scratch file**

```bash
rm docs/superpowers/scratch-refactor-audit-notes.md
```

Verify it is gone:
```bash
ls docs/superpowers/scratch-refactor-audit-notes.md 2>&1
```
Expected output: `ls: cannot access '...': No such file or directory`.

- [ ] **Step 3: Verify the working tree is clean**

```bash
git status --short
```
Expected output: empty (no uncommitted changes; no untracked files —
the scratch file is gone, and its gitignore entry remains in place).

- [ ] **Step 4: Final summary**

Print a one-paragraph summary of:
- Total findings count by category
- Estimated LOC savings (non-comment non-doc) and percent of baseline
- Hypothesis verdict (supported / partial / not)
- Number of dependency clusters
- Next step: writing impl plans, one per cluster (separate
  brainstorming sessions)

---

## Self-Review Checklist (for the executor before declaring done)

After Task 3.5, run this checklist:

- [ ] **Coverage**: every src file has a `### src/<file>` subsection in
  the findings doc, even if the file produced zero findings (note
  "no findings" under the subsection rather than omitting it).
- [ ] **Hypothesis numbers**: the executive-summary percent matches
  what the per-finding LOC savings actually sum to.
- [ ] **No scratch leaks**: `grep -r "scratch-refactor-audit-notes"
  docs/superpowers/2026-05-30-refactor-audit-findings.md` returns
  empty (no references to the throwaway scratch).
- [ ] **Hard constraints honored**: no finding's recommendation
  weakens `rate_equation` performance or net-reduces test coverage.
- [ ] **Doc-category findings do not inflate the LOC savings**: doc
  hygiene LOC saving is recorded as `0; doc category` per finding,
  and the §1 summary lists doc-hygiene count separately.
