# Rule: Write Minimal, Well-Documented Code

Code in this project should be as simple and small as possible while remaining clear and well-documented. Follow these principles:

## Minimize Code Volume (Not Line Count)

- Fewer functions, fewer types, fewer intermediate data structures.
- Do NOT compress multiple operations onto single lines to reduce line count — that produces hard-to-read code.
- Each line should do one clear thing. Readability is non-negotiable.

## Prefer Simpler Algorithms

- Before writing a multi-phase algorithm, check if a single unified approach works (e.g., Cartesian product over per-element options instead of separate enumeration passes with dedup).
- Inline helper functions that are called exactly once and are short — don't create a function just to name 5 lines of code.
- Avoid intermediate data structures that are built, returned, and immediately consumed. Inline the computation instead.

## DSL and Macro Design

- When extending DSL syntax, prefer leveraging existing parsing structures over adding new parsing functions.
- Example: `A[CX, 2]` reuses the existing `:ref` parser vs. `max_sites: A => 2` which requires 4 new AST-handling functions.

## Remove, Don't Toggle

- If a feature isn't needed, remove it entirely. Don't add a parameter to disable it.
- Accessor functions that exist only for user convenience but aren't used internally should be omitted until actually needed.

## Documentation and Tests

- Keep docstrings and comments — they don't count as code bloat.
- Keep tests — they don't count as code bloat.
- Only remove docs/tests when the functions they document/test are removed.

## Review Discipline

- After any refactor, re-read every changed file looking for dead code, unnecessary abstractions, and further simplification opportunities.
- Repeat review passes until no more improvements are found.
