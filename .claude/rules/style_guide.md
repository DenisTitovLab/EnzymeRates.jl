# Julia Style Guide (Merged from Julia Manual + SciML Style)

When Julia and SciML style guides conflict, prefer the Julia manual convention.

---

## Naming

| Element | Convention | Examples |
|---------|-----------|----------|
| Modules, Types | `CamelCase` | `SparseArrays`, `UnitRange` |
| Abstract types | `Abstract` prefix, `CamelCase` | `AbstractSolver`, `AbstractArray` |
| Functions | lowercase, words concatenated; underscores only when needed for readability | `isequal`, `haskey`, `maximum` |
| Variables | `snake_case` | `local_cache`, `step_count` |
| Constants | `UPPER_SNAKE_CASE` | `DEFAULT_TOLERANCE`, `MAX_ITER` |
| Type variables | Single capital letter | `T`, `N`, `S` |
| Mutating functions | Append `!` | `sort!`, `push!`, `fit!` |

- **No abbreviations** unless universally understood. Prefer `polynomial` over `poly`.
- **Unicode** is acceptable where it improves readability in internal code, but never in public API names (terminal compatibility).
- IO and RNG functions are exceptions to the `!` convention — `read(io)` mutates the stream, `rand(rng)` mutates the RNG.

## Formatting

### Indentation and Line Length

- **4 spaces** per indentation level. No tabs.
- **92 character** line length limit.
- Do not align assignment operators (`=`) across lines.
- Unix line endings (`\n`) only.
- Files end with a single newline.
- Trim trailing whitespace.

### Whitespace

- No spaces inside parentheses, brackets, or braces: `f(a, b)`, `[1, 2]`, `{T}`.
- No spaces before commas or semicolons.
- No spaces around `:` in ranges: `1:10`, `a[1:end]`.
- Single space around binary operators: `=`, `+=`, `==`, `!=`, `->`, `&&`, `||`.
- **No extra spaces** around `:`, `//`, `^`, or `=` in keyword arguments: `f(x; tol=1e-8)`, `a:b`, `x^2`.
- No space between unary operators and operands: `-x`, `!flag`.

### Blank Lines

- One blank line between function definitions.
- Group related single-line definitions together without blanks.
- Separate multi-line blocks with blank lines.
- No empty line immediately after a function signature or before `end`.
- Use a blank line between control flow blocks and return statements.

### Long Lines

Function calls exceeding the line limit break with arguments indented one level:

```julia
result = long_function_name(
    first_argument,
    second_argument;
    keyword_one=value,
    keyword_two=value,
)
```

### NamedTuples

```julia
xy = (x = 1, y = 2)    # space around = in named tuples
x = (x = 1,)           # trailing comma for single-element
x = (; kwargs...)       # semicolon for splatting
```

## Functions

### Argument Order

Follow this precedence:

1. Function argument (e.g., `f` in `map(f, x)`)
2. I/O stream
3. Input being mutated
4. Type
5. Input not being mutated
6. Key
7. Value
8. Everything else
9. Varargs
10. Keyword arguments

### Design Rules

- Single-line definitions only if the entire function fits on one line.
- Use keyword arguments with semicolons: `function solve(prob; tol=1e-8)`.
- Inputs should be required unless applicable to >95% of use cases.
- Prefer function instances over types as arguments: `Tsit5()` not `Tsit5`.
- Constructors must return an instance of their own type. `T(x)` returns a `T`.
- Prefer out-of-place operations and immutable structures. Reserve mutation for performance-critical non-allocating code only.
- A function should either fully embrace mutation (with `!` suffix) or treat all inputs as immutable — never mix.

### Functions Over Scripts

Divide programs into functions as soon as possible. Functions improve reusability, testability, and performance. Avoid relying on global variables except constants.

## Types and Generics

### Generic Programming

Write generic code. Prefer abstract types in function signatures to allow Julia's specialization:

```julia
# Good — accepts any AbstractArray and Integer
splicer(arr::AbstractArray, step::Integer) = arr[begin:step:end]

# Bad — limited to specific types
splicer(arr::Array{Int}, step::Int) = arr[begin:step:end]
```

- Use `complex(float(x))` instead of `Complex{Float64}(x)`.
- Use `similar(A)` instead of `Array(undef, size(A))`.
- Avoid hard-coded indexing; prefer `begin`, `end`, and broadcast operations.
- Avoid strange unions like `Union{Function, AbstractString}`.
- Use `Vector{Any}` over `Vector{Union{Int, String, Tuple}}`.

### Struct Field Types

Use concrete types in struct fields when possible. Parametric types for generality:

```julia
# Concrete — preferred when type is known
struct MySubString <: AbstractString
    string::String
    offset::Int
end

# Parametric — when generality is needed
struct MySubString{T<:Integer} <: AbstractString
    string::String
    offset::T
end
```

Avoid abstract field types like `AbstractString` — they prevent compiler inference. If the type truly varies, annotate with `::Any` explicitly.

### Type Stability

Prioritize type-stable code. Functions should return consistent types regardless of input values (not just input types).

### Unnecessary Type Parameters

Don't introduce type parameters that aren't used in the body:

```julia
# Bad — T is unused
foo(x::T) where {T<:Real} = 2x

# Good
foo(x::Real) = 2x
```

### Type Testing

Use `isa` and `<:` for type checks, not `==`:

```julia
# Good
x isa Float64
T <: AbstractArray

# Bad
typeof(x) == Float64
```

### Type Piracy

Never extend or redefine methods on types you don't own. This causes unexpected breakage in unrelated code.

## Numbers

- Floating-point literals must include leading and trailing zeros: `0.1` not `.1`, `2.0` not `2.`.
- Prefer `Int` over `Int32`/`Int64` unless specifically needed.
- Use `Cint`, `Clong`, etc. for C interop.
- In generic code, prefer integer literals (`2 * x`) over float literals (`2.0 * x`) to preserve argument types. Use rationals (`2//1`) for non-integer constants.

## Control Flow

### For Loops

Use `in` exclusively, never `=` or `∈`:

```julia
for i in 1:10
    # ...
end

[f(x) for x in xs]
```

### Conditionals

- No parentheses around conditions: `if a == b`, not `if (a == b)`.
- Keep ternary operators on a single line. Never chain them:

```julia
# Good
result = x > 0 ? x : -x

# Bad — use if-elseif-else instead
result = x > 0 ? x : x == 0 ? zero(x) : -x
```

### Error Handling

- Avoid excessive `try/catch`. Prevent errors through proper design.
- Define and throw custom exception types rather than using `error("string")`.
- Validate inputs early with domain-specific error messages.
- Include suggestions for correction when possible.
- Contextualize errors using user-facing terminology, not internal implementation details.

## Code Patterns

### Prefer

```julia
map(f, a)                    # not map(x -> f(x), a)
[a; b]                       # not [a..., b...]
map(Base.Fix2(getindex, i), vecs)  # not map(v -> v[i], vecs)
```

### Avoid

- Splicing abuse: `[a..., b...]` when `[a; b]` suffices.
- Trivial anonymous functions wrapping named functions.
- Mutable globals. Use immutable `const` at file top, after imports.
- Direct struct field access from outside the module. Prefer exported accessor methods.
- `unsafe_*` operations in public interfaces. Validate safety or include "unsafe" in the name.

## Macros

- Don't overuse macros. If a macro calls `eval`, it's often better as a function.
- Limit macros to syntactic sugar where generated code remains obvious.
- Acceptable macros: `@inbounds`, `@muladd`, `@view`, `@.`, `@simd`, `@threads`.

## Comments

- Use `TODO` for incomplete work, `XXX` for broken code.
- Quote code with backticks in comments: `` `variable_name` ``.
- Prefer code clarity over comments.
- Place comments above the code they describe, not inline (unless short).
- Include GitHub issue/PR URLs in relevant comments.

## Documentation

- Use triple-quoted `"""` docstrings.
- Only exported/public functions require docstrings.
- Wrap docstrings at 92 characters.
- Follow this template for functions:

```julia
"""
    mysearch(array::MyArray{T}, val::T; verbose=true) where {T} -> Int

Searches `array` for `val`.

# Arguments
- `array::MyArray{T}`: the array to search.
- `val::T`: the value to search for.

# Keywords
- `verbose::Bool = true`: print progress details.

# Returns
- `Int`: index where `val` is located.

# Throws
- `NotFoundError`: if `val` is not found.
"""
```

Follow this template for types:

```julia
"""
    MyArray{T, N}

My super awesome array wrapper!

# Fields
- `data::AbstractArray{T, N}`: the wrapped array.
- `metadata::Dict`: metadata about the array.
"""
```

## Imports and Modules

- Imports at file top, separated into `import` and `using` blocks.
- Large import sets on single lines, comma-separated.
- Exported names form the public API — changing them is a breaking change.
- Module files should contain only the module definition.

## Testing

- Tests must be reproducible in isolation.
- Group tests by category.
- Cover diverse input types where applicable: `Float64`, `Float32`, `Complex`, `BigFloat`.

## Dependencies and Versioning

- Use Semantic Versioning.
- All dependencies require upper bounds in compat.
- Prefer interface packages (e.g., `ChainRulesCore.jl`, `RecipesBase.jl`) over conditional loading.

## Security

- Avoid `unsafe_*` operations; minimize `@inbounds` and `ccall`.
- Validate all user inputs.
- Use `RandomDevice()` for cryptographic security, not `rand()`.
- Never use `eval` with user-influenced input.
- Initialize arrays explicitly (`zeros`, `fill`) rather than `undef` when values matter.
