# ABOUTME: PR #86's two-free-form fixtures (test/test_types.jl) as DSL step lines, shared by
# ABOUTME: the fixture scripts: (label, substrates, products, inhibitors, [(step, tag)]).
fixtures = [
    ("R1: F <--> E SS NonequalAI", ["S"], ["P"], String[],
     [("E + S ⇌ E(S)", "EqualAI"), ("E(S) <--> E(P)", "NonequalAI"),
      ("F + P ⇌ E(P)", "EqualAI"), ("F <--> E", "NonequalAI")]),
    ("R2: every step at E OnlyA", ["S"], ["P"], String[],
     [("E + S ⇌ E(S)", "OnlyA"), ("E(S) <--> E(P)", "EqualAI"),
      ("F + P ⇌ E(P)", "EqualAI"), ("F <--> E", "OnlyA")]),
    ("R3: E <--> F OnlyA, F binds I", ["S"], ["P"], ["I"],
     [("E + S ⇌ E(S)", "EqualAI"), ("E(S) <--> E(P)", "EqualAI"),
      ("E + P ⇌ E(P)", "EqualAI"), ("E <--> F", "OnlyA"), ("F + I ⇌ F(I)", "EqualAI")]),
    ("R4: F ⇌ E OnlyA (RE)", ["S"], ["P"], String[],
     [("E + S ⇌ E(S)", "EqualAI"), ("E(S) <--> E(P)", "NonequalAI"),
      ("E + P ⇌ E(P)", "EqualAI"), ("F ⇌ E", "OnlyA")]),
    ("A5: F ⇌ E NonequalAI (accepted by PR #86)", ["S"], ["P"], String[],
     [("E + S ⇌ E(S)", "EqualAI"), ("E(S) <--> E(P)", "NonequalAI"),
      ("E + P ⇌ E(P)", "EqualAI"), ("F ⇌ E", "NonequalAI")]),
    ("S1: stray free form F + B <--> E(B)", ["A", "B"], ["P"], ["A"],
     [("E + A <--> E(A)", "OnlyA"), ("E(A) + B <--> E(A, B)", "EqualAI"),
      ("E(A, B) <--> E(P)", "OnlyA"), ("E(P) <--> E + P", "EqualAI"),
      ("E + A::Inh ⇌ E(A::Inh)", "OnlyA"), ("E(A::Inh) + B ⇌ E(A::Inh, B)", "EqualAI"),
      ("E(B) + A::Inh ⇌ E(A::Inh, B)", "EqualAI"), ("F + B <--> E(B)", "EqualAI")]),
]
# PR #86's multiplicity: 2 for R1–R4 and A5, 1 for S1 (no multiplicity line).
written_n(label) = startswith(label, "S1") ? 1 : 2
