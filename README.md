# How far can the "ETH is false" method go?

An exploration of the limits of the algorithm in OpenAI's
[*The Exponential Time Hypothesis is False*](https://openai.com/index/building-scientific-understanding-with-academic-communities/) (October 2026),
which counts the solutions of a k-CNF formula in time about `2^(n / √log n)`.

The question we asked: can that method be pushed further — towards a much faster algorithm, or even
`P = NP`? The answer we found: **not within its framework**. We prove ceilings for the method,
show the paper's own design is close to optimal for fixed designs, and record the experiments and
dead ends along the way.

> **Status — please read.** This is an informal research note, not peer reviewed. The arguments were
> developed by Mehrshad together with an AI assistant (Claude). Earlier versions of Theorems A and C
> were each checked by a separate AI reviewer; the note as a whole has not been independently
> reviewed. One lemma is machine-checked in Lean. The main algorithmic result belongs to OpenAI;
> this repo only studies its limits. Corrections are very welcome — please open an issue.

## Results

The paper's method splits variables into small groups, lets groups exchange a few linear
"summary" measurements, deforms the problem with a parameter `T`, and recovers the count from
Taylor data at `Q` evaluation points with replacement degree `R`. Full statements and proofs are in
[`NOTE.md`](NOTE.md).

| # | Result | Status |
|---|---|---|
| Lemma | Triangular fooling set: if coordinates can be ordered so none reads an earlier one or itself, then `rank(I − M0) ≥` their number | **Machine-checked in Lean** ([`lean/Ceiling/FoolingSet.lean`](lean/Ceiling/FoolingSet.lean)) |
| A | The one-level routine (Sections 2–3, no recursion) needs time `2^Ω(√n)` on expander formulas, for any design | Proved (paper's Appendix D recursion escapes it) |
| B | Two-room designs: odd cycles of near-antipodal hybrids have length `≥ Q/(R−1)`, giving `N ≥ (c/τ)^(Q/2R)` | Proved (warm-up; conditions in note) |
| C | Multi-room designs with a fixed `M0`: `rank(I − M0) ≥ (φn/32)(4|X|)^(−8R/(μQ))` via random exponential-shift partitions | Proved, constants unoptimised |
| D | The "bell curve" is unavoidable: `R ≳ √(Q log(1/τ))` for every fixed design | Proved |
| Cor. | Every fixed design costs `2^Ω(n (log log n)² / log n)`; with `Q = O(log n)` (the paper's regime) `2^Ω(n/√log n)`, so the paper's design is optimal up to a `(log log n)²` factor | Proved for fixed designs |
| Escapes | Non-linear summaries, sum-product over measurements, other recovery schemes: no asymptotic gain | Argued in note (some parts sketched) |

**Open gap.** Theorems B–D assume one measurement matrix for all evaluation points. Designs whose
measured subspaces change from point to point are covered only by Theorem A.

## Experiments

All computational, none of it proved. Python 3 with NumPy/SciPy (MILP via HiGHS).

- [`experiments/codes/`](experiments/codes/) — exact MILP search for "hybrid" designs: the paper's
  cube, linear codes, cut codes of `K_{a,b}` (symmetry-reduced up to dimension 25), random sparse
  codes, Reed–Muller, geometry line codes, cycle blow-ups. Gains over the cube are at most ~1.5× at
  small sizes and fade by dimension ≈ 19; pull-back and information-set bounds explain why.
- [`experiments/escape/`](experiments/escape/) — cut-rank of satisfiability indicators for expander
  XOR systems (fixed groups cost `2^Ω(n)`).
- [`experiments/cancellation/`](experiments/cancellation/) — "correct only in the sum"
  representations: no gap with pointwise methods under rank measures; determinant and permanent
  have identical cut-rank, so the remaining door is an elimination-style algorithm (VP vs VNP).
- [`experiments/holographic/`](experiments/holographic/) — search for matchgate bases for
  planar read-twice monotone k-CNF mod p. It rediscovers that the moduli divide `2^k − 1`
  (Cai & Lu, ECCC TR06-145, generalising Valiant's mod-7 result); the extra "2 is a square mod p"
  condition we first saw is an artifact of searching only over `F_p`.

Run any script with `python3 <script>`; results and logs from our runs are committed next to them.

## Lean

```bash
cd lean
lake exe cache get
lake env lean Ceiling/FoolingSet.lean
```

Uses Mathlib (`v4.36.0-rc2`); the proof depends only on Lean's standard axioms
(`propext`, `Classical.choice`, `Quot.sound`).

## How this came about

This started as a casual "let's prove P = NP" and turned into a day of exploring with an AI
assistant, by someone with no background in complexity theory. Several steps came from everyday
intuitions — dinner-party seating, "hybrid" guests between two rooms, a messy room whose order only
its owner knows — that were translated into mathematics and then tested or proved.

## License

MIT for code and text written here. Papers referenced belong to their authors.
