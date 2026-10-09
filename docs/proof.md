# x⁵ + y⁵ = z¹¹: an independent unconditional proof along the Dahmen–Siksek route

## Statement

**Theorem.** There are no integers x, y, z with xyz ≠ 0 and gcd(x, y, z) = 1 such that

  x⁵ + y⁵ = z¹¹.

Since all three exponents are odd, the signatures (5, 5, 11), (5, 11, 5) and (11, 5, 5) are equivalent: move a term across and change its sign. This was first proven unconditionally by Stoll (2017); see below.

## Status before this work, and what is new here

**Correction (2026-10-09): this theorem was already known unconditionally.**

- M. Stoll, *Chabauty without the Mordell–Weil group*, in: Algorithmic and Experimental Methods in Algebra, Geometry, and Number Theory, Springer (2017), 623–663, [arXiv:1506.04286](https://arxiv.org/abs/1506.04286). Theorem 8.8 proves that x⁵ + y⁵ = zᵖ has only trivial primitive solutions for every prime 7 ≤ p ≤ 53, unconditionally for p ≤ 19 and assuming GRH for p ≥ 23.
- Stoll works with a different curve, C′ₚ : 5y² = 4xᵖ + 1, over ℚ, and uses "Selmer group Chabauty".
- The 2025 survey [arXiv:2412.11933](https://arxiv.org/abs/2412.11933) does not list this result. We missed it the first time.

**What this repository adds.** An *independent second proof* for p = 11 along Dahmen–Siksek's original route: the curve C₁₁,₀ over ℚ(√5), with the class group of a degree-22 field certified via Zimmert's bound. It is not a new case of Beal's conjecture.

The method itself is reusable: Zimmert's bound instead of Minkowski's at degree > 20, a range-split `bnfcertify`, and an exact-generator certificate. It applies to Stoll's GRH-conditional cases; see "Next case" below.

The original reference for the GRH-conditional route:

> S. R. Dahmen, S. Siksek, *Perfect powers expressible as sums of two fifth or seventh powers*,
> Acta Arith. 164 (2014), 65–100, doi:[10.4064/aa164-1-5](https://doi.org/10.4064/aa164-1-5),
> [arXiv:1309.4030](https://arxiv.org/abs/1309.4030). Referred to below as [DS].

## Where GRH enters [DS]

The case 5 | z is settled unconditionally for every exponent ≥ 2 ([DS, Prop. 3.3], by the modular method).

For 5 ∤ z the argument runs as follows:

1. **Modular method.** This step is unconditional: [DS, Prop. 5.3] for l = 11. It gives x + ζ₅y = β¹¹ with β ∈ ℤ[ζ₅].
2. **Reduction to a curve.** Taking norms down to K = ℚ(√5) = ℚ(θ), with θ² + θ − 1 = 0, every solution gives a point in C₁₁,₀(K)′ on the genus-5 curve

   C₁₁,₀ : 5Y² = (4θ + 12)X¹¹ − (4θ + 7).

   By [DS, Lemma 5.8] it suffices to show C₁₁,₀(K)′ = {∞, (1, ±1)}.
3. **Rank bound.** The rank of J₁₁,₀(K) is 1, because
   - [(1,1) − ∞] has infinite order, and
   - dim Sel⁽²⁾(K, J₁₁,₀) = 1 ([DS, Table 8], computed by 2-descent in Magma **under GRH**).
4. **Chabauty–Coleman** at the prime above 3. This step is unconditional ([DS, Lemma 5.9]).

[DS, §5.3] pins down exactly what GRH is used for:

> "we 'only' need to obtain certain class and unit group information unconditionally in order to carry out the 2-descent … it suffices to have available the class and unit group information of the number field L := K[x]/f(x) … L₁₁ = ℚ[t]/g₁₁(t) with g₁₁(t) := t²² + 2t¹¹ − 4 … It suffices in fact to know that our conditional unit group is a finite index 2-saturated subgroup of the (unconditional) unit group. This … reduces the problem to verifying that the class groups of the fields L₁₁ and L₁₃ are trivial."

So GRH is removed for (5,5,11) by two facts:

- **(A)** Cl(L₁₁) = 1, unconditionally.
- **(B)** The unit group used by the 2-descent is 2-saturated.

## What we prove

All computations use PARI/GP 2.19.0 on an Apple M3 (8 cores). The only change to PARI is a patch of about 60 lines that adds a range-restricted copy of PARI's own `bnftestprimes`, which is Phase 1 of `bnfcertify`. Magma checks were run on the public Magma calculator (V2.29-10).

### 0. The field

L₁₁ = ℚ[t]/(t²² + 2t¹¹ − 4):
- degree 22, signature (r₁, r₂) = (2, 10);
- d = 2²⁰·5¹¹·11²² ≈ 4.17·10³⁶.

`nfisisom` confirms that the étale algebra K[X]/((4θ+12)X¹¹ − (4θ+7)) of C₁₁,₀ is this field.

### 1. An unconditional bound for generators of Cl(L₁₁): Zimmert

PARI's `bnfcertify` uses tabulated Zimmert constants only up to degree 20. For degree 22 it falls back to the Minkowski bound, 75,252,913,349. That is about 3·10⁹ primes, or roughly 1100 core-hours.

Instead we use Zimmert's theorem directly:

> R. Zimmert, *Ideale kleiner Norm in Idealklassen und eine Regulatorabschätzung*, Invent. Math. 62 (1981), 367–380, **Satz 2**.

It states that every ideal class contains an integral ideal 𝔞 with

  log(|d|^{1/2} / N𝔞) ≥ Z(r₁, r₂; γ, α)

for all γ > α > 0, where Z is an explicit combination of Γ and ψ values.

We implement Z from the a,b-form in Zimmert's proof (`certify/zimmert.gp`) and check it two ways:
- it agrees with the printed form of Satz 2 to 10⁻³⁰;
- it reproduces every row of Zimmert's Tabelle 1 that we tested (8 rows, n = 2 to 100) to the 4 printed digits (`tests/test_certify.py`).

With α chosen by Zimmert's Bemerkung 1 and γ = 0.425, we get Z(2, 10) = 23.3182 and:

  **Zimmert bound = 152,405,135**, which is 494× smaller than Minkowski (8,571,784 primes).

### 2. (A): Cl(L₁₁) is trivial

`nfcertify` confirms that PARI's integral basis is correct, so the discriminant and the ring of integers are proven.

**Primary proof: an explicit generator for every small prime.** For *every* prime ideal P of norm ≤ 152,405,135, we compute a generator a and verify exactly that (a) = P, by comparing HNF matrices. This uses only exact integer arithmetic. The GRH-conditional data only helps *find* a; it plays no part in the verification. Every class contains such a P (Zimmert), and every such P is principal, so **Cl(L₁₁) = 1**.
- Scripts: `certify/generators.gp` and `certify/run_chunks.sh generators`.
- Results: `results/generators/`. All 64 chunks pass, covering **8,570,075 prime ideals**, about 8.8 core-hours on an M3.

**Supporting proof: PARI's `bnfcertify` Phase 1.** Let B be PARI's GRH-conditional `bnfinit` data, with factor base FB: 68 prime ideals, all of norm ≤ 389.

1. Every FB prime is principal. This is checked with explicit generators and exact HNF equality (`certify/setup.gp`, `results/setup.out`).
2. Every prime ideal of norm ≤ 152,405,135 lies in ⟨FB⟩. This uses PARI's own `bnftestprimes`, the Phase 1 code of `bnfcertify`, split into 64 contiguous ranges of p (`certify/phase1.gp`, `results/phase1/`). The sum of the per-chunk prime counts equals π(152,405,135) = 8,571,784.
   - For each P it finds α ∈ P such that (α)P⁻¹ factors over FB; the valuations are computed exactly.
   - When all primes above p are listed and the last one is unramified, the last one is implied by the others.
   - Caveat: `factorgen` accepts a candidate after rounding a *floating-point* norm to an integer (`buch2.c`, `embed_norm`). A false pass would need a relative error of at least 50%, but the check is not rigorous in principle. Stock `bnfcertify` has the same caveat. This is why the explicit-generator check above is the primary proof.

The factor base is the same in both runs: `bnfinit` is deterministic in its factor base, as verified by comparing `B[5]` across runs.

### 3. (B): the units are 2-saturated

**PARI's units.** B gives −1 and 11 fundamental units u₁…u₁₁. Each uᵢ is checked to be an algebraic integer of norm ±1.

For degree-1 primes 𝔭 = (p, t − r) with p > 1000, we form the matrix over F₂ of quadratic characters u ↦ (u mod 𝔭 / p). It reaches **rank 12** (p ≤ 1061).

What this shows:
- No nontrivial product of the 12 generators is a square in L₁₁.
- So they are multiplicatively independent, and generate a finite-index subgroup U of the unit group E. (r₁ + r₂ − 1 = 11, and the torsion of L₁₁ is {±1}.)
- The index is odd: if u ∈ E∖U had u² ∈ U, then u² would be a square in L₁₁ lying in U∖U².

Script: `certify/units.gp`; result in `results/units.out`.

**Magma's data (what the 2-descent actually uses).** [DS]'s Selmer bound was computed by Magma with Magma's own GRH-conditional data. An unsaturated unit group would make the computed Selmer group a *subgroup* of the true one, which would invalidate the rank bound. We therefore check Magma's object directly. Transcripts are in `results/magma/`, from Magma V2.29-10 on the public calculator.

- `selv`: run with verbose output, `TwoSelmerGroup(J₁₁,₀)` over K reports the bad primes of K, and that S consists of the primes of the algebra above them. Those are exactly the primes of L₁₁ above 2, 5 and 11. It then reports "Computing A(2,S) … A(2,S)=(Z/2Z)^18", and the final Selmer rank is 1, matching [DS, Table 8].
- `psel2`: in one session, `pSelmerGroup(2, S)` for L₁₁ (S = the 6 primes above 2, 5, 11) returns 18 representatives.
  - Each lies in L(S,2), checked by an exact factorization of (r): every prime outside S occurs to an even power.
  - They are independent modulo squares: F₂-rank 18 of their quadratic characters at degree-1 primes, p ≤ 1069.

  Since Cl(L₁₁) = 1, the true F₂-dimension of L(S,2) is |S| + r₁ + r₂ = 18. So Magma's L(S,2) for this field and S is correct and complete. The verbose output above shows the 2-descent builds a group of exactly this dimension for exactly this field and S.

### Conclusion

(A) and (B) are the GRH-dependent inputs identified in [DS, §5.3]. The rest of [DS]'s proof for l = 11 is unconditional:
- the modular method (Prop. 5.3, Lemmas 5.5–5.6);
- the case 5 | z (Prop. 3.3);
- J₁₁,₀(K)[2] = 0, since f is irreducible of odd degree;
- Chabauty–Coleman (Lemma 5.9).

Therefore x⁵ + y⁵ = z¹¹ has no non-trivial primitive solutions, unconditionally. ∎

## What must be trusted

- **Zimmert's Satz 2**, as published (refereed, Invent. Math.). The implementation is checked against his table.
- **PARI/GP 2.19.0**: exact HNF ideal arithmetic, `nfcertify`, and residue-field arithmetic for the characters. These are standard, heavily used routines.
- **Magma's `TwoSelmerGroup`**, given correct class and unit data. In particular, its internal "A(2,S)" step is assumed to use the same L(S,2) machinery as `pSelmerGroup`, which we verified for this field. We could not inspect Magma's source or run both in one 60-second calculator session. [DS] already trusted Magma's 2-descent; we remove only the GRH dependence.
- **[DS]'s published unconditional steps.**

Independent re-runs of the Magma step with a full Magma licence, ideally by the original authors, would remove the third item.

## Reproducing

```sh
certify/build_pari.sh                      # PARI 2.19.0 + patch -> .tools/pari  (~1 min)
.tools/pari/bin/gp -q -D parisizemax=2000000000 certify/setup.gp < /dev/null
.tools/pari/bin/gp -q -D parisizemax=1000000000 certify/units.gp < /dev/null
certify/run_chunks.sh phase1               # supporting proof, ~35 min on 8 cores
certify/run_chunks.sh generators           # primary proof of (A), ~75 min on 8 cores (8.8 core-hours)
python3 certify/check_results.py           # coverage and prime-count checks
.venv/bin/pytest tests
```

## Next case: done

x⁵ + y⁵ = z²³ has since been proven unconditionally with this method together with a free-software reimplementation of Stoll's criterion. See [5-5-23.md](5-5-23.md).
