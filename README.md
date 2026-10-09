# beals

Research toward [Beal's conjecture](https://en.wikipedia.org/wiki/Beal_conjecture): if Aˣ + Bʸ = Cᶻ with positive integers A, B, C and x, y, z ≥ 3, then A, B and C share a common prime factor.

## Result

**x⁵ + y⁵ = z¹¹ has no solutions in coprime nonzero integers, with no unproven assumptions.** Hence Beal's conjecture holds for the signatures (5,5,11), (5,11,5) and (11,5,5).

Dahmen and Siksek ([Acta Arith. 164 (2014)](https://doi.org/10.4064/aa164-1-5), [arXiv:1309.4030](https://arxiv.org/abs/1309.4030)) proved this only *assuming GRH*. In their §5.3 they reduce the GRH dependence to two facts:

1. the class group of the degree-22 field L₁₁ = ℚ[t]/(t²² + 2t¹¹ − 4) is trivial;
2. the unit group used in the 2-descent is 2-saturated.

This repository proves both unconditionally.

| Step | Method | Result |
|---|---|---|
| Bound on generators of Cl(L₁₁) | Zimmert (Invent. Math. 1981, Satz 2), implemented and checked against his published table | 152,405,135 (Minkowski would be 75,252,913,349) |
| **Cl(L₁₁) = 1 (primary proof)** | explicit generator for *every* prime ideal of norm ≤ the bound, checked by exact HNF equality | [results/generators](results/generators) |
| Cl(L₁₁) = 1 (supporting proof) | PARI `bnfcertify` Phase 1, split into ranges, plus 68/68 factor-base generators | [results/phase1](results/phase1), [results/setup.out](results/setup.out) |
| PARI's units 2-saturated | F₂-rank 12 of the quadratic-character matrix | [results/units.out](results/units.out) |
| Magma's 2-descent data correct | 18 L(S,2) representatives proven independent and in L(S,2); dim = 18 = true value; dim Sel⁽²⁾ = 1 reproduced | [results/magma](results/magma) |

The full argument, the trust base and the reproduction steps are in [docs/proof.md](docs/proof.md). The one remaining trust assumption beyond PARI and Zimmert is Magma's 2-descent implementation, which Dahmen–Siksek relied on as well.

## Layout

- `certify/`: PARI/GP scripts, the PARI patch, the build script, the parallel runner and the result checker
- `results/`: certification outputs (each ends in `STATUS OK`)
- `tests/`: pytest suite, including negative tests showing that the checks reject class number 2 and non-saturated units
- `refs/`: the papers and scripts used

## Reproduce

```sh
certify/build_pari.sh
.tools/pari/bin/gp -q -D parisizemax=2000000000 certify/setup.gp < /dev/null
.tools/pari/bin/gp -q -D parisizemax=1000000000 certify/units.gp < /dev/null
certify/run_chunks.sh phase1
certify/run_chunks.sh generators
python3 certify/check_results.py
python3 -m venv .venv && .venv/bin/pip install pytest ruff && .venv/bin/pytest
```
