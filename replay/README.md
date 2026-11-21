# Supplement: exact certificate for the fourth-power pair-clustering ceiling

This directory contains one certificate program, `replay.py`, that recomputes, in exact
rational and symbolic arithmetic, every finite witness quoted in the preprint
"A fourth-power strong-coupling ceiling for electronic pair clustering on hypercubic
lattices in three and more dimensions" (Davis Rens, SuperMatics Inc., 2026).

The program certifies finite statements only: exact rational enclosures, exact
identity checks, and exact counts of connected supports. It does not certify the
analytical many-body theorems; those are proved in the preprint, and the program
checks the arithmetic constants that appear in their statements.

## Requirements

* `uv` 0.9.5 (any recent `uv` works)
* Python 3.13.9 (fetched and used by `uv`)
* `sympy` 1.14.0
* `mpmath` 1.3.0

Both packages are public on the Python Package Index and are installed by `uv` into a
throwaway environment. No other dependency is used, and no network access is needed
after the first installation.

## Commands

Default run: every exact witness must certify.

```
uv run --no-project --with sympy --with mpmath python replay.py
```

Guarded mutant run: fifteen deliberately corrupted variants of the certificate must
each be rejected, and the program must exit zero only if all fifteen are rejected.

```
uv run --no-project --with sympy --with mpmath python replay.py --mutants
```

Observed exit codes, in a clean copy of this directory on Windows 11 with `uv` 0.9.5:

| command | exit code | output summary |
| --- | --- | --- |
| default | 0 | 17 `CERTIFIED` lines, then a one-line scope note |
| `--mutants` | 0 | 15 `MUTANT ... REJECTED` lines, then `ALL 15 MUTANTS REJECTED` |

A nonzero exit code means a certificate failed; in the mutant run it also means a
corrupted certificate was accepted, which is the failure mode that run is designed to
detect.

## What is certified

Default run, in order:

1. Rational enclosures of the exponential, `1957/720 < 11743/4320` with the upper
   enclosure `e_hi = 1175/432 < 3`, and the three power bounds `(4/3)^4 > e_hi`,
   `(11/8)^16 > e_hi^5`, `(6/5)^6 > e_hi` used in the geometric-series estimates.
2. The single-bond form on the sixteen-state four-orbital Fock space, in every
   particle-number sector from zero to four, including the exact sum of squares on the
   particle-number one and three block as an operator identity, its scalar form on a
   single doublet for real coefficients, and the nonnegativity of the multi-bond
   diagonal bound for dimensions two to five, with the exact minimum residual, which
   is zero in every case.
3. The bivariate polynomial identity `5 z_a - 4 z_b > 0` and the pair-source scalar
   penalty, which equals `1/8` and does not exceed `1/4`.
4. Exact counts of lattice animals against the Bethe-tree counts for dimensions two to
   five at sizes one to four, including the four-cycle deficit.
5. The theorem constants, symbolically in the dimension and then numerically for
   dimensions two to five: the scale `L(d) = 72(d+1)`, the ceiling
   `C(d) = 442368 d^2 (d+1)^2`, the rate `q(d)`, the prefactor `A(d)`, the physical and
   source penalties, the Bethe-tree evaluation point inside the convergence region,
   the two Kotecky-Preiss budgets below one, the decay hierarchy
   `q_exact < q_sharp < q(d) < 5/12`, the signed-deletion prefactor below fifteen,
   the finiteness of the shell row sum, the exact improvement factors
   `C_B2(d)/C(d) = 200/9` for every dimension and `C_ET(d)/C(d) = 8(d+1)^2/(3d)`,
   which equals `128/9` at dimension three, and the three-dimensional comparison
   `C(3) = 63700992 < 2^26`.
6. Each baseline route at dimensions two to five: the conditional thresholds, the
   constants, the source budgets and the exact ceiling comparisons; then the
   obstruction of the earlier canonical Euler-tour route in three dimensions, where
   the squared source activity exceeds `675/8192 > (4/19)^2` for every scale at
   `beta t^2 / U = 2^-27`; the minimum ceiling of that route at dimension three is
   `124761600 e`, whose coefficient is checked exactly and whose value is bracketed by
   the rational enclosures of the exponential between `2.52` and `2.53` times `2^27`.
7. The asymptotic limits `C(d)/d^4 -> 442368`, `q(d) -> 5/12`,
   `q_sharp(d) -> 6903125/20155392`, `A(d) -> 180/7` and `C(d)/(8d) / d^3 -> 55296`
   used in the discussion of scaling.

The mutant run corrupts one guarded claim at a time, fifteen corruptions in all, and
requires each of them to be rejected.

## Manifest

`MANIFEST.sha256` lists the SHA-256 digest of every file in this directory except
itself. Verify with

```
python -c "import hashlib,pathlib;[print(hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest(),p) for p in ['replay.py','README.md']]"
```

or with any standard sha256 checker that understands the `<digest>  <name>` format.
