# A fourth-power strong-coupling ceiling for electronic pair clustering on hypercubic lattices in three and more dimensions
Research record: [A fourth-power strong-coupling ceiling for electronic pair clustering on hypercubic lattices in three and more dimensions](https://supermatics.io/research/cubic-pair-clustering-ceiling).

Davis Rens, SuperMatics Inc. ([ORCID 0000-0001-5839-5197](https://orcid.org/0000-0001-5839-5197)).
Version 1.0.0.

Paper DOI: <https://doi.org/10.5281/zenodo.24412193>
Replay supplement DOI: <https://doi.org/10.5281/zenodo.24412195>
Source: <https://github.com/SuperMatics-Inc/cubic-pair-clustering-ceiling>

## What this deposit is

This repository carries the preprint `paper/main.pdf` (23 pages) and the complete replay
supplement for it. The manuscript proves exponential clustering of equal-time on-site
electronic pair correlations for the spinful Holstein model on the hypercubic torus
`(Z/MZ)^d` in every dimension `d >= 2`, with explicit constants that track the dimension.
With `U = F^2/K` the only hypothesis is `beta t^2 / U <= 1/C(d)` where
`C(d) = 442368 d^2 (d+1)^2`; no restriction on `t/U` is imposed, and the bound is uniform
in the phonon frequency, the torus side and the chemical potential. At distinct sites

    |<Delta_x^dag Delta_y>| <= A(d) q(d)^(|x-y|_1),
    A(d) = 360(d+1)/(14d+29), q(d) = 5(2d-1)/(24(d+1)),

with `Delta_x = c_{x down} c_{x up}`. Hence a temperature below which pair correlations
need not cluster obeys `k_B T_pair(d) <= C(d) t^2/U`.

In dimensions two and three the constants are explicit: `C(2) = 15 925 248 < 2^27`,
which improves the released square-lattice constant `2^27` by a factor `8.43`, and
`C(3) = 63 700 992 < 2^26`, the first bound of this kind in dimension three. The
mechanism has two independent parts: the hopping majorant is controlled edge by edge by
an exact sum of squares on the sixteen-state bond Fock space, valid in every
particle-number sector and for arbitrarily large hopping, so no smallness condition on
`t/U` arises; the combinatorial cost of connected bond supports is then paid by injecting
site animals into rooted subtrees of a Bethe lattice, which replaces doubled-edge Euler
tours and lowers the ceiling from order `d^5` to order `d^4`.

## Scope and labels

The bound is on `T_pair`, a clustering temperature of a two-point correlation function,
and it is not a critical temperature: no superconducting transition is claimed, no
material is discussed and no measurement is compared. The theorem covers the spinful
Holstein model with nearest-neighbour hopping, an on-site electron-phonon coupling and no
other term; the oscillator problem is treated at every `Omega > 0` and no statement is
made at `Omega = 0`; the torus side satisfies `M >= 4`.

| Result | Label | Evidence |
|---|---|---|
| `k_B T_pair <= 442368 d^2 (d+1)^2 t^2/U` with the clustering rates `A(d)` and `q(d)`, for every `d >= 2` | PROVED | `paper/main.tex`, Theorem 1 and its proofs in Sections 5 to 11 |
| `C(2) = 15 925 248 < 2^27` and `C(3) = 63 700 992 < 2^26` | CERTIFIED here | `replay/replay.py`, the statements `Theorem 1 (d=2)` and `Theorem 1 (d=3)`, which print `C=15925248` and `C=63700992` with the exact budgets |
| Every finite witness quoted in the paper: the rational enclosures of `e`, the single-bond and multi-bond forms, the animal counts, the theorem constants, the two Kotecky-Preiss budgets, the baseline routes and the asymptotic limits | CERTIFIED here | `replay/replay.py`, seventeen certified statements |
| The analytic cluster-expansion arguments that apply those witnesses | PROVED | `paper/main.tex`; the replay checks finite arithmetic only and does not formalize them |

## Layout

| Path | Contents |
|---|---|
| `paper/` | `main.tex`, `refs.bib` and `main.pdf`, the released bytes. The only change from the reviewed source is the publication metadata block: the `\PAPERDOI` and `\SUPPLEMENTDOI` macros and one archiving sentence at the start of the supporting-material section |
| `replay/` | the reviewed supplementary archive, byte identical to the reviewed bytes: the certificate program `replay.py`, its `README.md` and `MANIFEST.sha256`, which pins the two files |
| `MANIFEST.sha256` | SHA-256 of every tracked file of this repository except itself and `replay/MANIFEST.sha256` |
| `CITATION.cff` | citation metadata for the preprint and for this supplement |

The reviewed supplementary archive addresses the manuscript as the preprint and describes
its own commands relative to the archive root. In this repository the manuscript source is
`paper/` and the archive content is `replay/`.

## Replay

The supplement imports only `sympy` and `mpmath` beyond the Python standard library. The
deposited runs used Python 3.13.9 with `sympy==1.14.0` and `mpmath==1.3.0`, installed into
a throwaway environment by `uv`; the workflow here installs the same two pins with `pip`
and runs the same two commands.

```sh
sha256sum -c MANIFEST.sha256
(cd replay && sha256sum -c MANIFEST.sha256)
python -m pip install sympy==1.14.0 mpmath==1.3.0
(cd replay && python -B replay.py)
(cd replay && python -B replay.py --mutants)
```

The default run prints seventeen `CERTIFIED` lines and one scope line and exits zero only
if every exact statement holds. The `--mutants` run corrupts one guarded claim at a time,
fifteen times, and exits zero only if every corrupted claim is rejected; its last line is
`ALL 15 MUTANTS REJECTED`. Exit code zero in the first run means the constants are
reproduced, not that the analysis is formalized.

## Build

The manuscript is built from `paper/` with Tectonic 0.17.0 and a fixed timestamp:

```sh
cd paper && SOURCE_DATE_EPOCH=1795132800 tectonic -X compile main.tex --keep-logs
```

The build reports no overfull box, no underfull box, no undefined reference, no undefined
citation and no LaTeX warning of any kind, and 23 pages. Tectonic output is not byte
reproducible, so `paper/main.pdf` is the local build that was deposited, and the workflow
checks its own rebuild log and text rather than the hash of a rebuild.

## Licenses and provenance

The replay code and the documentation in this repository are licensed under Apache-2.0. The
preprint is distributed under CC-BY-4.0 through its Zenodo record.

AI agents generated proofs, certificates, code, and manuscript drafts and edits under
D. Rens' direction and review. Separate AI agents performed the independent reviews of the
cited records. This preprint has not been peer reviewed.
