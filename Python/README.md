# Python ports of the Mathematica codes

Two files, one-to-one ports of the scripts in `../Mathematica/`:

| file | port of | purpose |
|---|---|---|
| `sud_lookup.py` | `sud_generator_general_d.wl` | S.A.D. generators, support map, commutator lookup tables for general *d* |
| `two_site_effective_vs_exact_check_d3.py` | `two_site_effective_vs_exact_check_d3.wl` | two-site *d* = 3 bond-error check of H_eff vs. the exact driven propagator |

## Requirements

Python ≥ 3.10 with numpy (`sud_lookup.py`) and scipy
(`two_site_effective_vs_exact_check_d3.py`):

```
pip install -r requirements.txt
```

---

# `sud_lookup.py` — su(*d*) S.A.D. generators, support map, and commutator lookup tables

Python port of `sud_generator_general_d.wl`. Builds the S.A.D. generator basis for
arbitrary local dimension *d*, the support map, and the fully expanded commutator
lookup tables in the convention

```
[[T^a, T^b]]_u = 2^{-lambda} phi^{(p)} 2^{-(u-1) nu} T^{kappa_p(a,b)},   p = u mod 2
```

with `p = o` for odd *u* and `p = e` for even *u*. Pairs that commute are not listed;
by convention `phi = 0` and `lambda = nu = 0` for those.

All generator entries are dyadic rationals (±1/2, ±i/2), so `complex128` arithmetic
is exact here — the tolerances in the code mirror the Wolfram script rather than
compensating for round-off.

## Quick start

```bash
# reproduce the Wolfram script's default run (d = 4): prints the generator
# matrices and support map, and writes the four output files
python3 sud_lookup.py -d 4 -o .

# same, plus the copy-paste LaTeX on stdout
python3 sud_lookup.py -d 3 --print-latex

# just the table summary, no matrix dump
python3 sud_lookup.py -d 6 -q

# regression test against a Mathematica export
python3 sud_lookup.py -d 4 --no-write --verify sud_lookup_rows_d4.m

# check the geometric law directly for u = 1..8 rather than inferring it from u = 1,2
python3 sud_lookup.py -d 5 --no-write --selftest
```

## Command-line options

| flag | effect |
|---|---|
| `-d`, `--dim` | local Hilbert-space dimension, *d* ≥ 2 (default 4) |
| `-o`, `--outdir` | directory for written files (default: current) |
| `--print-generators` / `--no-print-generators` | print the generator matrices, part 1 of the Wolfram script (default: on) |
| `--print-support` / `--no-print-support` | print the support map, part 2 (default: on) |
| `-q`, `--quiet` | suppress parts 1 and 2; print only the table summary |
| `--print-latex` | print the copy-paste LaTeX to stdout |
| `--no-write` | suppress all file output |
| `--json` | additionally write `sud_lookup_rows_d<d>.json` |
| `--verify FILE` | compare row-for-row against a Mathematica `.m` export |
| `--selftest` | verify the nested-commutator law for *u* = 1…8 |

Exit status is 1 if `--verify` or `--selftest` fails, 0 otherwise, so both are usable
in CI.

## Output files

Identical names and contents to the Wolfram script, plus two additions:

- `sud_generators_d<d>.tex` — the S.A.D. generator matrices as `align`
  environments grouped by sector, each written as (1/2) times a matrix over
  {0, ±1, ±i}; matches the hand-written *d* = 3 list in Appendix A.
  (The Wolfram script prints these to the notebook but does not export them.)
- `sud_support_map_d<d>.tex` — the support map as an explicit `align` environment,
  one entry per input generator for the chosen *d* (not a generic cases-form rule)
- `sud_lookup_tables_d<d>.tex` — the four lookup tables as `table` environments,
  with the six columns *α* | *β* | *λ* | *ν* | (*κ*ₒ, *φ*⁽ᵒ⁾) | (*κ*ₑ, *φ*⁽ᵉ⁾)
- `sud_lookup_rows_d<d>.m` — machine-readable rows in Wolfram list syntax
- `sud_lookup_rows_d<d>.json` — same rows as JSON (`--json`; not in the Wolfram version)

The `.m` writer reproduces Mathematica's `Export` formatting, so `diff` against a
Mathematica-generated file works directly (up to trailing whitespace).

## Library use

```python
from sud_lookup import build_all_rows, partition_tables, latex_tables, gen_matrix

rows = build_all_rows(4)                    # list of Row tuples
tables = partition_tables(rows)             # dict of the four table classes
tex = latex_tables(4, tables)               # LaTeX string
M = gen_matrix(("S", 1, 3), 4)              # 4x4 numpy array for S^{13}
```

A `Row` is

```python
(alpha, beta, lam, nu, (kappa_o, phi_o), (kappa_e, phi_e))
```

where each label is a tuple `("S"|"A"|"H", m, n)` with 1-based level indices
`1 <= m < n <= d`, and each phase is one of the strings `"+1"`, `"-1"`, `"+i"`, `"-i"`.
The even-branch term is `None` only if the second nested commutator vanishes (it never
does for compact su(*d*), but the branch is kept for parity with the Wolfram code).

Useful entry points:

| function | purpose |
|---|---|
| `input_labels(d)` / `output_labels(d)` | 𝒜ᵢₙ (size *d*²−1) and 𝒜ₒᵤₜ |
| `gen_matrix(label, d)` | S.A.D. matrix for a label |
| `supp(label)` | support map |
| `decompose(M, d)` | write `M` as `c · T^γ` for a single output generator |
| `build_row(a, b, d)` | one lookup row, or `None` if the pair commutes |
| `build_all_rows(d)` | all non-commuting ordered input pairs |
| `partition_tables(rows)` | split into `same_support`, `diag_alpha`, `diag_beta`, `offdiag` |
| `latex_generators(d)`, `latex_support_map(d)`, `latex_tables(d, tables)` | LaTeX generation |
| `rows_to_wl(rows)`, `parse_wl_rows(path)` | Wolfram `.m` write / read |
| `compare_rows(got, want)` | order-independent row comparison |
| `check_nested_commutator_law(d, rows, umax)` | direct check of the geometric law |

## Basis convention

The input set is

- 𝒜_S = { S^{mn} : 1 ≤ m < n ≤ d }, with S^{mn} = ½(|m⟩⟨n| + |n⟩⟨m|)
- 𝒜_A = { A^{mn} : 1 ≤ m < n ≤ d }, with A^{mn} = (i/2)(|m⟩⟨n| − |n⟩⟨m|)
- 𝒜_D = { H^{1n} : n = 2…d }, with H^{mn} = ½(|m⟩⟨m| − |n⟩⟨n|)

so |𝒜ᵢₙ| = *d*² − 1. The diagonal inputs are pair-difference operators anchored at
level 1 — a non-orthogonal Cartan basis. For *d* = 3 this reproduces λ³ = H^{12} and
λ⁸ = H^{13} exactly, matching the manuscript.

The output label set enlarges this to all pair-difference diagonals,
𝒜ₒᵤₜ = 𝒜_S ∪ 𝒜_A ∪ { H^{mn} : 1 ≤ m < n ≤ d }. Diagonal commutator outputs that are
not themselves inputs — e.g. λ^h = H^{23} = λ⁸ − λ³ for *d* = 3 — appear only there,
as output-only bookkeeping generators.

Because the H^{mn} are not mutually orthogonal, `decompose` matches by exact
proportionality `M == c·G` rather than by trace projection, returning the first match
in 𝒜ₒᵤₜ order. This reproduces the Wolfram routine's behaviour exactly.

## Verification performed

| check | result |
|---|---|
| *d* = 4 rows vs. `sud_lookup_rows_d4.m` from Mathematica | 168/168 rows agree exactly |
| `.m` file byte-diff vs. the Mathematica export | identical up to trailing whitespace |
| `.m` write → read round trip, *d* = 2, 3, 4, 5 | 0 discrepancies |
| nested-commutator law, *u* = 1…8 | 4800 checks across *d* = 2–5, all pass |
| generated `.tex` compiles under `pdflatex` | *d* = 2, 3, 4, 5, 6 all clean |
| support-map `.tex` vs. the Wolfram script's | byte-identical, *d* = 2…7 |

Row counts: *d* = 2 → 6; *d* = 3 → 54 (14 / 8 / 8 / 24); *d* = 4 → 168 (24 / 24 / 24 / 96);
*d* = 5 → 372 (36 / 48 / 48 / 240), in the order same-support / diagonal-α /
diagonal-β / off-diagonal.

## Note on the self-test

The Wolfram script infers *λ*, *ν*, *κ*, and *φ* from *u* = 1 and *u* = 2 only, and the
geometric form for higher *u* is assumed. `--selftest` checks it directly: that *κ*
depends on the parity of *u* and not on *u* itself, and that the magnitudes are
geometric with ratio 2^{−ν}. This is worth quoting in the Data Availability Statement,
since it is exactly the assumption a referee would probe.

## Reproducibility

`numpy` 2.4.4, Python 3.12, no random number use anywhere — the output is
deterministic and depends only on *d*.

---

# `two_site_effective_vs_exact_check_d3.py` — two-site bond-error check

Line-by-line Python port of `../Mathematica/two_site_effective_vs_exact_check_d3.wl`,
with the same structure, parameter blocks, defaults, and output files. The
drive is piecewise constant, so the exact time-ordered propagator factorizes
into ordinary matrix exponentials over the constant-drive segments — no
truncation parameter is involved.

## Quick start

```
python3 two_site_effective_vs_exact_check_d3.py
```

This runs, for the active parameter block:

1. sanity checks — H_eff Hermitian and traceless, off-pulse limit
   H_eff → H₀, kick operator vanishing at cycle boundaries, validation of the
   exact propagator against a midpoint-sampled uniform grid and the
   cycle-power identity;
2. **scan 1** — bond error
   ε = ‖e^{−iK(t)} e^{−iH_eff t} e^{iK(0)} − U_exact‖ versus ω at fixed t;
   the printed log-log slope is −1 to seven decimal places;
3. **scan 2** — ε versus t at fixed ω;

and writes `epsilon_vs_omega.dat`, `epsilon_vs_time.dat`, and
`pulse_profile.dat` (tab-separated, `#`-commented headers; load with
`numpy.loadtxt`).

## Choosing the model

As in the `.wl`, two parameter blocks are provided in the source; uncomment
exactly one `aAmpVals`/`fFracVals` pair (and the matching `JVals` block):

- **single-pulse nematic model** (active as shipped): pulse on λ⁴, f₄ = 1;
- **two-pulse SU(3)-symmetric model**: pulses on λ¹, λ², a = 4, f = 1/2.

## Verification against the Mathematica version

| check | result |
|---|---|
| u_α, v_α, and H_eff vs. exact symbolic arithmetic (sympy reproduction of the Wolfram exact path), both parameter blocks | agree to one ulp (‖ΔH_eff‖ ≈ 6 × 10⁻¹⁷) |
| H_eff closed form vs. direct average of e^{iK(t)} H₀ e^{−iK(t)} over one cycle (20 000-point quadrature), both blocks | ≤ 8 × 10⁻¹⁴ |
| exact propagator vs. independent ODE integration (`scipy.solve_ivp`, rtol 10⁻¹¹) | 2 × 10⁻¹⁰ at t = 2.1672 T |
| scan quantities vs. 50-digit `mpmath` ground truth | agree to 11–12 significant figures |
| all inline "expect ~X" magnitudes of the `.wl` | reproduced, both blocks |

## Reproducibility

Developed with Python 3.12, numpy 2.4, scipy 1.17; deterministic, no random
number use.
