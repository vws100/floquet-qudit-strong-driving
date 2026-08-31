# Mathematica codes

Two standalone Wolfram Language scripts. Both run in the notebook front end
(open, then Evaluation → Evaluate Notebook) or from the command line with
`wolframscript`. Tested with Mathematica 13/14; any version ≥ 12 should work.

## 1. `sud_generator_general_d.wl` — generators, support map, lookup tables

Builds, for a single input dimension *d* ≥ 2:

1. the su(*d*) S.A.D. generator matrices (inputs S^{mn}, A^{mn}, H^{1n},
   then the output-only diagonals H^{mn});
2. the support map supp(α);
3. the fully expanded nested-commutator lookup tables in the convention

   ```
   [[T^a, T^b]]_u = 2^{-lambda} phi^{(p)} 2^{-(u-1) nu} T^{kappa_p(a,b)},  p = u mod 2,
   ```

   with the four table classes same-support / diagonal-α / diagonal-β /
   off-diagonal.

**How to run.** Set the dimension in the `USER INPUT` block near the top
(`d = 5;` as shipped), then

```
wolframscript -file sud_generator_general_d.wl
```

**Output.** Generator matrices and support map are printed; the script writes

- `sud_support_map_d<d>.tex` — support map as an `align` environment
- `sud_lookup_tables_d<d>.tex` — the four lookup tables as `table` environments
- `sud_lookup_rows_d<d>.m` — machine-readable rows in Wolfram list syntax

to the notebook/script directory. The `.m`, and `.tex` files shipped here for
*d* = 3, 4, 5 are the outputs of this script and serve as reference data for
the Python port (`../Python/sud_lookup.py --verify`).

## 2. `two_site_effective_vs_exact_check_d3.wl` — two-site bond-error check

Two-site *d* = 3 check of the effective description against the exact driven
propagator for the square-pulse drive. Because the drive is piecewise
constant, the exact time-ordered propagator factorizes into ordinary matrix
exponentials over the constant-drive segments — no truncation parameter is
involved. The script:

- builds H₀ and the effective Hamiltonian H_eff from the closed-form
  effective couplings J^eff (u_α, v_α) of the paper;
- runs sanity checks (Hermiticity, tracelessness, off-pulse limit
  H_eff → H₀, kick operator vanishing at cycle boundaries, validation of the
  exact propagator against a midpoint-sampled uniform grid and the
  cycle-power identity);
- scans the bond error ε = ‖e^{−iK(t)} e^{−iH_eff t} e^{iK(0)} − U_exact‖
  versus ω at fixed t (log-log slope → −1) and versus t at fixed ω.

**How to run.**

```
wolframscript -file two_site_effective_vs_exact_check_d3.wl
```

Two parameter blocks are provided in the script; exactly one should be
uncommented at a time (the same applies to the corresponding `JVals` blocks):

- **single-pulse nematic model** (active as shipped): pulse on λ⁴, f₄ = 1;
- **two-pulse SU(3)-symmetric model**: pulses on λ¹, λ², a = 4, f = 1/2.

**Output.** Parameter summary, check results, and both scans are printed; the
script writes `epsilon_vs_omega.dat`, `epsilon_vs_time.dat`, and
`pulse_profile.dat` (tab-separated, `#`-commented headers) and, in the front
end, displays the pulse-profile and error plots.

A line-by-line Python port with identical structure, defaults, and output
files is at `../Python/two_site_effective_vs_exact_check_d3.py`; the two
agree to 11–12 significant figures on the scan quantities (see the Python
README for the verification summary).
