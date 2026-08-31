# qutrit-floquet-ed

Exact diagonalization of strongly driven qutrit chains, comparing exact
driven dynamics with the lowest-order Floquet effective theory of

> R. Scott and V. W. Scarola, *Formalism for Floquet Engineering d-level Systems
> using Strong Driving*. arxivXXX

The code propagates an N-site qutrit chain under a global square-pulse drive
(exactly, as an ordered product of segment exponentials) and under the
effective Hamiltonian $\overline{\hat H^{(0)}_{\rm eff}}$ with kick operator
$\hat K^{(0)}(t)$ derived in the paper, and records two-site spin-1
dipole/quadrupole correlators versus time.

## Install

```
pip install -r requirements.txt
```
Python ≥ 3.9, numpy, scipy, matplotlib. No other dependencies.

## Quick start

```
python3 qutrit_chain_floquet_ed.py --selftest        # verifies couplings against the paper
python3 qutrit_chain_floquet_ed.py --N 6 --T 0.05 --ncycles 600 --ppc 2 --out run.dat
python3 plot_overlay_single_pulse.py run.dat --cols SzSz QxyQxy
./run_all.sh                                          # regenerates all manuscript figures
```

Full option list: `python3 qutrit_chain_floquet_ed.py -h`, or the docstring at
the top of each file.

## Files

| file | purpose |
|---|---|
| `qutrit_chain_floquet_ed.py` | ED engine; `--selftest` checks Eqs. a4_only_eff, symmetry_a4_eff, su3_single_block_BLBQ, and the $\omega^{-1}$ propagator error |
| `plot_compare_a4.py` | several runs on common axes (e.g. different $a_4$) |
| `plot_omega_panels.py` | stroboscopic exact points vs $H_{\rm eff}$, one panel per $\omega$ |
| `plot_micromotion.py` | exact vs $e^{-i\hat K^{(0)}(t)}e^{-i H_{\rm eff} t}$ at all $t$ |
| `plot_overlay_single_pulse.py` | single-run overlay |
| `run_all.sh` | reproduces the manuscript data and figures |

## Models

* `--model nematic` (Sec. "single pulse"): $H_0=\tfrac{J_3}{2}\sum V_{ij}\hat\lambda^3_i\hat\lambda^3_j$, one pulse on $\hat\lambda^4$, $f_4=1$; `--a4`, `--J3` tunable.
* `--model su3` (Sec. "dipole pulse"): anisotropic dipolar $H_0$ driven to the SU(3) point with two pulses on $\hat\lambda^{1,2}$, $a=4$, $f=1/2$.

Lattices: nearest-neighbour ring (`nn_pbc`, default), open chain (`nn`), dipolar $|i-j|^{-3}$ (`dipolar`).

## Output

Plain-text `.dat`, `#`-commented header with all parameters and effective
couplings, then columns `t`, `<O_iO_j>` for exact / effective / kick-dressed
evolution, `state_err`, `strobo` flag. Load with `numpy.loadtxt`.

## Notes

* $\overline{\hat H^{(0)}_{\rm eff}}$ alone is meaningful only at $t=nT$; use the `_kick` columns off-stroboscopically.
* $\langle\hat Q^0_i\hat Q^0_j\rangle$ is exactly conserved in the nematic model (each $\hat Q^0_j$ commutes with $H_0$, $\hat\lambda^4$, and $H_{\rm eff}$); it is recorded as a check.
* Deterministic; no random numbers. Runtime ≈ 10 s per N=6 run.

## Citation

If you use this code, please cite the paper above.

## License

MIT — see [../LICENSE](../LICENSE).
