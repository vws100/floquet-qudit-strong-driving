# floquet-qudit-strong-driving

Code accompanying

> R. Scott and V. W. Scarola, *Formalism for Floquet Engineering d-level Systems
> using Strong Driving*. arxivXXX

The repository has three parts, mirroring the structure of the paper's
numerical checks:

| directory | contents | language |
|---|---|---|
| [`Mathematica/`](Mathematica/) | su(*d*) S.A.D. generators, support map, and commutator lookup tables for general *d*; two-site *d* = 3 bond-error check of the effective Hamiltonian against the exact driven propagator | Wolfram Language |
| [`Python/`](Python/) | independent Python ports of both Mathematica codes, cross-validated against the Mathematica output | Python 3 |
| [`qutrit-floquet-ed/`](qutrit-floquet-ed/) | exact matrix representation of driven N-site qutrit chains: exact propagation vs. the lowest-order Floquet effective theory, with the scripts that regenerate the manuscript figures | Python 3 |

Each directory has its own README with explicit run instructions.

## Quick start

Mathematica (Wolfram Language ≥ 12; either open a `.wl` file in the front end
and evaluate all cells, or run from the command line):

```
cd Mathematica
wolframscript -file sud_generator_general_d.wl
wolframscript -file two_site_effective_vs_exact_check_d3.wl
```

Python (≥ 3.10; numpy and scipy suffice for `Python/`, matplotlib is needed
only for the plotting scripts in `qutrit-floquet-ed/`):

```
cd Python
pip install -r requirements.txt
python3 sud_lookup.py -d 3 --print-latex
python3 two_site_effective_vs_exact_check_d3.py
```

Exact diagonalization and manuscript figures:

```
cd qutrit-floquet-ed
pip install -r requirements.txt
python3 qutrit_chain_floquet_ed.py --selftest
./run_all.sh
```

All codes are deterministic; no random numbers are used anywhere.

## Citation

If you use this code, please cite the paper above and the Zenodo record for
this repository (DOI badge will appear here once minted).

## License

MIT — see [LICENSE](LICENSE).
