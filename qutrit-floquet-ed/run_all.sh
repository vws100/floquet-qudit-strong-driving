#!/bin/bash
# Regenerates all data and figures used in the manuscript from scratch (~2 min, N=6).
set -e
PY=${PYTHON:-python3}
$PY qutrit_chain_floquet_ed.py --selftest
for a in 2 1 0.5; do
  $PY qutrit_chain_floquet_ed.py --N 6 --T 0.05 --ncycles 600 --ppc 2 --a4 $a --out nem_a4_$a.dat
done
$PY plot_compare_a4.py nem_a4_2.dat nem_a4_1.dat nem_a4_0.5.dat --cols SzSz QxyQxy --out fig_compare_a4
for T in 0.5 2 5 10; do
  $PY qutrit_chain_floquet_ed.py --N 6 --T $T --ncycles $(python3 -c "print(max(3,int(60/$T)))") --ppc 4 --out nem_T$T.dat
done
$PY plot_omega_panels.py nem_T0.5.dat nem_T2.dat nem_T5.dat nem_T10.dat --col QxyQxy --out fig_omega_panels
$PY qutrit_chain_floquet_ed.py --N 6 --T 2 --ncycles 10 --ppc 40 --out nem_micromotion_T2.dat
$PY plot_micromotion.py nem_micromotion_T2.dat --cols SzSz QxyQxy --bare --out fig_micromotion_T2
$PY qutrit_chain_floquet_ed.py --N 6 --T 5 --ncycles 10 --ppc 40 --out nem_micromotion_T5.dat
$PY plot_micromotion.py nem_micromotion_T5.dat --cols SzSz QxyQxy --bare --out fig_micromotion_T5
echo "done"
