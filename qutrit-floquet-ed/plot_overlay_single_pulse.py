#!/usr/bin/env python3
"""plot_overlay_single_pulse.py -- exact vs effective correlators from a single run,
one panel per correlator.

Usage:
  python3 plot_overlay_single_pulse.py RUN.dat [--cols SzSz QxyQxy Q0Q0] [--kick] [--strobo] [--out NAME]
  --kick    add the kick-dressed effective curve
  --strobo  stroboscopic points only
Writes NAME.pdf and NAME.png (default NAME=RUN_overlay).
"""
import argparse, re, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt

LABELS = {'SzSz':   r'$\langle \hat S^z_i\hat S^z_{j}\rangle$',
          'QxyQxy': r'$\langle \hat Q^{xy}_i\hat Q^{xy}_{j}\rangle$',
          'Q0Q0':   r'$\langle \hat Q^{0}_i\hat Q^{0}_{j}\rangle$',
          'Qx2y2':  r'$\langle \hat Q^{x^2-y^2}_i\hat Q^{x^2-y^2}_{j}\rangle$'}

p = argparse.ArgumentParser()
p.add_argument('data', help='.dat file from qutrit_chain_floquet_ed.py')
p.add_argument('--cols', nargs='+', default=['SzSz', 'QxyQxy', 'Q0Q0'], choices=list(LABELS))
p.add_argument('--out', default=None, help='output basename (default: derived from data file)')
p.add_argument('--kick', action='store_true', help='also plot kick-dressed effective curve where available')
p.add_argument('--strobo', action='store_true', help='plot only stroboscopic points t=nT')
p.add_argument('--width', type=float, default=4.2)
p.add_argument('--height', type=float, default=2.1, help='per panel')
A = p.parse_args()

lines = open(A.data).read().splitlines()
hdr = lines[0]
cols = [l for l in lines if l.startswith('#')][-1].lstrip('# ').split()
D = np.loadtxt(A.data); d = {c: D[:, k] for k, c in enumerate(cols)}
g = lambda key: re.search(key + r'=([^\s]+)', hdr).group(1)
N, T, a4 = int(g('N')), float(g('T')), None
m = d['strobo'] == 1 if A.strobo else np.ones(len(d['t']), bool)
t = d['t'][m]

fig, ax = plt.subplots(len(A.cols), 1, figsize=(A.width, A.height*len(A.cols)), sharex=True, squeeze=False)
ax = ax[:, 0]
for a, c in zip(ax, A.cols):
    a.plot(t, d[c+'_exact'][m], '-', color='C0', lw=1.3, label='exact driven')
    a.plot(t, d[c+'_eff'][m], '--', color='C3', lw=1.3, label=r'$\overline{\hat H^{(0)}_{\rm eff}}$')
    if A.kick and c+'_kick' in d:
        a.plot(t, d[c+'_kick'][m], ':', color='C2', lw=1.3, label=r'$e^{-i\hat K^{(0)}}$-dressed')
    a.set_ylabel(LABELS[c])
ax[0].legend(frameon=False, fontsize=8)
ax[-1].set_xlabel(r'$tJ$')
ax[0].set_title(rf"$N={N}$, {g('V')}, {g('model')}, $\omega/J={2*np.pi/T:.0f}$", fontsize=9)
fig.tight_layout()
out = A.out or A.data.rsplit('.', 1)[0] + '_overlay'
fig.savefig(out + '.pdf'); fig.savefig(out + '.png', dpi=130)
print('wrote', out + '.pdf', out + '.png')
