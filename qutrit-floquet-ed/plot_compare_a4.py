#!/usr/bin/env python3
"""plot_compare_a4.py -- overlay exact vs effective correlators from several runs
(e.g. different a_4 or omega) on common axes.  Solid = exact driven, dashed = Heff.

Usage:
  python3 plot_compare_a4.py RUN1.dat RUN2.dat ... [--cols SzSz QxyQxy] [--labels L1 L2 ...]
                             [--strobo] [--out NAME] [--width W] [--height H]
  --cols    any of SzSz QxyQxy Q0Q0 Qx2y2 (one panel each)
  --labels  legend label per file (default: a_4 read from each header)
  --strobo  plot stroboscopic points only
Writes NAME.pdf and NAME.png (default NAME=compare_a4).
"""
import argparse, re, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt

LABELS = {'SzSz':   r'$\langle \hat S^z_i\hat S^z_{j}\rangle$',
          'QxyQxy': r'$\langle \hat Q^{xy}_i\hat Q^{xy}_{j}\rangle$',
          'Q0Q0':   r'$\langle \hat Q^{0}_i\hat Q^{0}_{j}\rangle$',
          'Qx2y2':  r'$\langle \hat Q^{x^2-y^2}_i\hat Q^{x^2-y^2}_{j}\rangle$'}

p = argparse.ArgumentParser()
p.add_argument('data', nargs='+')
p.add_argument('--cols', nargs='+', default=['SzSz', 'QxyQxy'], choices=list(LABELS))
p.add_argument('--labels', nargs='+', default=None, help='legend label per file (default: a_4 from header)')
p.add_argument('--out', default='compare_a4')
p.add_argument('--strobo', action='store_true')
p.add_argument('--width', type=float, default=4.2)
p.add_argument('--height', type=float, default=2.3)
A = p.parse_args()

def load(fn):
    lines = open(fn).read().splitlines()
    cols = [l for l in lines if l.startswith('#')][-1].lstrip('# ').split()
    D = np.loadtxt(fn); d = {c: D[:, k] for k, c in enumerate(cols)}
    d['hdr'] = lines[0]; return d
g = lambda h, k: (re.search(k + r'=([^\s]+)', h) or [None, '?'])[1]

runs = [load(f) for f in A.data]
labels = A.labels or [rf"$a_4={g(d['hdr'],'a4')}$" for d in runs]
fig, ax = plt.subplots(len(A.cols), 1, figsize=(A.width, A.height*len(A.cols)), sharex=True, squeeze=False); ax = ax[:, 0]
for k, (d, lab) in enumerate(zip(runs, labels)):
    m = d['strobo'] == 1 if A.strobo else np.ones(len(d['t']), bool)
    for a, c in zip(ax, A.cols):
        a.plot(d['t'][m], d[c+'_exact'][m], '-',  color=f'C{k}', lw=1.3, label=f'{lab} exact')
        a.plot(d['t'][m], d[c+'_eff'][m],   '--', color=f'C{k}', lw=1.3, label=f'{lab} eff.')
for a, c in zip(ax, A.cols): a.set_ylabel(LABELS[c])
ax[0].legend(frameon=False, fontsize=7, ncol=2)
ax[-1].set_xlabel(r'$tJ_3$')
h = runs[0]['hdr']; T = float(g(h, 'T'))
ax[0].set_title(rf"$N={g(h,'N')}$, {g(h,'V')}, $\omega/J_3={2*np.pi/T:.0f}$", fontsize=9)
fig.tight_layout(); fig.savefig(A.out+'.pdf'); fig.savefig(A.out+'.png', dpi=130); print('wrote', A.out+'.pdf/.png')
