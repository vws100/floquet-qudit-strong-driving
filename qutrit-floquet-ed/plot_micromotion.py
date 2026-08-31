#!/usr/bin/env python3
"""plot_micromotion.py -- exact driven dynamics vs the kick-dressed effective theory
exp(-i K0(t)) exp(-i Heff t) |psi0> at ALL sampled times, including inside the period.
Run the ED code with --ppc large (e.g. 40) so micromotion is resolved.

Usage:
  python3 plot_micromotion.py RUN.dat [--cols SzSz QxyQxy] [--bare] [--tmax T] [--out NAME]
  --bare   also draw the undressed Heff curve (meaningful only at t = nT)
  --tmax   truncate the time axis
Writes NAME.pdf and NAME.png (default NAME=RUN_micromotion).
"""
import argparse, re, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
LABELS={'SzSz':r'$\langle \hat S^z_i\hat S^z_{j}\rangle$','QxyQxy':r'$\langle \hat Q^{xy}_i\hat Q^{xy}_{j}\rangle$',
        'Q0Q0':r'$\langle \hat Q^{0}_i\hat Q^{0}_{j}\rangle$','Qx2y2':r'$\langle \hat Q^{x^2-y^2}_i\hat Q^{x^2-y^2}_{j}\rangle$'}
p=argparse.ArgumentParser(); p.add_argument('data'); p.add_argument('--cols',nargs='+',default=['SzSz','QxyQxy'],choices=list(LABELS))
p.add_argument('--bare',action='store_true',help='also draw undressed Heff curve'); p.add_argument('--tmax',type=float,default=None)
p.add_argument('--out',default=None); p.add_argument('--width',type=float,default=4.2); p.add_argument('--height',type=float,default=2.3)
A=p.parse_args()
L=open(A.data).read().splitlines(); cols=[l for l in L if l.startswith('#')][-1].lstrip('# ').split()
D=np.loadtxt(A.data); d={c:D[:,k] for k,c in enumerate(cols)}; g=lambda k: re.search(k+r'=([^\s]+)',L[0]).group(1)
T=float(g('T')); m=np.ones(len(d['t']),bool) if A.tmax is None else d['t']<=A.tmax; t=d['t'][m]
fig,ax=plt.subplots(len(A.cols),1,figsize=(A.width,A.height*len(A.cols)),sharex=True,squeeze=False); ax=ax[:,0]
for a,c in zip(ax,A.cols):
    a.plot(t,d[c+'_exact'][m],'-',color='C0',lw=1.4,label='exact driven')
    a.plot(t,d[c+'_kick'][m],'--',color='C3',lw=1.2,label=r'$e^{-i\hat K^{(0)}(t)}e^{-i\overline{\hat H^{(0)}_{\rm eff}}t}$')
    if A.bare: a.plot(t,d[c+'_eff'][m],':',color='C2',lw=1.2,label=r'$\overline{\hat H^{(0)}_{\rm eff}}$ only')
    a.set_ylabel(LABELS[c])
ax[0].legend(frameon=False,fontsize=7); ax[-1].set_xlabel(r'$tJ_3$')
ax[0].set_title(rf"$N={g('N')}$, {g('V')}, $a_4={g('a4')}$, $\omega/J_3={2*np.pi/T:.2g}$, {int(g('ppc'))} pts/period",fontsize=9)
fig.tight_layout(); out=A.out or A.data.rsplit('.',1)[0]+'_micromotion'; fig.savefig(out+'.pdf'); fig.savefig(out+'.png',dpi=130); print('wrote',out)
