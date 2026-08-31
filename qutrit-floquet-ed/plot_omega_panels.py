#!/usr/bin/env python3
"""plot_omega_panels.py -- one panel per run: exact dynamics at stroboscopic times
t = nT (markers) vs the continuous Heff curve.  Used to show the approach to
breakdown as omega -> Lambda.

Usage:
  python3 plot_omega_panels.py RUN1.dat RUN2.dat ... [--col QxyQxy] [--out NAME] [--width W] [--height H]
  --col   one of SzSz QxyQxy Q0Q0 Qx2y2
Writes NAME.pdf and NAME.png (default NAME=omega_panels).
"""
import argparse, re, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
LABELS={'SzSz':r'$\langle \hat S^z_i\hat S^z_{j}\rangle$','QxyQxy':r'$\langle \hat Q^{xy}_i\hat Q^{xy}_{j}\rangle$',
        'Q0Q0':r'$\langle \hat Q^{0}_i\hat Q^{0}_{j}\rangle$','Qx2y2':r'$\langle \hat Q^{x^2-y^2}_i\hat Q^{x^2-y^2}_{j}\rangle$'}
p=argparse.ArgumentParser(); p.add_argument('data',nargs='+'); p.add_argument('--col',default='QxyQxy',choices=list(LABELS))
p.add_argument('--out',default='omega_panels'); p.add_argument('--width',type=float,default=4.2); p.add_argument('--height',type=float,default=1.9)
A=p.parse_args()
def load(fn):
    L=open(fn).read().splitlines(); cols=[l for l in L if l.startswith('#')][-1].lstrip('# ').split()
    D=np.loadtxt(fn); d={c:D[:,k] for k,c in enumerate(cols)}; d['hdr']=L[0]; return d
g=lambda h,k: re.search(k+r'=([^\s]+)',h).group(1)
runs=[load(f) for f in A.data]
fig,ax=plt.subplots(len(runs),1,figsize=(A.width,A.height*len(runs)),sharex=True,squeeze=False); ax=ax[:,0]
for a,d in zip(ax,runs):
    T=float(g(d['hdr'],'T')); s=d['strobo']==1
    a.plot(d['t'],d[A.col+'_eff'],'-',color='C3',lw=1.2,label=r'$\overline{\hat H^{(0)}_{\rm eff}}$')
    a.plot(d['t'][s],d[A.col+'_exact'][s],'o',color='C0',ms=3,label='exact, $t=nT$')
    a.text(0.98,0.9,rf'$\omega/J_3={2*np.pi/T:.2g}$',transform=a.transAxes,ha='right',va='top',fontsize=9)
    a.set_ylabel(LABELS[A.col])
ax[0].legend(frameon=False,fontsize=8,loc='lower left'); ax[-1].set_xlabel(r'$tJ_3$')
ax[0].set_title(rf"$N={g(runs[0]['hdr'],'N')}$, {g(runs[0]['hdr'],'V')}, $a_4={g(runs[0]['hdr'],'a4')}$",fontsize=9)
fig.tight_layout(); fig.savefig(A.out+'.pdf'); fig.savefig(A.out+'.png',dpi=130); print('wrote',A.out)
