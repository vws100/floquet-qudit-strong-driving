#!/usr/bin/env python3
"""
qutrit_chain_floquet_ed.py
==========================

Exact driven dynamics vs. lowest-order Floquet effective theory for an
N-site qutrit (d=3) chain under strong, global, square-pulse driving.

Companion code for
  R. Scott and V. W. Scarola, "Formalism for Floquet Engineering d-level
  Systems using Strong Driving" arxivXXX.
Section and equation labels below refer to that manuscript.

Physics
-------
  H(t) = H_0 + V(t),      V(t) = omega * sum_alpha g_alpha(t) * sum_j lambda^alpha_j
  H_0  = (1/2) sum_{i!=j} V_ij sum_beta J_beta lambda^beta_i lambda^beta_j
  lambda^1..8 : su(3) generators in the S.A.D. basis, tr(lambda^a lambda^b)=delta_ab/2
                (App. "spin mapping"; lambda^3, lambda^8 are NOT orthogonal)
  g_alpha(t)  : square pulse, +a, -a, +a on the first f/4, middle f/2, last f/4
                of channel alpha's window, channels sequential, sum_alpha f_alpha = 1
                (Eq. assumed_square_pulse_g with r_alpha = 1)
  omega = 2 pi / T

Three propagations of a product initial state are compared:
  exact : U(t) = T-exp[-i int_0^t H(s) ds].  H is piecewise constant, so this is an
          ordered product of matrix exponentials -- no truncation parameter.
  eff   : exp(-i Heff t),  Heff = lowest-order time-averaged effective Hamiltonian,
          Eq. d_3_effective_H with couplings from App. "jeff_definition", valid at t = nT
  kick  : exp(-i K0(t)) exp(-i Heff t),  K0(t) = sum_alpha G_alpha(t) sum_j lambda^alpha_j,
          G_alpha = omega int_0^t g_alpha  (Eq. kick_zero_order) -- valid at all t

Observables (spin-1 basis |+1>,|0>,|-1> = qutrit |1>,|2>,|3>):
  <Sz_i Sz_j>, <Q0_i Q0_j>, <Qxy_i Qxy_j>, <Qx2y2_i Qx2y2_j>
  Q0 = Sz^2 - 2/3,  Qxy = {Sx,Sy},  Qx2y2 = Sx^2 - Sy^2

Models
------
  --model nematic : Sec. d_3_single_pulse.  H_0 = (J3/2) sum V_ij lambda^3 lambda^3,
                    one pulse on lambda^4 with f_4 = 1.  Tunables: --J3, --a4
                    (Table tab_single_pulse_assumptions: J3 = 1, a4 = 2).
  --model su3     : Sec. d_3_dipole_pulse.  Anisotropic dipolar H_0 (J_ex = 1,
                    J3 = 2/3, J8 = 4/3, J38 = 0), pulses on lambda^1, lambda^2 with
                    a = 4, f = 1/2 (Table tab_two_pulse_assumptions).

Usage
-----
  python3 qutrit_chain_floquet_ed.py [options]

  --N INT            sites (default 4).  dim = 3^N; N=6 ~ 10 s, N=7 ~ minutes.
  --model            nematic | su3               (default nematic)
  --V                nn_pbc | nn | dipolar       (default nn_pbc)
                       nn_pbc : V_ij = 1 for |i-j| = 1 mod N (ring)
                       nn     : open chain
                       dipolar: |i-j|^-3, open chain, theta_ij = pi/2
  --T FLOAT          drive period in units of 1/J; omega = 2 pi / T   (default 0.05)
  --ncycles INT      periods propagated; total time = ncycles * T      (default 60)
  --ppc INT          samples per period; 1 = stroboscopic only         (default 8)
  --pair i j         sites in <O_i O_j>, 0-indexed                     (default 0 1)
  --init             zz_neel | p0_neel | all0 | xnem | xnem_stag | tilted  (default zz_neel)
                       zz_neel   |+1,-1,+1,-1,...>
                       p0_neel   |+1, 0,+1, 0,...>
                       all0      |0,0,...>
                       xnem      (|+1>+|-1>)/sqrt2 on every site
                       xnem_stag (|+1> +/- |-1>)/sqrt2, sign alternating
                       tilted    (|+1> +/- |0> + |-1>)/sqrt3, sign alternating
  --a4 FLOAT         pulse height a_4          (nematic only, default 2)
  --J3 FLOAT         J_3                       (nematic only, default 1)
  --out FILE         output .dat               (default corr_vs_time.dat)

Examples
--------
  # Fig.: exact vs effective, high frequency
  python3 qutrit_chain_floquet_ed.py --N 6 --T 0.05 --ncycles 600 --ppc 2 --out nem_a4_2.dat
  # a_4 dependence
  python3 qutrit_chain_floquet_ed.py --N 6 --T 0.05 --ncycles 600 --ppc 2 --a4 1 --out nem_a4_1.dat
  # approach to breakdown, omega ~ Lambda
  python3 qutrit_chain_floquet_ed.py --N 6 --T 5 --ncycles 12 --ppc 4 --out nem_T5.dat
  # micromotion resolved
  python3 qutrit_chain_floquet_ed.py --N 6 --T 2 --ncycles 10 --ppc 40 --out mm_T2.dat

Output
------
Whitespace-delimited text.  Header lines (prefixed '#') record all run parameters,
u_alpha, v_alpha, J^eff_beta, J^eff_38, then column names.  Columns:
  t
  Q0Q0_exact   Q0Q0_eff   Q0Q0_kick
  QxyQxy_exact QxyQxy_eff QxyQxy_kick
  state_err            ||psi_exact(t) - psi_kick(t)||
  strobo               1 if t = nT else 0
  SzSz_exact   SzSz_eff
  Qx2y2_exact  Qx2y2_eff
  SzSz_kick    Qx2y2_kick
Read with numpy.loadtxt(FILE).  Plot scripts: plot_compare_a4.py, plot_omega_panels.py,
plot_micromotion.py, plot_overlay_single_pulse.py.

Notes
-----
* <Q0_i Q0_j> is exactly conserved for the nematic model at all t (each Q0_j commutes
  with H_0, with lambda^4, and with Heff).  It is included as a check, not as dynamics.
* Off-stroboscopic comparison must use the *_kick columns; *_eff alone is meaningful
  only at t = nT.
* Deterministic (no random numbers).  Requires numpy, scipy.

Self-test
---------
  python3 qutrit_chain_floquet_ed.py --selftest
reproduces: (i) Eq. a4_only_eff couplings at a4 = 0.7, 1.3, 2 and Eq. symmetry_a4_eff,
(ii) the SU(3) BLBQ form Eq. su3_single_block_BLBQ, (iii) two-site propagator error
scaling eps ~ omega^-1 for both models.

License: MIT.  Copyright (c) 2026 Scarola Research Group, Virginia Tech.
"""
import argparse, sys, time
from functools import reduce
import numpy as np
from scipy.linalg import expm, eigh

__version__ = "1.0.0"

# ----------------------------------------------------------------------------- su(3)
def sad_generators():
    """S.A.D. basis of su(3), normalised tr(l^a l^b) = delta_ab / 2 (App. spin mapping)."""
    L = np.zeros((8, 3, 3), dtype=complex)
    L[0] = 0.5*np.array([[0,1,0],[1,0,0],[0,0,0]])
    L[1] = 0.5*np.array([[0,1j,0],[-1j,0,0],[0,0,0]])
    L[2] = 0.5*np.array([[1,0,0],[0,-1,0],[0,0,0]])
    L[3] = 0.5*np.array([[0,0,1],[0,0,0],[1,0,0]])
    L[4] = 0.5*np.array([[0,0,1j],[0,0,0],[-1j,0,0]])
    L[5] = 0.5*np.array([[0,0,0],[0,0,1],[0,1,0]])
    L[6] = 0.5*np.array([[0,0,0],[0,0,1j],[0,-1j,0]])
    L[7] = 0.5*np.array([[1,0,0],[0,0,0],[0,0,-1]])
    return L
LAM = sad_generators()

def sc1(x):
    """sc_1(x) = sin(x)/x - 1, with the removable singularity at 0."""
    x = np.asarray(x, dtype=float)
    return np.where(x == 0, 0.0, np.sin(x)/np.where(x == 0, 1.0, x) - 1.0)

def u_coef(a, f): return (f/2)*sc1(np.pi*a*f)      # Eq. square_pulse_uvw
def v_coef(a, f): return (f/2)*sc1(np.pi*a*f/2)

# ------------------------------------------------------- effective couplings (d=3)
def effective_couplings(J, u, v):
    """J, u, v: length-8 sequences, index 0..7 <-> beta = 1..8.
    Transcription of App. jeff_definition.  Returns (Jeff[8], J38eff)."""
    J1,J2,J3,J4,J5,J6,J7,J8 = J
    u1,u2,u3,u4,u5,u6,u7,u8 = u
    v1,v2,v3,v4,v5,v6,v7,v8 = v
    J1e = (J1*(1+u2+u3+v4+v5+v6+v7+v8) - (u3+v8)*J2 - u2*J3 - v7*J4 - v6*J5
           - v5*J6 - v4*J7 - 0.25*u2*J8)
    J2e = (J2*(1+u1+u3+v4+v5+v6+v7+v8) - (u3+v8)*J1 - u1*J3 - v6*J4 - v7*J5
           - v4*J6 - v5*J7 - 0.25*u1*J8)
    J3e = (J3*(1+u1+u2+0.25*u6+0.25*u7+v6+v7) - u2*J1 - u1*J2 - u7*J6 - u6*J7
           + (0.25*u1+0.25*u2+0.25*u6+0.25*u7-v1-v2-v6-v7)*J8)
    J4e = (J4*(1+u5+u8+v1+v2+v3+v6+v7) - v7*J1 - v6*J2 - 0.25*u5*J3 - (u8+v3)*J5
           - v2*J6 - v1*J7 - u5*J8)
    J5e = (J5*(1+u4+u8+v1+v2+v3+v6+v7) - v6*J1 - v7*J2 - 0.25*u4*J3 - (u8+v3)*J4
           - v1*J6 - v2*J7 - u4*J8)
    J6e = (J6*(1+u7+v1+v2+v3+v4+v5+v8) - v5*J1 - v4*J2 - 0.25*u7*J3 - v2*J4 - v1*J5
           - (v3+v8)*J7 - 0.25*u7*J8)
    J7e = (J7*(1+u6+v1+v2+v3+v4+v5+v8) - v4*J1 - v5*J2 - 0.25*u6*J3 - v1*J4 - v2*J5
           - (v3+v8)*J6 - 0.25*u6*J8)
    J8e = (J8*(1+u4+u5+0.25*u6+0.25*u7+v6+v7)
           + (0.25*u4+0.25*u5+0.25*u6+0.25*u7-v4-v5-v6-v7)*J3 - u5*J4 - u4*J5 - u7*J6 - u6*J7)
    J38e = ((-0.25*u6-0.25*u7+v4+v5)*J3 + u7*J6 + u6*J7 + (-0.25*u6-0.25*u7+v1+v2)*J8)
    return np.array([J1e,J2e,J3e,J4e,J5e,J6e,J7e,J8e]), J38e

# ------------------------------------------------------------------ lattice operators
def site_op(op, i, N, d=3):
    ops = [np.eye(d)]*N
    ops[i] = op
    return reduce(np.kron, ops)

def two_body(Jbeta, J38, V, N, lam=LAM):
    """(1/2) sum_{i!=j} V_ij [ sum_b Jb l^b_i l^b_j + J38 (l^3_i l^8_j + l^8_i l^3_j) ].
    Returns (H, Ls) with Ls[i][b] = lambda^b on site i."""
    D = 3**N
    H = np.zeros((D, D), dtype=complex)
    Ls = [[site_op(lam[b], i, N) for b in range(8)] for i in range(N)]
    for i in range(N):
        for j in range(i+1, N):
            for b in range(8):
                if Jbeta[b] != 0:
                    H += V[i,j]*Jbeta[b]*(Ls[i][b] @ Ls[j][b])
            if J38 != 0:
                H += V[i,j]*J38*(Ls[i][2] @ Ls[j][7] + Ls[i][7] @ Ls[j][2])
    return H, Ls

def coupling_matrix(N, kind="nn_pbc"):
    V = np.zeros((N, N))
    for i in range(N):
        for j in range(N):
            if i == j: continue
            r = abs(i-j)
            if kind == "dipolar":  V[i,j] = r**-3
            elif kind == "nn":     V[i,j] = 1.0 if r == 1 else 0.0
            elif kind == "nn_pbc": V[i,j] = 1.0 if (r == 1 or r == N-1) else 0.0
            else: raise ValueError(kind)
    return V

# ---------------------------------------------------------------------- pulse train
class PulseTrain:
    """Sequential square-pulse channels.  Channel alpha occupies (F_{alpha-1} T, F_alpha T]."""
    def __init__(self, a, f, T):
        self.a, self.f, self.T = np.asarray(a, float), np.asarray(f, float), float(T)
        if abs(self.f.sum() - 1) > 1e-12:
            raise ValueError("sum_alpha f_alpha must equal 1")
        self.omega = 2*np.pi/self.T
        self.Fcum = np.concatenate([[0.0], np.cumsum(self.f)])
        segs = []                                   # (start/T, end/T, g-vector)
        for al in range(8):
            if self.f[al] == 0: continue
            t0, fa, A = self.Fcum[al], self.f[al], self.a[al]
            for (s, e, sgn) in [(0, .25, 1), (.25, .75, -1), (.75, 1, 1)]:
                g = np.zeros(8); g[al] = sgn*A
                segs.append((t0+s*fa, t0+e*fa, g))
        self.segs = segs
    def G(self, t):
        """G_alpha(t) = omega int_0^t g_alpha ds  (vector over alpha); vanishes at t = nT."""
        tau = np.mod(t, self.T)/self.T
        G = np.zeros(8)
        for al in range(8):
            if self.f[al] == 0: continue
            t0, fa, A = self.Fcum[al], self.f[al], self.a[al]
            x = tau - t0
            if x <= 0: continue
            if x <= fa/4:      G[al] = A*x
            elif x <= 3*fa/4:  G[al] = A*(fa/2 - x)
            elif x <= fa:      G[al] = A*(x - fa)
        return self.omega*self.T*G

# ---------------------------------------------------------------------- propagators
class Propagators:
    """Exact propagator within one period as an ordered product of segment exponentials."""
    def __init__(self, H0, drive_gens, pulse):
        self.H0, self.D, self.pulse = H0, drive_gens, pulse
        self._cache = {}
    def Hseg(self, g):
        return self.H0 + self.pulse.omega*sum(g[a]*self.D[a] for a in range(8) if g[a] != 0)
    def U_within_cycle(self, tau_frac):
        key = round(float(tau_frac), 12)
        if key in self._cache: return self._cache[key]
        T = self.pulse.T
        U = np.eye(self.H0.shape[0], dtype=complex)
        for (s, e, g) in self.pulse.segs:
            if tau_frac <= s + 1e-14: break
            dt = (min(tau_frac, e) - s)*T
            U = expm(-1j*self.Hseg(g)*dt) @ U
            if tau_frac <= e + 1e-14: break
        self._cache[key] = U
        return U

def onsite_kick_unitary(Gvec, N, lam=LAM):
    """exp(-i K0(t)) = tensor product of identical on-site unitaries."""
    K1 = sum(Gvec[a]*lam[a] for a in range(8))
    w = expm(-1j*K1)
    return reduce(np.kron, [w]*N)

# ----------------------------------------------------------------------- observables
def spin1_ops():
    s = 1/np.sqrt(2)
    Sx = s*np.array([[0,1,0],[1,0,1],[0,1,0]], complex)
    Sy = s*np.array([[0,-1j,0],[1j,0,-1j],[0,1j,0]], complex)
    Sz = np.diag([1.,0.,-1.]).astype(complex)
    return Sx, Sy, Sz

def quadrupoles():
    Sx, Sy, Sz = spin1_ops()
    return Sz@Sz - (2/3)*np.eye(3), Sx@Sy + Sy@Sx, Sx@Sx - Sy@Sy   # Q0, Qxy, Qx2y2

# ---------------------------------------------------------------------- model setup
def model_params(model, a4=2.0, J3=1.0):
    Jex = 1.0
    if model == "nematic":
        J = np.array([0,0,J3,0,0,0,0,0.]); a = np.array([0,0,0,a4,0,0,0,0.]); f = np.array([0,0,0,1,0,0,0,0.])
    elif model == "su3":
        J = np.array([Jex,Jex,2*Jex/3,Jex,Jex,Jex,Jex,4*Jex/3]); a = np.array([4,4,0,0,0,0,0,0.]); f = np.array([.5,.5,0,0,0,0,0,0])
    else:
        raise ValueError(model)
    return J, a, f

def initial_state(init, N):
    e = np.eye(3)
    if   init == "zz_neel":   kets = [e[0] if i%2==0 else e[2] for i in range(N)]
    elif init == "p0_neel":   kets = [e[0] if i%2==0 else e[1] for i in range(N)]
    elif init == "all0":      kets = [e[1]]*N
    elif init == "xnem":      kets = [(e[0]+e[2])/np.sqrt(2)]*N
    elif init == "xnem_stag": kets = [(e[0]+(-1)**i*e[2])/np.sqrt(2) for i in range(N)]
    elif init == "tilted":    kets = [(e[0]+(-1)**i*e[1]+e[2])/np.sqrt(3) for i in range(N)]
    else: raise ValueError(init)
    return reduce(np.kron, kets).astype(complex)

# ------------------------------------------------------------------------------ run
def run(N=4, model="nematic", Vkind="nn_pbc", T=0.05, ncycles=60, pts_per_cycle=8,
        pair=(0,1), init="zz_neel", out="corr_vs_time.dat", verbose=True, a4=2.0, J3=1.0):
    J, a, f = model_params(model, a4, J3)
    V = coupling_matrix(N, Vkind)
    u, v = u_coef(a, f), v_coef(a, f)
    Jeff, J38 = effective_couplings(J, u, v)
    H0, Ls   = two_body(J, 0.0, V, N)
    Heff, _  = two_body(Jeff, J38, V, N)
    D = [sum(Ls[i][al] for i in range(N)) for al in range(8)]
    pulse = PulseTrain(a, f, T)
    prop  = Propagators(H0, D, pulse)

    assert np.allclose(Heff, Heff.conj().T) and np.allclose(H0, H0.conj().T)
    if verbose:
        print(f"N={N} model={model} V={Vkind} T={T} omega={2*np.pi/T:.4f} a4={a4} J3={J3}")
        print("u    =", np.round(u, 6)); print("v    =", np.round(v, 6))
        print("Jeff =", np.round(Jeff, 6), " J38eff =", round(J38, 6))
        print("||K0(T)|| =", np.linalg.norm(pulse.G(T)), "(should be 0)")

    psi0 = initial_state(init, N)
    Q0, Qxy, Qx2y2 = quadrupoles(); Sx, Sy, Sz = spin1_ops()
    i, j = pair
    obs = {k: site_op(O, i, N) @ site_op(O, j, N) for k, O in
           [("Q0Q0", Q0), ("QxyQxy", Qxy), ("SzSz", Sz), ("Qx2y2", Qx2y2)]}
    ev = lambda psi, O: float(np.real(np.vdot(psi, O @ psi)))

    w, Wv = eigh(Heff); c0 = Wv.conj().T @ psi0
    psi_eff = lambda t: Wv @ (np.exp(-1j*w*t)*c0)

    t0 = time.time()
    UT = prop.U_within_cycle(1.0)
    fracs = np.arange(pts_per_cycle)/pts_per_cycle
    Ufr = [prop.U_within_cycle(fr) for fr in fracs]
    assert np.allclose(UT.conj().T @ UT, np.eye(3**N), atol=1e-10)
    if verbose: print(f"propagators built in {time.time()-t0:.1f}s (dim {3**N})")

    ts, rows = [], []
    psi_n = psi0.copy()
    for n in range(ncycles+1):
        for k, fr in enumerate(fracs):
            if n == ncycles and k > 0: break
            t = (n+fr)*T
            pe = Ufr[k] @ psi_n
            pf = psi_eff(t)
            pk = onsite_kick_unitary(pulse.G(t), N) @ pf
            rows.append([ev(pe,obs["Q0Q0"]), ev(pf,obs["Q0Q0"]), ev(pk,obs["Q0Q0"]),
                         ev(pe,obs["QxyQxy"]), ev(pf,obs["QxyQxy"]), ev(pk,obs["QxyQxy"]),
                         np.linalg.norm(pe-pk), int(k==0),
                         ev(pe,obs["SzSz"]), ev(pf,obs["SzSz"]),
                         ev(pe,obs["Qx2y2"]), ev(pf,obs["Qx2y2"]),
                         ev(pk,obs["SzSz"]), ev(pk,obs["Qx2y2"])])
            ts.append(t)
        psi_n = UT @ psi_n
    ts, rows = np.array(ts), np.array(rows)

    cols = ["t", "Q0Q0_exact", "Q0Q0_eff", "Q0Q0_kick", "QxyQxy_exact", "QxyQxy_eff", "QxyQxy_kick",
            "state_err", "strobo", "SzSz_exact", "SzSz_eff", "Qx2y2_exact", "Qx2y2_eff",
            "SzSz_kick", "Qx2y2_kick"]
    hdr = (f"N={N} model={model} V={Vkind} T={T} omega={2*np.pi/T:.10g} ncycles={ncycles} ppc={pts_per_cycle} "
           f"pair={pair} init={init} a4={a4} J3={J3} version={__version__}\n"
           f"u={np.array2string(u, precision=10, separator=',')}\n"
           f"v={np.array2string(v, precision=10, separator=',')}\n"
           f"Jeff={np.array2string(Jeff, precision=10, separator=',')}  J38eff={J38:.10g}\n"
           "kick = effective state dressed with exp(-i K0(t)); strobo=1 at t=nT\n"
           + "  ".join(f"{c:>14s}" for c in cols))
    np.savetxt(out, np.column_stack([ts, rows]), fmt="%16.10e", header=hdr, delimiter="  ")
    if verbose:
        sm = rows[:,7] == 1
        print("max |exact-eff|  (strobo) Q0Q0: %.3e  QxyQxy: %.3e  SzSz: %.3e" %
              (np.abs(rows[sm,0]-rows[sm,1]).max(), np.abs(rows[sm,3]-rows[sm,4]).max(), np.abs(rows[sm,8]-rows[sm,9]).max()))
        print("max |exact-kick| (all t)  Q0Q0: %.3e  QxyQxy: %.3e  SzSz: %.3e" %
              (np.abs(rows[:,0]-rows[:,2]).max(), np.abs(rows[:,3]-rows[:,5]).max(), np.abs(rows[:,8]-rows[:,12]).max()))
        print("max ||psi_exact - psi_kick|| : %.3e" % rows[:,6].max())
        print("wrote", out)
    return ts, rows, cols

# -------------------------------------------------------------------------- selftest
def selftest():
    ok = True
    def chk(name, val, tol):
        nonlocal ok
        flag = val < tol; ok &= flag
        print(f"  [{'PASS' if flag else 'FAIL'}] {name}: {val:.2e}")
    Sx, Sy, Sz = spin1_ops(); Q0, Qxy, Qx2y2 = quadrupoles()
    print("(i) nematic: code Heff vs Eq. a4_only_eff / symmetry_a4_eff, N=4 ring")
    N = 4; V = coupling_matrix(N, "nn_pbc"); so = lambda O, k: site_op(O, k, N)
    for a4 in [0.7, 1.3, 2.0]:
        J, a, f = model_params("nematic", a4, 1.0)
        Je, J38 = effective_couplings(J, u_coef(a, f), v_coef(a, f)); Hc, _ = two_body(Je, J38, V, N)
        Jzz = 1/32*(1+np.sin(np.pi*a4)/(np.pi*a4)); Jz0 = 3/8*np.sin(np.pi*a4/2)/(np.pi*a4)
        J00 = 9/16; Jxy = 1/32*(1-np.sin(np.pi*a4)/(np.pi*a4))
        Hp = np.zeros_like(Hc)
        for i in range(N):
            for j in range(N):
                if i == j: continue
                Hp += 0.5*V[i,j]*(Jzz*so(Sz,i)@so(Sz,j) + Jz0*(so(Sz,i)@so(Q0,j)+so(Q0,i)@so(Sz,j))
                                  + J00*so(Q0,i)@so(Q0,j) + Jxy*so(Qxy,i)@so(Qxy,j))
        chk(f"a4={a4}", np.linalg.norm(Hc-Hp), 1e-12)
    print("(ii) su3: code Heff vs BLBQ Eq. su3_single_block_BLBQ and SU(3) invariance, N=3")
    N = 3; V = coupling_matrix(N, "nn"); so = lambda O, k: site_op(O, k, N)
    J, a, f = model_params("su3"); Je, J38 = effective_couplings(J, u_coef(a, f), v_coef(a, f))
    Hc, _ = two_body(Je, J38, V, N); Hb = np.zeros_like(Hc)
    for i in range(N):
        for j in range(N):
            if i == j: continue
            SS = sum(so(S,i)@so(S,j) for S in (Sx,Sy,Sz))
            Hb += 0.25*V[i,j]*(SS + SS@SS - 4/3*np.eye(3**N))
    chk("||Heff - BLBQ||", np.linalg.norm(Hc-Hb), 1e-12)
    chk("max ||[Heff, Lambda^b]||", max(np.linalg.norm(Hc@L-L@Hc) for L in
        [sum(so(LAM[b],k) for k in range(N)) for b in range(8)]), 1e-12)
    print("(iii) two-site propagator error eps(omega) ~ omega^-1")
    for model in ["nematic", "su3"]:
        J, a, f = model_params(model); N = 2; V = coupling_matrix(N, "nn")
        Je, J38 = effective_couplings(J, u_coef(a, f), v_coef(a, f))
        H0, Ls = two_body(J, 0, V, N); Heff, _ = two_body(Je, J38, V, N)
        D = [sum(Ls[i][al] for i in range(N)) for al in range(8)]
        tF, ph = 0.01, 0.1672; eps = []
        for nc in [2, 8, 32, 128]:
            T = tF/(nc+ph); pulse = PulseTrain(a, f, T); prop = Propagators(H0, D, pulse)
            Uex = prop.U_within_cycle(ph) @ np.linalg.matrix_power(prop.U_within_cycle(1.0), nc)
            Uk = onsite_kick_unitary(pulse.G(tF), N) @ expm(-1j*Heff*tF)
            eps.append((2*np.pi/T, np.linalg.norm(Uex-Uk, 2)))
        eps = np.array(eps); slope = np.polyfit(np.log(eps[:,0]), np.log(eps[:,1]), 1)[0]
        chk(f"{model}: |slope+1|", abs(slope+1), 1e-3)
    print("ALL PASS" if ok else "SOME CHECKS FAILED"); return ok

if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__.split("Usage")[0], formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--N", type=int, default=4)
    p.add_argument("--model", default="nematic", choices=["nematic", "su3"])
    p.add_argument("--V", default="nn_pbc", choices=["nn_pbc", "nn", "dipolar"])
    p.add_argument("--T", type=float, default=0.05)
    p.add_argument("--ncycles", type=int, default=60)
    p.add_argument("--ppc", type=int, default=8)
    p.add_argument("--pair", type=int, nargs=2, default=[0, 1])
    p.add_argument("--init", default="zz_neel", choices=["zz_neel", "p0_neel", "all0", "xnem", "xnem_stag", "tilted"])
    p.add_argument("--a4", type=float, default=2.0)
    p.add_argument("--J3", type=float, default=1.0)
    p.add_argument("--out", default="corr_vs_time.dat")
    p.add_argument("--selftest", action="store_true", help="run verification checks and exit")
    A = p.parse_args()
    if A.selftest:
        sys.exit(0 if selftest() else 1)
    run(A.N, A.model, A.V, A.T, A.ncycles, A.ppc, tuple(A.pair), A.init, A.out, a4=A.a4, J3=A.J3)
