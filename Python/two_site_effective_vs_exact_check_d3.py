#!/usr/bin/env python3
# ================================================================
# Two-site d=3 check of Ueff vs the exact propagator for the SU(d)
# square-pulse drive.
#
# Python port of two_site_effective_vs_exact_check_d3.wl.
# Requires: numpy, scipy.
#
# The exact propagator exploits the piecewise-constant drive: the
# time-ordered exponential factorizes exactly into ordinary matrix
# exponentials over the constant-drive segments, so no truncation
# parameter is involved.
# ================================================================

import numpy as np
from functools import lru_cache
from scipy.linalg import expm

NDIM = 3


# ---------- safe helper for two-site tensor couplings ----------
def tensor_coupling(a, M1, M2):
    M1 = np.asarray(M1)
    M2 = np.asarray(M2)
    if M1.ndim != 2 or M1.shape[0] != M1.shape[1]:
        raise ValueError("M1 is not square")
    if M2.ndim != 2 or M2.shape[0] != M2.shape[1]:
        raise ValueError("M2 is not square")
    if M1.shape != M2.shape:
        raise ValueError("M1 and M2 have different dimensions")
    return a * np.kron(M1, M2)


# ---------- SU(3) lambda matrices (S.A.D. basis) ----------
_LAMBDA_LIST = [
    0.5 * np.array([[0, 1, 0], [1, 0, 0], [0, 0, 0]], dtype=complex),
    0.5 * np.array([[0, 1j, 0], [-1j, 0, 0], [0, 0, 0]], dtype=complex),
    0.5 * np.array([[1, 0, 0], [0, -1, 0], [0, 0, 0]], dtype=complex),
    0.5 * np.array([[0, 0, 1], [0, 0, 0], [1, 0, 0]], dtype=complex),
    0.5 * np.array([[0, 0, 1j], [0, 0, 0], [-1j, 0, 0]], dtype=complex),
    0.5 * np.array([[0, 0, 0], [0, 0, 1], [0, 1, 0]], dtype=complex),
    0.5 * np.array([[0, 0, 0], [0, 0, 1j], [0, -1j, 0]], dtype=complex),
    0.5 * np.array([[1, 0, 0], [0, 0, 0], [0, 0, -1]], dtype=complex),
]


def lam(k):
    """SU(3) generator lambda^k, k = 1..8 (S.A.D. basis, same as the .wl)."""
    if not 1 <= k <= 8:
        raise ValueError("k out of range")
    return _LAMBDA_LIST[k - 1]


# ---------- removable-singularity-safe sc1 ----------
def sc1(x):
    return np.sinc(x / np.pi) - 1.0  # np.sinc(y) = sin(pi y)/(pi y), sinc(0) = 1


# ================================================================
# Bare couplings J_beta of H0
#   H0 = sum_{beta=1}^{8} J_beta lambda^beta (x) lambda^beta
# Two sites with V_ij = 1, so the 1/2 in (1/2) sum_{i!=j} cancels
# against the two ordered pairs.  Heff is built the same way.
# ================================================================

# H_0 parameters for the 2-pulse SU(3) symmetric model
Jex = 1.0
JVals = np.array([Jex, Jex, (2.0 / 3.0) * Jex, Jex, Jex, Jex, Jex,
                  (4.0 / 3.0) * Jex])

# H_0 parameters for the single-pulse nematic model
# JVals = np.array([0, 0, 1, 0, 0, 0, 0, 0], dtype=float)

J1, J2, J3, J4, J5, J6, J7, J8 = JVals

H0 = sum(tensor_coupling(JVals[b - 1], lam(b), lam(b)) for b in range(1, 9))

# ================================================================
# Pulse-dependent coefficients u_alpha, v_alpha
# NOTE: these depend on a_alpha and f_alpha ONLY -- not on omega.
# Heff is therefore omega-independent, which is what makes the
# omega scan below a clean test of the truncation order.
# ================================================================


def u_coef(a, f):
    return (f / 2.0) * sc1(np.pi * a * f)


def v_coef(a, f):
    return (f / 2.0) * sc1(np.pi * a * f / 2.0)


# ================================================================
# Effective couplings J^eff (d=3 S.A.D. basis), verbatim from LaTeX
# ================================================================


def Hij(u, v):
    """Effective two-site Hamiltonian from pulse coefficients.

    u, v: sequences of length 8, indexed u[0] = u_1, ..., u[7] = u_8.
    """
    u1, u2, u3, u4, u5, u6, u7, u8 = u
    v1, v2, v3, v4, v5, v6, v7, v8 = v

    J1e = (J1 * (1 + u2 + u3 + v4 + v5 + v6 + v7 + v8)
           - (u3 + v8) * J2 - u2 * J3
           - v7 * J4 - v6 * J5
           - v5 * J6 - v4 * J7 - 0.25 * u2 * J8)

    J2e = (J2 * (1 + u1 + u3 + v4 + v5 + v6 + v7 + v8)
           - (u3 + v8) * J1 - u1 * J3
           - v6 * J4 - v7 * J5
           - v4 * J6 - v5 * J7 - 0.25 * u1 * J8)

    J3e = (J3 * (1 + u1 + u2 + 0.25 * u6 + 0.25 * u7 + v6 + v7)
           - u2 * J1 - u1 * J2 - u7 * J6 - u6 * J7
           + (0.25 * u1 + 0.25 * u2 + 0.25 * u6 + 0.25 * u7
              - v1 - v2 - v6 - v7) * J8)

    J4e = (J4 * (1 + u5 + u8 + v1 + v2 + v3 + v6 + v7)
           - v7 * J1 - v6 * J2 - 0.25 * u5 * J3
           - (u8 + v3) * J5
           - v2 * J6 - v1 * J7 - u5 * J8)

    J5e = (J5 * (1 + u4 + u8 + v1 + v2 + v3 + v6 + v7)
           - v6 * J1 - v7 * J2 - 0.25 * u4 * J3
           - (u8 + v3) * J4
           - v1 * J6 - v2 * J7 - u4 * J8)

    J6e = (J6 * (1 + u7 + v1 + v2 + v3 + v4 + v5 + v8)
           - v5 * J1 - v4 * J2 - 0.25 * u7 * J3
           - v2 * J4 - v1 * J5
           - (v3 + v8) * J7 - 0.25 * u7 * J8)

    J7e = (J7 * (1 + u6 + v1 + v2 + v3 + v4 + v5 + v8)
           - v4 * J1 - v5 * J2 - 0.25 * u6 * J3
           - v1 * J4 - v2 * J5
           - (v3 + v8) * J6 - 0.25 * u6 * J8)

    J8e = (J8 * (1 + u4 + u5 + 0.25 * u6 + 0.25 * u7 + v6 + v7)
           + (0.25 * u4 + 0.25 * u5 + 0.25 * u6 + 0.25 * u7
              - v4 - v5 - v6 - v7) * J3
           - u5 * J4 - u4 * J5 - u7 * J6 - u6 * J7)

    J38e = ((-0.25 * u6 - 0.25 * u7 + v4 + v5) * J3
            + u7 * J6 + u6 * J7
            + (-0.25 * u6 - 0.25 * u7 + v1 + v2) * J8)

    Je = [J1e, J2e, J3e, J4e, J5e, J6e, J7e, J8e]
    H = sum(tensor_coupling(Je[b - 1], lam(b), lam(b)) for b in range(1, 9))
    H = H + tensor_coupling(J38e, lam(3), lam(8))
    H = H + tensor_coupling(J38e, lam(8), lam(3))
    return H


# ================================================================
# Pulse amplitudes a_alpha and time fractions f_alpha
# ================================================================

# Two-pulse parameters for driving into the SU(3) symmetric point
# aAmpVals  = np.array([4, 4, 0, 0, 0, 0, 0, 0], dtype=float)
# fFracVals = np.array([0.5, 0.5, 0, 0, 0, 0, 0, 0], dtype=float)

# One-pulse parameters for driving into the nematic model
aAmpVals = np.array([0, 0, 0, 2, 0, 0, 0, 0], dtype=float)
fFracVals = np.array([0, 0, 0, 1, 0, 0, 0, 0], dtype=float)

MAX_PULSES = len(aAmpVals)

if abs(fFracVals.sum() - 1.0) > 1e-12:
    print(f"WARNING: sum_alpha f_alpha = {fFracVals.sum()} (should be 1)")

Fcum = np.concatenate(([0.0], np.cumsum(fFracVals)))  # Fcum[al] for al = 0..8

# ---------- build Heff from the pulse coefficients ----------
uVals = np.array([u_coef(aAmpVals[al], fFracVals[al]) for al in range(MAX_PULSES)])
vVals = np.array([v_coef(aAmpVals[al], fFracVals[al]) for al in range(MAX_PULSES)])

Heff = Hij(uVals, vVals)


def Ueff(time):
    return expm(-1j * Heff * time)


# ================================================================
# Square pulse g_alpha(t) and its integral G_alpha(t)
#
# Segment membership uses right-closed intervals (t0a, t0a + fa/4],
# etc., matching the Piecewise tests in the .wl.
# ================================================================


def gsq(al, T, t):
    """g_alpha within one cycle; al = 1..8, t in [0, T)."""
    t0a = Fcum[al - 1] * T
    fa = fFracVals[al - 1] * T
    A = aAmpVals[al - 1]
    if t0a < t <= t0a + fa / 4:
        return A
    if t0a + fa / 4 < t <= t0a + 3 * fa / 4:
        return -A
    if t0a + 3 * fa / 4 < t <= t0a + fa:
        return A
    return 0.0


def Gsq(al, T, t, omega):
    t0a = Fcum[al - 1] * T
    fa = fFracVals[al - 1] * T
    A = aAmpVals[al - 1]
    if t <= t0a:
        return 0.0
    if t0a < t <= t0a + fa / 4:
        return omega * A * (t - t0a)
    if t0a + fa / 4 < t <= t0a + 3 * fa / 4:
        return omega * A * (fa / 2 - (t - t0a))
    if t0a + 3 * fa / 4 < t <= t0a + fa:
        return omega * A * (t - (t0a + fa))
    return 0.0


def total_g_pulse(al, T, t, omega):
    return omega * gsq(al, T, t % T)


def total_G_pulse(al, T, t, omega):
    return Gsq(al, T, t % T, omega)


# ================================================================
# Driven Hamiltonian and memoized per-segment propagator
# ================================================================

ID3 = np.eye(NDIM, dtype=complex)

drive_generators = [
    tensor_coupling(1.0, lam(al), ID3) + tensor_coupling(1.0, ID3, lam(al))
    for al in range(1, MAX_PULSES + 1)
]


def drive_tuple(t, T):
    om = 2 * np.pi / T
    return np.array([total_g_pulse(al, T, t, om) for al in range(1, MAX_PULSES + 1)])


@lru_cache(maxsize=None)
def _step_prop_rounded(gkey, dt):
    H = H0.copy()
    for al in range(MAX_PULSES):
        if gkey[al] != 0.0:
            H = H + gkey[al] * drive_generators[al]
    return expm(-1j * H * dt)


def unitary_step(t, deltat, T):
    gkey = tuple(np.round(drive_tuple(t, T), 9))
    return _step_prop_rounded(gkey, float(deltat))


# retained only to validate UExact
def U_trotter_TO(t, num_trotter_steps, t0, T):
    dt = (t - t0) / num_trotter_steps
    U = np.eye(NDIM**2, dtype=complex)
    for k in range(num_trotter_steps):
        tk = t0 + k * dt
        U = unitary_step(tk, dt, T) @ U
    return U


# ================================================================
#  EXACT PROPAGATOR
#
#   U^exact_{t0->t} = T-exp[ -i int_{t0}^{t} ds H(s) ]
#                   = prod_{k=m}^{1} exp[ -i H(t_k) delta t_k ]
#
# where t_0 < t_1 < ... < t_m = t are t0, t, and every pulse-edge
# time in between, and delta t_k = t_k - t_{k-1}.  H is constant on
# each (t_{k-1}, t_k], so every factor is an ordinary matrix
# exponential and the time ordering across intervals is exact.
#
# H is sampled at the segment MIDPOINT rather than at the right
# endpoint t_k.  The two are identical, because H is constant on
# the interval, but the midpoint avoids floating-point ties at the
# interval endpoints in the segment tests of gsq.
# ================================================================

# pulse-edge times within one cycle, as fractions of T
edge_fracs = sorted(set(
    f
    for al in range(1, MAX_PULSES + 1)
    for f in (Fcum[al - 1],
              Fcum[al - 1] + fFracVals[al - 1] / 4,
              Fcum[al - 1] + 3 * fFracVals[al - 1] / 4)
))


def break_times(t0, t, Tc):
    tol = 1e-10 * Tc
    nlo = int(np.floor(t0 / Tc)) - 1
    nhi = int(np.ceil(t / Tc)) + 1
    cand = sorted(
        (n + ef) * Tc
        for n in range(nlo, nhi + 1)
        for ef in edge_fracs
        if t0 + tol < (n + ef) * Tc < t - tol
    )
    # merge candidates closer than tol (mirrors Union[..., SameTest -> ...])
    merged = []
    for c in cand:
        if not merged or abs(c - merged[-1]) >= tol:
            merged.append(c)
    return [t0] + merged + [t]


def U_exact(t, t0, Tc):
    pts = break_times(t0, t, Tc)
    U = np.eye(NDIM**2, dtype=complex)
    for k in range(1, len(pts)):
        U = unitary_step((pts[k - 1] + pts[k]) / 2, pts[k] - pts[k - 1], Tc) @ U
    return U


# ----------------------------------------------------------------
# Validation helper: uniform grid sampled at the step MIDPOINT.
#
# WARNING about U_trotter_TO above: it samples H at the LEFT endpoint
# t0 + k dt, which is inconsistent with gsq's right-closed intervals --
# the left endpoint of a segment belongs to the PREVIOUS segment.
# That biases the result by one sub-step at every pulse edge, an O(dt)
# error that persists even on a grid whose points coincide with the
# edges.  Right-endpoint sampling fails too, for the mirror-image
# reason plus floating-point ties exactly at the edges.  Midpoint
# sampling is unbiased, and on a commensurate grid it reproduces
# U_exact to machine precision for ANY M.
# ----------------------------------------------------------------


def U_uniform_mid(t, t0, Tc, M):
    dt = (t - t0) / M
    U = np.eye(NDIM**2, dtype=complex)
    for k in range(M):
        U = unitary_step(t0 + (k + 0.5) * dt, dt, Tc) @ U
    return U


# ================================================================
# Kick operator K^(0)(t) and effective propagator
# ================================================================


def kick(t, T):
    om = 2 * np.pi / T
    K = np.zeros((NDIM**2, NDIM**2), dtype=complex)
    for al in range(1, MAX_PULSES + 1):
        G = total_G_pulse(al, T, t, om)
        if G != 0.0:
            K = K + G * drive_generators[al - 1]
    return K


def Ueff_kick(t, t0, T):
    return expm(-1j * kick(t, T)) @ Ueff(t - t0) @ expm(1j * kick(t0, T))


# ---------- error metric: no truncation parameter ----------
def epsilon(t, t0, Tc):
    return np.linalg.norm(Ueff_kick(t, t0, Tc) - U_exact(t, t0, Tc), 2)


# also useful: the same with the kick omitted
def epsilon_no_kick(t, t0, Tc):
    return np.linalg.norm(Ueff(t - t0) - U_exact(t, t0, Tc), 2)


# ================================================================
# Parameters
# ================================================================

def main():
    t0 = 0.0

    # --- scan 1: fixed physical time, omega varied ---
    # omega is varied by varying Tcyc.  To keep the sampling PHASE
    # within the final cycle fixed as omega changes, Tcyc is chosen so
    # that tFixed/Tcyc = Ncyc + phaseFrac with Ncyc an integer.  With a
    # uniform omega grid instead, frac(t/Tcyc) wanders and the curve
    # acquires large scatter from the micromotion phase.
    tFixed = 0.01
    phaseFrac = 0.1672
    NcycList = [2, 4, 8, 16, 32, 64, 128]

    # --- scan 2: fixed omega, time varied ---
    TcFixed = 0.005
    ptsPerCycle = 32
    nCyclesScan = 2

    print("=== Parameter summary ===")
    print("JVals     =", JVals)
    print("aAmpVals  =", aAmpVals)
    print("fFracVals =", fFracVals)
    print("u_alpha   =", np.array2string(uVals, precision=8))
    print("v_alpha   =", np.array2string(vVals, precision=8))
    print(f"edgeFracs = {edge_fracs}  -> {len(edge_fracs)} segments per cycle")
    print(f"scan 1: tFixed = {tFixed}, phaseFrac = {phaseFrac}, "
          f"NcycList = {NcycList}")
    print(f"scan 2: TcFixed = {TcFixed} (omega = {2 * np.pi / TcFixed}), "
          f"{ptsPerCycle} points/cycle over {nCyclesScan} cycles")

    # ---------- sanity checks ----------
    print("=== Sanity checks ===")
    print("Heff Hermitian? ",
          np.linalg.norm(Heff - Heff.conj().T, 2) < 1e-10)
    print("Heff traceless? ", abs(np.trace(Heff)) < 1e-10)
    print("||K(Tcyc)|| ~ 0 at a cycle boundary:",
          np.linalg.norm(kick(TcFixed, TcFixed), 2))
    z = np.zeros(8)
    print("Off-pulse limit Heff -> H0:",
          np.linalg.norm(Hij(z, z) - H0, 2) < 1e-12)

    # ---------- U_exact validation ----------
    # (a) against a midpoint-sampled uniform grid at t = 2T.  Every edge
    #     fraction is a multiple of 1/8 (SU(3) case) or 1/4 (nematic),
    #     so any M divisible by 16 puts grid points on all edges and the
    #     two constructions must agree to machine precision.
    # (b) the cycle-power identity.
    print("UExact vs UUniformMid, t = 2T, M = 4096  (expect ~1e-13):",
          np.linalg.norm(U_exact(2 * TcFixed, 0, TcFixed)
                         - U_uniform_mid(2 * TcFixed, 0, TcFixed, 4096), 2))
    print("UExact vs UUniformMid, t = 2T, M = 65536 (expect ~1e-12):",
          np.linalg.norm(U_exact(2 * TcFixed, 0, TcFixed)
                         - U_uniform_mid(2 * TcFixed, 0, TcFixed, 65536), 2))
    print("||UExact[0->5T] - UExact[0->T]^5|| (expect ~1e-14):",
          np.linalg.norm(U_exact(5 * TcFixed, 0, TcFixed)
                         - np.linalg.matrix_power(
                             U_exact(TcFixed, 0, TcFixed), 5), 2))
    print("||Ueff(no kick) - UExact|| at t = 2.1672 T (expect O(1) failure):",
          epsilon_no_kick(2.1672 * TcFixed, 0, TcFixed))
    print("epsilon WITH kick at t = 2.1672 T:",
          epsilon(2.1672 * TcFixed, 0, TcFixed))
    # For contrast, the original left-endpoint Trotter reference at the
    # same commensurate M is biased at O(dt):
    print("UExact vs UTrotterTO (left endpoint), t = 2T, M = 4096:",
          np.linalg.norm(U_exact(2 * TcFixed, 0, TcFixed)
                         - U_trotter_TO(2 * TcFixed, 4096, 0, TcFixed), 2))

    # ================================================================
    # SCAN 1: epsilon versus omega at fixed t
    # ================================================================

    print("=== SCAN 1: epsilon vs omega at fixed t ===")
    omega_scan = []
    for nc in NcycList:
        Tn = tFixed / (nc + phaseFrac)
        omega_scan.append((2 * np.pi / Tn, epsilon(tFixed, t0, Tn)))
    for row in omega_scan:
        print(row)
    logw = np.log([r[0] for r in omega_scan])
    loge = np.log([r[1] for r in omega_scan])
    slope, intercept = np.polyfit(logw, loge, 1)
    print("log-log slope d(log eps)/d(log omega) =", slope)

    # ================================================================
    # SCAN 2: epsilon versus t at fixed omega
    # ================================================================

    print("=== SCAN 2: epsilon vs t at fixed omega ===")
    time_scan = []
    for n in range(1, nCyclesScan * ptsPerCycle + 1):
        time_scan.append((n / ptsPerCycle,
                          epsilon(n * TcFixed / ptsPerCycle, t0, TcFixed)))
    for row in time_scan:
        print(row)

    # ================================================================
    # Export plot data
    # ================================================================

    print("=== Writing .dat files ===")

    def export_dat(fname, data, col_header):
        with open(fname, "w") as fh:
            fh.write("#\t" + col_header + "\n")
            for x, y in data:
                fh.write(f"{x!r}\t{y!r}\n")
        print("  wrote", fname)

    export_dat("epsilon_vs_omega.dat", omega_scan, "omega\tnorm_diff")
    export_dat("epsilon_vs_time.dat", time_scan, "t_over_T\tnorm_diff")

    om = 2 * np.pi / TcFixed
    with open("pulse_profile.dat", "w") as fh:
        fh.write("#t\t" + "\t".join(f"g{al}" for al in range(1, MAX_PULSES + 1))
                 + "\n")
        for n in range(401):
            tt = n * (2 * TcFixed) / 400
            gvals = [total_g_pulse(al, TcFixed, tt, om)
                     for al in range(1, MAX_PULSES + 1)]
            fh.write("\t".join(repr(float(x)) for x in [tt] + gvals) + "\n")
    print("  wrote pulse_profile.dat")


if __name__ == "__main__":
    main()
