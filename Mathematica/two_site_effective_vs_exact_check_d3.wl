(* ::Package:: *)

ClearAll["Global`*"];

(* ================================================================ *)
(* Two-site d=3 check of Ueff vs the exact propagator for the SU(d)  *)
(* square-pulse drive.                                              *)
(*                                                                  *)
(* ================================================================ *)

ndim = 3;

(* ---------- safe helper for two-site tensor couplings ---------- *)
ClearAll[TensorCoupling];
TensorCoupling[a_?NumericQ | a_Symbol, M1_?MatrixQ, M2_?MatrixQ] := Module[
  {d1 = Dimensions[M1], d2 = Dimensions[M2]},
  If[Length[d1] != 2 || d1[[1]] != d1[[2]], Return[$Failed]];
  If[Length[d2] != 2 || d2[[1]] != d2[[2]], Return[$Failed]];
  If[d1 =!= d2, Return[$Failed]];
  a KroneckerProduct[M1, M2]
];

(* ---------- SU(3) lambda matrices (S.A.D. basis) ---------- *)
ClearAll[lambda];
$LambdaList = {
  (1/2) {{0, 1, 0}, {1, 0, 0}, {0, 0, 0}},
  (1/2) {{0, I, 0}, {-I, 0, 0}, {0, 0, 0}},
  (1/2) {{1, 0, 0}, {0, -1, 0}, {0, 0, 0}},
  (1/2) {{0, 0, 1}, {0, 0, 0}, {1, 0, 0}},
  (1/2) {{0, 0, I}, {0, 0, 0}, {-I, 0, 0}},
  (1/2) {{0, 0, 0}, {0, 0, 1}, {0, 1, 0}},
  (1/2) {{0, 0, 0}, {0, 0, I}, {0, -I, 0}},
  (1/2) {{1, 0, 0}, {0, 0, 0}, {0, 0, -1}}
};
lambda[k_Integer] /; 1 <= k <= 8 := $LambdaList[[k]];

(* ---------- removable-singularity-safe sc1 ---------- *)
ClearAll[sc1];
sc1[x_] := Piecewise[{{0, x == 0}}, Sin[x]/x - 1];

(* ================================================================ *)
(* Bare couplings J_beta of H0                                      *)
(*   H0 = sum_{beta=1}^{8} J_beta lambda^beta (x) lambda^beta        *)
(* Two sites with V_ij = 1, so the 1/2 in (1/2) sum_{i!=j} cancels   *)
(* against the two ordered pairs.  Heff is built the same way.       *)
(* ================================================================ *)

(* H_0 parameters for the 2-pulse SU(3) symmetric model *)
Jex = 1;
JVals = {Jex, Jex, (2/3) Jex, Jex, Jex, Jex, Jex, (4/3) Jex};
{J1, J2, J3, J4, J5, J6, J7, J8} = JVals;

(* H_0 parameters for the single-pulse nematic model *)
(*
JVals = {0, 0, 1, 0, 0, 0, 0, 0};
{J1, J2, J3, J4, J5, J6, J7, J8} = JVals;
*)

H0 = Sum[TensorCoupling[JVals[[b]], lambda[b], lambda[b]], {b, 1, 8}];

(* ================================================================ *)
(* Pulse-dependent coefficients u_alpha, v_alpha                    *)
(* NOTE: these depend on a_alpha and f_alpha ONLY -- not on omega.   *)
(* Heff is therefore omega-independent, which is what makes the      *)
(* omega scan below a clean test of the truncation order.            *)
(* ================================================================ *)

ClearAll[uCoef, vCoef];
uCoef[a_, f_] := (f/2) sc1[Pi a f];
vCoef[a_, f_] := (f/2) sc1[Pi a f/2];

(* ================================================================ *)
(* Effective couplings J^eff (d=3 S.A.D. basis), verbatim from LaTeX *)
(* ================================================================ *)

ClearAll[Hij];
Hij[
   u1_, u2_, u3_, u4_, u5_, u6_, u7_, u8_,
   v1_, v2_, v3_, v4_, v5_, v6_, v7_, v8_
  ] := Module[
  {J1e, J2e, J3e, J4e, J5e, J6e, J7e, J8e, J38e},

  J1e = J1 (1 + u2 + u3 + v4 + v5 + v6 + v7 + v8)
        - (u3 + v8) J2 - u2 J3
        - v7 J4 - v6 J5
        - v5 J6 - v4 J7 - (1/4) u2 J8;

  J2e = J2 (1 + u1 + u3 + v4 + v5 + v6 + v7 + v8)
        - (u3 + v8) J1 - u1 J3
        - v6 J4 - v7 J5
        - v4 J6 - v5 J7 - (1/4) u1 J8;

  J3e = J3 (1 + u1 + u2 + (1/4) u6 + (1/4) u7 + v6 + v7)
        - u2 J1 - u1 J2 - u7 J6 - u6 J7
        + ((1/4) u1 + (1/4) u2 + (1/4) u6 + (1/4) u7 - v1 - v2 - v6 - v7) J8;

  J4e = J4 (1 + u5 + u8 + v1 + v2 + v3 + v6 + v7)
        - v7 J1 - v6 J2 - (1/4) u5 J3
        - (u8 + v3) J5
        - v2 J6 - v1 J7 - u5 J8;

  J5e = J5 (1 + u4 + u8 + v1 + v2 + v3 + v6 + v7)
        - v6 J1 - v7 J2 - (1/4) u4 J3
        - (u8 + v3) J4
        - v1 J6 - v2 J7 - u4 J8;

  J6e = J6 (1 + u7 + v1 + v2 + v3 + v4 + v5 + v8)
        - v5 J1 - v4 J2 - (1/4) u7 J3
        - v2 J4 - v1 J5
        - (v3 + v8) J7 - (1/4) u7 J8;

  J7e = J7 (1 + u6 + v1 + v2 + v3 + v4 + v5 + v8)
        - v4 J1 - v5 J2 - (1/4) u6 J3
        - v1 J4 - v2 J5
        - (v3 + v8) J6 - (1/4) u6 J8;

  J8e = J8 (1 + u4 + u5 + (1/4) u6 + (1/4) u7 + v6 + v7)
        + ((1/4) u4 + (1/4) u5 + (1/4) u6 + (1/4) u7 - v4 - v5 - v6 - v7) J3
        - u5 J4 - u4 J5 - u7 J6 - u6 J7;

  J38e = (-(1/4) u6 - (1/4) u7 + v4 + v5) J3
        + u7 J6 + u6 J7
        + (-(1/4) u6 - (1/4) u7 + v1 + v2) J8;

  TensorCoupling[J1e, lambda[1], lambda[1]] +
  TensorCoupling[J2e, lambda[2], lambda[2]] +
  TensorCoupling[J3e, lambda[3], lambda[3]] +
  TensorCoupling[J4e, lambda[4], lambda[4]] +
  TensorCoupling[J5e, lambda[5], lambda[5]] +
  TensorCoupling[J6e, lambda[6], lambda[6]] +
  TensorCoupling[J7e, lambda[7], lambda[7]] +
  TensorCoupling[J8e, lambda[8], lambda[8]] +
  TensorCoupling[J38e, lambda[3], lambda[8]] +
  TensorCoupling[J38e, lambda[8], lambda[3]]
];

(* ================================================================ *)
(* Pulse amplitudes a_alpha and time fractions f_alpha              *)
(* ================================================================ *)

(* Two-pulse parameters for driving into the SU(3) symmetric point *)
(*
aAmpVals  = {4, 4, 0, 0, 0, 0, 0, 0};
fFracVals = {1/2, 1/2, 0, 0, 0, 0, 0, 0};
*)

(* One-pulse parameters for driving into the nematic model *)

aAmpVals  = {0, 0, 0, 2, 0, 0, 0, 0};
fFracVals = {0, 0, 0, 1, 0, 0, 0, 0};


maxnumberpulses = Length[aAmpVals];
aAmp[al_Integer]  /; 1 <= al <= maxnumberpulses := aAmpVals[[al]];
fFrac[al_Integer] /; 1 <= al <= maxnumberpulses := fFracVals[[al]];

If[Abs[Total[fFracVals] - 1] > 10^-12,
  Print["WARNING: sum_alpha f_alpha = ", Total[fFracVals], " (should be 1)"]];

Fcum[0] = 0;
Fcum[al_Integer] /; 1 <= al <= maxnumberpulses := Sum[fFracVals[[k]], {k, 1, al}];

(* ---------- build Heff from the pulse coefficients ---------- *)
uVals = Table[uCoef[aAmpVals[[al]], fFracVals[[al]]], {al, 1, 8}];
vVals = Table[vCoef[aAmpVals[[al]], fFracVals[[al]]], {al, 1, 8}];

HeffBuild = Hij @@ Join[uVals, vVals];
Heff = N[HeffBuild];   (* machine precision, to match UExact below *)

ClearAll[Ueff];
Ueff[time_] := MatrixExp[-I Heff time];

(* ================================================================ *)
(* Square pulse g_alpha(t) and its integral G_alpha(t)              *)
(* ================================================================ *)

ClearAll[gsq, Gsq, TotalgPulse, TotalGPulse];

gsq[al_Integer, T_?NumericQ, t_] := Module[
  {t0a = Fcum[al - 1] T, fa = fFrac[al] T, A = aAmp[al]},
  Piecewise[
    {
     {A,  t0a          < t <= t0a + fa/4},
     {-A, t0a + fa/4   < t <= t0a + 3 fa/4},
     {A,  t0a + 3 fa/4 < t <= t0a + fa}
    },
    0]
];

Gsq[al_Integer, T_?NumericQ, t_, omega_?NumericQ] := Module[
  {t0a = Fcum[al - 1] T, fa = fFrac[al] T, A = aAmp[al]},
  Piecewise[
    {
     {0, t <= t0a},
     {omega A (t - t0a),            t0a          < t <= t0a + fa/4},
     {omega A (fa/2 - (t - t0a)),   t0a + fa/4   < t <= t0a + 3 fa/4},
     {omega A (t - (t0a + fa)),     t0a + 3 fa/4 < t <= t0a + fa}
    },
    0]
];

TotalgPulse[al_Integer, T_?NumericQ, t_, omega_?NumericQ] := omega gsq[al, T, Mod[t, T, 0]];
TotalGPulse[al_Integer, T_?NumericQ, t_, omega_?NumericQ] := Gsq[al, T, Mod[t, T, 0], omega];

(* ================================================================ *)
(* Driven Hamiltonian and memoized per-segment propagator           *)
(* ================================================================ *)

ClearAll[DrivingValpha, driveTuple, driveGenerators, Hexact,
         stepPropRounded, Unitarystep, UTrotterTO];

driveGenerators = Table[
   TensorCoupling[1.0, lambda[al], IdentityMatrix[{ndim, ndim}]] +
   TensorCoupling[1.0, IdentityMatrix[{ndim, ndim}], lambda[al]],
   {al, 1, maxnumberpulses}
];

DrivingValpha[al_Integer, t_, T_?NumericQ, omega_?NumericQ] :=
  TotalgPulse[al, T, t, omega] driveGenerators[[al]];

Hexact[t_, T_?NumericQ] := H0 + Sum[DrivingValpha[nn, t, T, 2 Pi/T], {nn, 1, maxnumberpulses}];

driveTuple[t_, T_?NumericQ] :=
  Table[TotalgPulse[al, T, t, 2 Pi/T], {al, 1, maxnumberpulses}];

stepPropRounded[gkey_List, dt_?NumericQ] := stepPropRounded[gkey, dt] =
  MatrixExp[-I (H0 + Sum[gkey[[al]] driveGenerators[[al]], {al, 1, maxnumberpulses}]) dt];

Unitarystep[t_, deltat_, T_?NumericQ] :=
  stepPropRounded[Round[N[driveTuple[t, T]], 10^-9], N[deltat]];

(* retained only to validate UExact *)
UTrotterTO[t_?NumericQ, NumTrottersteps_Integer?Positive, t0_?NumericQ, T_?NumericQ] := Module[
  {dt, U, k, tk},
  dt = N[(t - t0)/NumTrottersteps];
  U = IdentityMatrix[ndim^2];
  Do[
    tk = N[t0 + k dt];
    U = Unitarystep[tk, dt, T] . U,
    {k, 0, NumTrottersteps - 1}
  ];
  U
];

(* ================================================================ *)
(*  EXACT PROPAGATOR                                                *)
(*                                                                  *)
(*   U^exact_{t0->t} = T-exp[ -i int_{t0}^{t} ds H(s) ]              *)
(*                   = prod_{k=m}^{1} exp[ -i H(t_k) delta t_k ]     *)
(*                                                                  *)
(* where t_0 < t_1 < ... < t_m = t are t0, t, and every pulse-edge   *)
(* time in between, and delta t_k = t_k - t_{k-1}.  H is constant on *)
(* each (t_{k-1}, t_k], so every factor is an ordinary matrix        *)
(* exponential and the time ordering across intervals is exact.      *)
(*                                                                  *)
(* In code H is sampled at the segment MIDPOINT rather than at the   *)
(* right endpoint t_k.  The two are identical, because H is constant *)
(* on the interval, but the midpoint avoids floating-point ties at   *)
(* the interval endpoints in the Piecewise tests of gsq.             *)
(* ================================================================ *)

(* pulse-edge times within one cycle, as exact fractions of T *)
edgeFracs = Union@Flatten@Table[
   {Fcum[al - 1], Fcum[al - 1] + fFrac[al]/4, Fcum[al - 1] + 3 fFrac[al]/4},
   {al, 1, maxnumberpulses}];

ClearAll[breakTimes, UExact, epsilon];

breakTimes[t0_?NumericQ, t_?NumericQ, Tc_?NumericQ] := Module[
  {tol, nlo, nhi, cand},
  tol = 10^-10 Tc;
  nlo = Floor[t0/Tc] - 1;
  nhi = Ceiling[t/Tc] + 1;
  cand = Flatten@Table[N[(n + ef) Tc], {n, nlo, nhi}, {ef, edgeFracs}];
  cand = Select[cand, t0 + tol < # < t - tol &];
  Join[{N[t0]}, Union[cand, SameTest -> (Abs[#1 - #2] < tol &)], {N[t]}]
];

UExact[t_?NumericQ, t0_?NumericQ, Tc_?NumericQ] := Module[
  {pts, U, k},
  pts = breakTimes[t0, t, Tc];
  U = IdentityMatrix[ndim^2];
  Do[
    U = Unitarystep[(pts[[k - 1]] + pts[[k]])/2, pts[[k]] - pts[[k - 1]], Tc] . U,
    {k, 2, Length[pts]}
  ];
  U
];

(* ---------------------------------------------------------------- *)
(* Validation helper: uniform grid sampled at the step MIDPOINT.     *)
(*                                                                  *)
(* WARNING about UTrotterTO above: it samples H at the LEFT endpoint *)
(* t0 + k dt, which is inconsistent with gsq's right-closed          *)
(* intervals -- the left endpoint of a segment belongs to the        *)
(* PREVIOUS segment.  That biases the result by one sub-step at every *)
(* pulse edge, an O(dt) error that persists even on a grid whose      *)
(* points coincide with the edges.  Right-endpoint sampling fails too, *)
(* for the mirror-image reason plus floating-point ties exactly at    *)
(* the edges.  Midpoint sampling is unbiased, and on a commensurate   *)
(* grid it reproduces UExact to machine precision for ANY M.          *)
(* ---------------------------------------------------------------- *)
ClearAll[UUniformMid];
UUniformMid[t_?NumericQ, t0_?NumericQ, Tc_?NumericQ, M_Integer?Positive] := Module[
  {dt, U, k},
  dt = N[(t - t0)/M];
  U = IdentityMatrix[ndim^2];
  Do[U = Unitarystep[N[t0 + (k + 1/2) dt], dt, Tc] . U, {k, 0, M - 1}];
  U
];

(* ================================================================ *)
(* Kick operator K^(0)(t) and effective propagator                  *)
(* ================================================================ *)

ClearAll[Kick, UeffKick];
Kick[t_, T_?NumericQ] := Sum[
  TotalGPulse[al, T, t, 2 Pi/T] (
    TensorCoupling[1.0, lambda[al], IdentityMatrix[{ndim, ndim}]] +
    TensorCoupling[1.0, IdentityMatrix[{ndim, ndim}], lambda[al]]
  ),
  {al, 1, maxnumberpulses}
];

UeffKick[t_, t0_, T_?NumericQ] :=
  MatrixExp[-I Kick[t, T]] . Ueff[t - t0] . MatrixExp[I Kick[t0, T]];

(* ---------- error metric: no truncation parameter ---------- *)
epsilon[t_?NumericQ, t0_?NumericQ, Tc_?NumericQ] :=
  Norm[UeffKick[t, t0, Tc] - UExact[t, t0, Tc]];

(* also useful: the same with the kick omitted *)
ClearAll[epsilonNoKick];
epsilonNoKick[t_?NumericQ, t0_?NumericQ, Tc_?NumericQ] :=
  Norm[Ueff[t - t0] - UExact[t, t0, Tc]];

(* ================================================================ *)
(* Parameters                                                       *)
(* ================================================================ *)

t0 = 0;

(* --- scan 1: fixed physical time, omega varied --- *)
(* omega is varied by varying Tcyc.  To keep the sampling PHASE      *)
(* within the final cycle fixed as omega changes, Tcyc is chosen so  *)
(* that tFixed/Tcyc = Ncyc + phaseFrac with Ncyc an integer.  With a *)
(* uniform omega grid instead, frac(t/Tcyc) wanders and the curve    *)
(* acquires large scatter from the micromotion phase.                *)
tFixed    = 0.01;
phaseFrac = 0.1672;
NcycList  = {2, 4, 8, 16, 32, 64, 128};

(* --- scan 2: fixed omega, time varied --- *)
TcFixed      = 0.005;
ptsPerCycle  = 32;
nCyclesScan  = 2;

Print["=== Parameter summary ==="];
Print["JVals     = ", JVals];
Print["aAmpVals  = ", aAmpVals];
Print["fFracVals = ", fFracVals];
Print["u_alpha   = ", N[uVals, 8]];
Print["v_alpha   = ", N[vVals, 8]];
Print["edgeFracs = ", edgeFracs, "  -> ", Length[edgeFracs], " segments per cycle"];
Print["scan 1: tFixed = ", tFixed, ", phaseFrac = ", phaseFrac, ", NcycList = ", NcycList];
Print["scan 2: TcFixed = ", TcFixed, " (omega = ", N[2 Pi/TcFixed], "), ",
      ptsPerCycle, " points/cycle over ", nCyclesScan, " cycles"];

(* ---------- sanity checks ---------- *)
Print["=== Sanity checks ==="];
Print["Heff Hermitian?  ", Norm[Heff - ConjugateTranspose[Heff]] < 10^-10];
Print["Heff traceless?  ", Abs[Tr[Heff]] < 10^-10];
Print["||K(Tcyc)|| ~ 0 at a cycle boundary: ", N[Norm[Kick[TcFixed, TcFixed]], 12]];
Print["Off-pulse limit Heff -> H0: ",
  Module[{u0 = ConstantArray[0, 8]},
   Norm[N[Hij @@ Join[u0, u0]] - N[H0]] < 10^-12]];

(* ---------- UExact validation ---------- *)
(* (a) against a midpoint-sampled uniform grid at t = 2T.  Every edge  *)
(*     fraction is a multiple of 1/8 (SU(3) case) or 1/4 (nematic),    *)
(*     so any M divisible by 16 puts grid points on all edges and the  *)
(*     two constructions must agree to machine precision.              *)
(* (b) the cycle-power identity.                                       *)
Print["UExact vs UUniformMid, t = 2T, M = 4096  (expect ~1e-13): ",
  Norm[UExact[2 TcFixed, 0, TcFixed] - UUniformMid[2 TcFixed, 0, TcFixed, 4096]]];
Print["UExact vs UUniformMid, t = 2T, M = 65536 (expect ~1e-12): ",
  Norm[UExact[2 TcFixed, 0, TcFixed] - UUniformMid[2 TcFixed, 0, TcFixed, 65536]]];
Print["||UExact[0->5T] - UExact[0->T]^5|| (expect ~1e-14): ",
  Norm[UExact[5 TcFixed, 0, TcFixed] - MatrixPower[UExact[TcFixed, 0, TcFixed], 5]]];
Print["||Ueff(no kick) - UExact|| at t = 2.1672 T (expect O(1) failure): ",
  epsilonNoKick[2.1672 TcFixed, 0, TcFixed]];
Print["epsilon WITH kick at t = 2.1672 T: ", epsilon[2.1672 TcFixed, 0, TcFixed]];
(* For contrast, the original left-endpoint Trotter reference at the    *)
(* same commensurate M is biased at O(dt): *)
Print["UExact vs UTrotterTO (left endpoint), t = 2T, M = 4096: ",
  Norm[UExact[2 TcFixed, 0, TcFixed] - UTrotterTO[2 TcFixed, 4096, 0, TcFixed]]];

(* ================================================================ *)
(* SCAN 1: epsilon versus omega at fixed t                          *)
(* ================================================================ *)

Print["=== SCAN 1: epsilon vs omega at fixed t ==="];
omegaScan = Table[
   Module[{Tn = tFixed/(nc + phaseFrac)},
     {N[2 Pi/Tn], N[epsilon[tFixed, t0, Tn], 20]}],
   {nc, NcycList}];
Print /@ omegaScan;
omegaFit = Fit[Log@omegaScan, {1, x}, x];
Print["log-log slope d(log eps)/d(log omega) = ", Coefficient[omegaFit, x]];

(* ================================================================ *)
(* SCAN 2: epsilon versus t at fixed omega                          *)
(* ================================================================ *)

Print["=== SCAN 2: epsilon vs t at fixed omega ==="];
timeScan = Table[
   {N[n/ptsPerCycle], N[epsilon[n TcFixed/ptsPerCycle, t0, TcFixed], 20]},
   {n, 1, nCyclesScan ptsPerCycle}];
Print /@ timeScan;

(* ================================================================ *)
(* Export plot data                                                 *)
(* ================================================================ *)

outDir = Quiet@Check[NotebookDirectory[], Directory[]] /. $Failed :> Directory[];

exportDat[fname_, data_, colHeader_] := Module[{path, rows},
  path = FileNameJoin[{outDir, fname}];
  rows = Prepend[Map[{CForm[N[#[[1]]]], CForm[N[#[[2]]]]} &, data], {"#", colHeader}];
  Export[path, rows, "Table", "FieldSeparators" -> "\t"];
  Print["  wrote ", path];
  path
];

Print["=== Writing .dat files ==="];
exportDat["epsilon_vs_omega.dat", omegaScan, "omega\tnorm_diff"];
exportDat["epsilon_vs_time.dat", timeScan, "t_over_T\tnorm_diff"];

pulseGrid = Table[n (2 TcFixed)/400, {n, 0, 400}];
pulseTable = Table[
   Prepend[Table[N[TotalgPulse[al, TcFixed, tt, 2 Pi/TcFixed]], {al, 1, maxnumberpulses}], N[tt]],
   {tt, pulseGrid}];
With[{pulsePath = FileNameJoin[{outDir, "pulse_profile.dat"}]},
  Export[pulsePath,
    Prepend[Map[Map[CForm, #] &, pulseTable],
      Prepend[Table["g" <> ToString[al], {al, 1, maxnumberpulses}], "#t"]],
    "Table", "FieldSeparators" -> "\t"];
  Print["  wrote ", pulsePath];
];

(* ================================================================ *)
(* Plots                                                            *)
(* ================================================================ *)

Plot[
  Evaluate@Table[TotalgPulse[al, TcFixed, t, 2 Pi/TcFixed], {al, 1, maxnumberpulses}],
  {t, 0, 2 TcFixed},
  PlotRange -> All, AxesLabel -> {"t", "\[Omega] g_\[Alpha](t)"},
  PlotLegends -> Table[Subscript["g", al], {al, 1, maxnumberpulses}],
  PlotLabel -> "Pulse channels from 0 to 2T"
]

Show[
  ListLogLogPlot[omegaScan, PlotMarkers -> Automatic, Joined -> True,
    Frame -> True, FrameLabel -> {"\[Omega]", "\[Epsilon]"},
    PlotLabel -> "\[Epsilon] vs \[Omega] at fixed t (slope \[Rule] -1)",
    PlotRange -> All],
  LogLogPlot[Exp[omegaFit /. x -> Log[om]], {om, omegaScan[[1, 1]], omegaScan[[-1, 1]]},
    PlotStyle -> {Dashed, Gray}]
]

ListLogPlot[timeScan,
  PlotMarkers -> Automatic, Joined -> True,
  Frame -> True, FrameLabel -> {"t/T", "\[Epsilon]"},
  PlotLabel -> "\[Epsilon] vs t at fixed \[Omega]",
  PlotRange -> All
]





