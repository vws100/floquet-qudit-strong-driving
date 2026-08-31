(* ::Package:: *)

(* ================================================================ *)
(*  General-d S.A.D. generator / support-map / lookup-table builder *)
(*                                                                  *)
(*  Set the single input d below.  The script then:                 *)
(*   (1) builds and prints the su(d) S.A.D. generator matrices,     *)
(*       INPUT matrices first, then OUTPUT-only matrices, in the    *)
(*       same format used for the d=3 lambda^a list;                *)
(*   (2) builds the support map supp(alpha);                        *)
(*   (3) builds the fully-expanded commutator lookup tables in the  *)
(*       exact convention of the d=3 example.                       *)
(*                                                                  *)
(*  DIAGONAL BASIS CONVENTION (matches the d=3 example):            *)
(*    The input set A_in is                                         *)
(*      A_S = { S^{mn} : 1<=m<n<=d }        (d(d-1)/2 generators)   *)
(*      A_A = { A^{mn} : 1<=m<n<=d }        (d(d-1)/2 generators)   *)
(*      A_D = { H^{1n} : n=2,...,d }        (d-1 diagonal gens)     *)
(*    so |A_in| = d^2 - 1.  The diagonal inputs are the pair-       *)
(*    difference operators anchored at level 1, H^{1n}, which form  *)
(*    a (non-orthogonal) Cartan basis.  For d=3 this reproduces     *)
(*    lambda^3 = H^{12} and lambda^8 = H^{13} exactly.              *)
(*                                                                  *)
(*    The OUTPUT label set A_out enlarges A_in to include ALL       *)
(*    pair-difference diagonals:                                    *)
(*      A_out = A_S U A_A U { H^{mn} : 1<=m<n<=d }.                 *)
(*    Diagonal commutator outputs that are not themselves inputs    *)
(*    (e.g. lambda^h = H^{23} = lambda^8 - lambda^3 for d=3) appear *)
(*    only here, as output-only bookkeeping generators.            *)
(*                                                                  *)
(*  Matrix conventions:                                             *)
(*    S^{mn} = (1/2)(|m><n| + |n><m|)                               *)
(*    A^{mn} = (i/2)(|m><n| - |n><m|)                               *)
(*    H^{mn} = (1/2)(|m><m| - |n><n|)                               *)
(*                                                                  *)
(*  Lookup convention (as in the d=3 tables):                       *)
(*    [[T^a, T^b]]_u = 2^{-lambda} phi^{(p)} 2^{-(u-1)nu}           *)
(*                     T^{kappa_p(a,b)},   p = u mod 2.             *)
(*    lambda in {0,1}: first commutator magnitude 1 (0) or 1/2 (1). *)
(*    nu     in {0,1}: successive-magnitude ratio 1 (0) or 1/2 (1). *)
(*    phi in {+1,-1,+i,-i}: unit phase of the u=1 (odd), u=2 (even) *)
(*                     outputs; magnitude lives entirely in lam,nu. *)
(* ================================================================ *)

ClearAll["Global`*"];

(* ----------------------- USER INPUT ----------------------------- *)
d = 5;    (* set the local Hilbert-space dimension here (d >= 2) *)
(* ---------------------------------------------------------------- *)

If[! (IntegerQ[d] && d >= 2),
  Print["ERROR: d must be an integer >= 2."]; Abort[]];

(* ================================================================ *)
(* (1)  S.A.D. generator matrices  (built AND printed)             *)
(*      inputs first, then output-only diagonals                   *)
(* ================================================================ *)

ClearAll[EU, Smat, Amat, Hmat, genMat, supp, labelString];
EU[m_, n_] := SparseArray[{{m, n} -> 1}, {d, d}];

Smat[m_, n_] := (1/2) (EU[m, n] + EU[n, m]);
Amat[m_, n_] := (I/2) (EU[m, n] - EU[n, m]);
Hmat[m_, n_] := (1/2) (EU[m, m] - EU[n, n]);

(* off-diagonal pair labels (m<n) *)
pairList = Flatten[Table[{m, n}, {m, 1, d}, {n, m + 1, d}], 1];

(* diagonal INPUT pairs: (1,n) for n=2..d  -> H^{1n} (Cartan basis) *)
diagInPairs = Table[{1, n}, {n, 2, d}];

(* diagonal OUTPUT-ONLY pairs: all H^{mn} that are NOT inputs *)
diagOutOnlyPairs = Complement[pairList, diagInPairs];

(* INPUT generator labels A_in = A_S U A_A U A_D  (d^2 - 1 total) *)
inputLabels = Join[
   {"S", #[[1]], #[[2]]} & /@ pairList,
   {"A", #[[1]], #[[2]]} & /@ pairList,
   {"H", #[[1]], #[[2]]} & /@ diagInPairs
];

(* full OUTPUT basis A_out = A_S U A_A U {all H^{mn}} *)
outputLabels = Join[
   {"S", #[[1]], #[[2]]} & /@ pairList,
   {"A", #[[1]], #[[2]]} & /@ pairList,
   {"H", #[[1]], #[[2]]} & /@ pairList
];

(* label -> matrix *)
genMat[{"S", m_, n_}] := Smat[m, n];
genMat[{"A", m_, n_}] := Amat[m, n];
genMat[{"H", m_, n_}] := Hmat[m, n];

(* label -> printable string *)
labelString[{t_String, m_, n_}] := "(" <> t <> "," <> ToString[m] <> "," <> ToString[n] <> ")";

Print["=== (1) S.A.D. generators for d = ", d, " ==="];
Print["off-diagonal pairs (m<n): ", pairList];
Print["diagonal INPUT pairs H^{1n}: ", diagInPairs];
Print["diagonal OUTPUT-ONLY pairs: ", diagOutOnlyPairs];
Print["|A_S| = |A_A| = ", Length[pairList],
      ",  |A_D input| = ", Length[diagInPairs],
      ",  total input = ", Length[inputLabels], " = d^2-1 = ", d^2 - 1];

hermOK  = AllTrue[inputLabels, Norm[Normal[genMat[#] - ConjugateTranspose[genMat[#]]]] < 10^-12 &];
traceOK = AllTrue[inputLabels, Abs[Tr[genMat[#]]] < 10^-12 &];
Print["all inputs Hermitian: ", hermOK, ",  traceless: ", traceOK,
      ",  count correct: ", Length[inputLabels] == d^2 - 1];

(* ---- PRINT INPUT matrices, grouped by sector ---- *)
Print["\n========== INPUT matrices (A_in) =========="];

Print["\n--- Symmetric sector  S^{mn} ---"];
Do[Print[labelString[{"S", p[[1]], p[[2]]}], " = ",
         MatrixForm[Normal[Smat[p[[1]], p[[2]]]]]], {p, pairList}];

Print["\n--- Antisymmetric sector  A^{mn} ---"];
Do[Print[labelString[{"A", p[[1]], p[[2]]}], " = ",
         MatrixForm[Normal[Amat[p[[1]], p[[2]]]]]], {p, pairList}];

Print["\n--- Diagonal input sector  H^{1n} (Cartan basis) ---"];
Do[Print[labelString[{"H", p[[1]], p[[2]]}], " = ",
         MatrixForm[Normal[Hmat[p[[1]], p[[2]]]]]], {p, diagInPairs}];

(* ---- PRINT OUTPUT-ONLY matrices ---- *)
Print["\n========== OUTPUT-ONLY matrices (A_out \\ A_in) =========="];
If[Length[diagOutOnlyPairs] == 0,
  Print["(none: every output generator is also an input for this d)"],
  Print["\n--- Diagonal output-only sector  H^{mn} ---"];
  Do[Print[labelString[{"H", p[[1]], p[[2]]}], " = ",
           MatrixForm[Normal[Hmat[p[[1]], p[[2]]]]]], {p, diagOutOnlyPairs}]
];

(* ================================================================ *)
(* (2)  Support map  supp(alpha)                                    *)
(* ================================================================ *)

supp[{"S", m_, n_}] := {m, n};
supp[{"A", m_, n_}] := {m, n};
supp[{"H", m_, n_}] := {m, n};

Print["\n=== (2) Support map supp(alpha) ==="];
Do[Print["  supp", labelString[lbl], " = ", supp[lbl]], {lbl, inputLabels}];

(* Explicit, d-specific support map: one entry per input generator, in the *)
(* same order as the printed listing above.  The generic cases-form rule   *)
(* is retained as a LaTeX comment for reference.                          *)
ClearAll[suppEntryTeX];
suppEntryTeX[lbl_] := StringJoin[
   "\\mathrm{supp}", labelString[lbl], " &= \\{",
   ToString[supp[lbl][[1]]], ",", ToString[supp[lbl][[2]]], "\\}"
];

suppMapTeX = StringJoin[
  "% d=" <> ToString[d] <> " support map, auto-generated.\n",
  "% Generic rule: supp(alpha) = {m,n} for alpha=(tau,m,n), tau in {S,A};\n",
  "%               supp(alpha) = {1,n} for alpha=(H,1,n).\n",
  "\\begin{align}\n",
  StringJoin[Riffle[suppEntryTeX /@ inputLabels, ",\n\\nonumber \\\\\n"]],
  ".\n\\label{eq_support_map_d", ToString[d], "}\n\\end{align}\n"
];

(* ================================================================ *)
(* (3)  Commutator lookup tables (correct lambda/nu/phi convention) *)
(* ================================================================ *)

comm[X_, Y_] := X . Y - Y . X;

(* decompose M as a single OUTPUT generator (times a clean scalar).  *)
(* The H^{mn} are not mutually orthogonal, so match by exact         *)
(* proportionality M == c*G rather than trace projection.            *)
ClearAll[decompose];
decompose[M_] := Module[{res = {None, 0}, Mn, G, idx, i, j, c},
  Mn = Normal[M];
  If[Norm[Mn] < 10^-10, Return[{None, 0}]];
  Do[
    G = Normal[genMat[lbl]];
    idx = Position[G, x_ /; Abs[x] > 10^-9, {2}, 1];
    If[idx =!= {},
      {i, j} = idx[[1]];
      c = Mn[[i, j]]/G[[i, j]];
      If[Abs[c] > 10^-9 && Norm[Mn - c G] < 10^-9 && res[[1]] === None,
        res = {lbl, c}]
    ],
    {lbl, outputLabels}];
  res
];

(* unit phase in {+1,-1,+i,-i} as a display string *)
ClearAll[unitPhase];
unitPhase[c_] := Module[{p = c/Abs[c]},
  Which[
    Abs[p - 1] < 10^-6,  "+1",
    Abs[p + 1] < 10^-6,  "-1",
    Abs[p - I] < 10^-6,  "+i",
    Abs[p + I] < 10^-6,  "-i",
    True, ToString[Chop[p]]]
];

(* build one row: {alpha, beta, lambda, nu, {ko,pho}, {ke,phe}} *)
ClearAll[buildRow];
buildRow[a_, b_] := Module[
  {Ga, c1, c2, k1, co1, k2, co2, m1, m2, lam, nu, oddT, evenT},
  Ga = genMat[a];
  c1 = comm[Ga, genMat[b]];
  If[Norm[Normal[c1]] < 10^-9, Return[Null]];   (* commute: not listed *)
  c2 = comm[Ga, c1];
  {k1, co1} = decompose[c1];
  {k2, co2} = decompose[c2];
  m1 = Abs[co1];
  m2 = If[k2 === None, 0, Abs[co2]];
  lam = If[Abs[m1 - 1] < 10^-6, 0, 1];               (* 1 -> 0, 1/2 -> 1 *)
  nu  = If[m2 == 0 || Abs[m2/m1 - 1] < 10^-6, 0, 1]; (* ratio 1 -> 0, 1/2 -> 1 *)
  oddT  = {k1, unitPhase[co1]};                      (* u=1 output (p=o) *)
  evenT = If[k2 === None, None, {k2, unitPhase[co2]}]; (* u=2 output (p=e) *)
  {a, b, lam, nu, oddT, evenT}
];

Print["\n=== (3) Building lookup tables (all ordered input pairs) ==="];
allRows = DeleteCases[
   Flatten[Table[buildRow[a, b], {a, inputLabels}, {b, inputLabels}], 1], Null];
Print["non-commuting ordered pairs: ", Length[allRows],
      " (of ", Length[inputLabels] (Length[inputLabels] - 1), " total)"];

(* partition into four table classes.  A "diagonal" input/target here *)
(* is an H^{1n} (the diagonal inputs of A_D).                         *)
isDiagInput[lbl_] := (lbl[[1]] == "H");
lambdaOf[row_] := row[[3]];
table1 = Select[allRows, lambdaOf[#] == 0 &];                       (* same-support *)
table2 = Select[allRows, lambdaOf[#] == 1 && isDiagInput[#[[1]]] &];
table3 = Select[allRows, lambdaOf[#] == 1 && isDiagInput[#[[2]]] && ! isDiagInput[#[[1]]] &];
table4 = Select[allRows, lambdaOf[#] == 1 && ! isDiagInput[#[[1]]] && ! isDiagInput[#[[2]]] &];
Print["Table 1 (same-support): ", Length[table1]];
Print["Table 2 (diagonal input alpha): ", Length[table2]];
Print["Table 3 (diagonal target beta): ", Length[table3]];
Print["Table 4 (off-diagonal one-overlap): ", Length[table4]];
Print["sum: ", Length[table1] + Length[table2] + Length[table3] + Length[table4],
      " of ", Length[allRows]];

(* ordering for readable output *)
typeOrder = <|"S" -> 0, "A" -> 1, "H" -> 2|>;
sortKey[row_] := {typeOrder[row[[1, 1]]], typeOrder[row[[2, 1]]],
                  Rest[row[[1]]], Rest[row[[2]]]};

(* LaTeX builders: six columns alpha|beta|lambda|nu|(ko,pho)|(ke,phe) *)
ClearAll[termString, rowTeX, tableTeX];
termString[None] := "--";
termString[{lbl_, ph_}] := labelString[lbl] <> ",\\," <> ph;

rowTeX[{a_, b_, lam_, nu_, oddT_, evenT_}] := StringJoin[
   labelString[a], " & ", labelString[b], " & $", ToString[lam], "$ & $",
   ToString[nu], "$ & $", termString[oddT], "$ & $", termString[evenT], "$ \\\\"
];

tableTeX[caption_, lbl_, rows_] := StringJoin[Riffle[Join[
   {"\\begin{table}[h]", "\\centering", "\\caption{" <> caption <> "}",
    "\\label{" <> lbl <> "}", "\\scriptsize", "\\setlength{\\tabcolsep}{3pt}",
    "\\begin{tabular}{|l|l|c|c|l|l|}", "\\hline",
    "$\\alpha$ & $\\beta$ & $\\lambda$ & $\\nu$ & $(\\kappa_{\\mathrm o},\\phi^{(\\mathrm o)})$ & $(\\kappa_{\\mathrm e},\\phi^{(\\mathrm e)})$ \\\\",
    "\\hline"},
   Flatten[{rowTeX[#], "\\hline"} & /@ SortBy[rows, sortKey]],
   {"\\end{tabular}", "\\end{table}"}
], "\n"]];

ds = ToString[d];
latexAll = StringJoin[Riffle[{
   "% d=" <> ds <> " commutator lookup tables (S.A.D. basis), auto-generated.",
   "% Diagonal inputs H^{1n}; convention [[T^a,T^b]]_u = 2^{-lambda} phi^{(p)} 2^{-(u-1)nu} T^{kappa_p}.",
   tableTeX["$d=" <> ds <> "$ same-support channels.",
            "tab_d" <> ds <> "_same_support", table1],
   tableTeX["$d=" <> ds <> "$ one-overlap channels with diagonal input $\\alpha$.",
            "tab_d" <> ds <> "_diag_alpha", table2],
   tableTeX["$d=" <> ds <> "$ one-overlap channels with diagonal target $\\beta$.",
            "tab_d" <> ds <> "_diag_beta", table3],
   tableTeX["$d=" <> ds <> "$ one-overlap off-diagonal channels.",
            "tab_d" <> ds <> "_offdiag", table4],
   "% All pairs not listed commute: phi=0, lambda=nu=0 by convention."
}, "\n\n"]];

(* ================================================================ *)
(* Output: print LaTeX and write files                             *)
(* ================================================================ *)

Print["\n=== copy-paste LaTeX: support map ==="];
Print[suppMapTeX];
Print["=== copy-paste LaTeX: lookup tables ==="];
Print[latexAll];

outDir = Quiet@Check[NotebookDirectory[], Directory[]] /. $Failed :> Directory[];
Export[FileNameJoin[{outDir, "sud_support_map_d" <> ds <> ".tex"}], suppMapTeX, "Text"];
Export[FileNameJoin[{outDir, "sud_lookup_tables_d" <> ds <> ".tex"}], latexAll, "Text"];

machineRows = {#[[1]], #[[2]], #[[3]], #[[4]], #[[5]], #[[6]]} & /@ allRows;
Export[FileNameJoin[{outDir, "sud_lookup_rows_d" <> ds <> ".m"}], machineRows];

Print["\nWrote:"];
Print["  ", FileNameJoin[{outDir, "sud_support_map_d" <> ds <> ".tex"}]];
Print["  ", FileNameJoin[{outDir, "sud_lookup_tables_d" <> ds <> ".tex"}]];
Print["  ", FileNameJoin[{outDir, "sud_lookup_rows_d" <> ds <> ".m"}]];





