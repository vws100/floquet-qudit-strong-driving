#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
sud_lookup.py -- general-d S.A.D. generator / support-map / lookup-table builder.

Python port of sud_generator_general_d.wl.  All functionality of the Wolfram
script is retained:

  (1) builds (and optionally prints) the su(d) S.A.D. generator matrices,
      INPUT matrices first, then OUTPUT-only matrices;
  (2) builds the support map supp(alpha);
  (3) builds the fully-expanded commutator lookup tables, partitioned into the
      same four classes, and writes the same LaTeX and machine-readable files.

Two capabilities are added that the Wolfram script does not have:

  --verify FILE   parse a Mathematica-exported sud_lookup_rows_d<d>.m and check
                  it row-for-row against the Python result (regression test);
  --selftest      check the geometric law of Eq. (nested commutator) directly,
                  [[T^a,T^b]]_u = 2^{-lambda} phi^{(p)} 2^{-(u-1)nu} T^{kappa_p},
                  for u = 1..UMAX, rather than assuming it from u = 1,2.

DIAGONAL BASIS CONVENTION (matches the d=3 example in the manuscript):
    The input set A_in is
      A_S = { S^{mn} : 1<=m<n<=d }        (d(d-1)/2 generators)
      A_A = { A^{mn} : 1<=m<n<=d }        (d(d-1)/2 generators)
      A_D = { H^{1n} : n=2,...,d }        (d-1 diagonal generators)
    so |A_in| = d^2 - 1.  The diagonal inputs are the pair-difference operators
    anchored at level 1, H^{1n}, which form a (non-orthogonal) Cartan basis.
    For d=3 this reproduces lambda^3 = H^{12} and lambda^8 = H^{13} exactly.

    The OUTPUT label set A_out enlarges A_in to include ALL pair-difference
    diagonals,  A_out = A_S U A_A U { H^{mn} : 1<=m<n<=d }.  Diagonal
    commutator outputs that are not themselves inputs (e.g. lambda^h = H^{23}
    = lambda^8 - lambda^3 for d=3) appear only here, as output-only
    bookkeeping generators.

Matrix conventions:
    S^{mn} = (1/2)(|m><n| + |n><m|)
    A^{mn} = (i/2)(|m><n| - |n><m|)
    H^{mn} = (1/2)(|m><m| - |n><n|)

Lookup convention:
    [[T^a, T^b]]_u = 2^{-lambda} phi^{(p)} 2^{-(u-1)nu} T^{kappa_p(a,b)},
    p = u mod 2  ("o" for odd u, "e" for even u).
    lambda in {0,1}: first commutator magnitude 1 (0) or 1/2 (1).
    nu     in {0,1}: successive-magnitude ratio 1 (0) or 1/2 (1).
    phi in {+1,-1,+i,-i}: unit phase of the u=1 (odd) and u=2 (even) outputs;
                          the magnitude lives entirely in lambda and nu.
    Pairs that commute are not listed: phi = 0 and lambda = nu = 0 by
    convention.

All generator entries are dyadic rationals (+-1/2, +-i/2), so complex128
arithmetic here is exact; the tolerances below mirror the Wolfram script
rather than compensating for round-off.

Requires: numpy.  Tested with numpy 2.x, Python 3.10+.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import sys
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np

__all__ = [
    "Label",
    "Row",
    "pair_list",
    "diag_input_pairs",
    "diag_output_only_pairs",
    "input_labels",
    "output_labels",
    "gen_matrix",
    "supp",
    "label_string",
    "decompose",
    "unit_phase",
    "build_row",
    "build_all_rows",
    "partition_tables",
    "latex_generators",
    "latex_support_map",
    "latex_tables",
    "rows_to_wl",
    "parse_wl_rows",
    "compare_rows",
    "check_nested_commutator_law",
]

# Tolerances, matching the Wolfram script.
TOL_ZERO = 1e-10      # Norm[Mn] < 10^-10  in decompose
TOL_ENTRY = 1e-9      # Abs[x] > 10^-9     entry selection / proportionality
TOL_UNIT = 1e-6       # magnitude and phase classification

Label = Tuple[str, int, int]
# Row = (alpha, beta, lambda, nu, (kappa_o, phi_o) | None, (kappa_e, phi_e) | None)
Row = Tuple[Label, Label, int, int, Optional[Tuple[Label, str]], Optional[Tuple[Label, str]]]

TYPE_ORDER = {"S": 0, "A": 1, "H": 2}


# ----------------------------------------------------------------------
# (1) S.A.D. generator matrices
# ----------------------------------------------------------------------

def _EU(m: int, n: int, d: int) -> np.ndarray:
    """|m><n| with 1-based level indices."""
    E = np.zeros((d, d), dtype=complex)
    E[m - 1, n - 1] = 1.0
    return E


def gen_matrix(label: Label, d: int) -> np.ndarray:
    """Matrix for an S.A.D. label ('S'|'A'|'H', m, n) with 1 <= m < n <= d."""
    t, m, n = label
    if t == "S":
        return 0.5 * (_EU(m, n, d) + _EU(n, m, d))
    if t == "A":
        return 0.5j * (_EU(m, n, d) - _EU(n, m, d))
    if t == "H":
        return 0.5 * (_EU(m, m, d) - _EU(n, n, d))
    raise ValueError(f"unknown generator type {t!r}")


def pair_list(d: int) -> List[Tuple[int, int]]:
    """Off-diagonal pair labels (m, n) with m < n, in Wolfram Table order."""
    return [(m, n) for m in range(1, d + 1) for n in range(m + 1, d + 1)]


def diag_input_pairs(d: int) -> List[Tuple[int, int]]:
    """Diagonal INPUT pairs (1, n), n = 2..d  ->  H^{1n} (Cartan basis)."""
    return [(1, n) for n in range(2, d + 1)]


def diag_output_only_pairs(d: int) -> List[Tuple[int, int]]:
    """Diagonal OUTPUT-ONLY pairs: all H^{mn} that are not inputs."""
    din = set(diag_input_pairs(d))
    return [p for p in pair_list(d) if p not in din]


def input_labels(d: int) -> List[Label]:
    """A_in = A_S U A_A U A_D, of size d^2 - 1, in the Wolfram Join order."""
    pl = pair_list(d)
    return (
        [("S", m, n) for (m, n) in pl]
        + [("A", m, n) for (m, n) in pl]
        + [("H", m, n) for (m, n) in diag_input_pairs(d)]
    )


def output_labels(d: int) -> List[Label]:
    """A_out = A_S U A_A U {all H^{mn}}, in the Wolfram Join order.

    The iteration order matters: decompose() returns the first proportional
    match, exactly as the Wolfram version does.
    """
    pl = pair_list(d)
    return (
        [("S", m, n) for (m, n) in pl]
        + [("A", m, n) for (m, n) in pl]
        + [("H", m, n) for (m, n) in pl]
    )


def label_string(label: Label) -> str:
    """Printable/LaTeX form '(T,m,n)'."""
    t, m, n = label
    return f"({t},{m},{n})"


# ----------------------------------------------------------------------
# (2) Support map
# ----------------------------------------------------------------------

def supp(label: Label) -> Tuple[int, int]:
    """supp(alpha): the basis labels on which T^alpha acts non-trivially."""
    _, m, n = label
    return (m, n)


# ----------------------------------------------------------------------
# (3) Commutator lookup tables
# ----------------------------------------------------------------------

def comm(X: np.ndarray, Y: np.ndarray) -> np.ndarray:
    return X @ Y - Y @ X


def decompose(M: np.ndarray, d: int,
              out_labels: Optional[Sequence[Label]] = None
              ) -> Tuple[Optional[Label], complex]:
    """Write M as c * T^gamma for a single OUTPUT generator gamma.

    The H^{mn} are not mutually orthogonal, so this matches by exact
    proportionality M == c*G rather than by trace projection, and returns the
    first match in A_out order (identical to the Wolfram routine).
    Returns (None, 0) if M is zero or is not proportional to any single
    output generator.
    """
    if out_labels is None:
        out_labels = output_labels(d)
    if np.linalg.norm(M) < TOL_ZERO:
        return None, 0.0
    for lbl in out_labels:
        G = gen_matrix(lbl, d)
        idx = np.argwhere(np.abs(G) > TOL_ENTRY)
        if idx.size == 0:
            continue
        i, j = idx[0]
        c = M[i, j] / G[i, j]
        if abs(c) > TOL_ENTRY and np.linalg.norm(M - c * G) < TOL_ENTRY:
            return lbl, complex(c)
    return None, 0.0


def unit_phase(c: complex) -> str:
    """Unit phase of c as one of '+1', '-1', '+i', '-i' (else a fallback str)."""
    p = c / abs(c)
    if abs(p - 1) < TOL_UNIT:
        return "+1"
    if abs(p + 1) < TOL_UNIT:
        return "-1"
    if abs(p - 1j) < TOL_UNIT:
        return "+i"
    if abs(p + 1j) < TOL_UNIT:
        return "-i"
    return f"{p:.6g}"


def build_row(a: Label, b: Label, d: int,
              out_labels: Optional[Sequence[Label]] = None) -> Optional[Row]:
    """One lookup row, or None if T^a and T^b commute (unlisted by convention)."""
    if out_labels is None:
        out_labels = output_labels(d)
    Ga = gen_matrix(a, d)
    c1 = comm(Ga, gen_matrix(b, d))
    if np.linalg.norm(c1) < TOL_ENTRY:
        return None
    c2 = comm(Ga, c1)
    k1, co1 = decompose(c1, d, out_labels)
    k2, co2 = decompose(c2, d, out_labels)
    m1 = abs(co1)
    m2 = 0.0 if k2 is None else abs(co2)
    lam = 0 if abs(m1 - 1) < TOL_UNIT else 1          # magnitude 1 -> 0, 1/2 -> 1
    nu = 0 if (m2 == 0 or abs(m2 / m1 - 1) < TOL_UNIT) else 1   # ratio 1 -> 0, 1/2 -> 1
    odd_t = None if k1 is None else (k1, unit_phase(co1))
    even_t = None if k2 is None else (k2, unit_phase(co2))
    return (a, b, lam, nu, odd_t, even_t)


def build_all_rows(d: int) -> List[Row]:
    """All non-commuting ordered input pairs, in Wolfram Table order."""
    labs = input_labels(d)
    outs = output_labels(d)
    rows: List[Row] = []
    for a in labs:
        for b in labs:
            r = build_row(a, b, d, outs)
            if r is not None:
                rows.append(r)
    return rows


def _is_diag_input(label: Label) -> bool:
    return label[0] == "H"


def partition_tables(rows: Sequence[Row]) -> Dict[str, List[Row]]:
    """Split rows into the four table classes used in the manuscript."""
    t1 = [r for r in rows if r[2] == 0]
    t2 = [r for r in rows if r[2] == 1 and _is_diag_input(r[0])]
    t3 = [r for r in rows if r[2] == 1 and _is_diag_input(r[1]) and not _is_diag_input(r[0])]
    t4 = [r for r in rows if r[2] == 1 and not _is_diag_input(r[0]) and not _is_diag_input(r[1])]
    return {"same_support": t1, "diag_alpha": t2, "diag_beta": t3, "offdiag": t4}


def _sort_key(row: Row):
    a, b = row[0], row[1]
    return (TYPE_ORDER[a[0]], TYPE_ORDER[b[0]], (a[1], a[2]), (b[1], b[2]))


# ----------------------------------------------------------------------
# LaTeX output
# ----------------------------------------------------------------------

_TEX_NAME = {"S": "S", "A": "A", "H": "H"}


def _tex_entry_halved(z: complex) -> str:
    """Entry of 2*T, which is always in {0, +-1, +-i} for S, A and H."""
    two = 2.0 * z
    re_, im_ = two.real, two.imag
    if abs(re_) < 1e-12 and abs(im_) < 1e-12:
        return "0"
    if abs(im_) < 1e-12:
        return f"{int(round(re_)):+d}".lstrip("+") if re_ > 0 else f"{int(round(re_))}"
    if abs(re_) < 1e-12:
        return "i" if im_ > 0 else "-i"
    raise ValueError(f"unexpected entry {z!r}")


def _tex_pmatrix(label: Label, d: int) -> str:
    """'\\hat{S}_j^{mn} &= \\frac{1}{2}\\begin{pmatrix}...\\end{pmatrix}'."""
    t, m, n = label
    M = gen_matrix(label, d)
    rows = ["&".join(_tex_entry_halved(M[i, j]) for j in range(d)) for i in range(d)]
    body = "\\\\\n".join(rows)
    return (f"\\hat{{{_TEX_NAME[t]}}}_j^{{{m}{n}}}\n&=\n\\frac{{1}}{{2}}\n"
            f"\\begin{{pmatrix}}\n{body}\n\\end{{pmatrix}}")


def latex_generators(d: int) -> str:
    """Copy-paste LaTeX for the S.A.D. generator matrices, grouped by sector.

    Not produced by the Wolfram script, which prints the matrices to the
    notebook only.  Entries are written as (1/2) times a matrix over
    {0, +-1, +-i}, matching the hand-written d=3 list in the appendix.
    """
    ds = str(d)
    pl, din, dout = pair_list(d), diag_input_pairs(d), diag_output_only_pairs(d)
    out = [f"% d={ds} S.A.D. generator matrices, auto-generated.",
           "% Conventions: S^{mn}=(1/2)(|m><n|+|n><m|), "
           "A^{mn}=(i/2)(|m><n|-|n><m|), H^{mn}=(1/2)(|m><m|-|n><n|).",
           "", "% --- Symmetric sector S^{mn} (input) ---", "\\begin{align}"]
    out.append(",\n\\nonumber \\\\\n".join(_tex_pmatrix(("S", m, n), d) for (m, n) in pl) + ".")
    out += [f"\\label{{eq_gen_S_d{ds}}}", "\\end{align}", "",
            "% --- Antisymmetric sector A^{mn} (input) ---", "\\begin{align}"]
    out.append(",\n\\nonumber \\\\\n".join(_tex_pmatrix(("A", m, n), d) for (m, n) in pl) + ".")
    out += [f"\\label{{eq_gen_A_d{ds}}}", "\\end{align}", "",
            "% --- Diagonal input sector H^{1n} (Cartan basis) ---", "\\begin{align}"]
    out.append(",\n\\nonumber \\\\\n".join(_tex_pmatrix(("H", m, n), d) for (m, n) in din) + ".")
    out += [f"\\label{{eq_gen_H_in_d{ds}}}", "\\end{align}", ""]
    if dout:
        out += ["% --- Diagonal output-only sector H^{mn} (bookkeeping) ---",
                "\\begin{align}"]
        out.append(",\n\\nonumber \\\\\n".join(_tex_pmatrix(("H", m, n), d) for (m, n) in dout) + ".")
        out += [f"\\label{{eq_gen_H_out_d{ds}}}", "\\end{align}", ""]
    else:
        out += ["% (no output-only generators for this d)", ""]
    return "\n".join(out)


def latex_support_map(d: int) -> str:
    """Explicit, d-specific support map: one entry per input generator.

    Matches the printed listing of part (2) and the .tex written by
    sud_generator_general_d.wl.  The generic cases-form rule is kept as a
    LaTeX comment for reference.
    """
    entries = [
        f"\\mathrm{{supp}}{label_string(l)} &= \\{{{supp(l)[0]},{supp(l)[1]}\\}}"
        for l in input_labels(d)
    ]
    return (
        f"% d={d} support map, auto-generated.\n"
        "% Generic rule: supp(alpha) = {m,n} for alpha=(tau,m,n), tau in {S,A};\n"
        "%               supp(alpha) = {1,n} for alpha=(H,1,n).\n"
        "\\begin{align}\n"
        + ",\n\\nonumber \\\\\n".join(entries)
        + f".\n\\label{{eq_support_map_d{d}}}\n\\end{{align}}\n"
    )


def _term_string(term: Optional[Tuple[Label, str]]) -> str:
    if term is None:
        return "--"
    lbl, ph = term
    return label_string(lbl) + ",\\," + ph


def _row_tex(row: Row) -> str:
    a, b, lam, nu, odd_t, even_t = row
    return (f"{label_string(a)} & {label_string(b)} & ${lam}$ & ${nu}$ & "
            f"${_term_string(odd_t)}$ & ${_term_string(even_t)}$ \\\\")


def _table_tex(caption: str, label: str, rows: Sequence[Row]) -> str:
    head = [
        "\\begin{table}[h]", "\\centering", "\\caption{" + caption + "}",
        "\\label{" + label + "}", "\\scriptsize", "\\setlength{\\tabcolsep}{3pt}",
        "\\begin{tabular}{|l|l|c|c|l|l|}", "\\hline",
        "$\\alpha$ & $\\beta$ & $\\lambda$ & $\\nu$ & "
        "$(\\kappa_{\\mathrm o},\\phi^{(\\mathrm o)})$ & "
        "$(\\kappa_{\\mathrm e},\\phi^{(\\mathrm e)})$ \\\\",
        "\\hline",
    ]
    body: List[str] = []
    for r in sorted(rows, key=_sort_key):
        body.append(_row_tex(r))
        body.append("\\hline")
    tail = ["\\end{tabular}", "\\end{table}"]
    return "\n".join(head + body + tail)


def latex_tables(d: int, tables: Dict[str, List[Row]]) -> str:
    ds = str(d)
    blocks = [
        f"% d={ds} commutator lookup tables (S.A.D. basis), auto-generated.",
        "% Diagonal inputs H^{1n}; convention "
        "[[T^a,T^b]]_u = 2^{-lambda} phi^{(p)} 2^{-(u-1)nu} T^{kappa_p}.",
        _table_tex(f"$d={ds}$ same-support channels.",
                   f"tab_d{ds}_same_support", tables["same_support"]),
        _table_tex(f"$d={ds}$ one-overlap channels with diagonal input $\\alpha$.",
                   f"tab_d{ds}_diag_alpha", tables["diag_alpha"]),
        _table_tex(f"$d={ds}$ one-overlap channels with diagonal target $\\beta$.",
                   f"tab_d{ds}_diag_beta", tables["diag_beta"]),
        _table_tex(f"$d={ds}$ one-overlap off-diagonal channels.",
                   f"tab_d{ds}_offdiag", tables["offdiag"]),
        "% All pairs not listed commute: phi=0, lambda=nu=0 by convention.",
    ]
    return "\n\n".join(blocks)


# ----------------------------------------------------------------------
# Machine-readable output / round-trip with Mathematica
# ----------------------------------------------------------------------

def _wl_label(lbl: Label) -> str:
    t, m, n = lbl
    return f'{{"{t}", {m}, {n}}}'


def _wl_term(term: Optional[Tuple[Label, str]]) -> str:
    if term is None:
        return "None"
    lbl, ph = term
    return f'{{{_wl_label(lbl)}, "{ph}"}}'


def rows_to_wl(rows: Sequence[Row]) -> str:
    """Serialise rows in the same Wolfram list syntax as Export[..., '.m']."""
    items = []
    for a, b, lam, nu, odd_t, even_t in rows:
        items.append(
            f"{{{_wl_label(a)}, {_wl_label(b)}, {lam}, {nu}, "
            f"{_wl_term(odd_t)}, {_wl_term(even_t)}}}"
        )
    body = ", \n ".join(items)
    return ("(* Created with the Wolfram Language : www.wolfram.com *)\n"
            "{" + body + "}\n")


_TOK = re.compile(r'\{|\}|"[^"]*"|-?\d+|None|,|\s+')


def _tokenize(text: str) -> List[str]:
    out = []
    for mo in _TOK.finditer(text):
        s = mo.group()
        if s.strip() and s != ",":
            out.append(s)
    return out


def _parse_expr(tokens: List[str], i: int):
    """Minimal recursive-descent parser for nested Wolfram lists."""
    tok = tokens[i]
    if tok == "{":
        items = []
        i += 1
        while tokens[i] != "}":
            item, i = _parse_expr(tokens, i)
            items.append(item)
        return items, i + 1
    if tok == "None":
        return None, i + 1
    if tok.startswith('"'):
        return tok[1:-1], i + 1
    return int(tok), i + 1


def parse_wl_rows(path: str) -> List[Row]:
    """Read a Mathematica-exported sud_lookup_rows_d<d>.m back into Row form."""
    with open(path, "r") as fh:
        text = fh.read()
    text = re.sub(r"\(\*.*?\*\)", " ", text, flags=re.S)   # strip comments
    tokens = _tokenize(text)
    parsed, _ = _parse_expr(tokens, 0)
    rows: List[Row] = []
    for r in parsed:
        a = (r[0][0], r[0][1], r[0][2])
        b = (r[1][0], r[1][1], r[1][2])
        lam, nu = r[2], r[3]
        odd_t = None if r[4] is None else ((r[4][0][0], r[4][0][1], r[4][0][2]), r[4][1])
        even_t = None if r[5] is None else ((r[5][0][0], r[5][0][1], r[5][0][2]), r[5][1])
        rows.append((a, b, lam, nu, odd_t, even_t))
    return rows


def compare_rows(got: Sequence[Row], want: Sequence[Row]) -> List[str]:
    """Row-for-row comparison, order-independent.  Returns a list of problems."""
    problems: List[str] = []
    gmap = {(r[0], r[1]): r for r in got}
    wmap = {(r[0], r[1]): r for r in want}
    if len(gmap) != len(got):
        problems.append("duplicate (alpha,beta) keys in computed rows")
    only_g = set(gmap) - set(wmap)
    only_w = set(wmap) - set(gmap)
    for k in sorted(only_g):
        problems.append(f"extra row in Python output: {k}")
    for k in sorted(only_w):
        problems.append(f"missing row in Python output: {k}")
    for k in sorted(set(gmap) & set(wmap)):
        if gmap[k] != wmap[k]:
            problems.append(f"mismatch at {k}:\n    python = {gmap[k]}\n    ref    = {wmap[k]}")
    return problems


# ----------------------------------------------------------------------
# Self-test: verify the geometric law directly for u = 1..UMAX
# ----------------------------------------------------------------------

_PHASE = {"+1": 1.0 + 0j, "-1": -1.0 + 0j, "+i": 1j, "-i": -1j}


def check_nested_commutator_law(d: int, rows: Optional[Sequence[Row]] = None,
                                umax: int = 8) -> Tuple[int, List[str]]:
    """Check  [[T^a,T^b]]_u == 2^{-lam} phi^{(p)} 2^{-(u-1)nu} T^{kappa_p}
    for every listed pair and u = 1..umax.

    The Wolfram script infers lambda, nu, kappa and phi from u = 1 and u = 2
    only; this confirms the law actually holds at higher u, i.e. that kappa
    depends on the parity of u and not on u itself, and that the magnitudes
    are geometric.  Returns (n_checks, problems).
    """
    if rows is None:
        rows = build_all_rows(d)
    problems: List[str] = []
    n = 0
    for a, b, lam, nu, odd_t, even_t in rows:
        Ga = gen_matrix(a, d)
        X = gen_matrix(b, d)
        for u in range(1, umax + 1):
            X = comm(Ga, X)
            term = odd_t if u % 2 == 1 else even_t
            if term is None:
                if np.linalg.norm(X) > TOL_ENTRY:
                    problems.append(f"{a},{b},u={u}: expected zero, got nonzero")
                continue
            kap, ph = term
            pred = (2.0 ** (-lam)) * _PHASE[ph] * (2.0 ** (-(u - 1) * nu)) * gen_matrix(kap, d)
            err = np.linalg.norm(X - pred)
            n += 1
            if err > 1e-9:
                problems.append(f"{a},{b},u={u}: law violated, ||diff|| = {err:.3e}")
    return n, problems


# ----------------------------------------------------------------------
# Printing (part 1 and 2 of the Wolfram script)
# ----------------------------------------------------------------------

def _fmt_entry(z: complex) -> str:
    re_, im_ = z.real, z.imag
    if abs(im_) < 1e-12:
        if abs(re_) < 1e-12:
            return "0"
        if abs(re_ - round(re_)) < 1e-12:
            return f"{int(round(re_))}"
        return f"{re_:g}"
    if abs(re_) < 1e-12:
        if abs(abs(im_) - 1) < 1e-12:
            return "i" if im_ > 0 else "-i"
        return f"{im_:g}i"
    return f"{re_:g}{im_:+g}i"


def _matrix_form(M: np.ndarray, indent: str = "    ") -> str:
    cells = [[_fmt_entry(M[i, j]) for j in range(M.shape[1])] for i in range(M.shape[0])]
    w = max(len(c) for row in cells for c in row)
    return "\n".join(indent + "[ " + "  ".join(c.rjust(w) for c in row) + " ]" for row in cells)


def print_generators(d: int, stream=sys.stdout) -> None:
    p = lambda *a: print(*a, file=stream)
    pl, din, dout = pair_list(d), diag_input_pairs(d), diag_output_only_pairs(d)
    labs = input_labels(d)
    p(f"=== (1) S.A.D. generators for d = {d} ===")
    p(f"off-diagonal pairs (m<n): {pl}")
    p(f"diagonal INPUT pairs H^{{1n}}: {din}")
    p(f"diagonal OUTPUT-ONLY pairs: {dout}")
    p(f"|A_S| = |A_A| = {len(pl)},  |A_D input| = {len(din)},  "
      f"total input = {len(labs)} = d^2-1 = {d * d - 1}")
    herm = all(np.linalg.norm(gen_matrix(l, d) - gen_matrix(l, d).conj().T) < 1e-12 for l in labs)
    trless = all(abs(np.trace(gen_matrix(l, d))) < 1e-12 for l in labs)
    p(f"all inputs Hermitian: {herm},  traceless: {trless},  "
      f"count correct: {len(labs) == d * d - 1}")

    p("\n========== INPUT matrices (A_in) ==========")
    p("\n--- Symmetric sector  S^{mn} ---")
    for (m, n) in pl:
        p(f'{label_string(("S", m, n))} =')
        p(_matrix_form(gen_matrix(("S", m, n), d)))
    p("\n--- Antisymmetric sector  A^{mn} ---")
    for (m, n) in pl:
        p(f'{label_string(("A", m, n))} =')
        p(_matrix_form(gen_matrix(("A", m, n), d)))
    p("\n--- Diagonal input sector  H^{1n} (Cartan basis) ---")
    for (m, n) in din:
        p(f'{label_string(("H", m, n))} =')
        p(_matrix_form(gen_matrix(("H", m, n), d)))

    p("\n========== OUTPUT-ONLY matrices (A_out \\ A_in) ==========")
    if not dout:
        p("(none: every output generator is also an input for this d)")
    else:
        p("\n--- Diagonal output-only sector  H^{mn} ---")
        for (m, n) in dout:
            p(f'{label_string(("H", m, n))} =')
            p(_matrix_form(gen_matrix(("H", m, n), d)))


def print_support_map(d: int, stream=sys.stdout) -> None:
    print("\n=== (2) Support map supp(alpha) ===", file=stream)
    for lbl in input_labels(d):
        print(f"  supp{label_string(lbl)} = {list(supp(lbl))}", file=stream)


# ----------------------------------------------------------------------
# CLI
# ----------------------------------------------------------------------

def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(
        description="Build su(d) S.A.D. generators, support map, and commutator "
                    "lookup tables (Python port of sud_generator_general_d.wl).")
    ap.add_argument("-d", "--dim", type=int, default=4,
                    help="local Hilbert-space dimension d >= 2 (default: 4)")
    ap.add_argument("-o", "--outdir", default=".",
                    help="directory for written files (default: current)")
    ap.add_argument("--print-generators", action=argparse.BooleanOptionalAction,
                    default=True,
                    help="print the generator matrices, part 1 of the script "
                         "(default: on; use --no-print-generators to suppress)")
    ap.add_argument("--print-support", action=argparse.BooleanOptionalAction,
                    default=True,
                    help="print the support map, part 2 of the script "
                         "(default: on; use --no-print-support to suppress)")
    ap.add_argument("-q", "--quiet", action="store_true",
                    help="suppress parts 1 and 2; print only the table summary")
    ap.add_argument("--print-latex", action="store_true",
                    help="print the copy-paste LaTeX to stdout")
    ap.add_argument("--no-write", action="store_true",
                    help="do not write any files")
    ap.add_argument("--json", action="store_true",
                    help="also write sud_lookup_rows_d<d>.json")
    ap.add_argument("--verify", metavar="FILE",
                    help="compare against a Mathematica-exported .m row file")
    ap.add_argument("--selftest", action="store_true",
                    help="verify the geometric nested-commutator law for u=1..8")
    args = ap.parse_args(argv)

    d = args.dim
    if d < 2:
        print("ERROR: d must be an integer >= 2.", file=sys.stderr)
        return 2

    if args.print_generators and not args.quiet:
        print_generators(d)
    if args.print_support and not args.quiet:
        print_support_map(d)

    rows = build_all_rows(d)
    n_in = len(input_labels(d))
    print(f"\n=== (3) Building lookup tables (all ordered input pairs) ===")
    print(f"non-commuting ordered pairs: {len(rows)} (of {n_in * (n_in - 1)} total)")

    tables = partition_tables(rows)
    print(f"Table 1 (same-support): {len(tables['same_support'])}")
    print(f"Table 2 (diagonal input alpha): {len(tables['diag_alpha'])}")
    print(f"Table 3 (diagonal target beta): {len(tables['diag_beta'])}")
    print(f"Table 4 (off-diagonal one-overlap): {len(tables['offdiag'])}")
    total = sum(len(v) for v in tables.values())
    print(f"sum: {total} of {len(rows)}")

    gens_tex = latex_generators(d)
    supp_tex = latex_support_map(d)
    tables_tex = latex_tables(d, tables)

    if args.print_latex:
        print("\n=== copy-paste LaTeX: generator matrices ===")
        print(gens_tex)
        print("=== copy-paste LaTeX: support map ===")
        print(supp_tex)
        print("=== copy-paste LaTeX: lookup tables ===")
        print(tables_tex)

    if not args.no_write:
        os.makedirs(args.outdir, exist_ok=True)
        written = []
        p0 = os.path.join(args.outdir, f"sud_generators_d{d}.tex")
        p1 = os.path.join(args.outdir, f"sud_support_map_d{d}.tex")
        p2 = os.path.join(args.outdir, f"sud_lookup_tables_d{d}.tex")
        p3 = os.path.join(args.outdir, f"sud_lookup_rows_d{d}.m")
        with open(p0, "w") as fh:
            fh.write(gens_tex)
        with open(p1, "w") as fh:
            fh.write(supp_tex)
        with open(p2, "w") as fh:
            fh.write(tables_tex)
        with open(p3, "w") as fh:
            fh.write(rows_to_wl(rows))
        written += [p0, p1, p2, p3]
        if args.json:
            p4 = os.path.join(args.outdir, f"sud_lookup_rows_d{d}.json")
            with open(p4, "w") as fh:
                json.dump([{"alpha": list(r[0]), "beta": list(r[1]),
                            "lambda": r[2], "nu": r[3],
                            "kappa_o": None if r[4] is None else list(r[4][0]),
                            "phi_o": None if r[4] is None else r[4][1],
                            "kappa_e": None if r[5] is None else list(r[5][0]),
                            "phi_e": None if r[5] is None else r[5][1]}
                           for r in rows], fh, indent=1)
            written.append(p4)
        print("\nWrote:")
        for w in written:
            print("  ", w)

    status = 0

    if args.selftest:
        n, problems = check_nested_commutator_law(d, rows)
        print(f"\n=== self-test: nested-commutator law, u = 1..8 ===")
        print(f"checks performed: {n}")
        if problems:
            status = 1
            print(f"FAILED ({len(problems)} problems):")
            for pr in problems[:20]:
                print("  ", pr)
        else:
            print("PASSED: law holds for every listed pair.")

    if args.verify:
        ref = parse_wl_rows(args.verify)
        problems = compare_rows(rows, ref)
        print(f"\n=== verify against {args.verify} ===")
        print(f"reference rows: {len(ref)},  python rows: {len(rows)}")
        if problems:
            status = 1
            print(f"FAILED ({len(problems)} problems):")
            for pr in problems[:20]:
                print("  ", pr)
        else:
            print("PASSED: all rows agree exactly.")

    return status


if __name__ == "__main__":
    raise SystemExit(main())
