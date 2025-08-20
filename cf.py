#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse, itertools, re, sys
import numpy as np
from pathlib import Path
from typing import List, Tuple, Optional, Dict

HEADER_RE = re.compile(
    r"^Hopping\s*<a\|H\|b>\s*between\s+"
    r"(?P<a>\S+)\s*\(.*?\)\s*<-->\s*"
    r"(?P<b>\S+)\s*\(.*?\)\s*in\s*sphere\s*#\s*"
    r"(?P<sphere>\d+).*?radius\s*"
    r"(?P<radius>[+-]?(?:\d+\.\d*|\.\d+|\d+)(?:[eEdD][+-]?\d+)?)",
    re.IGNORECASE,
)
RADIUS_VEC_RE = re.compile(
    r"""^Radius\s+vector\s+is:\s*
        (?P<x>[+-]?(?:\d+\.\d*|\.\d+|\d+)(?:[eEdD][+-]?\d+)?)\s+
        (?P<y>[+-]?(?:\d+\.\d*|\.\d+|\d+)(?:[eEdD][+-]?\d+)?)\s+
        (?P<z>[+-]?(?:\d+\.\d*|\.\d+|\d+)(?:[eEdD][+-]?\d+)?)\s*$""",
    re.VERBOSE | re.IGNORECASE,
)
FLOAT_RE = re.compile(r"[+-]?(?:\d+\.\d*|\.\d+|\d+)(?:[eEdD][+-]?\d+)?")

def parse_floats_from_line(line: str) -> Optional[List[float]]:
    toks = FLOAT_RE.findall(line.strip())
    if not toks: return None
    vals = []
    for t in toks:
        t = t.replace('D','E').replace('d','E')
        try: vals.append(float(t))
        except ValueError: return None
    return vals

def nearly_zero_vec(v, tol: float = 1e-8) -> bool:
    return abs(v[0])<tol and abs(v[1])<tol and abs(v[2])<tol

def atoms_match(pair_str: str, a: str, b: str) -> bool:
    p = pair_str.replace(" ","").lower()
    return p in (f"{a}-{b}".replace(" ","").lower(), f"{b}-{a}".replace(" ","").lower())

def find_onsite_block(lines: List[str], pair: str) -> Optional[np.ndarray]:
    i, nlines = 0, len(lines)
    while i < nlines:
        m = HEADER_RE.match(lines[i])
        if m:
            a, b = m.group("a"), m.group("b")
            sphere = int(m.group("sphere"))
            if atoms_match(pair, a, b) and sphere == 0:
                j = i+1
                while j < nlines and not lines[j].strip(): j += 1
                if j >= nlines: break
                r = RADIUS_VEC_RE.match(lines[j])
                if not r: i += 1; continue
                R = tuple(float(r.group(k).replace('D','E').replace('d','E')) for k in ("x","y","z"))
                if not nearly_zero_vec(R): i += 1; continue

                k = j+1
                while k < nlines and not lines[k].strip(): k += 1
                if k >= nlines: break
                first = parse_floats_from_line(lines[k])
                if not first: i += 1; continue
                N = len(first)
                mat = [first]
                k += 1
                while k < nlines and len(mat) < N:
                    row = parse_floats_from_line(lines[k])
                    if row and len(row) == N: mat.append(row); k += 1
                    else: break
                if len(mat) == N: return np.array(mat, float)
        i += 1
    return None

def analyze_block(H: np.ndarray, basis_labels: List[str], symmetrize: bool=True) -> Dict:
    if H.shape[0] != H.shape[1]: raise ValueError("Block not square")
    if H.shape[0] != len(basis_labels): raise ValueError("Dim mismatch vs basis")
    H_used = (H + H.T)/2.0 if symmetrize else H.copy()
    anti = H - H.T
    off = H_used - np.diag(np.diag(H_used))
    evals, evecs = np.linalg.eigh(H_used)
    evecs = evecs / np.linalg.norm(evecs, axis=0, keepdims=True)
    chars = (evecs**2); chars /= np.sum(chars, axis=0, keepdims=True)
    return {
        "H": H, "H_used": H_used,
        "max_antiherm": float(np.max(np.abs(anti))) if H.size else 0.0,
        "max_offdiag": float(np.max(np.abs(off))) if H.size else 0.0,
        "evals": evals, "evecs": evecs, "chars": chars, "basis": basis_labels,
    }

# ---------- pretty tables ----------
def _hline(widths: List[int]) -> str:
    return "-" * (sum(widths) + 3*(len(widths)-1))
def _fmt_row(cells: List[str], widths: List[int]) -> str:
    return " | ".join(c.ljust(w) for c, w in zip(cells, widths))
def print_table(headers: List[str], rows: List[List[str]]):
    widths = [max(len(h), *(len(r[i]) for r in rows)) for i, h in enumerate(headers)]
    line = _hline(widths)
    print(line); print(_fmt_row(headers, widths)); print(line)
    for r in rows: print(_fmt_row(r, widths))
    print(line)
def dominant_desc(basis: List[str], perc: np.ndarray, thr: float=5.0) -> str:
    idx = np.where(perc >= thr)[0]
    if idx.size==0: idx = [int(np.argmax(perc))]
    return ", ".join(f"{basis[j]}~{perc[j]:.1f}%" for j in idx)

# ---------- pairing by maximum overlap ----------
def pair_by_overlap(evecs_up: np.ndarray, evecs_dn: np.ndarray) -> Tuple[List[int], np.ndarray]:
    S2 = np.abs(evecs_up.T @ evecs_dn)**2
    N = S2.shape[0]
    best_perm, best_score = None, -1.0
    for perm in itertools.permutations(range(N)):
        score = sum(S2[i, perm[i]] for i in range(N))
        if score > best_score: best_perm, best_score = list(perm), score
    return best_perm, S2

# ---------- reporting ----------
def per_spin_table(spin_label: str, result: Dict, pct_threshold: float=5.0, decimals: int=6):
    basis, evals, chars = result["basis"], result["evals"], 100.0*result["chars"]
    headers = ["Level", "Energy (eV)"] + [f"{b} (%)" for b in basis] + [f"Dominant (≥{int(pct_threshold)}%)"]
    rows=[]
    for i, E in enumerate(evals):
        perc = chars[:, i]
        rows.append(
            [f"{i+1:d}", f"{E:.{decimals}f}"] +
            [f"{p:.1f}" if p>=pct_threshold else "0.0" for p in perc] +
            [dominant_desc(basis, perc, thr=pct_threshold)]
        )
    print(f"\n=== {spin_label} spin ===")
    print(f"Max |off-diagonal| = {result['max_offdiag']:.3e} ; max |H - H^T| = {result['max_antiherm']:.3e}")
    print_table(headers, rows)

def combined_table_avg(res_up: Dict, res_dn: Dict, pct_threshold: float=5.0, decimals: int=6):
    """
    Combined table matching ↑/↓ by max overlap, then averaging:
      Energy = (E↑ + E↓)/2
      Characters = 0.5*(chars↑ + chars↓)
    Columns mirror per-spin tables.
    """
    basis = res_up["basis"]
    perm_dn, _S2 = pair_by_overlap(res_up["evecs"], res_dn["evecs"])
    headers = ["Level", "Energy (eV)"] + [f"{b} (%)" for b in basis] + [f"Dominant (≥{int(pct_threshold)}%)"]
    rows = []
    for i_up, j_dn in enumerate(perm_dn):
        E_avg = 0.5*(res_up["evals"][i_up] + res_dn["evals"][j_dn])
        chars_avg = 50.0*(res_up["chars"][:, i_up] + res_dn["chars"][:, j_dn])  # already ×100/2
        rows.append(
            [f"{i_up+1:d}", f"{E_avg:.{decimals}f}"] +
            [f"{p:.1f}" if p>=pct_threshold else "0.0" for p in chars_avg] +
            [dominant_desc(basis, chars_avg, thr=pct_threshold)]
        )
    print("\n=== Combined (spin-averaged over overlap-paired states) ===")
    print_table(headers, rows)
    print("Note: Energy = (E↑+E↓)/2 ; characters = average of ↑/↓ percentages for the paired states.")

def main():
    p = argparse.ArgumentParser(description="On-site d-levels & orbital characters from Wannier90 hopping (spinful, no SOC).")
    p.add_argument("--up", default="out1.dat", help="Spin-up file (default: out1.dat)")
    p.add_argument("--down", default="out2.dat", help="Spin-down file (default: out2.dat)")
    p.add_argument("--pair", required=True, help='Atom pair like "Cr1-Cr1"')
    p.add_argument("--basis", default="dz2,dxz,dyz,dx2,dxy", help="5 labels for d-block order")
    p.add_argument("--no-sym", action="store_true", help="Disable (H+H^T)/2 symmetrization")
    p.add_argument("--decimals", type=int, default=6)
    p.add_argument("--dom-thr", type=float, default=10.0)
    args = p.parse_args()

    basis = [b.strip() for b in re.split(r"[,\s]+", args.basis.strip()) if b.strip()]
    if len(basis) != 5: print("ERROR: need 5 basis labels", file=sys.stderr); sys.exit(1)

    up_path, dn_path = Path(args.up), Path(args.down)
    if not up_path.exists(): print(f"ERROR: {up_path} not found", file=sys.stderr); sys.exit(1)
    if not dn_path.exists(): print(f"ERROR: {dn_path} not found", file=sys.stderr); sys.exit(1)

    up_lines = up_path.read_text(encoding="utf-8", errors="ignore").splitlines()
    dn_lines = dn_path.read_text(encoding="utf-8", errors="ignore").splitlines()

    H_up = find_onsite_block(up_lines, args.pair)
    H_dn = find_onsite_block(dn_lines, args.pair)
    if H_up is None: print(f"ERROR: on-site block '{args.pair}' not found in {up_path}", file=sys.stderr); sys.exit(1)
    if H_dn is None: print(f"ERROR: on-site block '{args.pair}' not found in {dn_path}", file=sys.stderr); sys.exit(1)

    res_up = analyze_block(H_up, basis, symmetrize=(not args.no_sym))
    res_dn = analyze_block(H_dn, basis, symmetrize=(not args.no_sym))

    per_spin_table("↑ (up)", res_up, pct_threshold=args.dom_thr, decimals=args.decimals)
    per_spin_table("↓ (down)", res_dn, pct_threshold=args.dom_thr, decimals=args.decimals)
    combined_table_avg(res_up, res_dn, pct_threshold=args.dom_thr, decimals=args.decimals)

if __name__ == "__main__":
    main()
