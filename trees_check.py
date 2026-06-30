#!/usr/bin/env python3
"""
trees_check.py -- certifies the shared-tree-extremes proposition:

    Among all trees on n vertices, both DEU and DSO are minimized uniquely by
    the path P_n and maximized uniquely by the star S_n.

Generates every tree on n vertices with nauty's `gentreeg`, computes DEU and DSO
(deu_core), and confirms the path attains the minimum and the star the maximum,
for BOTH indices, uniquely, over all trees -- well past the n <= 10 enumeration
range, since trees are cheap.

Usage:
    python3 trees_check.py                 # n = 4..18
    python3 trees_check.py --nmax 22

Requires nauty's gentreeg on PATH (binary name `gentreeg` or `nauty-gentreeg`).
Exit status 0 iff path=min and star=max (uniquely) for both indices, all n.
"""
import argparse
import shutil
import subprocess
import sys
from collections import Counter
import networkx as nx
from deu_core import invariants, TOL


def find_gentreeg() -> str:
    for cand in ("gentreeg", "nauty-gentreeg"):
        p = shutil.which(cand)
        if p:
            return p
    sys.exit("ERROR: nauty's gentreeg not found (looked for 'gentreeg' and "
             "'nauty-gentreeg'). Install nauty (provides 'nauty-gentreeg').")


def trees(gentreeg: str, n: int):
    """Yield every tree on n vertices (gentreeg emits sparse6)."""
    proc = subprocess.Popen([gentreeg, "-q", str(n)], stdout=subprocess.PIPE, bufsize=1 << 20)
    for raw in proc.stdout:
        line = raw.strip()
        if line:
            yield nx.from_sparse6_bytes(line)
    proc.stdout.close()
    proc.wait()


def degree_multiset(G):
    return Counter(d for _, d in G.degree())


def is_path(G):
    n = G.number_of_nodes()
    return degree_multiset(G) == Counter({1: 2, 2: n - 2}) if n >= 3 else True


def is_star(G):
    n = G.number_of_nodes()
    return degree_multiset(G) == Counter({1: n - 1, n - 1: 1}) if n >= 2 else True


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nmin", type=int, default=4)   # n<=3: path and star coincide
    ap.add_argument("--nmax", type=int, default=18)
    args = ap.parse_args()

    gentreeg = find_gentreeg()
    print(f"Using gentreeg at: {gentreeg}")
    print(f"Shared-tree-extremes check (DEU & DSO), n = {args.nmin}..{args.nmax}\n")
    print(f"{'n':>3} {'#trees':>8}   DEU min=P, max=S    DSO min=P, max=S   unique?")

    all_ok = True
    for n in range(args.nmin, args.nmax + 1):
        n_trees = 0
        best = {  # index -> (min_val, min_G, count_at_min, max_val, max_G, count_at_max)
            "DEU": [float("inf"), None, 0, float("-inf"), None, 0],
            "DSO": [float("inf"), None, 0, float("-inf"), None, 0],
        }
        for G in trees(gentreeg, n):
            n_trees += 1
            I = invariants(G)
            for key, val in (("DEU", I.DEU), ("DSO", I.DSO)):
                b = best[key]
                if val < b[0] - TOL:
                    b[0], b[1], b[2] = val, G, 1
                elif abs(val - b[0]) <= TOL:
                    b[2] += 1
                if val > b[3] + TOL:
                    b[3], b[4], b[5] = val, G, 1
                elif abs(val - b[3]) <= TOL:
                    b[5] += 1

        row_ok = True
        cells = []
        for key in ("DEU", "DSO"):
            mn_val, mn_G, mn_cnt, mx_val, mx_G, mx_cnt = best[key]
            min_is_path = is_path(mn_G) and mn_cnt == 1
            max_is_star = is_star(mx_G) and mx_cnt == 1
            cells.append((min_is_path, max_is_star))
            row_ok &= min_is_path and max_is_star
        unique = all(c[0] and c[1] for c in cells)
        all_ok &= row_ok
        print(f"{n:>3} {n_trees:>8}      {str(cells[0][0]):>5} {str(cells[0][1]):>5}        "
              f"{str(cells[1][0]):>5} {str(cells[1][1]):>5}      {unique}")

    print(f"\nPath=min & star=max (uniquely) for BOTH DEU and DSO, all n: {all_ok}")
    sys.exit(0 if all_ok else 1)


if __name__ == "__main__":
    main()
