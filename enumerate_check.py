#!/usr/bin/env python3
"""
enumerate_check.py -- LAYER 1 of the numeric certification: EXHAUSTIVE.

Streams every connected graph on n vertices (1 <= n <= NMAX) from nauty's
`geng`, builds it in NetworkX, and runs every bound / identity / equality-class
biconditional in deu_core.check_graph().  Graphs with m < 1 (the single edgeless
K_1) are skipped, matching the paper's m >= 1 convention.

The full run (NMAX=10) covers all 11,989,763 connected graphs and is what the
paper's "exhaustive certification over n <= 10" refers to.  It is CPU-heavy
(n=10 alone has 11,716,571 graphs); use --nmax to test on a smaller range.

Usage:
    python3 enumerate_check.py                 # n = 1..10  (full, slow)
    python3 enumerate_check.py --nmax 8        # quick self-test
    python3 enumerate_check.py --nmin 7 --nmax 7

Requires nauty's geng on PATH (binary name `geng` or `nauty-geng`).
Exit status 0 iff zero violations.
"""
import argparse
import shutil
import subprocess
import sys
import time
import networkx as nx
from deu_core import check_graph

EXPECTED_TOTAL_1_TO_10 = 11_989_763   # sum A001349(1..10) minus the edgeless K_1


def find_geng() -> str:
    for cand in ("geng", "nauty-geng"):
        p = shutil.which(cand)
        if p:
            return p
    sys.exit("ERROR: nauty's geng not found (looked for 'geng' and 'nauty-geng').\n"
             "Install nauty, e.g.  apt-get install nauty  or build from "
             "https://pallini.di.uniroma1.it/ .")


def enumerate_n(geng: str, n: int):
    """Yield connected graphs on n vertices as NetworkX graphs (graph6 via geng)."""
    # -c connected, -q quiet
    proc = subprocess.Popen([geng, "-cq", str(n)],
                            stdout=subprocess.PIPE, bufsize=1 << 20)
    for raw in proc.stdout:
        line = raw.strip()
        if line:
            yield nx.from_graph6_bytes(line)
    proc.stdout.close()
    proc.wait()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nmin", type=int, default=1)
    ap.add_argument("--nmax", type=int, default=10)
    args = ap.parse_args()

    geng = find_geng()
    print(f"Using geng at: {geng}")
    print(f"Exhaustive certification over connected graphs, n = {args.nmin}..{args.nmax}\n")

    grand_total = 0           # graphs with m >= 1 actually checked
    grand_checks = 0
    grand_viol = 0
    worst = {}                # name -> smallest slack seen (margin diagnostics)
    t0 = time.time()

    for n in range(args.nmin, args.nmax + 1):
        n_graphs = n_checked = n_viol = 0
        for G in enumerate_n(geng, n):
            n_graphs += 1
            if G.number_of_edges() < 1:
                continue                       # skip edgeless K_1
            n_checked += 1
            for c in check_graph(G):
                grand_checks += 1
                if c.slack < worst.get(c.name, float("inf")):
                    worst[c.name] = c.slack
                if not c.ok:
                    n_viol += 1
                    grand_viol += 1
                    if grand_viol <= 20:       # print first few only
                        print(f"  VIOLATION n={n} {c.name}: slack={c.slack:.3e} "
                              f"{c.detail}  graph6={nx.to_graph6_bytes(G, header=False).strip().decode()}")
        grand_total += n_checked
        print(f"  n={n:2d}: {n_graphs:>9d} connected  "
              f"({n_checked:>9d} with m>=1)  violations={n_viol}")

    dt = time.time() - t0
    print(f"\nTotal graphs checked (m>=1): {grand_total:,}")
    print(f"Total individual checks:     {grand_checks:,}")
    print(f"Total violations:            {grand_viol}")
    print(f"Elapsed:                     {dt:.1f}s")

    if (args.nmin, args.nmax) == (1, 10):
        ok = grand_total == EXPECTED_TOTAL_1_TO_10
        print(f"Count vs paper (11,989,763): {'MATCH' if ok else 'MISMATCH'} ({grand_total:,})")

    # tightest margins: which bounds are actually approached/attained
    print("\nTightest margins observed (min slack per check; ~0 => bound is tight somewhere):")
    for name in sorted(worst, key=lambda k: worst[k]):
        print(f"  {worst[name]:+.3e}  {name}")

    sys.exit(0 if grand_viol == 0 else 1)


if __name__ == "__main__":
    main()
