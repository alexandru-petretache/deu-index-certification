#!/usr/bin/env python3
"""
random_check.py -- LAYER 2 of the numeric certification: RANDOM SAMPLING.

Exhaustive enumeration stops at n = 10 (enumerate_check.py); beyond that the
number of connected graphs explodes.  This script extends the empirical
certification to LARGER n by sampling random connected graphs and running the
same deu_core.check_graph() battery on each.

To exercise the bounds and equality cases across genuinely different degree
regimes -- not just one model's typical graph -- it draws from a mix:

  * random trees                (varied degree sequences, always connected)
  * random d-regular graphs     (hit the regular equality cases at large n)
  * Barabasi-Albert             (hub-heavy: large Delta, the near-star regime
                                 where the separation lower bound is approached)
  * Erdos-Renyi G(n,p)          (dense, resampled until connected)
  * connected Watts-Strogatz    (intermediate structure)
  * random complete bipartite   (the (delta,Delta)-biregular equality class)

Everything is seeded, so a run is reproducible.

Usage:
    python3 random_check.py                       # n up to 60, default depth
    python3 random_check.py --nmax 200 --per 40   # larger / deeper
    python3 random_check.py --seed 7

Exit status 0 iff zero violations.
"""
import argparse
import random
import sys
import networkx as nx
from deu_core import check_graph


def connected_er(n, p, rng):
    """G(n,p) resampled until connected (give up after a few tries)."""
    for _ in range(20):
        G = nx.gnp_random_graph(n, p, seed=rng.randint(0, 2**31 - 1))
        if G.number_of_edges() >= 1 and nx.is_connected(G):
            return G
    return None


def sample_graphs(n, per, rng):
    """Yield a variety of connected graphs on (about) n vertices."""
    seed = lambda: rng.randint(0, 2**31 - 1)

    for _ in range(per):                                   # trees
        yield nx.random_labeled_tree(n, seed=seed())

    for _ in range(per):                                   # d-regular
        d = rng.choice([3, 4, 5, 6])
        if (n * d) % 2 == 0 and d < n:
            try:
                yield nx.random_regular_graph(d, n, seed=seed())
            except nx.NetworkXError:
                pass

    for _ in range(per):                                   # Barabasi-Albert (hubs)
        mm = rng.choice([1, 2, 3])
        if mm < n:
            yield nx.barabasi_albert_graph(n, mm, seed=seed())

    for _ in range(per):                                   # Erdos-Renyi (dense)
        p = rng.uniform(1.5 * (1.0 / n), 0.5)              # above connectivity threshold-ish
        G = connected_er(n, p, rng)
        if G is not None:
            yield G

    for _ in range(per):                                   # Watts-Strogatz
        k = rng.choice([2, 4, 6])
        if k < n:
            yield nx.connected_watts_strogatz_graph(n, k, 0.3, seed=seed())

    for _ in range(per):                                   # complete bipartite (biregular)
        x = rng.randint(1, n - 1)
        yield nx.complete_bipartite_graph(x, n - x)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nmin", type=int, default=11)
    ap.add_argument("--nmax", type=int, default=60)
    ap.add_argument("--step", type=int, default=1)
    ap.add_argument("--per", type=int, default=15,
                    help="samples per model per n")
    ap.add_argument("--seed", type=int, default=20260629)
    args = ap.parse_args()

    rng = random.Random(args.seed)
    print(f"Random-sampling certification, n = {args.nmin}..{args.nmax} step {args.step}, "
          f"~{args.per} per model per n, seed={args.seed}\n")

    total = checks = viol = 0
    max_Delta = 0
    near_star_min_sep = float("inf")   # track how close the strict lower sep bound gets

    for n in range(args.nmin, args.nmax + 1, args.step):
        n_total = n_viol = 0
        for G in sample_graphs(n, args.per, rng):
            if G is None or G.number_of_edges() < 1 or not nx.is_connected(G):
                continue
            n_total += 1
            total += 1
            max_Delta = max(max_Delta, max(dict(G.degree()).values()))
            for c in check_graph(G):
                checks += 1
                if c.name == "Prop:sep (lower, strict)":
                    near_star_min_sep = min(near_star_min_sep, c.slack)
                if not c.ok:
                    viol += 1
                    n_viol += 1
                    if viol <= 20:
                        print(f"  VIOLATION n={n} {c.name}: slack={c.slack:.3e} {c.detail}")
        if n % 10 == 0 or n == args.nmax or n == args.nmin:
            print(f"  n={n:3d}: sampled {n_total:4d} connected graphs, violations={n_viol}")

    print(f"\nTotal random graphs sampled: {total:,}")
    print(f"Total individual checks:     {checks:,}")
    print(f"Total violations:            {viol}")
    print(f"Largest Delta encountered:   {max_Delta}")
    print(f"Closest approach to strict lower separation bound: slack={near_star_min_sep:.3e} "
          f"(stays > 0, as it must)")
    sys.exit(0 if viol == 0 else 1)


if __name__ == "__main__":
    main()
