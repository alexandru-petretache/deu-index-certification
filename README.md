# Certification scripts — Diminished Euler Sombor (DEU) index

Reproducibility artifact for the paper
**"A Computer-Aided Study of the Diminished Euler Sombor Index."**

The index is

```
DEU(G) = sum_{uv in E(G)}  sqrt(d_u^2 + d_u d_v + d_v^2) / (d_u + d_v),
```

with per-edge kernel `phi(d_u, d_v) = sqrt(d_u^2 + d_u d_v + d_v^2)/(d_u + d_v)`,
which depends only on the degree ratio and lies in `[sqrt(3)/2, 1)`.

These scripts independently re-derive and certify **every** identity, bound, and
equality condition stated in the paper. The certification has two complementary
layers.

## Layer 0 — Symbolic (exact, holds for all graphs)

`symbolic_checks.py` verifies, in a computer-algebra system (SymPy), every
closed-form claim as an exact symbolic identity in the degree variables — so they
hold for **all** graphs, not merely a finite range. Among the 37 checks:

- kernel identities `phi^2 = (3 + q^2)/4`, `phi_DSO^2 = (1 + q^2)/2`,
  `phi^2 - phi_DSO^2 = (1 - q^2)/4`, where `q = (d_u - d_v)/(d_u + d_v)`;
- the ratio-only kernel `g(t) = phi(t,1)`: symmetry `g(t) = g(1/t)`,
  `g(1) = sqrt(3)/2`, `g(t)^2 = 1 - t/(t+1)^2`, and limits `g -> 1`;
- the size-estimate worst case `7 - 4 sqrt(3)` (at **both** ends of the range);
- the quadratic-refinement identity `B^2 - (1/4) m(3m + sigma*) = (sigma*)^2/48`
  and its optimality limit `sqrt(3)/12`;
- convexity of `phi` in `|q|` (the secant bound) and the `sqrt(3)/48` constant;
- the separation constants `1/8` and `(sqrt3 - sqrt2)/2 = 1/(2(sqrt3 + sqrt2))`,
  and the denominator range `[2(sqrt3 + sqrt2), 8)`;
- the dominance threshold `Delta/delta <= 13 + 8 sqrt(3)`;
- the cosine-rule family relations (`theta = 2pi/3` is DEU, `pi/2` is DSO, `pi`
  gives the constant kernel summing to `m`);
- the example closed forms for `K_{a,b}`, `K_{a,a}`, and the star `S_n`.

## Layer 1 — Exhaustive enumeration (`n <= 10`)

`enumerate_check.py` streams **every** connected graph on `n` vertices from
nauty's `geng`, builds it in NetworkX, and runs the full per-graph battery
(`deu_core.check_graph`). The complete run covers all

```
11,989,763 connected graphs on n <= 10 vertices,
```

which is `sum_{n=1}^{10} A001349(n) = 11,989,764` minus the single edgeless `K_1`
(excluded by the paper's `m >= 1` convention). Every inequality is checked, and
every equality condition is checked as an **exact biconditional** (the bound is
tight **iff** the graph lies in the stated class). The script also reports the
tightest margin observed per bound, making visible which bounds are attained in
range and which (e.g. the strict lower separation bound) are only approached.

## Layer 2 — Random sampling (`n > 10`) — *the extension beyond exhaustive range*

Exhaustive enumeration is infeasible past `n = 10`. `random_check.py` extends the
empirical certification to **larger** `n` by sampling random connected graphs and
running the same battery on each. To exercise the bounds across genuinely
different degree regimes — not one model's typical graph — it mixes random trees,
random `d`-regular graphs (to hit the regular equality cases at large `n`),
Barabasi–Albert graphs (hub-heavy: large `Delta`, the near-star regime where the
lower separation bound is approached), dense Erdos–Renyi graphs, Watts–Strogatz
graphs, and random complete bipartite graphs (the `(delta, Delta)`-biregular
equality class). Runs are seeded and reproducible.

## Files

| file                  | role                                                            |
|-----------------------|-----------------------------------------------------------------|
| `deu_core.py`         | shared invariants + `check_graph` (every per-graph check)       |
| `symbolic_checks.py`  | Layer 0: exact symbolic identities (SymPy)                      |
| `enumerate_check.py`  | Layer 1: exhaustive enumeration over `n <= 10` (needs `geng`)   |
| `random_check.py`     | Layer 2: random sampling for `n > 10`                           |
| `requirements.txt`    | Python pins + nauty note                                        |

## Install

```bash
pip install -r requirements.txt
# nauty (system dependency, only for enumerate_check.py):
#   Debian/Ubuntu:  apt-get install nauty     # provides `nauty-geng`
#   macOS:          brew install nauty
#   source:         https://pallini.di.uniroma1.it/
```

## Run

```bash
# Layer 0 — exact identities (seconds)
python3 symbolic_checks.py

# Layer 1 — exhaustive. Quick self-test on small n:
python3 enumerate_check.py --nmax 8
# Full certification over all 11,989,763 graphs on n <= 10 (CPU-heavy: n=10
# alone is 11.7M graphs — expect a long run):
python3 enumerate_check.py

# Layer 2 — random sampling for larger n (reproducible via --seed)
python3 random_check.py --nmax 60
```

Each script exits `0` iff all its checks pass (symbolic: all identities verified;
numeric: zero violations).

## Check → paper statement

| check name (in `deu_core` / output)                  | paper statement                                  |
|------------------------------------------------------|--------------------------------------------------|
| `Prop:global (lower/upper/DEU<m)` + equality iffs    | global range `sqrt3/2 m <= DEU < m`              |
| `Cor:est size estimate`                              | edge-count size estimate, error `7 - 4 sqrt3`    |
| `Prop:CS second-moment identity (exact)`             | `sum phi^2 = (3m + sigma*)/4`                     |
| `Prop:CS bound`                                      | Cauchy–Schwarz upper bound                       |
| `Prop:quad refinement`                               | `DEU - sqrt3/2 m <= sqrt3/12 sigma*`             |
| `Cor:quad sigma-index bound`                         | `sqrt3/(48 delta^2) sigma` refinement            |
| `Prop:excess (...)`                                  | excess vs normalized / Albertson irregularity    |
| `Prop:normirr (irr*/sigma*)`                         | normalized-irregularity bounds                   |
| `Prop:interlock (lower/upper)`                       | `irr*^2/m <= sigma* <= ((D-d)/(D+d)) irr*`       |
| `Prop:sep (lower/upper)` + equality iff              | DEU–DSO separation, constants `1/8`, `(√3−√2)/2` |
| `Remark:dominance ...`                               | `sigma*`-bound sharper for `Delta/delta <= 13+8√3` |

All inequalities use an absolute tolerance of `1e-9`; observed float noise sits
at machine epsilon (`~1e-15`), comfortably inside it.
