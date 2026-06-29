"""
deu_core.py -- shared invariants and bound checks for the
Diminished Euler Sombor (DEU) index certification.

Every numeric check below corresponds to a labelled statement in the paper
"A Computer-Aided Study of the Diminished Euler Sombor Index".  The mapping
proposition <-> check is given in the `name` field of each Check.

Conventions (paper Section II): G is a simple connected graph with m >= 1 edges;
d_u is the degree of u; Delta, delta are the max/min degree.  For an edge uv,
    phi(d_u,d_v)     = sqrt(d_u^2 + d_u d_v + d_v^2) / (d_u + d_v)   (DEU kernel)
    phi_DSO(d_u,d_v) = sqrt(d_u^2 + d_v^2)         / (d_u + d_v)   (DSO kernel)
    q                = (d_u - d_v) / (d_u + d_v)                   (signed imbalance)
"""

from __future__ import annotations
from dataclasses import dataclass
from math import sqrt, isclose
import networkx as nx

# ----- exact constants used in the paper -----------------------------------
S3 = sqrt(3.0)
S2 = sqrt(2.0)
PHI_MIN = S3 / 2.0                       # sqrt(3)/2, the kernel floor
MID = (2.0 + S3) / 4.0                   # midpoint estimate coefficient, DEU_hat = MID * m
EPS_REL = 7.0 - 4.0 * S3                 # 7 - 4 sqrt(3) ~= 0.0718, size-estimate worst case
C_QUAD = S3 / 12.0                       # sqrt(3)/12, Prop. quadratic refinement
SEP_LO = 1.0 / 8.0                       # 1/8, separation lower constant (strict)
SEP_HI = (S3 - S2) / 2.0                 # (sqrt3 - sqrt2)/2, separation upper constant
DOM_RATIO = 13.0 + 8.0 * S3             # 13 + 8 sqrt(3) ~= 26.86, sigma*-bound dominance threshold

TOL = 1e-9                               # absolute float tolerance for (in)equalities


# ----- per-edge kernels ------------------------------------------------------
def phi(a: int, b: int) -> float:
    return sqrt(a * a + a * b + b * b) / (a + b)


def phi_dso(a: int, b: int) -> float:
    return sqrt(a * a + b * b) / (a + b)


# ----- graph invariants ------------------------------------------------------
@dataclass
class Invariants:
    n: int
    m: int
    Delta: int
    delta: int
    DEU: float
    DSO: float
    sigma_star: float        # sum q^2
    irr_star: float          # sum |q|
    irr: float               # Albertson  sum |d_u - d_v|
    sigma: float             # sigma index sum (d_u - d_v)^2
    sum_phi2: float          # sum phi^2  (for the second-moment identity)
    regular: bool
    biregular_class: bool    # every edge joins a delta-vertex to a Delta-vertex


def invariants(G: "nx.Graph") -> Invariants:
    deg = dict(G.degree())
    m = G.number_of_edges()
    if m < 1:
        raise ValueError("graph has no edges (m >= 1 required)")
    Delta = max(deg.values())
    delta = min(deg.values())
    DEU = DSO = sigma_star = irr_star = irr = sigma = sum_phi2 = 0.0
    biregular = True
    for u, v in G.edges():
        du, dv = deg[u], deg[v]
        DEU += phi(du, dv)
        DSO += phi_dso(du, dv)
        q = (du - dv) / (du + dv)
        sigma_star += q * q
        irr_star += abs(q)
        irr += abs(du - dv)
        sigma += (du - dv) ** 2
        p = phi(du, dv)
        sum_phi2 += p * p
        if {du, dv} != {delta, Delta}:
            biregular = False
    regular = (Delta == delta)
    return Invariants(G.number_of_nodes(), m, Delta, delta, DEU, DSO,
                      sigma_star, irr_star, irr, sigma, sum_phi2,
                      regular, biregular)


# ----- a single check result -------------------------------------------------
@dataclass
class Check:
    name: str            # paper reference, e.g. "Prop:global (upper)"
    ok: bool
    slack: float         # signed margin; >= -TOL means satisfied
    detail: str = ""


def _ge(lhs: float, rhs: float, name: str, detail: str = "") -> Check:
    """Assert lhs >= rhs (with tolerance).  slack = lhs - rhs."""
    s = lhs - rhs
    return Check(name, s >= -TOL, s, detail)


def _le(lhs: float, rhs: float, name: str, detail: str = "") -> Check:
    """Assert lhs <= rhs (with tolerance).  slack = rhs - lhs."""
    s = rhs - lhs
    return Check(name, s >= -TOL, s, detail)


def _iff(equality: bool, in_class: bool, name: str) -> Check:
    """Biconditional: bound is tight  <=>  graph lies in the stated class."""
    ok = (equality == in_class)
    return Check(name, ok, 0.0 if ok else -1.0,
                 f"equality={equality}, in_class={in_class}")


def check_graph(G: "nx.Graph") -> list[Check]:
    """Run every paper inequality / identity / equality-class biconditional on G."""
    I = invariants(G)
    m, D, d = I.m, I.Delta, I.delta
    out: list[Check] = []

    # -- Prop: global range  (sqrt3/2) m <= DEU <= sqrt(D^2+Dd+d^2)/(D+d) m < m
    mid_kernel = sqrt(D * D + D * d + d * d) / (D + d)
    out.append(_ge(I.DEU, PHI_MIN * m, "Prop:global (lower)"))
    out.append(_le(I.DEU, mid_kernel * m, "Prop:global (upper)"))
    out.append(_le(I.DEU, float(m), "Prop:global (DEU<m)"))
    out.append(_iff(isclose(I.DEU, PHI_MIN * m, abs_tol=TOL), I.regular,
                    "Prop:global lower equality <=> regular"))
    out.append(_iff(isclose(I.DEU, mid_kernel * m, abs_tol=TOL), I.biregular_class,
                    "Prop:global upper equality <=> regular-or-biregular"))

    # -- Cor: size estimate   |DEU - MID m| <= EPS_REL * MID m
    DEU_hat = MID * m
    out.append(_le(abs(I.DEU - DEU_hat), EPS_REL * DEU_hat, "Cor:est size estimate"))

    # -- Prop: second-moment identity   sum phi^2 = (3m + sigma*)/4   (exact)
    out.append(Check("Prop:CS second-moment identity (exact)",
                     isclose(I.sum_phi2, (3 * m + I.sigma_star) / 4.0, abs_tol=TOL),
                     I.sum_phi2 - (3 * m + I.sigma_star) / 4.0))
    # -- Prop: CS upper bound   DEU <= (1/2) sqrt(m (3m + sigma*))
    out.append(_le(I.DEU, 0.5 * sqrt(m * (3 * m + I.sigma_star)), "Prop:CS bound"))

    # -- Prop: quadratic refinement   DEU - (sqrt3/2) m <= (sqrt3/12) sigma*
    out.append(_le(I.DEU - PHI_MIN * m, C_QUAD * I.sigma_star, "Prop:quad refinement"))

    # -- Cor: sigma-index bound   DEU <= (sqrt3/2) m + sqrt3/(48 d^2) sigma
    out.append(_le(I.DEU, PHI_MIN * m + (S3 / (48.0 * d * d)) * I.sigma,
                   "Cor:quad sigma-index bound"))

    # -- Prop: excess via Albertson irregularity
    #     0 <= DEU - (sqrt3/2) m <= (1 - sqrt3/2) irr* <= (1 - sqrt3/2)/(2 delta) irr
    excess = I.DEU - PHI_MIN * m
    out.append(_ge(excess, 0.0, "Prop:excess (lower 0)"))
    out.append(_le(excess, (1 - PHI_MIN) * I.irr_star, "Prop:excess (irr* bound)"))
    out.append(_le((1 - PHI_MIN) * I.irr_star,
                   (1 - PHI_MIN) / (2.0 * d) * I.irr, "Prop:excess (Albertson relax)"))

    # -- Prop: normalized-irregularity bounds
    out.append(_le(I.irr_star, m * (D - d) / (D + d), "Prop:normirr (irr*)"))
    out.append(_le(I.sigma_star, m * ((D - d) / (D + d)) ** 2, "Prop:normirr (sigma*)"))

    # -- Prop: interlock   irr*^2/m <= sigma* <= ((D-d)/(D+d)) irr*
    out.append(_le(I.irr_star ** 2 / m, I.sigma_star, "Prop:interlock (lower, CS)"))
    rhs_int = (D - d) / (D + d) * I.irr_star if D != d else 0.0
    out.append(_le(I.sigma_star, rhs_int, "Prop:interlock (upper)"))

    # -- Prop: separation   (1/8)(m - sigma*) < DEU - DSO <= ((sqrt3-sqrt2)/2)(m - sigma*)
    gap = I.DEU - I.DSO
    mms = m - I.sigma_star
    out.append(_ge(gap, SEP_LO * mms, "Prop:sep (lower, strict)"))
    out.append(_le(gap, SEP_HI * mms, "Prop:sep (upper)"))
    out.append(_iff(isclose(gap, SEP_HI * mms, abs_tol=TOL), I.regular,
                    "Prop:sep upper equality <=> regular"))
    # lower bound is strict: for a graph with >=1 edge it is never attained
    out.append(Check("Prop:sep lower strict (gap > LO(m-sigma*))",
                     gap - SEP_LO * mms > TOL or m == 0,
                     gap - SEP_LO * mms))

    # -- Remark: sigma*-bound dominates irr*-bound when Delta/delta <= 13 + 8 sqrt3
    if d > 0 and (D / d) <= DOM_RATIO + TOL:
        out.append(_le(C_QUAD * I.sigma_star, (1 - PHI_MIN) * I.irr_star,
                       "Remark:dominance (sigma*-bound sharper, Delta/delta<=13+8sqrt3)"))

    return out
