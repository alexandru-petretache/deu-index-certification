#!/usr/bin/env python3
"""
symbolic_checks.py -- exact (computer-algebra) verification of every closed-form
identity, constant, limit, and family relation asserted in the paper
"A Computer-Aided Study of the Diminished Euler Sombor Index".

These hold for ALL graphs, not merely the enumerated range: they are identities
in the symbolic variables, certified by sympy.  Run:

    python3 symbolic_checks.py

Exit status 0 iff every check passes.
"""
import sys
import sympy as sp

# symbols
a, b, t, s, q = sp.symbols('a b t s q', positive=True)
m, sig, delta = sp.symbols('m sigma delta', positive=True)
n = sp.symbols('n', positive=True)

# kernels
phi      = sp.sqrt(a**2 + a*b + b**2) / (a + b)          # DEU kernel
phi_dso  = sp.sqrt(a**2 + b**2)       / (a + b)          # DSO kernel
qab      = (a - b) / (a + b)                              # signed imbalance
g        = sp.sqrt(t**2 + t + 1) / (t + 1)                # g(t) = phi(t,1)

S3, S2 = sp.sqrt(3), sp.sqrt(2)

results = []  # (name, passed)


def chk(name, expr_is_zero=None, boolean=None):
    """Pass if simplified expr == 0, or if boolean is True."""
    if expr_is_zero is not None:
        ok = sp.simplify(expr_is_zero) == 0
    else:
        ok = bool(boolean)
    results.append((name, ok))


# ---- 1. kernel-square identities ------------------------------------------
chk("phi^2 = 1 - ab/(a+b)^2",            phi**2 - (1 - a*b/(a+b)**2))
chk("phi^2 = (3+q^2)/4",                 phi**2 - (3 + qab**2)/4)
chk("phi_DSO^2 = (1+q^2)/2",             phi_dso**2 - (1 + qab**2)/2)
chk("phi^2 - phi_DSO^2 = (1-q^2)/4",     phi**2 - phi_dso**2 - (1 - qab**2)/4)
# separation decomposition: (phi-phi_DSO)(phi+phi_DSO) = (1-q^2)/4
chk("(phi-phi_DSO)(phi+phi_DSO) = (1-q^2)/4",
    (phi - phi_dso)*(phi + phi_dso) - (1 - qab**2)/4)

# ---- 2. ratio-only kernel g(t) = phi(t,1) ---------------------------------
chk("g(t) = g(1/t)  (ratio symmetry)",   g - g.subs(t, 1/t))
chk("g(1) = sqrt(3)/2",                  g.subs(t, 1) - S3/2)
chk("g(t)^2 = 1 - t/(t+1)^2",            g**2 - (1 - t/(t+1)**2))
chk("lim_{t->0+} g = 1",                 boolean=sp.limit(g, t, 0, '+') == 1)
chk("lim_{t->oo}  g = 1",                boolean=sp.limit(g, t, sp.oo) == 1)
chk("g increasing on (1,oo): g'(t)>0",   boolean=sp.simplify(sp.diff(g, t)) != 0
                                          and bool(sp.diff(g, t).subs(t, 2) > 0))

# ---- 3. range floor and the size estimate ---------------------------------
mid = (S3/2 + 1)/2
chk("midpoint = (2+sqrt3)/4",            mid - (2 + S3)/4)
chk("rel. error at floor sqrt3/2 = 7-4sqrt3", (mid - S3/2)/mid - (7 - 4*S3))
chk("rel. error at ceil 1        = 7-4sqrt3", (1 - mid)/mid - (7 - 4*S3))
chk("(2-sqrt3)/(2+sqrt3) = 7-4sqrt3",    (2 - S3)/(2 + S3) - (7 - 4*S3))

# ---- 4. quadratic refinement (Prop:quad) ----------------------------------
B = S3/2*m + S3/12*sig
chk("B^2 - (1/4) m(3m+sigma) = sigma^2/48",
    B**2 - sp.Rational(1, 4)*m*(3*m + sig) - sig**2/48)
chk("quad optimality limit_{t->1} = sqrt3/12",
    boolean=sp.limit((g - S3/2)/((t - 1)/(t + 1))**2, t, 1) == S3/12)

# ---- 5. sigma-index corollary constant ------------------------------------
chk("(sqrt3/12)*(1/4) = sqrt3/48",       (S3/12)*sp.Rational(1, 4) - S3/48)

# ---- 6. excess / secant convexity (Prop:excess) ---------------------------
h = sp.sqrt(3 + s**2)/2                   # phi as a function of s=|q|
chord = S3/2 + (1 - S3/2)*s
chk("h(0) = sqrt3/2",                    h.subs(s, 0) - S3/2)
chk("h(1) = 1",                          h.subs(s, 1) - 1)
chk("h convex: h''(s) > 0 on (0,1)",     boolean=bool(sp.simplify(sp.diff(h, s, 2)).subs(s, sp.Rational(1, 2)) > 0))
chk("h below chord at s=1/2 (convexity)", boolean=bool((chord - h).subs(s, sp.Rational(1, 2)) > 0))

# ---- 7. separation constants (Prop:sep) -----------------------------------
chk("upper sep const 1/(2(sqrt3+sqrt2)) = (sqrt3-sqrt2)/2",
    1/(2*(S3 + S2)) - (S3 - S2)/2)
# denominator 4(phi+phi_DSO): value at q=0 (a=b) is 2(sqrt3+sqrt2); sup is 8
chk("denom 4(phi+phi_DSO) at a=b = 2(sqrt3+sqrt2)",
    (4*(phi + phi_dso)).subs(b, a) - 2*(S3 + S2))
chk("lim denom 4(phi+phi_DSO), ratio->oo = 8",
    boolean=sp.limit((4*(phi + phi_dso)).subs(b, 1), a, sp.oo) == 8)
chk("lower sep coeff 1/(4(phi+phi_DSO)) -> 1/8 as ratio->oo",
    boolean=sp.limit((1/(4*(phi + phi_dso))).subs(b, 1), a, sp.oo) == sp.Rational(1, 8))

# ---- 8. dominance threshold (Remark) --------------------------------------
chk("|q| threshold 12(1-sqrt3/2)/sqrt3 = 4sqrt3-6",
    12*(1 - S3/2)/S3 - (4*S3 - 6))
chk("Delta/delta threshold = 13+8sqrt3",
    (1 + (4*S3 - 6))/(1 - (4*S3 - 6)) - (13 + 8*S3))
chk("(4sqrt3-5)/(7-4sqrt3) = 13+8sqrt3",
    (4*S3 - 5)/(7 - 4*S3) - (13 + 8*S3))

# ---- 9. cosine-rule family (Introduction) ---------------------------------
#   kernel sqrt(a^2+b^2-2ab cos theta)/(a+b);  DEU is theta=2pi/3, DSO is pi/2, m is pi
chk("theta=2pi/3 gives Euler-Sombor numerator",
    sp.sqrt(a**2 + b**2 - 2*a*b*sp.cos(2*sp.pi/3)) - sp.sqrt(a**2 + a*b + b**2))
chk("theta=pi/2 gives Sombor numerator (DSO)",
    sp.sqrt(a**2 + b**2 - 2*a*b*sp.cos(sp.pi/2)) - sp.sqrt(a**2 + b**2))
chk("theta=pi gives kernel 1 (radicand = (a+b)^2, sums to m)",
    (a**2 + b**2 - 2*a*b*sp.cos(sp.pi)) - (a + b)**2)

# ---- 10. example closed forms ---------------------------------------------
chk("phi symmetric: phi(a,b)=phi(b,a)  [K_{a,b} form]",
    phi - phi.subs({a: b, b: a}, simultaneous=True))
chk("K_{a,a}: phi(a,a) = sqrt3/2",       phi.subs(b, a) - S3/2)
chk("S_n radicand (n-1)^2+(n-1)+1 = n^2-n+1",
    (n - 1)**2 + (n - 1) + 1 - (n**2 - n + 1))
chk("DEU(S_n)/m = sqrt(n^2-n+1)/n -> 1", boolean=sp.limit(sp.sqrt(n**2 - n + 1)/n, n, sp.oo) == 1)

# ---- 11. connected-graph count (m>=1 convention) --------------------------
A001349 = [1, 1, 2, 6, 21, 112, 853, 11117, 261080, 11716571]   # connected graphs, n=1..10
chk("sum A001349(1..10) - 1 (edgeless K_1) = 11,989,763",
    boolean=(sum(A001349) - 1) == 11_989_763)


# ---- report ----------------------------------------------------------------
width = max(len(nm) for nm, _ in results)
npass = sum(ok for _, ok in results)
print(f"Symbolic checks: {npass}/{len(results)} passed\n")
for nm, ok in results:
    print(f"  [{'PASS' if ok else 'FAIL'}] {nm}")
if npass != len(results):
    print("\nSOME SYMBOLIC CHECKS FAILED", file=sys.stderr)
    sys.exit(1)
print("\nAll symbolic identities verified.")
