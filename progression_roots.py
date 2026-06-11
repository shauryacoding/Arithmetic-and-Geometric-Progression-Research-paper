"""
Coefficient Conditions for Polynomials with Roots in Arithmetic
and Geometric Progressions

Symbolic verification and extension of results from:
  Mehra, S. "Coefficient Conditions for Polynomials with Roots in
  Arithmetic and Geometric Progressions" (2024)

This script:
  1. Verifies all identities proven in the paper (cubic and quartic, AP and GP)
  2. Corrects the quartic AP second identity (the paper's version has a sign error)
  3. Derives and verifies NEW degree-5 AP coefficient relations
  4. Provides a general framework for computing elimination ideals
     via parameter substitution

Dependencies: sympy >= 1.9
Run:  python progression_roots.py
      python progression_roots.py --degree 5 --type AP
"""

from itertools import combinations
from sympy import (
    symbols, expand, simplify, Rational,
    solve, Poly
)
import sys
import argparse
import time


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _prod(iterable):
    result = 1
    for x in iterable:
        result *= x
    return result


def vieta_coefficients(roots):
    """
    Return monic polynomial coefficients [a_{n-1}, ..., a_0] via Vieta's
    formulas given a list of symbolic root expressions.

    For P(x) = x^n + a_{n-1} x^{n-1} + ... + a_0:
        a_{n-k} = (-1)^k * e_k(roots)
    """
    n = len(roots)
    e = [
        sum(
            _prod(roots[i] for i in combo)
            for combo in combinations(range(n), k)
        )
        for k in range(n + 1)
    ]
    e[0] = 1
    return [(-1)**k * expand(e[k]) for k in range(1, n + 1)]


def check(expr, label, expected=0):
    """Expand expr, compare to expected, print PASS/FAIL."""
    val = expand(expr)
    ok = (val == expected)
    status = "PASS" if ok else "FAIL"
    print(f"  [{status}] {label}")
    if not ok:
        print(f"         Expected {expected}, got {val}")
    return ok


def section(title):
    w = 62
    print()
    print("=" * w)
    print(f"  {title}")
    print("=" * w)


def coeffs_from_roots(root_vals):
    """Return [a_{n-1}, ..., a_0] as integers for a numeric root list."""
    x = symbols('x')
    p = Poly(_prod(x - r for r in root_vals), x)
    return p.all_coeffs()[1:]   # drop leading 1


# ---------------------------------------------------------------------------
# Section 1: Cubic AP
# ---------------------------------------------------------------------------

def cubic_ap():
    section("CUBIC AP  —  identity:  2a2^3 - 9*a2*a1 + 27*a0 = 0")

    al, d = symbols('alpha d')
    roots = [al - d, al, al + d]
    a2, a1, a0 = vieta_coefficients(roots)

    check(2*a2**3 - 9*a2*a1 + 27*a0,
          "symbolic  (roots: alpha-d, alpha, alpha+d)")

    # Paper example: roots 1, 2, 3
    a2v, a1v, a0v = coeffs_from_roots([1, 2, 3])
    check(2*a2v**3 - 9*a2v*a1v + 27*a0v,
          f"numeric: roots 1,2,3  (a2={a2v}, a1={a1v}, a0={a0v})")


# ---------------------------------------------------------------------------
# Section 2: Quartic AP
# ---------------------------------------------------------------------------

def quartic_ap():
    section("QUARTIC AP  —  two coefficient relations")

    al, be = symbols('alpha beta')
    roots = [al - 3*be, al - be, al + be, al + 3*be]
    a3, a2, a1, a0 = vieta_coefficients(roots)

    print("\n  Identity 1:  a3^3 - 4*a3*a2 + 8*a1 = 0")
    check(a3**3 - 4*a3*a2 + 8*a1, "symbolic")

    # The paper's second identity contains an error.
    # The correct relation, derived by Groebner basis elimination, is:
    #   400*a0 - 22*a1*a3 - 36*a2^2 + 13*a2*a3^2 = 0
    print("\n  Identity 2 (CORRECTED from paper):")
    print("    400*a0 - 22*a1*a3 - 36*a2^2 + 13*a2*a3^2 = 0")
    print("    (The paper's original second identity is incorrect.)")
    check(400*a0 - 22*a1*a3 - 36*a2**2 + 13*a2*a3**2,
          "symbolic")

    examples = [
        ((-3, -1,  1,  3), "roots -3,-1,1,3"),
        ((-1,  0,  1,  2), "roots -1, 0,1,2"),
        (( 1,  2,  3,  4), "roots  1, 2,3,4"),
    ]
    print()
    for root_vals, label in examples:
        a3v, a2v, a1v, a0v = coeffs_from_roots(root_vals)
        check(a3v**3 - 4*a3v*a2v + 8*a1v,               f"id1 | {label}")
        check(400*a0v - 22*a1v*a3v - 36*a2v**2
              + 13*a2v*a3v**2,                           f"id2 | {label}")


# ---------------------------------------------------------------------------
# Section 3: Cubic GP
# ---------------------------------------------------------------------------

def cubic_gp():
    section("CUBIC GP  —  identity:  a1^3 = a2^3 * a0")

    a, r = symbols('a r')
    roots = [a/r, a, a*r]
    a2, a1, a0 = vieta_coefficients(roots)

    check(simplify(a1**3 - a2**3 * a0),
          "symbolic  (roots: a/r, a, a*r)")

    for root_vals, label in [
        ([1, 2, 4], "roots 1,2,4 (paper example)"),
        ([1, 3, 9], "roots 1,3,9"),
        ([2, 6, 18], "roots 2,6,18"),
    ]:
        a2v, a1v, a0v = coeffs_from_roots(root_vals)
        check(a1v**3 - a2v**3 * a0v, f"numeric: {label}")


# ---------------------------------------------------------------------------
# Section 4: Quartic GP
# ---------------------------------------------------------------------------

def quartic_gp():
    section("QUARTIC GP  —  two coefficient relations")

    a, q = symbols('a q', positive=True)
    roots = [a*q**Rational(-3, 2),
             a*q**Rational(-1, 2),
             a*q**Rational( 1, 2),
             a*q**Rational( 3, 2)]
    a3, a2, a1, a0 = vieta_coefficients(roots)

    print("\n  Identity 1:  a1^2 = a3^2 * a0")
    check(simplify(a1**2 - a3**2 * a0), "symbolic")

    print("\n  Identity 2:  a1^2*a3 - 2*a1*a2^2 + a1*a2*a3^2 - a1*a3^4 + a2^3*a3 = 0")
    id2_sym = simplify(
        a1**2*a3 - 2*a1*a2**2 + a1*a2*a3**2 - a1*a3**4 + a2**3*a3
    )
    check(id2_sym, "symbolic")

    examples = [
        ([1, 2,  4,  8 ], "roots 1,2,4,8   (ratio 2)"),
        ([1, 3,  9,  27], "roots 1,3,9,27  (ratio 3)"),
        ([2, 6, 18, 54 ], "roots 2,6,18,54 (ratio 3)"),
    ]
    print()
    for root_vals, label in examples:
        a3v, a2v, a1v, a0v = coeffs_from_roots(root_vals)
        check(a1v**2 - a3v**2 * a0v,                          f"id1 | {label}")
        check(a1v**2*a3v - 2*a1v*a2v**2
              + a1v*a2v*a3v**2 - a1v*a3v**4 + a2v**3*a3v,    f"id2 | {label}")


# ---------------------------------------------------------------------------
# Section 5: Degree-5 AP  (NEW — not in paper)
# ---------------------------------------------------------------------------

def quintic_ap():
    """
    Derive and verify the THREE coefficient relations for degree-5 polynomials
    whose roots form an arithmetic progression.

    This extends the paper beyond its stated scope and provides computational
    evidence for the n-2 dimension-counting heuristic (Section 4.2 of paper).

    Method
    ------
    Use the centred parametrisation rk = alpha + k*delta, k in {-2,-1,0,1,2}.

    Step 1: From the sum of roots,  a4 = -5*alpha  =>  alpha = -a4/5.
    Step 2: Substitute alpha into the a3 expression, solve for delta^2:
               delta^2 = (2*a4^2 - 5*a3) / 25
    Step 3: Substitute both into a2, a1, a0 to eliminate both parameters.

    This yields three polynomial relations among the five coefficients:

      R1:  25*a2 - 15*a3*a4 + 4*a4^3              = 0
      R2:  625*a1 - 100*a3^2 + 5*a3*a4^2 + 9*a4^4 = 0
      R3:  3125*a0 - 100*a3^2*a4 + 55*a3*a4^3
             - 7*a4^5                              = 0

    Three relations for n=5 matches the n-2 = 3 heuristic exactly.
    """
    section("DEGREE-5 AP  —  three new coefficient relations  (original)")

    al, ds = symbols('alpha delta')
    a4s, a3s, a2s, a1s, a0s = symbols('a4 a3 a2 a1 a0')

    roots = [al + k*ds for k in [-2, -1, 0, 1, 2]]
    a4, a3, a2, a1, a0 = vieta_coefficients(roots)

    # Step 1: alpha = -a4s / 5
    al_val = Rational(-1, 5) * a4s

    # Step 2: delta^2 from a3
    a3_sub = expand(a3.subs(al, al_val))
    ds2_expr = solve(a3_sub - a3s, ds**2)[0]

    # Step 3: eliminate delta^2
    a2_elim = expand(a2.subs(al, al_val).subs(ds**2, ds2_expr))
    a1_elim = expand(a1.subs(al, al_val).subs(ds**2, ds2_expr))
    a0_elim = expand(a0.subs(al, al_val).subs(ds**2, ds2_expr))

    # Relations (multiplied through to clear denominators)
    R1 = expand(25   * (a2s - a2_elim))   # 25*a2 - 15*a3*a4 + 4*a4^3
    R2 = expand(625  * (a1s - a1_elim))   # 625*a1 - 100*a3^2 + 5*a3*a4^2 + 9*a4^4
    R3 = expand(3125 * (a0s - a0_elim))   # 3125*a0 - 100*a3^2*a4 + 55*a3*a4^3 - 7*a4^5

    print("\n  R1:  25*a2 - 15*a3*a4 + 4*a4^3 = 0")
    print("  R2:  625*a1 - 100*a3^2 + 5*a3*a4^2 + 9*a4^4 = 0")
    print("  R3:  3125*a0 - 100*a3^2*a4 + 55*a3*a4^3 - 7*a4^5 = 0")

    # Symbolic checks (substitute back into parametric expressions)
    def eval_R(R, a4_expr, a3_expr, a2_expr, a1_expr, a0_expr):
        return expand(R.subs([
            (a4s, a4_expr), (a3s, a3_expr),
            (a2s, a2_expr), (a1s, a1_expr), (a0s, a0_expr)
        ]))

    # Substitute the parametric expressions for all coefficients
    check(eval_R(R1, a4, a3, a2, a1, a0), "R1 symbolic")
    check(eval_R(R2, a4, a3, a2, a1, a0), "R2 symbolic")
    check(eval_R(R3, a4, a3, a2, a1, a0), "R3 symbolic")

    # Numeric verification
    examples = [
        ((-2, -1,  0,  1,  2), "roots -2,-1,0,1,2"),
        (( 1,  2,  3,  4,  5), "roots  1, 2,3,4,5"),
        (( 0,  3,  6,  9, 12), "roots  0, 3,6,9,12"),
        ((-5,  0,  5, 10, 15), "roots -5, 0,5,10,15"),
    ]
    print()
    for root_vals, label in examples:
        a4v, a3v, a2v, a1v, a0v = coeffs_from_roots(root_vals)
        r1v = 25*a2v - 15*a3v*a4v + 4*a4v**3
        r2v = 625*a1v - 100*a3v**2 + 5*a3v*a4v**2 + 9*a4v**4
        r3v = 3125*a0v - 100*a3v**2*a4v + 55*a3v*a4v**3 - 7*a4v**5
        check(r1v, f"R1 | {label}")
        check(r2v, f"R2 | {label}")
        check(r3v, f"R3 | {label}")

    print()
    print("  n=5 yields 3 = n-2 relations, consistent with the heuristic")
    print("  in Section 4.2 of the paper (Jacobian rank / dimension count).")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def run_all():
    print()
    print("=" * 62)
    print("  Progression Roots — Symbolic Verification Suite")
    print("  Companion code for Mehra (2024)")
    print("=" * 62)

    t0 = time.time()
    cubic_ap()
    quartic_ap()
    cubic_gp()
    quartic_gp()
    quintic_ap()
    elapsed = time.time() - t0

    print()
    print("=" * 62)
    print(f"  All checks complete  ({elapsed:.1f}s)")
    print("=" * 62)
    print()


def run_single(degree, prog_type):
    dispatch = {
        (3, "AP"): cubic_ap,
        (4, "AP"): quartic_ap,
        (3, "GP"): cubic_gp,
        (4, "GP"): quartic_gp,
        (5, "AP"): quintic_ap,
    }
    fn = dispatch.get((degree, prog_type.upper()))
    if fn is None:
        print(f"No result for degree={degree}, type={prog_type}.")
        print("Available: (3,AP), (4,AP), (3,GP), (4,GP), (5,AP)")
        sys.exit(1)
    fn()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Verify coefficient identities for AP/GP root polynomials"
    )
    parser.add_argument("--degree", type=int, help="Degree: 3, 4, or 5")
    parser.add_argument("--type",   type=str, help="AP or GP")
    args = parser.parse_args()

    if args.degree or args.type:
        if not (args.degree and args.type):
            print("Provide both --degree and --type, or neither to run all.")
            sys.exit(1)
        run_single(args.degree, args.type)
    else:
        run_all()
