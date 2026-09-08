"""
test_abox_independence.py

Check whether a_box only changes the internal normalized representation of the
nonlinear invariant, while the reconstructed physical polynomials Ix and Iy
remain unchanged.

The test intentionally uses:
    m = 8
    d = 1
"""

import numpy as np
import sympy as sp

import linear_lattice as lin
import nonlinear as nl


M = 8
D = 1

RTOL = 1e-8
ATOL = 1e-12

REFERENCE_BOX = np.array(
    [0.05e-2, 3.5e-3, 1.0e-3, 0.8e-3, 0.6e-3],
    dtype=float,
)

A_BOXES = [
    REFERENCE_BOX,
    0.5 * REFERENCE_BOX,
    2.0 * REFERENCE_BOX,
    REFERENCE_BOX * np.array([1.5, 0.8, 1.3, 0.7, 1.2]),
]


def coefficient_dict(poly):
    """Physical polynomial coefficients indexed by exponent tuple."""
    P = sp.Poly(sp.expand(poly), *nl.vars)

    return {
        monomial: float(coeff)
        for monomial, coeff in P.terms()
    }


def compare_polynomials(reference, candidate):
    """Compare two physical polynomials coefficient by coefficient."""
    ref = coefficient_dict(reference)
    new = coefficient_dict(candidate)

    keys = sorted(set(ref) | set(new))

    r = np.array([ref.get(k, 0.0) for k in keys], dtype=float)
    c = np.array([new.get(k, 0.0) for k in keys], dtype=float)

    diff = c - r
    abs_diff = np.abs(diff)

    ref_norm = np.linalg.norm(r)
    diff_norm = np.linalg.norm(diff)

    relative_l2 = diff_norm / max(ref_norm, 1e-300)

    scale = ATOL + RTOL * np.maximum(np.abs(r), np.abs(c))
    scaled_error = abs_diff / scale

    worst = int(np.argmax(scaled_error)) if len(keys) else 0

    return {
        "same": bool(np.all(scaled_error <= 1.0)),
        "max_abs": float(np.max(abs_diff)) if len(keys) else 0.0,
        "relative_l2": float(relative_l2),
        "max_scaled": float(np.max(scaled_error)) if len(keys) else 0.0,
        "worst_monomial": keys[worst] if len(keys) else None,
        "reference_coeff": float(r[worst]) if len(keys) else 0.0,
        "candidate_coeff": float(c[worst]) if len(keys) else 0.0,
    }


def print_result(name, result):
    status = "PASS" if result["same"] else "FAIL"

    print(f"  {name}: {status}")
    print(f"    max |Δc|       = {result['max_abs']:.6e}")
    print(f"    ||Δc||/||c||   = {result['relative_l2']:.6e}")
    print(f"    max scaled err = {result['max_scaled']:.6e}")
    print(f"    worst monomial = {result['worst_monomial']}")
    print(
        "    coefficient    = "
        f"{result['reference_coeff']:.16e} -> "
        f"{result['candidate_coeff']:.16e}"
    )


def main():

    # Prepare the lattice ONCE.  Only a_box changes between nonlinear runs.
    magnets, cell, data, correction, p = lin.prepare_lattice()

    print("=" * 72)
    print("a_box INDEPENDENCE TEST")
    print(f"m = {M}, d = {D}")
    print(f"RTOL = {RTOL:.1e}, ATOL = {ATOL:.1e}")
    print("=" * 72)

    polynomials = []

    for i, a_box in enumerate(A_BOXES):

        print(f"\nComputing box {i}:")
        print(a_box)

        Ix, Iy = nl.invariant_polynomial(
            cell,
            data,
            order=M,
            delta_order=D,
            a_box=a_box,
        )

        polynomials.append((Ix, Iy))

    Ix_ref, Iy_ref = polynomials[0]

    all_pass = True

    for i in range(1, len(A_BOXES)):

        print("\n" + "-" * 72)
        print(f"REFERENCE BOX vs BOX {i}")

        Ix_test, Iy_test = polynomials[i]

        rx = compare_polynomials(Ix_ref, Ix_test)
        ry = compare_polynomials(Iy_ref, Iy_test)

        print_result("Ix", rx)
        print_result("Iy", ry)

        all_pass &= rx["same"] and ry["same"]

    print("\n" + "=" * 72)

    if all_pass:
        print("PASS")
        print(
            "Within the chosen tolerance, changing a_box only changes the "
            "internal normalized representation; the reconstructed physical "
            "polynomials are unchanged."
        )
    else:
        print("FAIL")
        print(
            "The reconstructed physical polynomial changes with a_box. "
            "Therefore a_box is affecting more than an internal relabeling "
            "at the chosen tolerance."
        )

    print("=" * 72)


if __name__ == "__main__":
    main()
