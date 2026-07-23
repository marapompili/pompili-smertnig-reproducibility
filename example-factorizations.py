"""
Reproduction of the factorization examples from Pompili--Smertnig.

This script demonstrates the factorization algorithms implemented for
Banff cluster algebras in SageMath. It reproduces two examples:

1. a rank-eight Banff cluster algebra over the Gaussian rational field
   Q(i), where all factorizations of a specified element are computed;

2. a cluster algebra of type A3 over Q, where a divisor is assembled from
   selected prime divisors and converted to a principal generator whose
   factorization invariants are then computed.

For each example, the script prints the element, the number of
factorizations, the number of atoms dividing it, and its set of lengths.
The first example additionally prints every factorization in LaTeX form.

Run with SageMath:

    sage factorization_examples.py
"""

from sage.algebras.banff_cluster_algebra import BanffClusterAlgebra
from sage.matrix.constructor import matrix
from sage.misc.latex import latex 
from sage.rings.number_field.number_field import QuadraticField
from sage.rings.rational_field import QQ

def print_invariants(tag, CF):
    print("\n" + "="*70)
    print(f"[{tag}]")
    print("-"*70)
    print(f"Element f = {latex(CF.element())}")
    print(f"Total number of factorizations: {CF.num_factorizations()}")
    print(f"Number of atoms dividing f: {CF.num_atoms()}")
    print(f"Set of lengths L(f): {CF.set_of_lengths()}")


def factorization_rhs_latex(CF, fac):
    pieces = []

    for atom, exp in fac:
        atom_ltx = latex(atom)
        if exp == 1:
            pieces.append(rf"\left({atom_ltx}\right)")
        else:
            pieces.append(rf"\left({atom_ltx}\right)^{{{exp}}}")

    return r" \cdot ".join(pieces) if pieces else "1"


def print_selected_factorizations(CF, indices_1based):
    facs = CF.sort_all_factorizations()
    f_ltx = latex(CF.element())

    print("\nSelected factorizations (ordered with sort_all_factorizations):")
    for k in indices_1based:
        if k < 1 or k > len(facs):
            print(f"  - index {k} out of range (1..{len(facs)})")
            continue

        rhs = factorization_rhs_latex(CF, facs[k-1])

        print(f"\n  #{k}:")
        print(r"\[" + rf"{f_ltx} = {rhs}" + r"\]")


def print_all_factorizations(tag, CF):
    print("\nAll factorizations:")
    print(CF.latex(env="aligned", include_display_math=True))


# ----------------------------------------------------------------------
# Example (1)
# ----------------------------------------------------------------------

def run_example_1():

    B = matrix([
        [ 0, 0, 0, 2, 0, 0, 0, 0],
        [ 0, 0, 0,-1, 1, 0, 0, 0],
        [ 0, 0, 0, 0,-1, 0, 0, 0],
        [-2, 1, 0, 0, 0, 1,-1, 0],
        [ 0,-1, 1, 0, 0, 0, 1,-1],
        [ 0, 0, 0,-1, 0, 0, 0, 0],
        [ 0, 0, 0, 1,-1, 0, 0, 0],
        [ 0, 0, 0, 0, 1, 0, 0, 0],
    ])

    A = BanffClusterAlgebra(B, scalars=QuadraticField(-1))
    x = A.gens()  

    f_expr = (x[3]**3*x[4] + x[3]**2*x[4]**2 + x[3]**3 + x[3]**2*x[4] + x[3]*x[4] + x[4]**2 + x[3] + x[4]) / x[1]
    f = A(f_expr)

    CF = f.factor()

    print_invariants("Example 1", CF)
    print_all_factorizations("Example 1", CF)


# ----------------------------------------------------------------------
# Example (2)
# ----------------------------------------------------------------------

def run_example_2():
    B = matrix([
        [ 0, 1,  0],
        [ -1,  0, 1],
        [0,  -1,  0],
    ])
    A = BanffClusterAlgebra(B, scalars=QQ)
    x1, x2, x3 = A.gens()
    _, P2, _, P4 = A.extra_primes()
    P2 = A.divisor({P2: 1})
    P4 = A.divisor({P4: 1})
    g = A((x1*x2*x3)**2 + ((x1-x3)**2)*(x3+1))
    f = A((x1*x2*x3)**2 + (((x2+1)**2)/(x1*x3))*(x3+1))
    Q1 = g.divisor()
    Q2 = f.divisor()
    div = P2 + P2 + Q1 + Q2 + P4 + P4 
    h = A.principal_generator(div)
    

    

    CF = h.factor()

    print_invariants("Example 2 (Cluster Algebra of type A3 su Q)", CF)
    print(CF.set_of_lengths())


# ----------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------

if __name__ == "__main__":
    import time
    t0 = time.time()

    run_example_1()
    run_example_2()

    t1 = time.time()
    print("\n" + "="*70)
    print(f"Total execution time: {t1 - t0:.2f} seconds")



