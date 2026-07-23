"""
Performance comparison for class-group computations in Banff cluster algebras.

This script reproduces the benchmark family used in the paper of Pompili and
Smertnig.  For exchange matrices B(c, p) coming from triangulated punctured
discs, it compares two approaches to computing the rank of the divisor class
group:

1. the finite Laurent intersection ring (FLIR) method implemented by
   ``BanffClusterAlgebra.class_group()``;
2. the presentation-based method, which computes generators, constructs a
   finite presentation by elimination, and then determines the relevant
   height-one primes using associated-prime computations.

Each expensive step is run in a separate process with a fixed timeout.  The
script records the number of variables, generators and Laurent charts, checks
that both methods return the same rank when both terminate, and prints the
results as a LaTeX table.

Run this file with SageMath, for example::

    sage performance_comparison.py

The running times are machine-dependent.  Algebraic outputs such as ranks,
numbers of generators and numbers of charts should be reproducible for the
same SageMath implementation and input matrices.
"""

import time
from multiprocessing import Process, Queue
from queue import Empty  


from sage.algebras.banff_cluster_algebra import BanffClusterAlgebra
from sage.matrix.constructor import matrix
from sage.misc.misc_c import prod
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.rational_field import QQ



# ========================================================
# Disc Matrices
# ========================================================

class DiscMatrix:
    """
    Construct the exchange matrix ``B(c, p)`` used in the benchmark.

    The matrix comes from a triangulation of a disc with

    - ``c`` marked points on the boundary, and
    - ``p`` punctures in the interior.

    The construction is recursive.  For ``p = 1`` it starts from ``B(2, 1)``
    and repeatedly adds a boundary marked point.  For ``p >= 2`` it starts
    from ``B(2, 2)``, first adds boundary marked points, and then adds
    punctures.  Adding a boundary marked point increases the matrix size by
    one, while adding a puncture increases it by three.

    The resulting matrix is skew-symmetric and has size

    .. MATH::

        3p + c - 3.

    """
    def __init__(self, c, p):
        if c < 2 or p < 1:
            raise ValueError("Must have c >= 2 and p >= 1.")
        self.c = c
        self.p = p
        self.B = self._build_matrix(c, p)

    def _base_matrix(self, c, p):
        """
        Return one of the initial matrices used by the recursive construction.

        Base cases are provided for ``(c, p) = (2, 1)``, ``(2, 2)`` and
        ``(3, 1)``.  The main constructor currently builds all matrices from
        ``B(2, 1)`` or ``B(2, 2)``.
        """
        if c == 2 and p == 1:
            return matrix(2, 2, [[0, 0], [0, 0]])
        elif c == 2 and p == 2:
            return matrix([
                [ 0,  -1,  1, -1, 0],
                [ 1,  0, -1,  0,  1],
                [-1,  1,  0,  1, -1],
                [ 1,  0, -1,  0,  1],
                [ 0, -1,  1,  -1,  0]
            ])
        elif c == 3 and p == 1:
            return matrix([
                [ 0,  1,  -1],
                [ -1,  0, 1],
                [1,  -1,  0]
            ])
        else:
            raise ValueError("Base matrix defined only for (c,p) = (2,1) or (2,2) or (3,1).")

    def _extend_c(self, B):
        """
        Add one marked point on the boundary.

        If ``B`` is the current exchange matrix, return the matrix for the
        next value of ``c``.  The construction adds one row and one column,
        fills the local arrows near the new boundary vertex, restores
        skew-symmetry, and finally permutes the last two indices so that the
        indexing agrees with the recursive convention used in the examples.
        """
        n = B.nrows()
        B_new = matrix(n + 1)

        for i in range(n):
            for j in range(i + 1, n):
                B_new[i, j] = B[i, j]

        # new entries
        B_new[n - 2, n - 1] = 0
        B_new[n - 2, n]     = -1
        B_new[n - 1, n]     = 1

        # antisymmetrize
        for i in range(n + 1):
            for j in range(i + 1, n + 1):
                B_new[j, i] = -B_new[i, j]

        # Build permutation matrix that swaps (n-1) and n

        P = matrix.identity(n + 1)
        P[n-1, n-1], P[n, n] = 0, 0
        P[n-1,   n], P[n, n-1] = 1, 1

        # Conjugate: P * B_new * P
        B_new = P * B_new * P

        return B_new


    def _extend_p(self, B):
        """
        Add one puncture to the triangulated disc.

        This operation adds three rows and three columns.  The old matrix is
        embedded into the enlarged matrix and the final entries describe the
        local quiver configuration introduced by the new puncture.
        """
        n = B.nrows()
        B_new = matrix(n + 3)
        for i in range(n + 3):
            for j in range(i + 1, n + 3):
                if i <= n - 3 and j in {n, n + 1, n + 2}:
                    B_new[i, j] = 0 
                if i <= n and j <= n - 1:
                    B_new[i, j] = B[i, j]
        
        B_new[n - 2, n - 1] = 0
        B_new[n -2, n]      = -1
        B_new[n -2, n + 1]  = 1
        B_new[n -1, n]      = 1
        B_new[n - 1, n + 2] =-1
        B_new[n, n + 1]     = -1
        B_new[n, n + 2]     = 1
        B_new[n + 1, n + 2] = -1

        for i in range(n + 3):
            for j in range(i + 1, n + 3):
                B_new[j, i] = -B_new[i, j]
        
        return B_new

    def _build_matrix(self, c, p):
        """
        Build ``B(c, p)`` from the appropriate base matrix.

        Boundary extensions are performed before puncture extensions, in the
        same order as in the recursive definition of the benchmark family.
        """
        if p == 1:
            B = self._base_matrix(2, 1)
            for _ in range(c - 2):
                B = self._extend_c(B)
        else:
            B = self._base_matrix(2, 2)
            for _ in range(c - 2):
                B = self._extend_c(B)
            for _ in range(p - 2):
                B = self._extend_p(B)
        return B

    def matrix(self):
        return self.B

    def size(self):
        return self.B.nrows()


# ========================================================
# Compute associated primes
# ========================================================

def compute_associated_primes(T, R, I):
    primes = set()
    xvars = [g for g in R.gens() if str(g).startswith("X")]
    for x in xvars:
        J = R.ideal(x) + I
        primes_in_R = J.associated_primes()
        for P in primes_in_R:
            if all(g in P for g in I.gens()):
                gen_in_T = [T(g) for g in P.gens()]
                P_in_T = T.ideal(gen_in_T)
                primes.add(P_in_T)
    return primes

# ============================================================
# Compute presentation
# ============================================================

def presentation_with_given_generators(A, generators):

        L = A.ambient()
        scalars = L.base_ring()
        n = L.ngens()

        yvars = [ f"Y{i}" for i in range(n) ]
        xvars = [ f"X{i}" for i in range(n) ]

        R = PolynomialRing(scalars, xvars, order='degrevlex')
        I = R.ideal(0)

        numgens = len(generators)

        for j in range(numgens):
            f = generators[j]

            curgens  = list(R.gens())
            currels = I.gens()

            curgens = yvars + [f"T{j}"] + curgens

            S = PolynomialRing(scalars, curgens)

            Xsub = dict([ (L(f"x{i}"), S(f"X{i}")) for i in range(n) ])
            Ysub = dict([ (L(f"x{i}"), S(f"Y{i}")) for i in range(n) ])

            
            xyrels = [ S(f"X{i}*Y{i} - 1") for i in range (n) ]
            currels = [ S(r) for r in currels ]
            currels += xyrels
                
            # construct the denominator exponents of f
            dexp = [0]*(n+1)
            for e in f.exponents():
                for i in range(n):
                    if e[i] < 0:
                        dexp[i] = max(-e[i], dexp[i])

            d = prod([ L.gen(i) ** dexp[i] for i in range(n) ])
            g = f*d
            G=g.substitute(Xsub)
            D=d.substitute(Ysub)

            currels.append(G*D - S(f"T{j}"))
            
            Iprime = S.ideal(currels)
            
            J = Iprime.elimination_ideal([ S(f"Y{i}") for i in range(n) ])  

            R = PolynomialRing(scalars, list(R.gens()) + [f"T{j}"])
            I = R.ideal(J.gens())

        return R.quo(I), R, I

# ============================================================
# Utility: run function with timeout
# ============================================================


def run_with_timeout(func, args, timeout_seconds):
    q = Queue()
    p = Process(target=func, args=(*args, q))

    start = time.perf_counter()
    p.start()

    try:
        res = q.get(timeout=timeout_seconds)
        p.join()
        end = time.perf_counter()
        return res, end - start

    except Empty:
        if p.is_alive():
            p.terminate()
            p.join()
        end = time.perf_counter()
        return ("timeout",), end - start

    except Exception as e:
        if p.is_alive():
            p.terminate()
            p.join()
        end = time.perf_counter()
        return ("error", f"run_with_timeout: {type(e).__name__}: {e}"), end - start



# ============================================================
# Target functions 
# ============================================================

def rank_using_FLIRs_target(matrix, result_queue):
    try:
        A = BanffClusterAlgebra(matrix.B, scalars=QQ)
        G = A.class_group()
        rank = G.rank()
        number_charts = len(A.charts)
        result_queue.put(("success", rank, number_charts))
    except Exception as e:
        result_queue.put(("error", str(e)))


def banff_algorithm_target(A, result_queue):
    try:
        generators = A.generators()
        numgens = len(generators)
        result_queue.put(("success", generators, numgens))
    except Exception as e:
        result_queue.put(("error", str(e)))


def presentation_target(A, generators, result_queue):
    try:
        T, R, I = presentation_with_given_generators(A, generators)
        result_queue.put(("success", T, R, I))
    except Exception as e:
        result_queue.put(("error", str(e)))


def primary_decomposition_target(T, R, I, result_queue):
    try:
        primes = compute_associated_primes(T, R, I)
        result_queue.put(("success", primes))
    except Exception as e:
        result_queue.put(("error", str(e)))


# ============================================================
# Main comparison function
# ============================================================

def run_comparison(m, p, timeout_seconds=600):
    B_obj = DiscMatrix(m, p)

    assert B_obj.B.rank() != B_obj.B.nrows()

    result = {
        "status_FLIRs": None,
        "status_generators": None,
        "status_presentation": None,
        "status_primary_decomposition": None,
        "values": {},
        "times": {},
        "errors": {}
    }


    print(f"Computing FLIRs for (m,p)=({m},{p})")

    # ---------------- FLIRs ----------------
    res, t = run_with_timeout(
        rank_using_FLIRs_target,
        (B_obj,),
        timeout_seconds
    )
    result["times"]["FLIRs"] = t

    if res[0] == "timeout":
        result["status_FLIRs"] = "timeout"
        result["errors"]["FLIRs"] = "timeout"
        return result
    elif res[0] == "error":
        result["status_FLIRs"] = "error"
        result["errors"]["FLIRs"] = res[1] if len(res) > 1 else "unknown error"
        return result
    elif res[0] != "success":
        result["status_FLIRs"] = "error"
        result["errors"]["FLIRs"] = f"unexpected response: {res}"
        return result

    result["status_FLIRs"] = "success"
    result["values"]["rank"] = res[1]
    result["values"]["number_charts"] = res[2]

    # ---------------- Generators ----------------
    print(f"Computing generators for (m,p)=({m},{p})")

    res, t_gens = run_with_timeout(
        banff_algorithm_target,
        (BanffClusterAlgebra(B_obj.B, scalars = QQ),),
        timeout_seconds
    )

    if res[0] == "timeout":
        result["status_generators"] = "timeout"
        result["errors"]["generators"] = "timeout"
        return result
    elif res[0] == "error":
        result["status_generators"] = "error"
        result["errors"]["generators"] = res[1] if len(res) > 1 else "unknown error"
        return result
    elif res[0] != "success":
        result["status_generators"] = "error"
        result["errors"]["generators"] = f"unexpected response: {res}"
        return result


    result["status_generators"] = "success"
    generators = res[1]
    result["values"]["number_generators"] = res[2]

    # ---------------- Presentation ----------------
    print(f"Computing presentation for (m,p)=({m},{p})")
    res, t_pres = run_with_timeout(
        presentation_target,
        (BanffClusterAlgebra(B_obj.B, scalars = QQ), generators,),
        timeout_seconds
    )
    result["times"]["presentation"] = t_gens + t_pres


    if res[0] == "timeout":
        result["status_presentation"] = "timeout"
        result["errors"]["presentation"] = "timeout"
        return result
    elif res[0] == "error":
        result["status_presentation"] = "error"
        result["errors"]["presentation"] = res[1] if len(res) > 1 else "unknown error"
        return result
    elif res[0] != "success":
        result["status_presentation"] = "error"
        result["errors"]["presentation"] = f"unexpected response: {res}"
        return result


    result["status_presentation"] = "success"
    T, R, I = res[1], res[2], res[3]

    # ---------------- Primary decomposition ----------------
    print(f"Computing primary decomposition for (m,p)=({m},{p})")

    t = result["times"]["presentation"]

    res, t_pri = run_with_timeout(
        primary_decomposition_target,
        (T, R, I),
        timeout_seconds
    )

    if res[0] == "timeout":
        result["status_primary_decomposition"] = "timeout"
        result["errors"]["primary_decomposition"] = "timeout"
        return result
    elif res[0] == "error":
        result["status_primary_decomposition"] = "error"
        result["errors"]["primary_decomposition"] = res[1] if len(res) > 1 else "unknown error"
        return result
    elif res[0] != "success":
        result["status_primary_decomposition"] = "error"
        result["errors"]["primary_decomposition"] = f"unexpected response: {res}"
        return result
    t = t + t_pri

    #If we reach here, all primary decomposition of x1,...,xn succeeded
    result["status_primary_decomposition"] = "success"
    result["values"]["rank_presentation"] = len(res[1]) - 3*p-m + 3
    result["times"]["rank_presentation"] = t

    # Consistency check
    if result["values"]["rank"] != result["values"]["rank_presentation"]:
        raise ValueError(
            f"Rank mismatch: FLIR={result['values']['rank']}, "
            f"pres={result['values']['rank_presentation']}"
        )

    return result

def status_mark(status):
    if status == "success":
        return None
    if status == "timeout":
        return "*"
    if status == "error":
        return "ERR"
    return "?"



# ============================================================
# LaTeX table output
# ============================================================

def format_time_nicely(t):
    hours = int(t) // 3600
    minutes = (int(t) - hours * 3600) // 60
    seconds = t - hours * 3600 - minutes * 60

    if hours > 0:
        return f"{hours}h{minutes}m"
    elif minutes > 0:
        return f"{minutes}m{seconds:.0f}s"
    else:
        return f"{seconds:.1f}s"

def print_latex_table(results):

    # HEADER
    header = " "
    for (m, p, data) in results:
        header += f" & $B({m},{p})$"
    header += " \\\\"

    print("\\begin{tabular}{l" + " c"*len(results) + "}")
    print("\\toprule")
    print(header)
    print("\\midrule")

    # RANK
    row = "rank"
    for (_, _, data) in results:
        if data["status_FLIRs"] != "success":
            row += " & --"
        else:
            row += " & " + str(data['values'].get('rank', '--'))
    print(row + "\\\\")

    # VARIABLES
    row = "vars"
    for (m, p, _) in results:
            row += f" & {3*p+m-3}"
    print(row + "\\\\")

    # GENERATORS
    row = "gens"
    for (_, _, data) in results:
        if data["status_generators"] != "success":
            row += " & --"
        else:
            row += " & " + str(data['values'].get('number_generators', '--'))
    print(row + "\\\\")

    # CHARTS
    row = "charts"
    for (_, _, data) in results:
        if data["status_FLIRs"] != "success":
            row += " & --"
        else:
            row += " & " + str(data['values'].get('number_charts', '--'))
    print(row + "\\\\")

    # TIME FLIRs
    row = "\\textsf{L}"
    for (_, _, data) in results:
        if data["status_FLIRs"] != "success":
            row += " & $" + status_mark(data["status_FLIRs"]) + "$"
        else:
            t = data["times"].get("FLIRs")
            row += " & " + f"{format_time_nicely(t)}"
    print(row + "\\\\")



    # TIME presentation
    row = "\\textsf{P}"
    for (_, _, data) in results:
        if data["status_generators"] != "success":
            row += " & --"
        elif data["status_presentation"] != "success":
            row += " & $" + status_mark(data["status_presentation"]) + "$"
        else:
            t = data["times"].get("presentation")
            row += " & " + f"{format_time_nicely(t)}"
    print(row + "\\\\")

    # TIME primary decomposition
    row = "\\textsf{PD}"
    for (_, _, data) in results:
        if  data["status_generators"] != "success" or data["status_presentation"] != "success":
            row += " & --"
        elif data["status_primary_decomposition"] == "timeout":
            row += " & $*$"
        elif data["status_primary_decomposition"] == "error":
            row += " & ERR"
        else:
            t = data["times"].get("rank_presentation")
            row += " & " + f"{format_time_nicely(t)}"
    print(row + "\\\\")


    print("\\bottomrule")
    print("\\end{tabular}")


# ============================================================
# Main
# ============================================================

if __name__ == "__main__":
    MAX_TIME = 3600  # 1 hour

    cases = [
        (2,2),
        (2, 3), (3, 3), 
        (4, 3), (3, 4), 
        (4, 4), 
    ]

    results = []
    for m, p in cases:
        results.append((m, p, run_comparison(m, p, MAX_TIME)))

    print_latex_table(results)
