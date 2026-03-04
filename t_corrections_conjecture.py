from sage.all import *
import time

# To build Krawtchouk and Vandermonde matrices.
load("build_krawtchouk.sage")
load("build_vandermonde.sage")

"""
Return an iterator on all `n`-bits integers with hamming weight `k`.
"""
def kbits(n, k):
    limit = 1 << n
    val = (1 << k) - 1
    while val < limit:
        yield val
        minbit = val & -val
        fillbit = (val + minbit) & ~val  #rightmost 0 to the left of that bit
        val = val + minbit | (fillbit // (minbit << 1)) - 1

"""
Exhaustively try all combinations of vectors v with coefficients in {±1, 0} and
count how many columns of the n-th order krawtchouk matrix they cancel. The
maximum number of columns such that the dot product is null is the t-correction.

An actual output of this function ran up to n = 62 with execution time is
available in `t_corrections.out`.
"""
def compute_t_corrections(target):
    time_start = time.time()
    for n in range(1, target + 1):
        K = build_krawtchouk_matrix(n).T
        t = 0
        non0_offsets = []
        for k in range(n + 1):
            if binomial(n, k) % 2 == 0:
                t += 1
            else:
                non0_offsets.append(k)

        v = vector(ZZ, n + 1)
        max_T = 0

        for balanced_word in kbits(n + 1 - t, (n + 1 - t) // 2):
            for m, i in enumerate(non0_offsets):
                b = (balanced_word >> m) & 1
                if b == 0:
                    v[i] = 1
                else:
                    v[i] = -1

            T = 0
            for col in K.columns():
                if v * col == 0:
                    T += 1
                else:
                    break
            if T > max_T:
                max_T = T

        print(f"{n}: t-corrector {max_T - 1} (wt(n) - 1 = {vector(GF(2), bin(n)[2:]).hamming_weight() - 1})")
        print("Time elapsed:", time.time() - time_start)

"""
Same as `compute_t_corrections` but with the columns of the vandermonde matrix
as described in the paper.
"""
def compute_t_corrections_equivalence(target):
    time_start = time.time()
    for n in range(1, target + 1):
        V = build_vandermonde_matrix(n).T
        t = 0
        non0_offsets = []
        for k in range(n + 1):
            if binomial(n, k) % 2 == 0:
                t += 1
            else:
                non0_offsets.append(k)

        v = vector(ZZ, n + 1)
        max_T = 0

        for balanced_word in kbits(n + 1 - t, (n + 1 - t) // 2):
            for m, i in enumerate(non0_offsets):
                b = (balanced_word >> m) & 1
                if b == 0:
                    v[i] = 1
                else:
                    v[i] = -1

            T = 0
            for col in V.columns():
                if v * col == 0:
                    T += 1
                else:
                    break
            if T > max_T:
                max_T = T

        print(f"{n}: t-corrector {max_T - 1} (wt(n) - 1 = {vector(GF(2), bin(n)[2:]).hamming_weight() - 1})")
        print("Time elapsed:", time.time() - time_start)


compute_t_corrections(30)
compute_t_corrections_equivalence(30)
