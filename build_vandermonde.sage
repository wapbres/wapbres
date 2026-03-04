from itertools import combinations

load("build_krawtchouk.sage")

def hw(n):
    return bin(n).count("1")

"""
Build the `n`-order Vandermonde matrix as described in Theorem 2.
"""
def build_vandermonde_matrix(n):
    V = matrix(ZZ, n + 1, n + 1)
    for k in range(n + 1):
        for l in range(n + 1):
            V[l, k] = k^l
    return V

"""
Helper functions to compute the quantities described as alpha, beta and gamma in
Theorem 2.
"""
def a_k(k, n):
    return 0 if k < 0 else -k / 2

def b_k(k, n):
    return 0 if k < 0 else n / 2

def g_k(k, n):
    return 0 if k < 0 else -(n - k) / 2

"""
Build the `r`-order matrix H as described in Theorem 2.
"""
def Hr(r, n):
    assert r <= n, "Invalid argument: r > n"

    H = matrix(QQ, r, r + 1)
    for i in range(r):
        g_i = g_k(i - 1, n)
        b_i = b_k(i + 0, n)
        a_i = a_k(i + 1, n)

        H[i, i - 1] = g_i
        H[i, i + 0] = b_i
        H[i, i + 1] = a_i
    return H

def all_k_supported_by_n(n):
    t = hw(n)
    bits_set = []
    i = 0
    while n >> i > 0:
        if n >> i & 1:
            bits_set.append(i)
        i += 1

    all_k = []
    for i in range(t + 1):
        for bits in combinations(bits_set, i):
            k = 0
            for b in bits:
                k += 2^b
            all_k.append(k)
    return sorted(all_k, reverse=False)

def build_supp_vandermonde_matrix(n, nb_cols):
    t = hw(n)
    K = all_k_supported_by_n(n)
    V = matrix(ZZ, 2^t, nb_cols)
    # Build the matrix columns wise.
    for j in range(V.ncols()):
        for i in range(V.nrows()):
            V[i, j] = K[i]^j
    return V

"""
Verify the relation between the Krawtchouk matrix and the Vandermonde matrix
from theorem 2 for small spaces.
"""
def theorem2():
    max_n = 31
    for n in range(1, max_n):
        M = matrix(1, n + 1)
        M[0, 0] = 1
        M_row = 1
        for r in range(1, n + 1):
            M_row *= Hr(r, n)
            M = M.stack(M_row.augment(matrix(1, n - r)))

            K = build_krawtchouk_matrix(n)
            V = build_vandermonde_matrix(n)

        assert M * K == V
