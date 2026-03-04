"""
Build the n-th order Krawtchouk matrix.
"""
def build_krawtchouk_matrix(n):
    K = []
    for k in range(n + 1):
        row = []
        for l in range(n + 1):
            row.append(codes.bounds.krawtchouk(n, 2, k, l))
        K.append(row)
    return matrix(ZZ, K)

"""
Build the n-th order Krawtchouk matrix with rows only when n choose k is odd.
"""
def build_krawtchouk_matrix_odd_binomial(n):
    K = []
    for k in range(n + 1):
        row = []
        if binomial(n, k) % 2 == 0:
            continue
        for l in range(n + 1):
                row.append(codes.bounds.krawtchouk(n, 2, k, l))
        K.append(row)
    return matrix(ZZ, K)
