#

Sagemath implementation of Theorem 2 and Conjecture 1 of the paper [On the Resilience Order of W(A)PB Functions]().
It exposes the algebraic relation between Krawtchouk matrices and Vandermonde matrices.

The code in `t_corrections_conjecture.py` verifies Conjecture 1, incrementally for all positive $n$. It has been run until $n \leq 62$, since the Hamming weight of $63$ incudes a search space of approximately $2^{64}$ balanced vectors. As a guarantee of comparison, each $n$ with Hamming weight $5$ required approximately 2h30 of computations on an Intel-i5. In `t_corrections_conjecture.out` is  the output of the call to
```py
compute_t_corrections(62)
```
