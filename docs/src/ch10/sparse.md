### Sparse matrix tricks

If direct matrix-matrix multiplications are used when evaluating the local trace, the $$F$$-matrix must be very sparse. Thus, we can convert them into sparse matrices in compressed sparse row (CSR) format, and then the sparse matrix multiplication can be applied to obtain a significant speedup.