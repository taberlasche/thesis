"""
TODO:
+ checking positive semidefiniteness
+ partial traces
+ tensor product
"""
import numpy as np

# for eigenvalues we can use eigvalsh for hermitian matrices which doesnt compute the eigenvectors and should therefore be more efficient np.linalg.eigvalsh()

# pos. def. For symmetric matrices the cholesky decomposition has cubic complexity. Implementing a check for symmetry using np.array_equal(A, A.H) and else compute eigenvalues

def pos_def(A)






    # np.all(np.linalg.eigvalsh(A) > 0)
    # chol_A = np.linalg.cholesky(A)
