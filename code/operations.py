"""
TODO:
+ Checking positive semidefiniteness
+ partial traces
+ eigenvalues
"""
import numpy as np

# pos. def. For symmetric matrices the cholesky decomposition has cubic complexity. Implementing a check for symmetry using np.array_equal and else compute eigenvalues

def pos_def(A)






    # np.all(np.linalg.eigvals(A) > 0)
    # chol_A = np.linalg.cholesky(A)
