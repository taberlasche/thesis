"""
TODO:
+ [x] eigenvalues
+ [ ] checking positive semidefiniteness
+ [ ] partial traces
+ [x] tensor product
"""
import numpy as np

# for eigenvalues we can use eigvalsh for hermitian matrices which doesnt compute the eigenvectors and should therefore be more efficient np.linalg.eigvalsh()

# pos. def. For hermitian matrices the cholesky decomposition has cubic complexity. Implementing a check for hermicity using np.array_equal(A, A.H) and else compute eigenvalues (which isnt needed for bloch states but for completeness)

def pos_def(A)





    # np.all(np.linalg.eigvalsh(A) > 0)
    # chol_A = np.linalg.cholesky(A)

# for the tensor product and partial traces qutip could be used
