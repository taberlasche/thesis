"""
TODO:
+ [x] eigenvalues
+ [ ] checking positive semidefiniteness
+ [ ] partial traces
+ [x] tensor product
+ [x] maximal eigenvalue
"""
import numpy as np
from math import *
# for eigenvalues we can use eigvalsh for hermitian matrices which doesnt compute the eigenvectors and should therefore be more efficient np.linalg.eigvalsh()

# pos. def. For hermitian matrices the cholesky decomposition has cubic complexity. Implementing a check for hermicity using np.array_equal(A, A.H) and else compute eigenvalues (which isnt needed for bloch states but for completeness)

def pos_sdef(A):
     return np.all(np.linalg.eigvalsh(A) >= 0)
    # chol_A = np.linalg.cholesky(A)

# for the tensor product qutip could be used
# for partial traces np.reshape can be used alternatively to qutip

# the gell-mann matrices
gm = {
        1: np.array([[0,1,0],[1,0,0],[0,0,0]]),
        2: np.array([[0,-1j,0],[1j,0,0],[0,0,0]]),
        3: np.array([[1,0,0],[0,-1,0],[0,0,0]]),
        4: np.array([[0,0,1],[0,0,0],[1,0,0]]),
        5: np.array([[0,0,-1j],[0,0,0],[1j,0,0]]),
        6: np.array([[0,0,0],[0,0,1],[0,1,0]]),
        7: np.array([[0,0,0],[0,0,-1j],[0,1j,0]]),
        8: 1/sqrt(3)*np.array([[1,0,0],[0,1,0],[0,0,-2]])
        }
