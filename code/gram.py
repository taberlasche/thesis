import numpy as np


def findGenVecGram(gram):
    """Finds upon input of the Gram matrix gram a collection of vectors
    {v_1,...,v_n} such that gram_{ij} = <v_i, v_j>."""

    # Calculate the Decomposition such that gram = L*L^T
    L = choleskyPSDSym(gram)

    # Take the columns of L and write them into a list
    result = []
    for i in range(len(gram)):
        result += [L[:,i]]
    return result



def calcGram(v):
    """Returns upon input of the collection of vectors v = {v_1, ..., v_n}
    the associated Gram matrix gram, i.e., the matrix gram such that
    gram_{ij} = <v_i, v_j>."""
    n = len(v)
    result = np.zeros((n, n), dtype=float)
    for i in range(n):
        for j in range(n):
            result[i][j] = np.inner(v[i], v[j])
    return result



def choleskyPSDSym(A):
    """Returns a matrix L such that L*L^T = A, where A is a positive
    semidefinite symmetric matrix"""

    # Calculate spectral decomposition (eigenvalues + eigenvectors) of A
    w, v = np.linalg.eig(A)

    # Return the orthogonal matrix v multiplied with the square root of
    # the diagonal matrix whose diagonal consists of the eigenvalues
    return np.transpose(np.matmul(v, np.sqrt(np.diag(w))))




main()
