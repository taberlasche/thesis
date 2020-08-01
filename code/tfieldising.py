# This code returns the eigenvalues of a spinchain of length lattice_length through a jordan wigner transformation
import numpy as np

def nested_kronecker_product(a):
    if len(a) == 2:
        return np.kron(a[0],a[1])
    else:
        return np.kron(a[0], nested_kronecker_product(a[1:]))

def jordan_wigner_transform(j, lattice_length):
    sigma = np.array([[0, 1], [0, 0]])
    sigma_z = np.array([[1, 0], [0, -1]])
    I = np.eye(2)
    operators = []
    for k in range(j):
        operators.append(sigma_z)
    operators.append(sigma)
    for k in range(lattice_length-j-1):
        operators.append(I)
    return -nested_kronecker_product(operators)

lattice_length = 12
a = []
for i in range(lattice_length):
    a.append(jordan_wigner_transform(i, lattice_length))

def hamiltonian(gam, lam, a, lattice_length):
    H = 0
    for i in range(lattice_length - 1):
        H += a[i].T.dot(a[i+1]) - a[i].dot(a[i+1].T)
        H -= gam*(a[i].T.dot(a[i+1].T) - a[i].dot(a[i+1]))
    for i in range(lattice_length):
        H -= 2*lam*(a[i].dot(a[i].T))
    return H

gam, lam =1, 1
H = hamiltonian(gam, lam, a, lattice_length)
eigenvalues = np.linalg.eig(H)[0]
sorted(eigenvalues)
np.amax(eigenvalues)
