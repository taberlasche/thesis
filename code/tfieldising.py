import numpy as np
# implements the nested kronecker product
def nested_kronecker_product(a):
    if len(a) == 2:
        return np.kron(a[0],a[1])
    else:
        return np.kron(a[0], nested_kronecker_product(a[1:]))

# returns the jordan wigner transform
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

#builds the transverse field ising model hamiltonian, takes a list of jw trafos
def hamiltonian(gam, lam, a, lattice_length):
    H = 0
    for i in range(lattice_length - 1):
        H += a[i].T.dot(a[i+1]) - a[i].dot(a[i+1].T)
        H -= gam*(a[i].T.dot(a[i+1].T) - a[i].dot(a[i+1]))
    for i in range(lattice_length):
        H -= 2*lam*(a[i].dot(a[i].T))
    return H
# returns the maximal eigenvalue of such a chain with a certain length upon input of gamma and lambda
def tfisingmaxeig(gam, lam, lattice_length):
    a = []
    for i in range(lattice_length):
        a.append(jordan_wigner_transform(i, lattice_length))
    return np.amax(np.linalg.eig(hamiltonian(gam, lam, a, lattice_length)
)[0])
