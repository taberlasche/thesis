# TODO:
# + Merge two loops

import numpy as np

def round(v, c):

    #Pauli matrices
    X = np.array([[0, 1],[1, 0]])
    Y = np.array([[0, -1j],[1j, 0]])
    Z = np.array([[1, 0],[0, -1]])
    I = np.array([1, 0],[0, 1])
    #n should still be saved at this point from the previous function
    T = c*sqrt(log(n)) #c=O(1)
    r = np.random.random(3n) #vector of 3n i.d.d. N(0,1) random variables

    y = []
    rho = []
    for i in range(1,3n+1):
        z = np.inner(r, v[i])/T
            if np.linalg.norm(z, 2) > 1/sqrt(3): y.append(np.sign(z)/sqrt(3))
            else: y.append(z)
    for a in range(1, n+1):
        rho.append(0.5*(I + y[3a-2]*X + y[3a-1]*Y + y[3a]*Z))
    return rho
