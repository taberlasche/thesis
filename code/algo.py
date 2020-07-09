# TODO:
# + Merge two loops

import numpy as np
import math
# returns bloch vector and state upon input of C and c
def round(v, c):

    #Pauli matrices
    X = np.array([[0, 1],[1, 0]])
    Y = np.array([[0, -1j],[1j, 0]])
    Z = np.array([[1, 0],[0, -1]])
    I = np.array([[1, 0],[0, 1]])
    N=len(v)
    print(v)
    n=N/3
    T = c*math.sqrt(math.log(n)) #c=O(1)
    r = np.random.random(int(3*n)) #vector of 3n i.d.d. N(0,1) random variables

    y = []
    rho = []
    for i in range(0,int(n)):
        z = np.inner(r, v[i])/T
        if np.linalg.norm(z) > 1/math.sqrt(3): y.append(np.sign(z)/math.sqrt(3))
        else: y.append(z)
    for a in range(0, 3*int(n)):
        rho.append(0.5*(I + y[3*a-2]*X + y[3*a-1]*Y + y[3*a]*Z))
    return y
    return rho
