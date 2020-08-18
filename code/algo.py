import numpy as np
import math
# for a 2-local Hamiltonian on qubits returns bloch vector #and state# upon input of C and c
def rund(v, c, C):
    N=C.size[0]
    n=N/3
    T = c*math.sqrt(math.log(n)) #c=O(1)
    rd = {}
    y = []
    r = np.random.normal(0,1,N) #vector of 3n i.d.d. N(0,1) random variables
    for i in range(0,N):
        z = np.inner(r, v[i])/T
        if np.linalg.norm(z) > 1/math.sqrt(3): y.append(np.sign(z)/math.sqrt(3))
        else: y.append(z)
    rd["blochvec"] = y
    return rd
