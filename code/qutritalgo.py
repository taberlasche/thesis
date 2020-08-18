import numpy as np
import math
# for a 2-local hamiltonian on qutrits returns bloch vector  upon input of C and c
def trund(v, c, C):
    N=C.size[0]
    n=N/8
    T = c*math.sqrt(math.log(n)) #c=O(1)
    rd = {}
    y = []
    r = np.random.normal(0,1,N) #vector of 8n i.d.d. N(0,1) random variables
    for i in range(0,N):
        z = np.inner(r, v[i])/T
        if np.linalg.norm(z) > 1/2*math.sqrt(8): y.append(np.sign(z)/2*math.sqrt(8))
        else: y.append(z)
    rd["blochvec"] = y
    return rd
