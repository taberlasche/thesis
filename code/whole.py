import numpy as np
import picos as pic
import cvxopt as cvx
import cvxopt.lapack

from sdp import *
from gram import *
from algo import *
from operations import *


# returns upon input of H and rho a ratio
def ratio(C, y, maxeig):
    k = np.inner(np.inner(np.transpose(y),C),y)
    rat = maxeig/k

# returns rho upon input of C and c(O(1))
def solve(C, c):
    mc(C)
    V = mc(C)
    K = (np.conj(V)+V)/2
    v = findGenVecGram(K)
    round(v, c)


# returns approximation ratios upon input of C, c, o, maxeig, where o is the number of rounds and maxeig the maximum eigenvalue of the Hamiltonian corresponding to C
def approx(C, c, o, maxeig):
    for i in range(1, o+1):
        solve(C, c)
        ratio(C, y)
        #print(i ":" rat)
        ratlist=[]
        ratlist.append(rat)
