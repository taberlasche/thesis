import numpy as np
import picos as pic
import cvxopt as cvx
import cvxopt.lapack

from sdp import *
from gram import *
from algo import *
from operations import *


# returns upon input of H and rho a ratio
def ratio(C, y):
    k = np.inner(np.inner(np.transpose(y),C),y)
    H = hamil(C)
    rat = maxeig(H)/k

# returns rho upon input of C and c(O(1))
def solve(C, c):
    mc(C)
    Xr = (X+np.conjugate(X))/2
    findGenVecGram(Xr)
    result = v
    round(v, c)

# returns approximation ratios upon input of C, c and o, where o is the count of approximation rounds
def approx(C, c, o):
    for i in range(1, o+1):
        solve(C, c)
        ratio(C, y)
        print(i ':' rat)
