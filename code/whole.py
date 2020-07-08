import numpy as np
import picos as pic
import cvxopt as cvx
import cvxopt.lapack

from sdp import *
from gram import *
from algo import *
from operations import *

def solve(C, c):
    mc(C)
    Xn = (X+np.conjugate(X))/2
    findGenVecGram(Xn)
    result = v
    round(v, c)
