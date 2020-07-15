import numpy as np
import picos as pic
import cvxopt as cvx
import cvxopt.lapack

from sdp import *
from gram import *
from algo import *
from operations import *
from hamiltonians import *

# returns 1. the solution to the SDP 2. the gram vectors 3. dictionary with the state and bloch vector upon input of C and c(O(1))
def solve(C, c):
    V = np.array(mc(C))
    K = (np.conj(V)+V)/2
    v = findGenVecGram(K)
    sol = round(v,c,C)
    return sol

# returns upon input of C, y and the maximal eigenvalue corresponding to C the ratio
def ratio(C, y, maxeig):
    return np.inner(np.dot(y,C),np.transpose(y))/maxeig

# returns list of ratios upon input of number of qubits, c and number of interations
def sample(n, c, o):
    C = buildC(n)
    K = np.add(C, np.identity(3*n))
    v = findGenVecGram(K)
    ratlist=[]
    #for i in range(0,o+1):
    i=0
    while i < o+1:
        ratlist.append(ratio(C, round(v,c,cvx.matrix(C)), n))
        i=+1
    return ratlist

# returns a log plot
#def plotn(ni, nf, c, o)



# returns approximation ratios upon input of C, c, o, maxeig, where o is the number of rounds and maxeig the maximum eigenvalue of the Hamiltonian corresponding to C
# def approx(C, c, o, maxeig):
#    ratlist=[]
#    i=1
#    while i < o+1:
#        sol = solve(C, c)
#        y=sol["blochvec"]
#        ratlist.append(ratio(C,y,maxeig)
#    return ratlist

ratio(buildC(4),round(findGenVecGram(np.add(buildC(4),np.identity(12))), 2, cvx.matrix(buildC(4))), 4)
