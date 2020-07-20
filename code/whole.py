import numpy as np
import picos as pic
import cvxopt as cvx
import cvxopt.lapack
import matplotlib.pyplot as plt

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
    print(sol)
    return sol

# returns upon input of C, y and the maximal eigenvalue corresponding to C the ratio
def ratio(C, y, maxeig):
    return np.inner(np.dot(y,C),np.transpose(y))/maxeig

# returns list of ratios upon input of number of qubits, c and number of iterations
def sample(n, c, o, C):
    K = np.add(C, np.identity(3*n))
    v = findGenVecGram(K)
    ratlist=[]
    i=0
    while i < o:
        ratlist.append(ratio(C, round(v,c,cvx.matrix(C))["blochvec"], n))
        i=i+1
    return ratlist

# returns the highest eigenvalue of the n-qubit productstate after performing the rounding o times
def sample2(n,c,o,C):
    K = np.add(C, np.identity(3*n))
    v = findGenVecGram(K)
    eiglist=[]
    for i in range(o):
        y=round(v,c,cvx.matrix(C))["blochvec"]
        eiglist.append(np.inner(np.dot(y,C),np.transpose(y)))
    return eiglist



# returns a plot of samples in a range of qubitnumbers upon input of start and end of range, c and the number of rounding attemps per qubitnumber
def plotn(ni, nf, c, o):
    C = buildC(n)
    x=[]
    y=[]
    for i in range(ni,nf+1,4):
        x.extend([i for j in range(o)])
    for i in range(ni,nf+1,4):
        y=y+sample(i,c,o,C)
    plt.scatter(x,y,marker=".",s=2)
    plt.xlabel('n')
    plt.ylabel('log(ratio)')
    plt.yscale('log')
    plt.title('nplot')
    plt.show()

# returns a plot of samples for a range of c's to hopefully find out something about what the optimal c is
def plotc(ci, cf, n, o):
    C = buildC(n)
    x=[]
    y=[]
    for i in range(ci,cf+1):
        x.extend([i for j in range(o)])
    for i in range(ci,cf+1):
        y=y+sample(n,i,o,C)
    plt.scatter(x,y,marker=".",s=2)
    plt.xlabel('c')
    plt.ylabel('ratio')
    plt.title('cplot')
    plt.show()
