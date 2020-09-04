import numpy as np
import math
from whole import *

# For a traceless 2-local hamiltonian on qutrits returns bloch vector  upon input of the solution vectors, C and c.
def trund(v, c, C):
    N=C.size[0]
    n=N/8
    T = c*math.sqrt(math.log(n))
    rd = {}
    y = []
    r = np.random.normal(0,1,N)
    for i in range(0,N):
        z = np.inner(r, v[i])/T
        if np.linalg.norm(z) > 1/2*math.sqrt(8): y.append(np.sign(z)/2*math.sqrt(8))
        else: y.append(z)
    rd["blochvec"] = y
    return rd

# Builds the C for a simple closed chain of the form H_i = gm[1]_igm[1]_(i+1) where gm[i] is the ith Gell-Mann matrix.
def tC(n):
    a = []
    a.extend([1,0,0,0,0,0,0,0] for i in range(n-1))
    C = np.diagflat(a,8) + np.diagflat(a,-8)
    C[0][8*n-8]=1
    C[8*n-8][0]=1
    return C

# Returns the SDP solution for aforementioned C.
def tsdp(n):
    A = np.eye(8*n)
    for i in range(0,n):
        for j in range(0,n):
            A[i*8][j*8] = 1
    return A

# Computes the the ratio upon C, y (output of 'trund'), and the maximal eigenvalue of the model.
def tratio(C, y, maxeig):
    return np.inner(np.dot(y,C),np.transpose(y))/maxeig

# Outputs a list of o ratios for a given qutrit number, c, o and C.
def tsample(n, c, o, C):
    v = findGenVecGram(tsdp(n))
    ratlist=[]
    i=0
    while i < o:
        ratlist.append(tratio(C, trund(v,c,cvx.matrix(C))["blochvec"], n))
        i=i+1
    return ratlist

# Plots the average of aforementioned samplelists within a range of qutrit numbers, upon input of c, sample size o per step, and total steps s. The steps are made such that they are equidistant on a logarithmic scale.
def tavgplot(ni,nf,c,o,s):
    y=[]
    x=[]
    a = np.exp(np.linspace(np.log(ni), np.log(nf), s))
    for i in a:
        x.append(i)
    for i in range(len(x)):
        y.append(avg(tsample(x[i],c,o,tC(x[i]))))
        print(x[i])
    plt.scatter(x,y,marker="o",s=5)
    plt.xlabel('log(n)')
    plt.ylabel('ratio')
    plt.xscale('log')
    plt.title('Qutrit Approximation Algorithm')
    plt.show()
