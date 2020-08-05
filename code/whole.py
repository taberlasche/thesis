import numpy as np
import picos as pic
import cvxopt as cvx
import cvxopt.lapack
import matplotlib.pyplot as plt

from hamiltonians import *
from sdp import *
from gram import *
from algo import *
from operations import *

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
    v = findGenVecGram(sdp1(n))
    ratlist=[]
    i=0
    while i < o:
        ratlist.append(ratio(C, round(v,c,cvx.matrix(C))["blochvec"], n))
        i=i+1
    return ratlist

def chainsample(n, c, o, C):
    v = findGenVecGram(sdp2(n))
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

def avg(l):
    return sum(l)/len(l)


# for the buildC(n) returns a plot of samples in a range of qubitnumbers upon input of start and end of range, c and the number of rounding attemps per qubitnumber
def plotn(ni, nf, c, o, s):
    x=[]
    y=[]
    for i in range(ni,nf+1,s):
        x.extend([i for j in range(o)])
        y=y+sample(i,c,o,buildC(i))
    plt.scatter(x,y,marker=".",s=3)
    plt.xlabel('n')
    plt.ylabel('ratio')
    #plt.ylabel('log(ratio)')
    #plt.yscale('log')
    plt.title('nplot')
    plt.show()

# does the same as plotn but returns the averages over the lists per n, and the steps are choosable for efficiency
def avgplotn(ni,nf,c,o,s):
    x=[]
    y=[]
    for i in range(ni,nf+1,s):
        x.append(i)
        y.append(avg(sample(i,c,o,buildC(i))))
        print(i)
    plt.scatter(x,y,marker="o",s=5)
    plt.xlabel('log(n)')
    plt.ylabel('ratio')
    #plt.ylabel('log(ratio)')
    plt.xscale('log')
    plt.title('nplot')
    plt.show()

# for the buildC(n) returns a plot of samples for a range of c's to hopefully find out something about what the optimal c is
def plotc(ci, cf, n, o):
    C = buildC(n)
    x=[]
    y=[]
    for i in range(ci,cf+1):
        x.extend([i for j in range(o)])
        y=y+sample(n,i,o,C)
    plt.scatter(x,y,marker=".",s=3)
    plt.xlabel('c')
    plt.ylabel('ratio')
    plt.title('cplot')
    plt.show()

# for the chainC(n) returns a plot of samples in a range of qubitnumbers upon input of start and end of range, c and the number of rounding attemps per qubitnumber
def chainplotn(ni, nf, c, o,s):
    x=[]
    y=[]
    for i in range(ni,nf+1):
        x.append(i)
        y=y+avg(chainsample(i,c,o,chainC(i)))
    plt.scatter(x,y,marker=".",s=3)
    plt.xlabel('n')
    plt.ylabel('log(ratio)')
    plt.yscale('log')
    plt.title('nplot')
    plt.show()

# for the chainbuildC(n) returns a plot of samples for a range of c's to hopefully find out something about what the optimal c is
def chainplotc(ci, cf, n, o):
    C = chainC(n)
    x=[]
    y=[]
    for i in range(ci,cf+1):
        x.extend([i for j in range(o)])
    for i in range(ci,cf+1):
        y=y+chainsample(n,i,o,C)
    plt.scatter(x,y,marker=".",s=3)
    plt.xlabel('c')
    plt.ylabel('ratio')
    plt.title('cplot')
    plt.show()

def avgchainplotn(ni,nf,c,o,s):
    x=[]
    y=[]
    for i in range(ni,nf+1,s):
        x.append(i)
    for i in range(ni,nf+1,s):
        y.append(avg(chainsample(i,c,o,chainC(i))))
    plt.scatter(x,y,marker="o",s=5)
    plt.xlabel('log(n)')
    plt.ylabel('ratio')
    #plt.ylabel('log(ratio)')
    plt.xscale('log')
    plt.title('nplot')
    plt.show()
