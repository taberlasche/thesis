import numpy as np
import picos as pic
import cvxopt as cvx
import cvxopt.lapack
import matplotlib.pyplot as plt
import scipy
from scipy import special
import math
import matplotlib.scale as scale
import statistics

from hamiltonians import *
from sdp import *
from gram import *
from algo import *
from operations import *

# returns 1. the solution to the SDP 2. the gram vectors 3. dictionary with the state and bloch vector upon input of C and c(O(1)).
def solve(C, c):
    V = np.array(mc(C))
    K = (np.conj(V)+V)/2
    v = findGenVecGram(K)
    sol = rund(v,c,C)
    print(sol)
    return sol

# returns upon input of C, y and the maximal eigenvalue corresponding to C the ratio.
def ratio(C, y, maxeig):
    return np.inner(np.dot(y,C),np.transpose(y))/maxeig

# returns list of ratios upon input of number of qubits, c and number of iterations o, for buildC(n).
def sample(n, c, o, C):
    v = findGenVecGram(sdp1(n))
    ratlist=[]
    i=0
    while i < o:
        ratlist.append(ratio(C, rund(v,c,cvx.matrix(C))["blochvec"], n))
        i=i+1
    return ratlist

def chainsample(n, c, o, C):
    v = findGenVecGram(sdp2(n))
    ratlist=[]
    i=0
    while i < o:
        ratlist.append(ratio(C, rund(v,c,cvx.matrix(C))["blochvec"], n))
        i=i+1
    return ratlist

# returns the highest eigenvalue of the n-qubit productstate after performing the rounding o times
def maxsample(n, c, o, C):
    v = findGenVecGram(sdp1(n))
    ratlist=[]
    i=0
    while i < o:
        ratlist.append(ratio(C, rund(v,c,cvx.matrix(C))["blochvec"], n))
        i=i+1
    return max(ratlist)

def maxchainsample(n,c,o,C):
    v = findGenVecGram(sdp2(n))
    ratlist=[]
    i=0
    while i < o:
        ratlist.append(ratio(C, rund(v,c,cvx.matrix(C))["blochvec"], n))
        i=i+1
    return max(ratlist)

# Returns the average of the list.
def avg(l):
    return sum(l)/len(l)

def forward(x):
    return x**2
def inverse(x):
    return x**(1/2)

# For the buildC(n) returns a plot of samples in a range of qubitnumbers upon input of start and end of range, c and the number of rounding attemps per qubitnumber
def plotn(ni, nf, c, o, s):
    x=[]
    y=[]
    for i in range(ni,nf+1,s):
        x.extend([i for j in range(o)])
        y=y+sample(i,c,o,buildC(i))
    plt.scatter(x,y,marker=".",s=3)
    plt.xlabel('n')
    plt.ylabel('ratio')
    plt.title('nplot')
    plt.show()

# Does the same as plotn but returns the averages over the lists per n, and the steps are choosable for efficiency
def avgplot(ni,nf,c,o,s):
    y=[]
    x=[]
    var=[]
    a = np.exp(np.linspace(np.log(ni), np.log(nf), s))
    for i in a:
        i = int(4*round(i/4.))
        x.append(i)
    for i in range(len(x)):
        print(x[i])
        list=sample(x[i],c,o,buildC(x[i]))
        y.append(avg(list))
        var.append(statistics.variance(list))
    plt.errorbar(x,y,yerr=var,fmt='o',elinewidth=1,capsize=3,lw=0)
    plt.xlabel('n')
    plt.ylabel('ratio')
    plt.xscale('log')
    plt.title('Hamiltonian With Linear Terms')
    plt.show()

def maxplot(ni,nf,c,o,s):
    y=[]
    x=[]
    a = np.exp(np.linspace(np.log(ni), np.log(nf), s))
    for i in a:
        i = int(4*round(i/4.))
        x.append(i)
    for i in range(len(x)):
        print(x[i])
        y.append(maxsample(x[i],c,o,buildC(x[i])))
    plt.scatter(x,y,marker="o",s=5)
    plt.xlabel('log(n)')
    plt.ylabel('ratio')
    plt.xscale('log')
    plt.title('Hamiltonian With Linear Terms')
    plt.show()




# For the buildC(n) returns a plot of samples for a range of c's to hopefully find out something about what the optimal c is. Not used.
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

# For the chainC(n) returns a plot of samples in a range of qubitnumbers upon input of start and end of range, c and the number of rounding attemps per qubitnumber.
def chainplot(ni, nf, c, o,s):
    x=[]
    y=[]
    var=[]
    for i in range(ni,nf+1):
        x.append(i)
        y=y+avg(chainsample(i,c,o,chainC(i)))
    plt.scatter(x,y,marker=".",s=3)
    plt.xlabel('n')
    plt.ylabel('log(ratio)')
    plt.yscale('log')
    plt.title('nplot')
    plt.show()

def maxchainplot(ni,nf,c,o,s):
    y=[]
    x=[]
    a = np.exp(np.linspace(np.log(ni), np.log(nf), s))
    for i in a:
        i = int(4*round(i/4.))
        x.append(i)
    for i in range(len(x)):
        y.append(avg(maxsample(x[i],c,o,buildC(x[i]))))
        print(x[i])
    plt.scatter(x,y,marker="o",s=5)
    plt.xlabel('log(n)')
    plt.ylabel('ratio')
    plt.xscale('log')
    plt.title('s')
    plt.show()

# For the chainbuildC(n) returns a plot of samples for a range of c's to hopefully find out something about what the optimal c is.
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

# Plots the averages of the lists of ratios for chainC(n).
def avgchainplot(ni,nf,c,o,s):
    x=[]
    y=[]
    var=[]
    a = np.exp(np.linspace(np.log(ni), np.log(nf), s))
    for i in a:
        i = int(round(i))
        x.append(i)
    for i in range(len(x)):
        print(x[i])
        list=chainsample(x[i],c,o,chainC(x[i]))
        y.append(avg(list))
        var.append(statistics.variance(list))
    plt.errorbar(x,y,yerr=var,fmt='o',elinewidth=1,capsize=3,lw=0)
    plt.xlabel('n')
    plt.ylabel('ratio')
    plt.xscale('log')
    plt.title('One Dimensional Ising Model')
    plt.show()

# Returns the maximal eigenvalue of the transverse field Ising model upon input of a, b, n and a constant d.
def tfimaxeig(n,a,b,d):
    x = a+b/2+(2*a*n/math.pi)*(1+b/(2*a))*d
    return x

# Returns list of ratios for the transverse field Ising model.
def tfisample(n, c, o, C,a,b,d):
    v = findGenVecGram(tfisdp(n))
    ratlist=[]
    i=0
    while i < o:
        ratlist.append(ratio(C, rund(v,c,cvx.matrix(C))["blochvec"], tfimaxeig(n,a,b,d)))
        i=i+1
    return ratlist

# Plots the average of the list of ratios for the transverse field Ising model.
def tfiplot(ni,nf,c,iterations,steps,a,b):
    l = math.sqrt((2*b/a)/((1+b/2*a)**2))
    d=scipy.special.ellipeinc(math.pi/2,l)
    x=[]
    y=[]
    var=[]
    a = np.exp(np.linspace(np.log(ni), np.log(nf), steps))
    for i in a:
        i = int(2*round(i/2.))
        x.append(i)
    for i in range(len(x)):
        print(x[i])
        list=tfisample(x[i],c,iterations,tfiC(x[i],a,b),a,b,d)
        y.append(avg(list))
        plt.errorbar(x,y,yerr=var,fmt='o',elinewidth=1,capsize=3,lw=0)
        var.append(statistics.variance(list))
    plt.xlabel('n')
    plt.ylabel('ratio')
    plt.xscale('log')
    plt.title('Transverse Field Ising Model')
    plt.show()
