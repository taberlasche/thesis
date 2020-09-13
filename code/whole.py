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

from timeit import default_timer as timer
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
    print('finding gram vectors...')
    start = timer()
    v = findGenVecGram(sdp1(n))
    end = timer()
    print('elapsed time:',end-start)
    ratlist=[]
    i=0
    while i < o:
        ratlist.append(ratio(C, rund(v,c,cvx.matrix(C))["blochvec"], n))
        print('step',i)
        i=i+1
    return ratlist

# returns list of ratios upon input of number of qubits, c and number of iterations o, for chainC(n).
def chainsample(n, c, o, C):
    print('finding gram vector...')
    start = timer()
    v = findGenVecGram(sdp2(n))
    end = timer()
    print('elapsed time:',end-start)
    ratlist=[]
    i=0
    while i < o:
        ratlist.append(ratio(C, rund(v,c,cvx.matrix(C))["blochvec"], n))
        print('step',i)
        i=i+1
    return ratlist

# Returns the average of the list.
def avg(l):
    return sum(l)/len(l)

# Outputs the average, the variance and the maximum of o iterations for s equidistant steps on a log scale in the range [ni,nf].
def avgplot(ni,nf,c,o,s):
    y=[]
    x=[]
    var=[]
    maxlist=[]
    a = np.exp(np.linspace(np.log(ni), np.log(nf), s))
    for i in a:
        i = int(4*round(i/4.))
        x.append(i)
    for i in range(len(x)):
        print(x[i])
        list=np.real(sample(x[i],c,o,buildC(x[i])))
        y.append(avg(list))
        var.append(statistics.variance(list))
        maxlist.append(np.max(list))
    print('x:',x)
    print('y',y)
    print('var',var)
    f=open("avgplot().txt","a+")
    f.writelines(['Parameters:',str(ni),str(nf),str(c),str(o),str(s),'x:',str(x),'y:',str(y),'var:',str(var),'maxlist:',str(maxlist),"\n"])
    f.close()
    #plt.errorbar(x,y,yerr=var,fmt='o',elinewidth=1,capsize=3,lw=0)
    #plt.xlabel('n')
    #plt.ylabel('ratio')
    #plt.xscale('log')
    #plt.title('Hamiltonian With Linear Terms')
    #plt.show()

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


# For the chainbuildC(n) returns a plot of samples for a range of c's to hopefully find out something about what the optimal c is. Not used.
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

# Outputs the average, the variance and the maximum of o iterations for s equidistant steps on a log scale in the range [ni,nf].
def avgchainplot(ni,nf,c,o,s):
    x=[]
    y=[]
    var=[]
    maxlist=[]
    a = np.exp(np.linspace(np.log(ni), np.log(nf), s))
    for i in a:
        i = int(2*round(i/2.))
        x.append(i)
    for i in range(len(x)):
        print(x[i])
        list=np.real(chainsample(x[i],c,o,chainC(x[i])))
        y.append(avg(list))
        var.append(statistics.variance(list))
        maxlist.append(np.max(list))
    print('x:',x)
    print('y',y)
    print('var',var)
    f=open("avgchainplot().txt","a+")
    f.writelines(['Parameters:',str(ni),str(nf),str(c),str(o),str(s),'x:',str(x),'y:',str(y),'var:',str(var),'max:',str(maxlist),"\n"])
    f.close()
    #plt.errorbar(x,y,yerr=var,fmt='o',elinewidth=1,capsize=3,lw=0)
    #plt.xlabel('n')
    #plt.ylabel('ratio')
    #plt.xscale('log')
    #plt.title('One Dimensional Ising Model')
    #plt.show()

# Returns the maximal eigenvalue of the transverse field Ising model upon input of a, b, n and a constant d.
def tfimaxeig(n,a,b,d):
    x = a+b/2+(2*a*n/math.pi)*(1+b/(2*a))*d
    return x

# Returns list of ratios for the transverse field Ising model.
def tfisample(n, c, o, C,a,b,d):
    print('finding gram vector...')
    start = timer()
    v = findGenVecGram(tfisdp(n))
    end = timer()
    print('elapsed time:', end-start)
    ratlist=[]
    i=0
    while i < o:
        ratlist.append(ratio(C, rund(v,c,cvx.matrix(C))["blochvec"], tfimaxeig(n,a,b,d)))
        print('step',i)
        i=i+1
    return ratlist

# Outputs the average, the variance and the maximum of o iterations for s equidistant steps on a log scale in the range [ni,nf].
def tfiplot(ni,nf,c,iterations,steps,a,b):
    l = math.sqrt((2*b/a)/((1+b/2*a)**2))
    d=scipy.special.ellipeinc(math.pi/2,l)
    x=[]
    y=[]
    var=[]
    maxlist=[]
    k = np.exp(np.linspace(np.log(ni), np.log(nf), steps))
    for i in k:
        i = int(2*round(i/2.))
        x.append(i)
    for i in range(len(x)):
        print(x[i])
        alist=np.real(tfisample(x[i],c,iterations,tfiC(x[i],a,b),a,b,d))
        y.append(avg(alist))
        var.append(statistics.variance(alist))
        maxlist.append(np.max(alist))
    f=open("tfiplot().txt","a+")
    f.writelines(['Parameters:',str(ni),str(nf),str(c),str(iterations),str(steps),str(a),str(b),'x:',str(x),'y:',str(y),'var:',str(var),'max:',str(maxlist),"\n"])
    f.close()
    print('x:',x)
    print('y',y)
    print('var',var)
    #plt.errorbar(x,y,yerr=var,fmt='o',elinewidth=1,capsize=3,lw=0)
    #plt.xlabel('n')
    #plt.ylabel('ratio')
    #plt.xscale('log')
    #plt.title('Transverse Field Isisng Model')
    #plt.show()
