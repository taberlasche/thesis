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

def tC(n):
    a = []
    a.extend([1,0,0,0,0,0,0,0] for i in range(n-1))
    C = np.diagflat(a,8) + np.diagflat(a,-8)
    C[0][8*n-8]=1
    C[8*n-8][0]=1
    return C

def tsdp(n):
    return M

def ratio(C, y, maxeig):
    return np.inner(np.dot(y,C),np.transpose(y))/maxeig

def tsample(n, c, o, C):
    v = findGenVecGram(tsdp(n))
    ratlist=[]
    i=0
    while i < o:
        ratlist.append(ratio(C, trund(v,c,cvx.matrix(C))["blochvec"], n))
        i=i+1
    return ratlist

def avgplot(ni,nf,c,o,s):
    y=[]
    x=[]
    a = np.exp(np.linspace(np.log(ni), np.log(nf), s))
    for i in a:
        i = int(4*round(i/4.))
        x.append(i)
    for i in range(len(x)):
        y.append(avg(sample(x[i],c,o,buildC(x[i]))))
        print(x[i])
    plt.scatter(x,y,marker="o",s=5)
    plt.xlabel('log(n)')
    plt.ylabel('ratio')
    #plt.ylabel('log(ratio)')
    plt.xscale('log')
    plt.title('nplot')
    plt.show()
