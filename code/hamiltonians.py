import numpy as np
from scipy.linalg import block_diag

#This function assumes that n is a multiple of 4 it returns C for H' the form X1X2+Z1Z2#Z+X5X6+Z5Z6... upon input of the number of qubits
def buildC(n):
    o1 = np.zeros((3,3))
    a = np.array([[1,0,0], [0,0,0], [0,0,1]])
    b = np.block([[o1, a], [a, o1]])
    o2 = np.zeros((6,6))
    d = np.block([[b, o2], [o2, o2]])
    s = []
    for i in range(1,int(n/4)+1):
        s.append(d)
    s.append(o1)
    A = np.array(block_diag(*s))
    for i in range(0,int(n/4)):
        A[3*n+2][6+i*12] = 1
        A[6+i*12][3*n+2] = 1
        A[3*n+2][9+i*12] = 1
        A[9+i*12][3*n+2] = 1
    return A

def sdp1(n):
    A = np.add(buildC(n),np.eye(3*n+3))
    if int(n/4)>1:k=1
    else: k=0
    for i in range(0,int(n/4)):
        A[6+12*i][9+i*12] = 1
        A[9+i*12][6+12*i] = 1
        for j in range(0,int(n/4)-i-1):
            A[6+i*12][9+i*12+k*9+12*j] = 1
            A[9+i*12+k*9+12*j][6+i*12] = 1
            A[9+i*12][9+i*12+k*9+12*j] = 1
            A[9+i*12+k*9+12*j][9+i*12] = 1
            A[6+i*12][12+i*12+k*9+12*j] = 1
            A[12+i*12+k*9+12*j][6+i*12] = 1
            A[9+i*12][12+i*12+k*9+12*j] = 1
            A[12+i*12+k*9+12*j][9+i*12] = 1
    return A


# returns the C for a chain of spins
def chainC(n):
    a = []
    a.extend([1,0,0] for i in range(n-1))
    C = np.diagflat(a,3) + np.diagflat(a,-3)
    C[0][3*n-3]=1
    C[3*n-3][0]=1
    return C
# returns the sdp solution for chainC(n)
def sdp2(n):
    A = np.eye(3*n)
    for i in range(0,n):
        for j in range(0,n):
            A[i*3][j*3] = 1
    return A

def tfiC(n,a,b):
    x = []
    x.extend([1,0,0] for i in range(n-1))
    x.append([0,0,0])
    C = -1*a*(np.diagflat(x,3) + np.diagflat(x,-3))
    C[3*n-3][0]=-1*a
    C[0][3*n-3]=-1*a
    for i in range(2,3*n,3):
        C[3*n+2][i]=-1*b
        C[i][3*n+2]=-1*b
    return C

def tfisdp(n):
    a = np.eye(3*n+3)
    b,c=[],[]
    for i in range(int(n/2)):
        b.extend([1,0,0,-1,0,0])
        c.extend([-1,0,0,1,0,0])
    b.extend([0,0,0])
    c.extend([0,0,0])
    np.array(b).reshape(-1)
    np.array(c).reshape(-1)
    for i in range(0,3*n-5,6):
        a[i]=b
    for i in range(3,3*n-2,6):
        a[i]=c
    f=[]
    for i in range(n):
        f.extend([0,0,1])
    f.extend([0,0,0])
    np.array(f).reshape(-1)
    for i in range(2,3*n,3):
        a[i]=f
    for i in range(2,3*n,3):
        a[3*n+2][i]=-1
        a[i][3*n+2]=-1
    return a

def qutrith(n):
    x=[]
    for i in range():
    x.extend([1,0,0,0,0,0,0,0] for i in range(n-1))
    return np.diagflat(x,3) + np.diagflat(x,-3)
