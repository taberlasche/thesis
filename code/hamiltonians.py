import numpy as np
from scipy.linalg import block_diag

#This function assumes that n is a multiple of 4 it returns C for Hamiltonians of the form X1X2+Z1Z2+X5X6+Z5Z6... upon input of the number of qubits
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
        A[14+i*12][6+i*12] = 1
        A[6+i*12][14+i*12] = 1
        A[14+i*12][9+i*12] = 1
        A[9+i*12][14+i*12] = 1
    return A


def chainC(n):
    a = []
    a.extend([1,0,0] for i in range(n-1))
    C = np.diagflat(a,3) + np.diagflat(a,-3)
    C[0][3*n-3]=1
    C[3*n-3][0]=1
    return C

def chainsdp(n):
    O = np.zeros((3,3))
    I = np.eye(3*n)
    A = np.array([[1,0,0], [0,0,0], [0,0,0]])
    B = [O]
    for i in range(n-1):
        B.extend(A for j in range(n))
        B.append(O)
    return np.reshape(np.block(np.vsplit(np.stack(B),n)), (3*n,3*n)) + I
