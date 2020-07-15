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
    return block_diag(*s)
