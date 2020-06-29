"""
© 2018 Stefan Huber, Robert König, Marco Tomamichel

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""
from tree import *
from problem import *
import numpy as np
import picos as sdp
import cvxopt as cvx

# Pauli sigmas
sigmas = [[[1.,0.],[0.,1.]], [[0.,1.],[1.,0.]], [[0.,1.j],[-1.j,0.]], [[1.,0.],[0.,-1.]]]
# their transpose
sigmast = [[[1.,0.],[0.,1.]], [[0.,1.],[1.,0.]], [[0.,-1.j],[1.j,0.]], [[1.,0.],[0.,-1.]]]
# their trace
tracesigma = [2.0, 0.0, 0.0, 0.0]

def Pauli(j, nqubits, skip={}, transp={}):
    # returns element of Pauli basis for Hermitians acting on Hilbert space
    # of dimension 2^n normalized with respect to trace
    # j: index in [0, 4^n-1]
    # n: number of qubits
    assert(isinstance(j, int) and isinstance(nqubits, int))
    assert(j >= 0 and j < 4**nqubits)

    val = np.matrix([[1]])
    for idx in range(nqubits):
        if idx not in skip:
            # tensor corresponding sigmas[j % 4]
            if idx in transp:
                val = np.kron(val, sigmast[j % 4])
            else:
                val = np.kron(val, sigmas[j % 4])
        else:
            # instead just take the trace of sigmas[j % 4]
            val = val * tracesigma[j % 4]
        j = int(np.floor(j / 4))
    return val/sqrt(2**nqubits)


class QuantumBilinearProblem(BilinearProblem):

    def __init__(self):
        BilinearProblem.__init__(self)
        self.dims = [[]]

    def init_matrix_form(self, dims, J, A = [], B = [], maximize = False, verb = 1):
        # syntax for dims - [[a_1, a_2, a_3, ...], [b_1, b_2, b_3, ...]]
        # all these values currently have to be powers of 2
        # interprets first variable as a direct sum of a_i-dimensional subsystems,
        # and the second variable as a direct sum of b_i-dimensional subsystems
        # J is a array of the form of dims where every entry is a a_i*b_j matrix
        # defining the correlations between the i-th X-variable and the j-th Y-variable

        self.dims = dims

        # total dimension of the vector problem
        dimX = np.sum(np.array(dims[0])**2)
        dimY = np.sum(np.array(dims[1])**2)
        numX = len(dims[0])
        numY = len(dims[1])

        if verb >= 1:
            print('problem shape:', numX, 'x', numY)
        assert(len(J) == numX and len(J[0]) == numY)

        # expand default arguments
        if A == []:
            for j in range(numX):
                A.append(np.zeros((dims[0][j], dims[0][j])))
        if B == []:
            for j in range(numY):
                B.append(np.zeros((dims[1][j], dims[1][j])))

        # initialize Q matrix and vectors a and b
        Q = np.zeros((dimX, dimY))
        a = np.zeros(dimX)
        b = np.zeros(dimY)

        basei = 0
        for k in range(numX):
            nqX = int(np.log2(dims[0][k])) # dimension of X variable (in qubits)

            basej = 0
            for l in range(numY):
                nqY = int(np.log2(dims[1][l])) # dimension of Y variable (in qubits)

                if verb >= 2:
                    print('correlations matrix: ', nqX, 'x', nqY, 'qbits.')

                corr = np.matrix(J[k][l])  # correlation matrix
                # create entries for Q matrix
                for i in range(dims[0][k]**2):
                    for j in range(dims[1][l]**2):
                        Q[basei+i][basej+j] = \
                            np.real(np.trace( np.kron(Pauli(i, nqX), Pauli(j, nqY)) * corr ))

                # create vector a
                for i in range(dims[0][k]**2):
                    a[basei+i] = np.real(np.trace( Pauli(i, nqX) * np.matrix(A[k]) ))
               # create vector b
                for j in range(dims[1][l]**2):
                    b[basej+j] = np.real(np.trace( Pauli(j, nqY) * np.matrix(B[l]) ))

                basej += dims[1][l]**2
            basei += dims[0][k]**2

        if verb >= 1:
            print('convert to vector problem...')
            print('** Q =', Q)
            print('** a =', a)
            print('** b =', b)

        self.init_vector_form(Q=Q, a=a, b=b, maximize=maximize, verb=verb)

    # helper functions that create matrix variable and matrix solutions
    def __getvarparams(self, var, idx):
        dim = self.dims[var][idx]
        return (dim, int(np.log2(dim)), int(np.sum(np.array(self.dims[var][0:idx])**2)))

    def matvarX(self, idx=0, qubits=None, transp={}):
        (sqrtdim, nqubits, baseidx) = self.__getvarparams(0, idx)

        if qubits == None:
            skip = {}
        else:
            skip = set(range(nqubits)) - qubits

        X = self.varX()[baseidx]*cvx.matrix(Pauli(0, nqubits, skip))
        for i in range(1, sqrtdim**2):
            X += self.varX()[baseidx+i]*cvx.matrix(Pauli(i, nqubits, skip, set(transp)))
        return X

    def matvarY(self, idx=0, qubits=None, transp={}):
        (sqrtdim, nqubits, baseidx) = self.__getvarparams(1, idx)

        if qubits == None:
            skip = {}
        else:
            skip = set(range(nqubits)) - qubits

        Y = self.varY()[baseidx]*cvx.matrix(Pauli(0, nqubits, skip))
        for i in range(1, sqrtdim**2):
            Y += self.varY()[baseidx+i]*cvx.matrix(Pauli(i, nqubits, skip, set(transp)))
        return Y

    def matsolX(self, idx = 0):
        (sqrtdim, nqubits, baseidx) = self.__getvarparams(0, idx)

        X = np.zeros((sqrtdim, sqrtdim), dtype=np.complex_)
        for i in range(0, sqrtdim**2):
            X += self.solX[baseidx+i]*Pauli(i,nqubits)
        return np.array(X)

    def matsolY(self, idx = 0):
        (sqrtdim, nqubits, baseidx) = self.__getvarparams(1, idx)

        Y = np.zeros((sqrtdim, sqrtdim), dtype=np.complex_)
        for i in range(0, sqrtdim**2):
            Y += self.solY[baseidx+i]*Pauli(i,nqubits)
        return np.array(Y)
