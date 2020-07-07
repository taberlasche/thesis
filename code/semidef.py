import cvxopt as cvx
import numpy as np

# a cvxopt.matrix can be converted to numpy using the array() function B = array(A) and vice versa A = matrix(B)

# M: in herm(3n)
# C: C^T = C real n x n with C_ij=0 for i=j

# maximize Tr(CM) over M in herm(3n)
# subject to M >= 0  and M_i,i = 1 for all i
c = cvx.matrix()
n = int(np.sqrt(len(c)))
# the constraint M_ii =1 in the Form G*vec(M)=h
Gx = []
for i in range(1,n):
    Gx.append(1)
    Gx.extend([0 for i in range(n)])
Gx.append(1)

G = cvx.matrix(np.reshape(np.diagflat(Gx), (n**2,n**2)).astype(np.double))
h = cvx.matrix(np.asarray(Gx).flatten('F').astype(np.double))

sol = cvx.solvers.sdp(c, Gl=G, hl=h)

print(sol['x'])
