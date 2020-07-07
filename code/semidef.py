import cvxopt as cvx
import numpy as np


def solve(c):

    n = c.size[0]
    G = -1*cvx.matrix(np.identity(n**2).astype(np.double))
    h = cvx.matrix(np.zeros((n**2,1), dtype=np.double))
    Y = []
    for i in range(1,n):
        Y.append(1)
        Y.extend([0 for i in range(n)])
    Y.append(1)

    A = cvx.matrix(np.reshape(np.diagflat(Y), (n**2,n**2)).astype(np.double))
    b = cvx.matrix(np.asarray(Y).flatten('F').astype(np.double))

    sol = cvx.solvers.sdp(c, Gl=G, hl=h, A, b)
    print(sol['x'])

#ValueError: Rank(A) < p or Rank([G; A]) < n
