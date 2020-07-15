import picos as pic
import cvxopt as cvx
import cvxopt.lapack
import numpy as np
def mc(C):
    N = C.size[0]
    maxcut = pic.Problem()

    # Add the symmetric matrix variable.
    X=maxcut.add_variable('X', (N,N), 'hermitian')

    # Constrain X to have ones on the diagonal.
    maxcut.add_constraint(pic.maindiag(X) == 1)

    # Constrain X to be positive semidefinite.
    maxcut.add_constraint(X >> 0)

    # Set the objective.
    maxcut.set_objective('max', C|X)

    print(maxcut)

    # Solve the problem.
    maxcut.solve(solver='cvxopt')
    print(X)
    return X.value
