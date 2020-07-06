import picos as pic
import cvxopt as cvx
import cvxopt.lapack
import numpy as np

# Make G undirected.
G=nx.Graph(G)

# Allocate weights to the edges.
for (i,j) in G.edges():
  G[i][j]['weight']=c[i,j]+c[j,i]

maxcut = pic.Problem()

# Add the symmetric matrix variable.
X=maxcut.add_variable('X', (N,N), 'symmetric')

# Retrieve the Laplacian of the graph.
LL = 1/4.*nx.laplacian_matrix(G).todense()
L=pic.new_param('L', LL)

# Constrain X to have ones on the diagonal.
maxcut.add_constraint(pic.diag_vect(X) == 1)

# Constrain X to be positive semidefinite.
maxcut.add_constraint(X >> 0)

# Set the objective.
maxcut.set_objective('max', L|X)

#print(maxcut)

# Solve the problem.
maxcut.solve(solver='cvxopt')

#print('bound from the SDP relaxation: {0}'.format(maxcut.obj_value()))
