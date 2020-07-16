from whole import *

C = cvx.matrix(buildC(4))
K = np.add(buildC(4),np.identity(12))
v = findGenVecGram(K)
y=round(v, 2, C)["blochvec"]
print(ratio(C,y,4))
