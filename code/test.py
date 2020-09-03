from whole import *
from tfieldising import *
from qutritalgo import *
import sys
np.set_printoptions(threshold=sys.maxsize)
# returns upon input of n the unabridged prettified solution of the sdp
def test(n):
    a = np.real(np.array(mc(cvx.matrix(tC(n)))))
    tol = 1e-8
    a.real[abs(a.real) < tol] = 0.0
    print(n)
    return a

def roundeig(M):
    w = np.linalg.eigvals(M)
    for i in range(len(w)):
        if w[i]<0:
            if abs(w[i])<10e-8: w[i]=0
    return w

tavgplot(4,2000,2,20,20)
