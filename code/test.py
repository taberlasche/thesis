from whole import *
import sys
np.set_printoptions(threshold=sys.maxsize)
# returns upon input of n the unabridged prettified solution of the sdp
def test(n):
    a = np.real(np.array(mc(cvx.matrix(buildC(n)))))
    tol = 1e-8
    a.real[abs(a.real) < tol] = 0.0
    print(n)
    return a
print(sdp1(4))
print(sdp1(8))
print(sdp1(12))
