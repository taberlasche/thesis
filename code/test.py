from whole import *
import sys
np.set_printoptions(threshold=sys.maxsize)
def test(n):
    a = np.real(np.array(mc(cvx.matrix(buildC(n)))))
    tol = 1e-8
    a.real[abs(a.real) < tol] = 0.0
    print(np.nonzero(a))
    print(a)
    return a
test(4)
test(8)
test(12)
test(16)
