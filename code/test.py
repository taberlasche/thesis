from whole import *
import sys
np.set_printoptions(threshold=sys.maxsize)
# returns upon input of n the unabridged prettified solution of the sdp
def test(n):
    a = np.real(np.array(mc(cvx.matrix(chainC(n)))))
    tol = 1e-8
    a.real[abs(a.real) < tol] = 0.0
    print(n)
    return a
print(sdp2(4))
print(sdp2(8))
print(sdp2(12))
print(sdp2(16))



#print(test(4))
#print(test(8))
#print(test(12))
#print(test(16))
