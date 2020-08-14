from whole import *
from tfieldising import *
import sys
np.set_printoptions(threshold=sys.maxsize)
# returns upon input of n the unabridged prettified solution of the sdp
def test(n,a,b):
    a = np.real(np.array(mc(cvx.matrix(tfiC(n,a,b)))))
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
<<<<<<< HEAD

x=[]
y=[]
for i in range(2,20):
    x.append(i)
    y.append(tfisingmaxeig(1,1,i))
plt.scatter(x,y,marker=".",s=3)
plt.xlabel('n')
plt.ylabel('maxeig')
#plt.ylabel('log(ratio)')
#plt.yscale('log')
plt.title('nplot')
plt.show()
=======
print(tfisdp(4))
print(tfisdp(6))
print(tfisdp(8))
print(tfisdp(10))
print(tfisdp(12))
print(tfisdp(14))
>>>>>>> 6e8cda01c29b9e9c46601cc4f0c8b812c9fe58d1
