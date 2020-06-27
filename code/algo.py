import numpy as np

n #Integer
T = c*sqrt(log(n)) #c=O(1)
r = np.random.random(3n) #vector of 3n i.d.d. N(0,1) random variables
v = [] #array for the v_i
y = []
for i in range(1,3n+1):
z = np.inner(r, v[i])
    if abs(z) > 1/sqrt(3): y.append(np.sign(z)/sqrt(3))
    else: y.append(z)
return rho
