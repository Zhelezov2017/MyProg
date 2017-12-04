import math
import cmath
import matplotlib.pyplot as plt
import numpy as np
import scipy.special as special
from sympy import besselj, hankel2

from Xan import dispeq_gyrotropic_cylindersArrayFull
from Dispeq_girotropic_cylinder import dispeq_gyrotropic_cylinder

a_0 = 100
L = 5 * a_0
n_e = 1e6
e_0 = 4.8e-10
m_e = 9.1e-28
c = 3e10
H0 = 0.5
w_H = (e_0 * H0) / (m_e * c)
w_p = cmath.sqrt(4 * math.pi * e_0 ** 2 * n_e / m_e)
w_0 = 4.97 * w_H
k_0 = w_0 / c
EE = 1 + w_p ** 2 / (w_H ** 2 - w_0 ** 2)
GG = -w_p ** 2 * w_H / ((w_H ** 2 - w_0 ** 2) * w_0)
HH = 1 - w_p ** 2 / w_0 ** 2
p = [i for i in np.arange(1.1,10,0.1)]
cylXY1 = np.array([[-L, 0.001], [0, 0.001]])

X=[]
Y=[]
mMax = 1
for op in np.arange(1.1,10,0.1):
    zz_single = 1
    for mm in np.arange(-mMax,mMax+1,1):
        zz_single = (zz_single * (dispeq_gyrotropic_cylinder(op, mm, a_0, k_0, EE, GG, HH)))
    X.append(cmath.log(zz_single**2))
for yp in np.arange(1.1,10,0.1):
    result =cmath.log(dispeq_gyrotropic_cylindersArrayFull(1, w_0, a_0, 1, cylXY1, k_0, yp, EE, GG, HH, 3e10))
    Y.append(result)
plt.plot(p, X, p, Y)
plt.show()
print(X)