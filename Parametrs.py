import math
import cmath
import matplotlib.pyplot as plt
import numpy as np
import scipy.special as special
from sympy import besselj, hankel2

a_0 = 100
L = 10 * a_0
n_e = 1e6
e_0 = 4.8e-10
m_e = 9.1e-28
c = 3e10
H0 = 0.5
w_H = (e_0 * H0) / (m_e * c)
w_p = cmath.sqrt(4 * math.pi * e_0 ** 2 * n_e / m_e)
w_0 = 4.97 * w_H
k_0 = w_0 / c
    # w_0
    # N_max=mMax=...=3
    # a_0
    # ee=1
cylXY1 = np.array([[-L, 0.001], [0, 0.001]])
    # k=  w_0 / c
    # p
    # EE
    # GG
    # HH
    # c
    #
    # xmin =  25
    # xmax =  40
    # ymin =  -4
    # ymax =   4
    # Npntx = 100
    # Npnty = 100
cylXY = np.array([[0, 0.001]])
EE = 1 + w_p ** 2 / (w_H ** 2 - w_0 ** 2)
GG = -w_p ** 2 * w_H / ((w_H ** 2 - w_0 ** 2) * w_0)
HH = 1 - w_p ** 2 / w_0 ** 2


