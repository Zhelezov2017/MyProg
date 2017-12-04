import cmath
import matplotlib.pyplot as plt
import numpy as np
import scipy.special as special
from sympy import besselj, hankel2


def dispeq_gyrotropic_cylinder(p, m, a_0, k0, EE, GG, HH):

    # вычисляем элементы матрицы рассеяния для внутренней области цилиндра
    mainq = EE ** 2 - GG ** 2 + EE * HH - (HH + EE) * p ** 2
    #F=cmath.sqrt(EE*GG**2*HH*(GG**2-(HH-EE)**2))
    #Pb=cmath.sqrt(EE-(HH+EE)*GG**2/(HH-EE)**2-2*F/(HH-EE)**2)
    #Pc=cmath.sqrt(EE-(HH+EE)*GG**2/(HH-EE)**2+2*F/(HH-EE)**2)
    radq = cmath.sqrt((HH - EE) ** 2 * p **4 + 2 * ((GG ** 2)*(HH + EE) - EE * (HH - EE) **2) * p **2 +(EE ** 2 - GG ** 2 - EE * HH) ** 2)
    #radq = (HH-EE)*cmath.sqrt((p**2-Pb**2)*(p**2-Pc**2))

    q1 = cmath.sqrt(0.5 * (mainq - radq) / EE)
    q2 = cmath.sqrt(0.5 * (mainq + radq) / EE)

    n1 = -(EE * (p * GG)**(-1)) * (p ** 2 + q1 ** 2 + (GG ** 2)/EE - EE)
    n2 = -(EE * (p * GG)**(-1)) * (p ** 2 + q2 ** 2 + (GG ** 2)/EE - EE)

    alp1 = -1 + (p ** 2 + q1 ** 2 - EE) / GG
    alp2 = -1 + (p ** 2 + q2 ** 2 - EE) / GG

    bet1 = 1 + p / n1
    bet2 = 1 + p / n2
#что делать с Q она мнимая, брать действительную или мнимую часть?
    Q1 = (q1 * a_0) * k0
    Q2 = (q2 * a_0) * k0
    JM = np.zeros(((1, 22)), dtype=np.complex)

    JM1 = special.jv(m + 1, Q1) #JM1
    JM2 = special.jv(m + 1, Q2) #JM2
    Jm1 = special.jv(m, Q1)  #Jm1
    Jm2 = special.jv(m, Q2)  #Jm2
    Jm1_Q1 = Jm1 / Q1  #Jm1_Q1
    Jm2_Q2 = Jm2/ Q2  #Jm2_Q2

    Ez1 = (((1j / HH) * n1) * q1) * Jm1  #Ez1
    Ez2 = (((1j / HH) * n2) * q2) * Jm2 #Ez2
    Ephi1 = 1j * (JM1 + (alp1 * m) * Jm1_Q1)  #Ephi1
    Ephi2= 1j * (JM2 + (alp2 * m) * Jm2_Q2)  #Ephi2
    Hz1= - q1 * Jm1  #Hz1
    Hz2= - q2 * Jm2 #Hz2
    Hphi1 = - n1 * (JM1 - (bet1 * m) * Jm1_Q1)  #Hphi1
    Hphi2= - n2 * (JM2 - (bet2 * m) * Jm2_Q2)  #Hphi2

        # вычисляем элементы матрицы рассеяния для внешней области цилиндра

    q = cmath.sqrt(1 - p ** 2)
    q = q * (2 * ((q.imag) <= 0) - 1)
    Q = k0 * a_0 * q

    H2m = special.hankel2(m, Q)  #H2m
    dH2m = (H2m * m)/Q - special.hankel2(m + 1, Q)  #dH2m
        #Jm_out = besselj(m, Q)
        #dJm_out = Jm_out * m / Q - besselj(m + 1, Q)

    Ez_sct1 = -q * H2m #Ez_sct
    Hz_sct1 = -q * H2m  #Hz_sct

    A = 1 / (k0 * (1 - p **2))
    Ephi_sctE1 = q * (A * ((p * (m / a_0)) * H2m)) #Ephi_sctE1
    Ephi_sctH1 = -1j*q * (A * (k0 * q * dH2m)) #Ephi_sctH1
    Hphi_sctE1 = 1j*q * (A * k0 * q * dH2m) #Hphi_sctE1
    Hphi_sctH1 = q * (A * ((p * (m / a_0)) * H2m )) #Hphi_sctH1


    aa11 = Ephi1
    aa12 = Ephi2
    aa13 = Ephi_sctE1
    aa14 = Ephi_sctH1


    aa21 = Ez1
    aa22 = Ez2
    aa23 = Ez_sct1
    aa24 = Q1 - Q1

    aa31 = Hphi1
    aa32 = Hphi2
    aa33 = Hphi_sctE1
    aa34 = Hphi_sctH1

    aa41 = Hz1
    aa42 = Hz2
    aa43 = Q1 - Q1
    aa44 = Hz_sct1

    zz = abs(aa14 * aa23 * aa32 * aa41 - aa13 * aa24 * aa32 * aa41 -
             aa14 * aa22 * aa33 * aa41 + aa12 * aa24 * aa33 * aa41 +
             aa13 * aa22 * aa34 * aa41 - aa12 * aa23 * aa34 * aa41 -
             aa14 * aa23 * aa31 * aa42 + aa13 * aa24 * aa31 * aa42 +
             aa14 * aa21 * aa33 * aa42 - aa11 * aa24 * aa33 * aa42 -
             aa13 * aa21 * aa34 * aa42 + aa11 * aa23 * aa34 * aa42 +
             aa14 * aa22 * aa31 * aa43 - aa12 * aa24 * aa31 * aa43 -
             aa14 * aa21 * aa32 * aa43 + aa11 * aa24 * aa32 * aa43 +
             aa12 * aa21 * aa34 * aa43 - aa11 * aa22 * aa34 * aa43 -
             aa13 * aa22 * aa31 * aa44 + aa12 * aa23 * aa31 * aa44 +
             aa13 * aa21 * aa32 * aa44 - aa11 * aa23 * aa32 * aa44 -
             aa12 * aa21 * aa33 * aa44 + aa11 * aa22 * aa33 * aa44)

    return zz

a_0 = 100
L = 100 * a_0
n_e = 1e6
e_0 = 4.8e-10
m_e = 9.1e-28
c = 3e10
H0 = 0.5
w_H = (e_0 * H0) / (m_e * c)
w_p = cmath.sqrt(4 * cmath.pi * e_0 ** 2 * n_e / m_e)
w_0 = 4.97 * w_H
k_0 = w_0 / c
EE = 1 + w_p ** 2 / (w_H ** 2 - w_0 ** 2)
GG = -w_p ** 2 * w_H / ((w_H ** 2 - w_0 ** 2) * w_0)
HH = 1 - w_p ** 2 / w_0 ** 2
p = [i for i in np.arange(1.2,10,0.1)]
X=[]
mMax=1
#zz_single = 1
#for mm in np.arange(-mMax,mMax+1,1):
#    zz_single =(zz_single * dispeq_gyrotropic_cylinder(1.2, mm, a_0, k_0, EE, GG, HH))
#X.append(cmath.log(zz_single))
#print(X)

#plt.plot(p, X)
#plt.show()