import math
import cmath
import matplotlib.pyplot as plt
import numpy as np
from sympy import besselj, hankel2


def dispeq_gyrotropic_cylindersArrayFull(N_max, w_01, a_0, ee, cylXY, k, p, EE, GG, HH, c, ):
    m = list(range(-N_max, N_max + 1, 1))
    if N_max == 0:
        m = 0

    nu = m
    w = w_01
    k0 = w / c
    q = cmath.sqrt(1 - p **2)
    # q = q * (2 * ((q.imag) <= 0) - 1)
    Q = k0 * a_0 * q

    # вычисляем элементы матрицы рассеяния для внутренней области цилиндра
    mainq = EE ** 2 - GG ** 2 + EE * HH - (HH + EE) * p ** 2
    F=cmath.sqrt(EE*GG**2*HH*(GG**2-(HH-EE)**2))
    Pb=cmath.sqrt(EE-(HH+EE)*GG**2/(HH-EE)**2-2*F/(HH-EE)**2)
    Pc=cmath.sqrt(EE-(HH+EE)*GG**2/(HH-EE)**2+2*F/(HH-EE)**2)
    #radq = math.sqrt((HH - EE) ** 2 * p **4 + 2 * ((GG ** 2)*(HH + EE) - EE * (HH - EE) **2) * p **2 +(EE ** 2 - GG ** 2 - EE * HH) ** 2)
    radq = -(HH-EE)*cmath.sqrt((p**2-Pb**2)*(p**2-Pc**2))
    q1 = cmath.sqrt(0.5 * (mainq - radq) / EE)
    q2 = cmath.sqrt(0.5 * (mainq + radq) / EE)

    n1 = -(EE / (p * GG) ) * (p ** 2 + q1 ** 2 + (GG ** 2)/EE - EE)
    n2 = -(EE / (p * GG) ) * (p ** 2 + q2 ** 2 + (GG ** 2)/EE - EE)

    alp1 = -1 + (p ** 2 + q1 ** 2 - EE) / GG
    alp2 = -1 + (p ** 2 + q2 ** 2 - EE) / GG

    bet1 = 1 + p / n1
    bet2 = 1 + p / n2
#что делать с Q она мнимая, брать действительную или мнимую часть?
    Q1 = (q1 * a_0) * k0
    Q2 = (q2 * a_0) * k0
    mMax = 2 * N_max + 1
    JM = np.eye((31) , dtype=np.complex)
    for k in range(0, 31, 1):
        JM[k][k] = 0
    n=0
    for per in range(-N_max,N_max+1,1):
        JM[n][0] = besselj(per + 1, Q1) #JM1
        JM[n][1] = besselj(per + 1, Q2) #JM2
        JM[n][2]= besselj(per, Q1)  #Jm1
        JM[n][3]= besselj(per, Q2)  #Jm2
        JM[n][4]= JM[n][per + 2] / Q1  #Jm1_Q1
        JM[n][5]= JM[n][per + 3]/ Q2  #Jm2_Q2

        JM[n][6]= (((1j / HH) * n1) * q1) * JM[n][per + 2]  #Ez1
        JM[n][7] = (((1j / HH) * n2) * q2) * JM[n][per + 3] #Ez2
        JM[n][8] = 1j * (JM[n][per] + (alp1 * per) * JM[n][per + 4])  #Ephi1
        JM[n][9]= 1j * (JM[n][per+1] + (alp2 * per) *JM[n][per + 5] )  #Ephi2
        JM[n][10]= - q1 * JM[n][per + 2]  #Hz1
        JM[n][11]= - q2 * JM[n][per + 3]  #Hz2
        JM[n][12] = - n1 * (JM[n][per] - (bet1 * per) * JM[n][per + 4])  #Hphi1
        JM[n][13]= - n2 * (JM[n][per+1] - (bet2 * per) * JM[n][per + 5])  #Hphi2

        # вычисляем элементы матрицы рассеяния для внешней области цилиндра

        JM[n][14] = hankel2(per, Q)  #H2m
        JM[n][15]= -(JM[n][per + 14]/ Q)* per + hankel2(per + 1, Q)  #dH2m
        #Jm_out = besselj(m, Q)
        #dJm_out = Jm_out * m / Q - besselj(m + 1, Q)

        JM[n][16] =-1j * q * JM[n][per + 14]  #Ez_sct
        JM[n][17]= q * JM[n][per + 14]  #Hz_sct

        A = 1 / (k0 * (1 - p ^ 2))
        JM[n][18] = -1j*q * (A * ((p * (per / a_0)) * JM[n][per + 14])) #Ephi_sctE
        JM[n][19]= 1j*q * (A * (k0 * q * JM[n][per + 15])) #Ephi_sctH
        JM[n][20]= q * (A * ((p * (per / a_0)) * JM[n][per + 14])) #Hphi_sctH
        JM[n][21] = -q * (A * (k0 * q * JM[n][per + 15])) #Hphi_sctE
        n = n + 1
    matrixGS = np.eye((4 * mMax * 2) , dtype=np.complex)
    for k in range(0, 4 * mMax * 2, 1):
        matrixGS[k][k] = 0
    for jj in range(0, 2):
        for jm in range(-N_max, N_max + 1, 1):
            for ll in range(0, 2):
                for jnu in range(-N_max, N_max + 1, 1):
                    jmmm = (jm - 1) * 4 + 1
                    jnunu = (jnu - 1) * 4 + 1
                        # спросить про проход ведь будет совпадение номеров
                    if jj == ll:
                        if jm == jnu:
                                matrixGS[jmmm + jj * 4 * mMax, jnunu + ll * 4 * mMax] = JM[jm + 1][6]  # Ez1
                                matrixGS[jmmm + jj * 4 * mMax, jnunu + ll * 4 * mMax + 1] = JM[jm + 1][7]  # Ez2(jm)
                                matrixGS[jmmm + jj * 4 * mMax, jnunu + ll * 4 * mMax + 2] = JM[jm + 1][16]  # Ez_sct(jm)
                                matrixGS[jmmm + jj * 4 * mMax, jnunu + ll * 4 * mMax + 3] = 0

                                matrixGS[jmmm + jj * 4 * mMax + 1, jnunu + ll * 4 * mMax] = JM[jm + 1][10]  # Hz1
                                matrixGS[jmmm + jj * 4 * mMax + 1, jnunu + ll * 4 * mMax + 1] = JM[jm + 1][11]  # Hz2
                                matrixGS[jmmm + jj * 4 * mMax + 1, jnunu + ll * 4 * mMax + 2] = 0
                                matrixGS[jmmm + jj * 4 * mMax + 1, jnunu + ll * 4 * mMax + 3] = - JM[jm + 1][17]  # -Hz_sct

                                matrixGS[jmmm + jj * 4 * mMax + 2, jnunu + ll * 4 * mMax] = JM[jm + 1][8]  # Ephi1
                                matrixGS[jmmm + jj * 4 * mMax + 2, jnunu + ll * 4 * mMax + 1] = JM[jm + 1][9]  # Ephi2
                                matrixGS[jmmm + jj * 4 * mMax + 2, jnunu + ll * 4 * mMax + 2] = -JM[jm + 1][18]  # -Ephi_sctE
                                matrixGS[jmmm + jj * 4 * mMax + 2, jnunu + ll * 4 * mMax + 3] = -JM[jm + 1][19]  # -Ephi_sctH

                                matrixGS[jmmm + jj * 4 * mMax + 3, jnunu + ll * 4 * mMax] = JM[jm + 1][12]  # Hphi1
                                matrixGS[jmmm + jj * 4 * mMax + 3, jnunu + ll * 4 * mMax + 1] = JM[jm + 1][13]  # Hphi2
                                matrixGS[jmmm + jj * 4 * mMax + 3, jnunu + ll * 4 * mMax + 2] = -JM[jm + 1][21]  # -Hphi_sctE
                                matrixGS[jmmm + jj * 4 * mMax + 3, jnunu + ll * 4 * mMax + 3] = -JM[jm + 1][20]  # -Hphi_sctH
                    else:
                        Ljl = float(math.sqrt((cylXY[jj , 0] - cylXY[ll , 0]) ** 2 + (cylXY[jj , 1] - cylXY[ll , 1]) ** 2))

                        thetaIJ =float(math.atan(abs(cylXY[jj , 1] - cylXY[ll , 1]) / abs(cylXY[jj , 0] - cylXY[ll , 0])))

                        if ((cylXY[jj , 0] - cylXY[ll , 0]) <= 0 and (cylXY[jj , 1] - cylXY[ll , 1]) > 0):
                            thetaIJ = cmath.pi - thetaIJ
                        elif ((cylXY[jj , 0] - cylXY[ll , 0]) <= 0 and (cylXY[jj, 1] - cylXY[ll , 1]) <= 0):
                            thetaIJ = cmath.pi + thetaIJ
                        elif ((cylXY[jj , 1] - cylXY[ll , 1]) <= 0 and (cylXY[jj, 0] - cylXY[ll , 0]) > 0):
                            thetaIJ = 2 * cmath.pi - thetaIJ
                        r=0
                        for per in range(-N_max, N_max + 1, 1):
                            JM[r][22] = -cmath.exp(- 1j * (jnu - per) * thetaIJ) * q * hankel2((per - jnu),k * Ljl)  # HmInll

                            JM[r][23] = besselj(per, Q)  # Jm_out
                            JM[r][24] =  JM[r][23] * per / Q - besselj(per + 1, Q)  # dJm_out
#спросить почему они одинаковые?
                            JM[r][25] =  JM[r][22]*JM[r][23] # Ez_inc
                            JM[r][26] =  JM[r][22]* JM[r][23]  # Hz_inc

                            A = 1 / (k0 * (1 - p ^ 2))
                            JM[r][27] = - JM[r][22]* (A * ((p * (per / a_0)) *JM[r][23] ))  # Ephi_incE
                            JM[r][28] = - JM[r][22]* (A * (- 1j * k0 * q * JM[r][24]))  # Ephi_incH
                            JM[r][29] = - JM[r][22]* (A * (1j * k0 * q * JM[r][24]))  # Hphi_incE
                            JM[r][30] = -JM[r][22]* (A * ((p * (per / a_0)) * JM[r][23]))  # Hphi_incH
                            r = r + 1

                        matrixGS[jmmm + jj * 4 * mMax, jnunu + ll * 4 * mMax] = 0
                        matrixGS[jmmm + jj * 4 * mMax, jnunu + ll * 4 * mMax + 1] = 0
                        matrixGS[jmmm + jj * 4 * mMax, jnunu + ll * 4 * mMax + 2] = JM[jm + 1][25]  # Ez_inc
                        matrixGS[jmmm + jj * 4 * mMax, jnunu + ll * 4 * mMax + 3] = 0

                        matrixGS[jmmm + jj * 4 * mMax + 1, jnunu + ll * 4 * mMax] = 0
                        matrixGS[jmmm + jj * 4 * mMax + 1, jnunu + ll * 4 * mMax + 1] = 0
                        matrixGS[jmmm + jj * 4 * mMax + 1, jnunu + ll * 4 * mMax + 2] = 0
                        matrixGS[jmmm + jj * 4 * mMax + 1, jnunu + ll * 4 * mMax + 3] = JM[jm + 1][25]  # Hz_inc

                        matrixGS[jmmm + jj * 4 * mMax + 2, jnunu + ll * 4 * mMax] = 0
                        matrixGS[jmmm + jj * 4 * mMax + 2, jnunu + ll * 4 * mMax + 1] = 0
                        matrixGS[jmmm + jj * 4 * mMax + 2, jnunu + ll * 4 * mMax + 2] = JM[jm + 1][ 27]  # Ephi_incE
                        matrixGS[jmmm + jj * 4 * mMax + 2, jnunu + ll * 4 * mMax + 3] = JM[jm + 1][28]  # Ephi_incH

                        matrixGS[jmmm + jj * 4 * mMax + 3, jnunu + ll * 4 * mMax] = 0
                        matrixGS[jmmm + jj * 4 * mMax + 3, jnunu + ll * 4 * mMax + 1] = 0
                        matrixGS[jmmm + jj * 4 * mMax + 3, jnunu + ll * 4 * mMax + 2] = JM[jm + 1][29]  # Hphi_incE
                        matrixGS[jmmm + jj * 4 * mMax + 3, jnunu + ll * 4 * mMax + 3] = JM[jm + 1][30]  # Hphi_incH

    return abs(np.linalg.det(matrixGS))

a_0 = 100
L = 3 * a_0
n_e = 1e13
e_0 = 4.8e-10
m_e = 9.1e-28
c = 3e10
H0 = 800
w_H = (e_0 * H0) / (m_e * c)
w_p = math.sqrt(4 * math.pi * e_0 ** 2 * n_e / m_e)
w_0 = 0.5 * w_H
k_0 = w_0 / c
    # w_0
    # N_max=mMax=...=3
    # a_0
    # ee=1
cylXY = np.array([[-L, 0.001], [0, 0.001]])
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
EE = 1 + w_p ** 2 / (w_H ** 2 - w_0 ** 2)
GG = -w_p ** 2 * w_H / ((w_H ** 2 - w_0 ** 2) * w_0)
HH = 1 - w_p ** 2 / w_0 ** 2
p = [i for i in range(2,10)]
n=0
X=[]
#спросить про случай p=1 ,следовательно Q=0 и получаем деление на ноль при J[n][15]
#for оp in range(2, 10):
#    result = dispeq_gyrotropic_cylindersArrayFull(1, w_0, 100, 1, cylXY, k_0, оp, EE, GG, HH, 3e10)
#   X.append(result)
print(dispeq_gyrotropic_cylindersArrayFull(1, w_0, 100, 1, cylXY, k_0,3, EE, GG, HH, 3e10))
#plt.plot(p, X)
#plt.show()