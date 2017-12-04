import math
import cmath
import matplotlib.pyplot as plt
import numpy as np
import scipy.special as special
from sympy import besselj, hankel2


def dispeq_gyrotropic_cylindersArrayFull(N_max, w_0, a_0, ee, cylXY, k, p, EE, GG, HH, c, ):
    m = list(range(-N_max, N_max + 1, 1))
    if N_max == 0:
        m = 0

    nu = m
    w = w_0
    k0 = w / c
    q = cmath.sqrt(1 - p **2)
    q = q * (2 * ((q.imag) <= 0) - 1)
    Q = k0 * a_0 * q

    # вычисляем элементы матрицы рассеяния для внутренней области цилиндра
    mainq = EE ** 2 - GG ** 2 + EE * HH - (HH + EE) * p ** 2
    #F=cmath.sqrt(EE*GG**2*HH*(GG**2-(HH-EE)**2))
    #Pb=cmath.sqrt(EE-(HH+EE)*GG**2/(HH-EE)**2-2*F/(HH-EE)**2)
    #Pc=cmath.sqrt(EE-(HH+EE)*GG**2/(HH-EE)**2+2*F/(HH-EE)**2)
    radq = cmath.sqrt((HH - EE) ** 2 * p **4 + 2 * ((GG ** 2)*(HH + EE) - EE * (HH - EE) **2) * p **2 +(EE ** 2 - GG ** 2 - EE * HH) ** 2)
    #radq = -(HH-EE)*cmath.sqrt((p**2-Pb**2)*(p**2-Pc**2))
    q1 = cmath.sqrt(0.5 * (mainq - radq) / EE)
    q2 = cmath.sqrt(0.5 * (mainq + radq) / EE)

    n1 = -(EE / (p * GG)) * (p ** 2 + q1 ** 2 + (GG ** 2)/EE - EE)
    n2 = -(EE / (p * GG)) * (p ** 2 + q2 ** 2 + (GG ** 2)/EE - EE)

    alp1 = -1 + (p ** 2 + q1 ** 2 - EE) / GG
    alp2 = -1 + (p ** 2 + q2 ** 2 - EE) / GG

    bet1 = 1 + p / n1
    bet2 = 1 + p / n2
#что делать с Q она мнимая, брать действительную или мнимую часть?
    Q1 = (q1 * a_0) * k0
    Q2 = (q2 * a_0) * k0
    mMax = 2 * N_max + 1
    JM = np.zeros(((2 * N_max + 1, 31)), dtype=np.complex)
    n=0
    for per in range(-N_max,N_max+1,1):
        JM[n][0] = special.jv(per + 1, Q1) #JM1
        JM[n][1] = special.jv(per + 1, Q2) #JM2
        JM[n][2]= special.jv(per, Q1)  #Jm1
        JM[n][3]= special.jv(per, Q2)  #Jm2
        JM[n][4]= JM[n][2] / Q1  #Jm1_Q1
        JM[n][5]= JM[n][3]/ Q2  #Jm2_Q2

        JM[n][6]= (((1j / HH) * n1) * q1) * JM[n][2]  #Ez1
        JM[n][7] = (((1j / HH) * n2) * q2) * JM[n][3] #Ez2
        JM[n][8] = 1j * (JM[n][0] + (alp1 * per) * JM[n][4])  #Ephi1
        JM[n][9]= 1j * (JM[n][1] + (alp2 * per) *JM[n][5] )  #Ephi2
        JM[n][10]= - q1 * JM[n][2]  #Hz1
        JM[n][11]= - q2 * JM[n][3]  #Hz2
        JM[n][12] = - n1 * (JM[n][0] - (bet1 * per) * JM[n][4])  #Hphi1
        JM[n][13]= - n2 * (JM[n][1] - (bet2 * per) * JM[n][5])  #Hphi2

        # вычисляем элементы матрицы рассеяния для внешней области цилиндра

        JM[n][14] = special.hankel2(per, Q)  #H2m
        JM[n][15]= (JM[n][14]* per)/Q - special.hankel2(per + 1, Q)  #dH2m
        #Jm_out = besselj(m, Q)
        #dJm_out = Jm_out * m / Q - besselj(m + 1, Q)

        JM[n][16] = q * JM[n][14]  #Ez_sct
        JM[n][17]= q * JM[n][14]  #Hz_sct

        A = 1 / (k0 * q**2)
        JM[n][18] = -q * (A * ((p * (per / a_0)) * JM[n][14])) #Ephi_sctE
        JM[n][19]= 1j*q * (A * (k0 * q * JM[n][15])) #Ephi_sctH
        JM[n][20]= -q * (A * ((p * (per / a_0)) * JM[n][14])) #Hphi_sctH
        JM[n][21] = -1j*q * (A * (k0 * q * JM[n][15])) #Hphi_sctE
        n = n + 1
    N_cylinders = int((cylXY.size)/2)
    matrixGS = np.zeros(((mMax*4*N_cylinders, mMax*4*N_cylinders)), dtype=np.complex)

    for jj in range( 0, N_cylinders ):
        for jm in range(1,  mMax+1, 1):
            for ll in range( 0, N_cylinders ):
                for jnu in range(1, mMax+1, 1):
                    jmmm = (jm - 1) * 4
                    jnunu = (jnu - 1) * 4
                        # спросить про проход ведь будет совпадение номеров
                    if jj == ll:
                        if jm == jnu:
                                matrixGS[jmmm + jj * 4 * mMax, jnunu + ll * 4 * mMax] = JM[jm-1][6]  # Ez1
                                matrixGS[jmmm + jj * 4 * mMax, jnunu + ll * 4 * mMax + 1] = JM[jm-1][7]  # Ez2(jm)
                                matrixGS[jmmm + jj * 4 * mMax, jnunu + ll * 4 * mMax + 2] = -JM[jm-1][16]  # Ez_sct(jm)
                                matrixGS[jmmm + jj * 4 * mMax, jnunu + ll * 4 * mMax + 3] = 0

                                matrixGS[jmmm + jj * 4 * mMax + 1, jnunu + ll * 4 * mMax] = JM[jm - 1][10]  # Hz1
                                matrixGS[jmmm + jj * 4 * mMax + 1, jnunu + ll * 4 * mMax + 1] = JM[jm - 1][11]  # Hz2
                                matrixGS[jmmm + jj * 4 * mMax + 1, jnunu + ll * 4 * mMax + 2] = 0
                                matrixGS[jmmm + jj * 4 * mMax + 1, jnunu + ll * 4 * mMax + 3] = -JM[jm - 1][17]  # -Hz_sct

                                matrixGS[jmmm + jj * 4 * mMax + 2, jnunu + ll * 4 * mMax] = JM[jm - 1][8]  # Ephi1
                                matrixGS[jmmm + jj * 4 * mMax + 2, jnunu + ll * 4 * mMax + 1] = JM[jm - 1][9]  # Ephi2
                                matrixGS[jmmm + jj * 4 * mMax + 2, jnunu + ll * 4 * mMax + 2] = -JM[jm - 1][18]  # -Ephi_sctE
                                matrixGS[jmmm + jj * 4 * mMax + 2, jnunu + ll * 4 * mMax + 3] = -JM[jm - 1][19]  # -Ephi_sctH

                                matrixGS[jmmm + jj * 4 * mMax + 3, jnunu + ll * 4 * mMax] = JM[jm - 1][12]  # Hphi1
                                matrixGS[jmmm + jj * 4 * mMax + 3, jnunu + ll * 4 * mMax + 1] = JM[jm - 1][13]  # Hphi2
                                matrixGS[jmmm + jj * 4 * mMax + 3, jnunu + ll * 4 * mMax + 2] = -JM[jm - 1][21]  # -Hphi_sctE
                                matrixGS[jmmm + jj * 4 * mMax + 3, jnunu + ll * 4 * mMax + 3] = -JM[jm - 1][20]  # -Hphi_sctH
                    else:
                        Ljl = float(math.sqrt((cylXY[jj , 0] - cylXY[ll , 0]) ** 2 + (cylXY[jj , 1] - cylXY[ll , 1]) ** 2))

                        thetaIJ =float(math.atan(abs(cylXY[jj , 1] - cylXY[ll , 1]) / abs(cylXY[jj , 0] - cylXY[ll , 0])))

                        if ((cylXY[jj , 0] - cylXY[ll , 0]) <= 0 and (cylXY[jj , 1] - cylXY[ll , 1]) > 0):
                            thetaIJ = cmath.pi - thetaIJ
                        elif ((cylXY[jj , 0] - cylXY[ll , 0]) <= 0 and (cylXY[jj, 1] - cylXY[ll , 1]) <= 0):
                            thetaIJ = cmath.pi + thetaIJ
                        elif ((cylXY[jj , 1] - cylXY[ll , 1]) <= 0 and (cylXY[jj, 0] - cylXY[ll , 0]) > 0):
                            thetaIJ = 2 * cmath.pi - thetaIJ



                        JM[jnu-1][22] =  -cmath.exp(-1j * (nu[jnu-1] - m[jm-1]) * thetaIJ) * q *  special.hankel2((m[jm-1] - nu[jnu-1]), k0* q * Ljl)  # HmInll

                        JM[jnu-1][23] = special.jv(m[jnu-1], Q)  # Jm_out
                        JM[jnu-1][24] = JM[jnu-1][23] * m[jnu-1] / Q - special.jv(m[jnu-1] + 1, Q)  # dJm_out
# спросить почему они одинаковые?
                        JM[jnu-1][25] = JM[jnu-1][22] * JM[jnu-1][23]  # Ez_inc
                        JM[jnu-1][26] = JM[jnu-1][22] * JM[jnu-1][23]  # Hz_inc

                        A = 1 / (k0 * (1 - p ** 2))
                        JM[jnu-1][27] = - JM[jnu-1][22] * (A * ((p * (m[jnu-1] / a_0)) * JM[jnu-1][23]))  # Ephi_incE
                        JM[jnu-1][28] = - JM[jnu-1][22] * (A * (- 1j * k0 * q * JM[jnu-1][24]))  # Ephi_incH
                        JM[jnu-1][29] = - JM[jnu-1][22] * (A * (1j * k0 * q * JM[jnu-1][24]))  # Hphi_incE
                        JM[jnu-1][30] = -JM[jnu-1][22] * (A * ((p * (m[jnu-1] / a_0)) * JM[jnu-1][23]))  # Hphi_incH



                        matrixGS[jmmm + jj * 4 * mMax, jnunu + ll * 4 * mMax] = 0
                        matrixGS[jmmm + jj * 4 * mMax, jnunu + ll * 4 * mMax + 1] = 0
                        matrixGS[jmmm + jj * 4 * mMax, jnunu + ll * 4 * mMax + 2] = JM[jnu-1][25]  # Ez_inc
                        matrixGS[jmmm + jj * 4 * mMax, jnunu + ll * 4 * mMax + 3] = 0

                        matrixGS[jmmm + jj * 4 * mMax + 1, jnunu + ll * 4 * mMax] = 0
                        matrixGS[jmmm + jj * 4 * mMax + 1, jnunu + ll * 4 * mMax + 1] = 0
                        matrixGS[jmmm + jj * 4 * mMax + 1, jnunu + ll * 4 * mMax + 2] = 0
                        matrixGS[jmmm + jj * 4 * mMax + 1, jnunu + ll * 4 * mMax + 3] = JM[jnu-1][26]  # Hz_inc

                        matrixGS[jmmm + jj * 4 * mMax + 2, jnunu + ll * 4 * mMax] = 0
                        matrixGS[jmmm + jj * 4 * mMax + 2, jnunu + ll * 4 * mMax + 1] = 0
                        matrixGS[jmmm + jj * 4 * mMax + 2, jnunu + ll * 4 * mMax + 2] = JM[jnu-1][27]  # Ephi_incE
                        matrixGS[jmmm + jj * 4 * mMax + 2, jnunu + ll * 4 * mMax + 3] = JM[jnu-1][28]  # Ephi_incH

                        matrixGS[jmmm + jj * 4 * mMax + 3, jnunu + ll * 4 * mMax] = 0
                        matrixGS[jmmm + jj * 4 * mMax + 3, jnunu + ll * 4 * mMax + 1] = 0
                        matrixGS[jmmm + jj * 4 * mMax + 3, jnunu + ll * 4 * mMax + 2] = JM[jnu-1][29]  # Hphi_incE
                        matrixGS[jmmm + jj * 4 * mMax + 3, jnunu + ll * 4 * mMax + 3] = JM[jnu-1][30]  # Hphi_incH

    return abs(np.linalg.det(matrixGS))
#matrixGS
#abs(np.linalg.det(matrixGS))
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
p = [i for i in np.arange(1.1,10,0.1)]
n=0
X=[]
Y=[]

#for оp in np.arange(1.1,10,0.1):
#    result =math.log(dispeq_gyrotropic_cylindersArrayFull(1, w_0, a_0, 1, cylXY1, k_0, оp, EE, GG, HH, 3e10))
#    X.append(result)
#print(dispeq_gyrotropic_cylindersArrayFull(1, w_0, 100, 1, cylXY, k_0,1.1, EE, GG, HH, 3e10))
#for yp in np.arange(1.1,10,0.1):
#    result =(dispeq_gyrotropic_cylindersArrayFull(1, w_0, a_0, 1, cylXY, k_0, yp, EE, GG, HH, 3e10))
#    Y.append(math.log(result**2))
#plt.plot(p, X, p, Y)
#plt.show()