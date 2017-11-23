% /*===========================================================================
% 
% DESCRIPTION
%       Функция вычисляет значение "дисперсионной" функции для заданной
%       постоянной распространения для изотропного цилиндра.
% 
%   Copyright (c) 2005-2016 by Vasiliy Es'kin. All Rights Reserved.
% ===========================================================================*/
% 
%                       EDIT HISTORY FOR FILE
% 
%   This section contains comments describing changes made to the module.
%   Notice that changes are listed in reverse chronological order.
% 
% when       who              what, where, why
% --------   ---       ----------------------------------------------------------
% 27/01/20016 Vasiliy Es'kin   Create programma.
% ==========================================================================*/
function [z, Q1, Q2, Q] = dispeq_gyrotropic_cylinder(m, pp)

global R EE GG HH ee R_psi k_0 a_0

p = pp;

    mainq = EE.^2 - GG.^2 + EE.* HH - (HH + EE).* p.^2;
    radq = sqrt((HH - EE).^2 * p.^4 + 2 * ((GG.^2).* (HH + EE) - EE.* (HH - EE).^2) * p.^2 +...
        (EE.^2 - GG.^2 - EE.* HH).^2);
    q1 = sqrt(0.5 * (mainq - radq)./ EE);
    q2 = sqrt(0.5 * (mainq + radq)./ EE);

    n1 = -(EE.* (p.* GG).^(-1)).* (p.^2 + q1.^2 + (GG.^2)./ EE - EE);
    n2 = -(EE.* (p * GG).^(-1)).* (p.^2 + q2.^2 + (GG.^2)./ EE - EE);
    
    alp1 = -1 + (p.^2 + q1.^2 - EE)./ GG;
    alp2 = -1 + (p.^2 + q2.^2 - EE)./ GG;
    
    bet1 = 1 + p./ n1;
    bet2 = 1 + p./ n2;

    Q1 = (q1.* a_0).* k_0;
    Q2 = (q2.* a_0).* k_0;  
    
    JM1    = besselj(m+1, Q1);
    JM2    = besselj(m+1, Q2);
    Jm1    = besselj(m, Q1);
    Jm2    = besselj(m, Q2);
    Jm1_Q1 = Jm1./ Q1;
    Jm2_Q2 = Jm2./ Q2;
    
    Ez1    = (((1i./HH).* n1).* q1).* Jm1;
    Ez2    = (((1i./HH).* n2).* q2).* Jm2;
    Ephi1  = 1i * (JM1 + (alp1.* m).* Jm1_Q1);
    Ephi2  = 1i * (JM2 + (alp2.* m).* Jm2_Q2);
    Hz1    = - q1.* Jm1;
    Hz2    = - q2.* Jm2;
    Hphi1  = - n1.* (JM1 - (bet1 * m).* Jm1_Q1);
    Hphi2  = - n2.* (JM2 - (bet2 * m).* Jm2_Q2);
    
    
    %%% компоненты рассеянного поля
    q = sqrt(1 - p.^2);
    q = q.* (2*(imag(q) <= 0)-1);
    Q = k_0.* a_0 * q;
%     q = q.* (2*(imag(q) <= 0)-1); %%% мнимые части S и p должный быть разных знаков
%     %%% необходимо, чтобы при стремлении аргумента функции Ханкеля к бесконечности сама функция была ограниченной  
    H2m  = besselh(m, 2, Q);
    dH2m = (H2m * m)./ Q  - besselh(m + 1, 2, Q);

    %%% equation for E_rho
    A = 1./ (k_0 * (1 - p.^2));   
    Ez_sct1  = -q.* H2m;
    Hz_sct1  = -q.* H2m;
    
    Ephi_sctE1  =  q.* (A.* ((p.* (m./a_0)).* H2m));
    Ephi_sctH1  =  q.* (A.* (- 1i * k_0 * q.* dH2m));
    Hphi_sctE1  =  q.* (A.* (  1i * k_0 * q.* dH2m));
    Hphi_sctH1  =  q.* (A.* ((p.* (m./a_0)).* H2m));

 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%   %E_azimutal
    aa11 = Ephi1;
    aa12 = Ephi2;
    aa13 = Ephi_sctE1;
    aa14 = Ephi_sctH1;
% %  %E_z   
    aa21 = Ez1;
    aa22 = Ez2;
    aa23 = Ez_sct1;
    aa24 = Q1 - Q1;
% %  %H_azimutal     
    aa31 = Hphi1;
    aa32 = Hphi2;
    aa33 = Hphi_sctE1;
    aa34 = Hphi_sctH1;
% %  %H_z   
    aa41 = Hz1;
    aa42 = Hz2;
    aa43 = Q1 - Q1;
    aa44 = Hz_sct1;
% %%%%%%%%%%%s%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% product decomposition of matrix
zz = abs(aa14.* aa23.* aa32.* aa41 - aa13.* aa24.* aa32.* aa41 -...
         aa14.* aa22.* aa33.* aa41 + aa12.* aa24.* aa33.* aa41 +...
         aa13.* aa22.* aa34.* aa41 - aa12.* aa23.* aa34.* aa41 -...
         aa14.* aa23.* aa31.* aa42 + aa13.* aa24.* aa31.* aa42 +...
         aa14.* aa21.* aa33.* aa42 - aa11.* aa24.* aa33.* aa42 -...
         aa13.* aa21.* aa34.* aa42 + aa11.* aa23.* aa34.* aa42 +...
         aa14.* aa22.* aa31.* aa43 - aa12.* aa24.* aa31.* aa43 -...
         aa14.* aa21.* aa32.* aa43 + aa11.* aa24.* aa32.* aa43 +...
         aa12.* aa21.* aa34.* aa43 - aa11.* aa22.* aa34.* aa43 -...
         aa13.* aa22.* aa31.* aa44 + aa12.* aa23.* aa31.* aa44 +...
         aa13.* aa21.* aa32.* aa44 - aa11.* aa23.* aa32.* aa44 -...
         aa12.* aa21.* aa33.* aa44 + aa11.* aa22.* aa33.* aa44);

     
z = zz;


