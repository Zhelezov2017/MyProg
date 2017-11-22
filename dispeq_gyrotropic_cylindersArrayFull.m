%%%%%%%
%%%   N_max      - maximum number of azimuthal index m,
%%%   N_cylinder - number of cylinder
%%%   cylX       - x coordinates of cylinders
%%%   cylY       - y coordinates of cylinders


function zz = dispeq_gyrotropic_cylindersArrayFull(N_max, w_0, a_0, ee, cylXY, k, p, EE, GG, HH, c)

m  = [-N_max:N_max];
if(N_max == 0) m = 0; end
nu = m;
ee = 1.
w = w_0;
k0 = w / c;
q = sqrt(1 - p.^2);
q = q.* (2*(imag(q) <= 0)-1);
Q = k0.* a_0 * q;

%%%%%% ��������� �������� ������� ��������� ��� ���������� ������� ��������
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

    Q1 = (q1.* a_0).* k0;
    Q2 = (q2.* a_0).* k0;
    
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


%%%%%% ��������� �������� ������� ��������� ��� ������� ������� ��������    
    H2m  = besselh(m, 2, Q);
    dH2m = (H2m.* m)./ Q  - besselh(m + 1, 2, Q);
%     H2m  = m./ m;
%     dH2m = (m)./ Q  - besselh(m + 1, 2, Q)./ besselh(m, 2, Q);
    
    Jm_out  = besselj(m, Q);
    dJm_out = Jm_out.* m / Q  - besselj(m + 1, Q);

    Ez_sct  = q.* H2m;    
    Hz_sct  = q.* H2m; 
    
    A = 1./ (k0 * (1 - p.^2));
    Ephi_sctE  = -q.* (A.* ((p.* (m./a_0)).* H2m));
    Ephi_sctH  = -q.* (A.* (- 1i * k0 * q.* dH2m));    
    Hphi_sctH  = -q.* (A.* ((p.* (m./a_0)).* H2m));
    Hphi_sctE  = -q.* (A.* (  1i * k0 * q.* dH2m));


N_cylinders = size(cylXY,1);
     
mMax = (2 * N_max + 1); 
%%% ������ ������������ ������� ������. 
                           
       matrixGS = zeros(4 *mMax * N_cylinders, 4 *mMax * N_cylinders);
%%%��������� ����� �� 2 ��-�� ������� ���� �����������.

for jj = 0:N_cylinders-1
    for jm = 1:size(m,2)
        for ll = 0:N_cylinders-1
            for jnu = 1:size(nu,2)
                
                jmmm  = ((jm-1)*4 + 1);
                jnunu = ((jnu-1)*4 + 1);
                
                if(jj == ll)
                   if(jm == jnu)
                       %%% continuity of Ez
                       matrixGS(jmmm + jj * 4 *mMax    , jnunu + ll * 4 *mMax)     = Ez1(jm);
                       matrixGS(jmmm + jj * 4 *mMax    , jnunu + ll * 4 *mMax + 1) = Ez2(jm);
                       matrixGS(jmmm + jj * 4 *mMax    , jnunu + ll * 4 *mMax + 2) = -Ez_sct(jm);
                       matrixGS(jmmm + jj * 4 *mMax    , jnunu + ll * 4 *mMax + 3) = 0;
                       %%% continuity of Hz
                       matrixGS(jmmm + jj * 4 *mMax + 1, jnunu + ll * 4 *mMax)     = Hz1(jm);
                       matrixGS(jmmm + jj * 4 *mMax + 1, jnunu + ll * 4 *mMax + 1) = Hz2(jm);
                       matrixGS(jmmm + jj * 4 *mMax + 1, jnunu + ll * 4 *mMax + 2) = 0;
                       matrixGS(jmmm + jj * 4 *mMax + 1, jnunu + ll * 4 *mMax + 3) = -Hz_sct(jm);
                       %%% continuity of Ephi
                       matrixGS(jmmm + jj * 4 *mMax + 2, jnunu + ll * 4 *mMax)     = Ephi1(jm);
                       matrixGS(jmmm + jj * 4 *mMax + 2, jnunu + ll * 4 *mMax + 1) = Ephi2(jm);
                       matrixGS(jmmm + jj * 4 *mMax + 2, jnunu + ll * 4 *mMax + 2) = -Ephi_sctE(jm);
                       matrixGS(jmmm + jj * 4 *mMax + 2, jnunu + ll * 4 *mMax + 3) = -Ephi_sctH(jm);
                       %%% continuity of Hphi
                       matrixGS(jmmm + jj * 4 *mMax + 3, jnunu + ll * 4 *mMax)     = Hphi1(jm);
                       matrixGS(jmmm + jj * 4 *mMax + 3, jnunu + ll * 4 *mMax + 1) = Hphi2(jm);
                       matrixGS(jmmm + jj * 4 *mMax + 3, jnunu + ll * 4 *mMax + 2) = -Hphi_sctE(jm);
                       matrixGS(jmmm + jj * 4 *mMax + 3, jnunu + ll * 4 *mMax + 3) = -Hphi_sctH(jm);
                   end
                   
                   
                else
                    %%% distance from jth cylinder to lth
                    Ljl     = sqrt((cylXY(jj + 1, 1) - cylXY(ll + 1, 1))^2 + (cylXY(jj + 1, 2) - cylXY(ll + 1, 2))^2);
                    %%% angle of jth cylinder at the coordinate of lth
                    thetaIJ = atan(abs(cylXY(jj + 1,2) - cylXY(ll + 1,2)) / abs(cylXY(jj + 1,1) - cylXY(ll + 1,1)));

                    if((cylXY(jj + 1,1) - cylXY(ll + 1,1)) <=0 && (cylXY(jj + 1,2) - cylXY(ll + 1,2)) > 0)
                        thetaIJ = pi - thetaIJ;
                    elseif((cylXY(jj + 1,1) - cylXY(ll + 1,1)) <= 0 && (cylXY(jj + 1,2) - cylXY(ll + 1,2)) <= 0)
                        thetaIJ = pi + thetaIJ;
                    elseif((cylXY(jj + 1,2) - cylXY(ll + 1,2)) <= 0 && (cylXY(jj + 1,1) - cylXY(ll + 1,1)) > 0)
                        thetaIJ = 2 * pi - thetaIJ;
                    end
                    
                    
                    HmInll = -exp(- 1i * (nu(jnu) - m(jm)) * thetaIJ) * besselh((m(jm) - nu(jnu)),2, k* q* Ljl); 
                    
                    Jm_out  = besselj(m(jnu), Q);
                    dJm_out = Jm_out * m(jnu) / Q  - besselj(m(jnu) + 1, Q);
                    
                    Ez_inc  = HmInll.* Jm_out;
                    Hz_inc  = HmInll.* Jm_out;
                    
                    A = 1./ (k0 * (1 - p.^2));
                    Ephi_incE  = -HmInll.* (A.* ((p.* (m(jnu)./a_0)).* Jm_out));
                    Ephi_incH  = -HmInll.* (A.* (- 1i * k0 * q.* dJm_out));
                    Hphi_incE  = -HmInll.* (A.* (  1i * k0 * q.* dJm_out));
                    Hphi_incH  = -HmInll.* (A.* ((p.* (m(jnu)./a_0)).* Jm_out));
                    
                    %%% continuity of Ez
                    matrixGS(jmmm + jj * 4 *mMax    , jnunu + ll * 4 *mMax)     = 0;
                    matrixGS(jmmm + jj * 4 *mMax    , jnunu + ll * 4 *mMax + 1) = 0;
                    matrixGS(jmmm + jj * 4 *mMax    , jnunu + ll * 4 *mMax + 2) = Ez_inc;
                    matrixGS(jmmm + jj * 4 *mMax    , jnunu + ll * 4 *mMax + 3) = 0;
                    %%% continuity of Hz
                    matrixGS(jmmm + jj * 4 *mMax + 1, jnunu + ll * 4 *mMax)     = 0;
                    matrixGS(jmmm + jj * 4 *mMax + 1, jnunu + ll * 4 *mMax + 1) = 0;
                    matrixGS(jmmm + jj * 4 *mMax + 1, jnunu + ll * 4 *mMax + 2) = 0;
                    matrixGS(jmmm + jj * 4 *mMax + 1, jnunu + ll * 4 *mMax + 3) = Hz_inc;
                    %%% continuity of Ephi
                    matrixGS(jmmm + jj * 4 *mMax + 2, jnunu + ll * 4 *mMax)     = 0;
                    matrixGS(jmmm + jj * 4 *mMax + 2, jnunu + ll * 4 *mMax + 1) = 0;
                    matrixGS(jmmm + jj * 4 *mMax + 2, jnunu + ll * 4 *mMax + 2) = Ephi_incE;
                    matrixGS(jmmm + jj * 4 *mMax + 2, jnunu + ll * 4 *mMax + 3) = Ephi_incH;
                    %%% continuity of Hphi
                    matrixGS(jmmm + jj * 4 *mMax + 3, jnunu + ll * 4 *mMax)     = 0;
                    matrixGS(jmmm + jj * 4 *mMax + 3, jnunu + ll * 4 *mMax + 1) = 0;
                    matrixGS(jmmm + jj * 4 *mMax + 3, jnunu + ll * 4 *mMax + 2) = Hphi_incE;
                    matrixGS(jmmm + jj * 4 *mMax + 3, jnunu + ll * 4 *mMax + 3) = Hphi_incH;                    
                end
            end
        end
    end
end

zz = abs(det(matrixGS));



            



