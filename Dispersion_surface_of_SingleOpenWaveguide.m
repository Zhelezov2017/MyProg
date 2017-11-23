% /*===========================================================================
% 
% DESCRIPTION
%       ѕрограмма вычисл€ет зависимость "диспесионной" функции от
%       действительной и мнимой частей посто€нной распространени€ при
%       фиксированной частоте.
% 
%   Copyright (c) 2005-2013 by Vasiliy Es'kin. All Rights Reserved.
% ===========================================================================*/
% 
%                       EDIT HISTORY FOR FILE
% 
%   This section contains comments describing changes made to the module.
%   Notice that changes are listed in reverse chronological order.
% 
% when       who              what, where, why
% --------   ---       ----------------------------------------------------------
% 11/09/2007 Vasiliy Es'kin   Create programma.
% ==========================================================================*/

clear all
 GPC_systemParameters
global EE GG HH k_0 a_0 m
% EE = 10;
% HH = EE;
% HH = 5;
% GG = -1e-6;
% GG = 1;


eps = 2
y0  = 0
x0  = 0
xmin =  25;
xmax =  40;
ymin =  -4;
ymax =   4;
Npntx = 100;
Npnty = 100;

systemParameters
m = 1;
EE = EE1;
GG = GG1;
HH = HH1;
R   = k_0 * a_0;


[px,py] = meshgrid(xmin : (xmax - xmin) / Npntx :xmax, ymin : (ymax - ymin) / Npnty :ymax);
p = px + 1i * py;

zz1 = (dispeq_gyrotropic_cylinder(m,p));
zz2 = abs(zz1);


% zz2 = log(zz2);
figure
mesh(px,py, zz2)
title('zz11');
xlabel('Re p');
ylabel('Im p');
set(gca,'ZScale','log');
view(0,0)


% ff = fminsearch(@(x) dispeq_gyrotropic_cylinder(x(1) + 1i * x(2)), [6,   -96], optimset('TolX',1e-8))
% ff1= fminsearch(@(x) dispeq_gyrotropic_cylinder(x(1) + 1i * x(2)), [6.4, -51.2], optimset('TolX',1e-8))
% ff2= fminsearch(@(x) dispeq_gyrotropic_cylinder(x(1) + 1i * x(2)), [6.4, -73.6], optimset('TolX',1e-8))
% 
% x1 = [6.21567438022368, -436.356136009314];
% for in = 1:20
%     x1 = fminsearch(@(x) dispeq_gyrotropic_cylinder(x(1) + 1i * x(2)), [x1(1), x1(2)-22.5105], optimset('TolX',1e-8));
%     ppp(2*in-1,1) = x1(1) + 1i * x1(2);
%     ppp(2*in,1) = x1(1) - 1i * x1(2);
% %     ppp(in,1) = x1(1) + 1i * x1(2);
% end
% ppp


% figure
% mesh(px,py,zz10)
% title('zz10');
% xlabel('Re p');
% ylabel('Im p');
% set(gca,'ZScale','log');

% for n=[1:1:3]
%     ff = fminsearch(@(x) dispeq_gyrotropic_cylinder(x(1) + 1i * x(2)),...
%         [0.0428, 124.434+n*0.93], optimset('TolX',1e-8));
%     p_n(n) = ff(1) + 1i*ff(2)
% end
% p_n = p_n'
% 
% ff = fminsearch(@(x) dispeq_isotropic_cylinder(x(1) + 1i * x(2)), [9, -35], optimset('TolX',1e-8))
% ff1= fminsearch(@(x) dispeq_isotropic_cylinder(x(1) + 1i * x(2)), [11, -61], optimset('TolX',1e-8))
% ff2= fminsearch(@(x) dispeq_isotropic_cylinder(x(1) + 1i * x(2)), [11, -86], optimset('TolX',1e-8))
% 
% x1 = [13.0111710758198, -569.793726883961];
% for in = 1:10
%     x1 = fminsearch(@(x) dispeq_isotropic_cylinder(x(1) + 1i * x(2)), [x1(1), x1(2)-24.1], optimset('TolX',1e-8));
%     ppp(2*in-1,1) = x1(1) + 1i * x1(2);
%     ppp(2*in,1) = x1(1) - 1i * x1(2);
% end
% ppp


% ff3= fminsearch(@(x) dispeq_isotropic_cylinder(x(1) + 1i * x(2)), [2.76, 0], optimset('TolX',1e-8))
% ff4= fminsearch(@(x) dispeq_isotropic_cylinder(x(1) + 1i * x(2)), [2.9, 0], optimset('TolX',1e-8))
% ff5= fminsearch(@(x) dispeq_isotropic_cylinder(x(1) + 1i * x(2)), [2.92, 0], optimset('TolX',1e-8))
% ff6= fminsearch(@(x) dispeq_isotropic_cylinder(x(1) + 1i * x(2)), [3.03, 0], optimset('TolX',1e-8))
% ff7= fminsearch(@(x) dispeq_isotropic_cylinder(x(1) + 1i * x(2)), [3.11, 0], optimset('TolX',1e-8))
% ff8= fminsearch(@(x) dispeq_isotropic_cylinder(x(1) + 1i * x(2)), [2.92, 0], optimset('TolX',1e-8))

% ff = fminsearch(@(x) dispeq_isotropic_cylinder(x(1) + 1i * x(2)), [1.131, 0], optimset('TolX',1e-8))
% ff1= fminsearch(@(x) dispeq_isotropic_cylinder(x(1) + 1i * x(2)), [1.158, 0], optimset('TolX',1e-8))
% ff2= fminsearch(@(x) dispeq_isotropic_cylinder(x(1) + 1i * x(2)), [1.302, 0], optimset('TolX',1e-8))

% ff = fminsearch(@(x) dispeq_gyrotropic_cylinder(x(1) + 1i * x(2)), [2.9235, 0], optimset('TolX',1e-8))
% ff1= fminsearch(@(x) dispeq_gyrotropic_cylinder(x(1) + 1i * x(2)), [0.03636, 14.29283], optimset('TolX',1e-8))
% ff2= fminsearch(@(x) dispeq_gyrotropic_cylinder(x(1) + 1i * x(2)), [0.03636, 15.29283], optimset('TolX',1e-8))
% ff3= fminsearch(@(x) dispeq_gyrotropic_cylinder(x(1) + 1i * x(2)), [0.03636, 16.29283], optimset('TolX',1e-8))
% ff4= fminsearch(@(x) dispeq_gyrotropic_cylinder(x(1) + 1i * x(2)), [3.5, 0], optimset('TolX',1e-8))
% ff5= fminsearch(@(x) dispeq_gyrotropic_cylinder(x(1) + 1i * x(2)), [3.85, 0], optimset('TolX',1e-8))
% ff6= fminsearch(@(x) dispeq_gyrotropic_cylinder(x(1) + 1i * x(2)), [0.35, -1.2], optimset('TolX',1e-8))
% ff6= fminsearch(@(x) dispeq_gyrotropic_cylinder(x(1) + 1i * x(2)), [2.976, 0], optimset('TolX',1e-8))


% %%%%% search of mode constant
% %%% the search of transverse constant of surface modes
% p_min = 0;
% p_Max = 5;
%     p = [p_min:0.1:p_Max];
%     
% for in = 1:size(p,2)
%     p_n(in,:) = fminsearch(@(x) dispeq_gyrotropic_cylinder(x(1) + 1i * x(2)),...
%                 [p(in),0], optimset('TolX',1e-8))
% end
% 
% %%% sorting of transverse propagation constant
% sizep_n=size(p_n);
% for i=1:sizep_n(1)-1
% for j=i+1:sizep_n(1)
%     if p_n(i,1)>p_n(j,1)
%         p_ntemp=p_n(i,:);
%         p_n(i,:)=p_n(j,:);
%         p_n(j,:)=p_ntemp;
%     end
% end
% end
% 
% %%% forming no equal mode 
% zz = p_n(1,:);
% sizey = size(p_n(:,1));
% for t = 1:sizey(1)-1
%     sizezz = size(zz(:,1));
%           y1 = p_n(t,1);
%           y2 = p_n(t + 1,1);
%           if abs(y1-y2)>0.0001 
%              zz(sizezz(1) + 1,1) = p_n(t+1,1);
% %              zz(sizezz + 1,2) = q_n(t,2);
%           end
% end
% p_n = zz
% % p_n = p_n(:,p_n(:,1)<p_Max)



% %%% the search of transverse constant of complex modes
% p_min = 0;
% p_Max = 15;
%     p = [p_min:0.1:p_Max];
%     
% for in = 1:size(p,2)
%     p_n(in,:) = fminsearch(@(x) dispeq_gyrotropic_cylinder(x(1) + 1i * x(2)),...
%                 [1, p(in)], optimset('TolX',1e-8))
% end
% 
% %%% sorting of transverse propagation constant
% sizep_n=size(p_n);
% for i=1:sizep_n(1)-1
% for j=i+1:sizep_n(1)
%     if p_n(i,1)>p_n(j,1)
%         p_ntemp=p_n(i,:);
%         p_n(i,:)=p_n(j,:);
%         p_n(j,:)=p_ntemp;
%     end
% end
% end
% 
% %%% forming no equal mode 
% zz = p_n(1,:);
% sizey = size(p_n(:,1));
% for t = 1:sizey(1)-1
%     sizezz = size(zz(:,1));
%           y1 = p_n(t,1);
%           y2 = p_n(t + 1,1);
%           if abs(y1-y2)>0.00001 
%              zz(sizezz(1) + 1,:) = p_n(t+1,:);
% %              zz(sizezz + 1,2) = q_n(t,2);
%           end
% end
% p_n = zz



