%%% TDTR Testing Program
%% clear memory
clear all
close all

%% Material Properties to Simulate
lambda=[17 0.183 6.5]; %W/m-K
phi = ones(size(lambda)); %radial anisotropic ratio, usually = 1 for isotropic.
lambda_tensor = MakeRadial_ktensor(lambda,phi); %initialize tensor assuming radial symmetry in all layers
%lambda_tensor is of form lamba_tensor(layer_number,direction_index)
%direction index is 1 (kxx), 2 (kyy), 3 (kzz), 4 (kxy), 5 (kxz), 6 (kyz)
%where z is the through plane direction
lambda_tensor(3,2) = 14; %change to 3D anisotropic for some layers
C = [2.65 0.1 2.0]*1e6;
h = [36.6 1 1e6]*1e-9;

f=1.11e6; %pump Modulation frequency, Hz
wp_x = 1.0334e-6; %x-direction 1/e2 radius, pump
wp_y = 1.0264e-6; %y-direction 1/e2 radius, pump

ws_x = wp_x; %x-direction 1/e2 radius, probe
ws_y = wp_y; %x-direction 1/e2 radius, probe

tau_rep=1/80e6; %laser repetition rate
Qp=1e-3; %average pump power
TCR=1e-4; %thermoreflectance coefficient

%% which time delays?
%tdelay = logspace(log10(200e-12),log10(4e-9),30)';

tdelay = -200e-12;

lambda_perp_vect = linspace(3,12,10) %conditions to simulate

for i = 1:length(lambda_perp_vect)
    lambda_tensor(3,1) = lambda_perp_vect(i); % conditions to simulate, a-dir (low k), 

%% First calculate concentric case
xoffset = 0;
yoffset = 0;

[Vout_max,Vin_max,ratio] = TDTR_3DAni_TDTR_Sig(tdelay,lambda_tensor,C,h,f,tau_rep,wp_x,wp_y,Qp,ws_x,ws_y,TCR,xoffset,yoffset);

%% Now find FWHM of xoffset direction
yoffset =0;
Lpx = sqrt(lambda_tensor(3,1)/C(3)*1/(2*pi*f));
xmin_guess = sqrt(wp_x^2+Lpx^2)*0.8;
xmax_guess = sqrt(wp_x^2+Lpx^2)*1.3;
x_HWHM(i) = fzero(@(x) 0.5*Vout_max - TDTR_3DAni_TDTR_Sig(tdelay,lambda_tensor,C,h,f,tau_rep,wp_x,wp_y,Qp,ws_x,ws_y,TCR,x,yoffset),[xmin_guess,xmax_guess],optimset('Display','iter','TolX',1e-9))

%% Now find FWHM of yoffset direction
xoffset =0;
Lpy = sqrt(lambda_tensor(3,2)/C(3)*1/(2*pi*f));
ymin_guess = sqrt(wp_y^2+Lpy^2)*0.8;
ymax_guess = sqrt(wp_y^2+Lpy^2)*1.3;
y_HWHM(i) = fzero(@(x) 0.5*Vout_max - TDTR_3DAni_TDTR_Sig(tdelay,lambda_tensor,C,h,f,tau_rep,wp_x,wp_y,Qp,ws_x,ws_y,TCR,xoffset,x),[ymin_guess,ymax_guess],optimset('Display','iter','TolX',1e-9))

FWHM_ratio(i) = x_HWHM(i)/y_HWHM(i);
eta(i) =  lambda_tensor(3,1)/ lambda_tensor(3,2);

end

%% post-processing
x_FWHM = 2*x_HWHM; %in the program, it actually measures the half-width half-max
y_FWHM = 2*y_HWHM;

%% file outputs
save('FWHM_Map_XQuartz_varykxdir=y_fixkz=6p5_fixky=14_f=1p11Mhz.mat')
dlmwrite('FWHM_Map_XQuartz_varykxdir=y_fixkz=6p5_fixky=14_f=1p11Mhz.txt',[lambda_perp_vect',x_FWHM',y_FWHM'],'delimiter','\t');

%% plots
plot(lambda_perp_vect',x_FWHM',lambda_perp_vect',y_FWHM')