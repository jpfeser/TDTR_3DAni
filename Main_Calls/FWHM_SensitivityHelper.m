%%% TDTR Testing Program
%% clear memory
clear all
close all

%% More Complicated
lambda=[17 0.183 6.5]; %W/m-K
phi = ones(size(lambda)); %radial anisotropic ratio, usually = 1 for isotropic.
lambda_tensor = MakeRadial_ktensor(lambda,phi); %initialize tensor assuming radial symmetry in all layers
%lambda_tensor is of form lamba_tensor(layer_number,direction_index)
%direction index is 1 (kxx), 2 (kyy), 3 (kzz), 4 (kxy), 5 (kxz), 6 (kyz)
lambda_tensor(3,2) = 10.5; %change to 3D anisotropic for some layers
C = [2.65 0.1 2.1]*1e6;
h = [36.6 1 1e6]*1e-9;

f=1.11e6; %laser Modulation frequency, Hz
tdelay = -200e-12;

wp_x = 1.0334e-6; %a-direction
wp_y = 1.0264e-6; %c-direction

ws_x = wp_x; %pump 1/e^2 radius, m
ws_y = wp_y; %probe 1/e^2 radius, m

tau_rep=1/80e6;
Qp=1e-3;
TCR=1e-4;

%save base case
lambda_tensor_base = lambda_tensor;
h_base = h;
C_base = C;
wp_x_base = wp_x;
wp_y_base = wp_y;

i = 1 %base case
CalcFWHMs

fid = fopen('sensitivities.dat','w');
%%
i = 2 %vary kM
lambda_tensor(1,1:3)=lambda_tensor_base(1,1:3)*1.01;
CalcFWHMs
S_x(i) = ((x_FWHM(i)-x_FWHM(1))/x_FWHM(1))/0.01;
S_y(i) = ((y_FWHM(i)-y_FWHM(1))/y_FWHM(1))/0.01;
fprintf(fid,'kM\tx\t%6.3f\ty\t%6.3f\n',[S_x(i) S_y(i)])
lambda_tensor=lambda_tensor_base; %set the properties back when your done

i = 3 %vary kz
lambda_tensor(3,3)=lambda_tensor_base(3,3)*1.01;
CalcFWHMs
S_x(i) = ((x_FWHM(i)-x_FWHM(1))/x_FWHM(1))/0.01;
S_y(i) = ((y_FWHM(i)-y_FWHM(1))/y_FWHM(1))/0.01;
fprintf(fid,'kz\tx\t%6.3f\ty\t%6.3f\n',[S_x(i) S_y(i)])
lambda_tensor=lambda_tensor_base; %set the properties back when your done
%%
i = 4 %vary kx
lambda_tensor(3,1)=lambda_tensor_base(3,1)*1.05;
CalcFWHMs
S_x(i) = ((x_FWHM(i)-x_FWHM(1))/x_FWHM(1))/0.05;
S_y(i) = ((y_FWHM(i)-y_FWHM(1))/y_FWHM(1))/0.05;
fprintf(fid,'kx\tx\t%6.3f\ty\t%6.3f\n',[S_x(i) S_y(i)])
lambda_tensor=lambda_tensor_base; %set the properties back when your done

i = 5 %vary ky
lambda_tensor(3,2)=lambda_tensor_base(3,2)*1.05;
CalcFWHMs
S_x(i) = ((x_FWHM(i)-x_FWHM(1))/x_FWHM(1))/0.05;
S_y(i) = ((y_FWHM(i)-y_FWHM(1))/y_FWHM(1))/0.05;
fprintf(fid,'ky\tx\t%6.3f\ty\t%6.3f\n',[S_x(i) S_y(i)])
lambda_tensor=lambda_tensor_base; %set the properties back when your done

i = 6 %vary CM
C(1)=C_base(1)*1.05;
CalcFWHMs
S_x(i) = ((x_FWHM(i)-x_FWHM(1))/x_FWHM(1))/0.05;
S_y(i) = ((y_FWHM(i)-y_FWHM(1))/y_FWHM(1))/0.05;
fprintf(fid,'CM\tx\t%6.3f\ty\t%6.3f\n',[S_x(i) S_y(i)])
C=C_base; %set the properties back when your done

i = 7 %vary hM
h(1)=h_base(1)*1.05;
CalcFWHMs
S_x(i) = ((x_FWHM(i)-x_FWHM(1))/x_FWHM(1))/0.05;
S_y(i) = ((y_FWHM(i)-y_FWHM(1))/y_FWHM(1))/0.05;
fprintf(fid,'hM\tx\t%6.3f\ty\t%6.3f\n',[S_x(i) S_y(i)])
h=h_base; %set the properties back when your done
%%
i = 8 % vary wx_0
wp_x = wp_x_base*1.05;
ws_x = wp_x_base*1.05;
CalcFWHMs
S_x(i) = ((x_FWHM(i)-x_FWHM(1))/x_FWHM(1))/0.05;
S_y(i) = ((y_FWHM(i)-y_FWHM(1))/y_FWHM(1))/0.05;
fprintf(fid,'wx\tx\t%6.3f\ty\t%6.3f\n',[S_x(i) S_y(i)])
wp_x = wp_x_base;
ws_x = wp_x_base;

i = 9 % vary wy_0
wp_y = wp_y_base*1.05;
ws_y = wp_y_base*1.05;
CalcFWHMs
S_x(i) = ((x_FWHM(i)-x_FWHM(1))/x_FWHM(1))/0.05;
S_y(i) = ((y_FWHM(i)-y_FWHM(1))/y_FWHM(1))/0.05;
fprintf(fid,'wy\tx\t%6.3f\ty\t%6.3f\n',[S_x(i) S_y(i)])
wp_y = wp_y_base;
ws_y = wp_y_base;

i = 10 %vary G
lambda_tensor(2,3)=lambda_tensor_base(2,3)*1.05;
CalcFWHMs
S_x(i) = ((x_FWHM(i)-x_FWHM(1))/x_FWHM(1))/0.05;
S_y(i) = ((y_FWHM(i)-y_FWHM(1))/y_FWHM(1))/0.05;
fprintf(fid,'G\tx\t%6.3f\ty\t%6.3f\n',[S_x(i) S_y(i)])
lambda_tensor=lambda_tensor_base; %set the properties back when your done
%%
i = 11 %vary C substrate
C(3)=C_base(3)*1.05;
CalcFWHMs
S_x(i) = ((x_FWHM(i)-x_FWHM(1))/x_FWHM(1))/0.05;
S_y(i) = ((y_FWHM(i)-y_FWHM(1))/y_FWHM(1))/0.05;
fprintf(fid,'CS\tx\t%6.3f\ty\t%6.3f\n',[S_x(i) S_y(i)])
C(3)=C_base(3); %set the properties back when your done

fclose(fid);
save('sensitivities.mat')