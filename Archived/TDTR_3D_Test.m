%%% TDTR Testing Program
%% clear memory
clear all

%% More Complicated
lambda = [142 0.2 1.3];
phi = ones(size(lambda));
lambda_tensor = MakeRadial_ktensor(lambda,phi);
C = [2.42 0.1 1.62]*1e6;
h = [70 1 1e6]*1e-9;
%lambda_tensor(:,1:2)=0;
TCR = 1e-4;

%% Generate laser parameters
wp = 25e-6;
ws = 25e-6;
Qp = 30e-3*pi/2;
f = 9.8e6;
tau_rep = 1/80e6;

%% which time delays?
tdelay = logspace(log10(100e-12),log10(4e-9),30)';

%% Do simulation
tstart_3D = tic;
[ratio_3D,deltaR_3D] = TDTR_3DAni_getREFL(tdelay,lambda_tensor,C,h,f,tau_rep,wp,Qp,ws,TCR)
telapsed_3D = toc(tstart_3D)

%% Now for radial case
tstart_radial = tic;
[deltaR_radial,ratio_radial]=TDTR_REFL_VV3(tdelay,TCR,tau_rep,f,lambda,C,h,phi,wp,ws,Qp)
telapsed_radial = toc(tstart_radial)

%% plot a comparison of the answers
semilogx(tdelay,ratio_3D,'o',tdelay,ratio_radial,'-')
axis([100e-12 4e-9 0 2.0])
axis square
