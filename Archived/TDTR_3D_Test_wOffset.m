%%% TDTR Testing Program
%% clear memory
clear all

%% More Complicated
lambda = [142 0.2 5.5];
phi = ones(size(lambda));
phi(3)=2000/lambda(3);
lambda_tensor = MakeRadial_ktensor(lambda,phi);
C = [2.42 0.1 1.62]*1e6;
h = [70 1 1e6]*1e-9;
%lambda_tensor(:,1:2)=0;
TCR = 1e-4;

%% Generate laser parameters
wp = 2.5e-6;
ws = 2.5e-6;
Qp = 30e-3*pi/2;
f = 0.625e6;
tau_rep = 1/80e6;

%% which time delays?
tdelay = logspace(log10(100e-12),log10(4e-9),30)';
xoffset = 0;

%% Do simulation
tstart_3D = tic;
[ratio_3D,deltaR_3D] = TDTR_3DAni_getREFL_wOffset(tdelay,lambda_tensor,C,h,f,tau_rep,wp,Qp,ws,TCR,xoffset,0)
telapsed_3D = toc(tstart_3D)

%% Now for radial case
tstart_radial = tic;
SysParam.lambda=lambda;
SysParam.C=C;
SysParam.t=h;
SysParam.r_probe=ws;
SysParam.r_pump=wp;
SysParam.eta=phi;
SysParam.tau_rep=tau_rep;
SysParam.A_pump=Qp;
SysParam.TCR=TCR;
SysParam.f=f;;
[deltaR_radial,ratio_radial]=TDTR_REFL_DOUGHNUT_V2(tdelay,SysParam,xoffset)
telapsed_radial = toc(tstart_radial)

%% plot a comparison of the answers
semilogx(tdelay,ratio_3D,'o',tdelay,ratio_radial,'-')
%axis([100e-12 4e-9 0 2.0])
axis square
