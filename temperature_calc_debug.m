%Check the frequency response against Cahill's
clear all

%% More Complicated
lambda = [150 0.2 141];
phi = ones(size(lambda));
lambda_tensor = MakeRadial_ktensor(lambda,phi);
C = [2.42 0.1 1.62]*1e6;
h = [100 1 1e6]*1e-9;
%lambda_tensor(:,1:2)=0;
TCR = 1e-4;

%% Generate laser parameters
wp = 5e-6;
ws = 5e-6;
Qp = 20e-3;

%% Generate integration parameters
kmax=1/sqrt(wp^2+ws^2)*2;
a = -2*pi*kmax;
b = 2*pi*kmax;
c = -2*pi*kmax;
d = 2*pi*kmax;

%% Choose the range of frequencies in Hz
%fvect = logspace(log10(100e3),log10(100e9),200);
fvect = 1e6 + (-1250:1250)*80e6;

%% Do calculation for radial code
dT_radial=rombint_VV3(@(kvect) TDTR_TEMP_VV3(kvect,fvect,lambda,C,h,phi,wp,ws,Qp),0,kmax,length(fvect)); %f must be a row vector

total = tic;
%% Test the 3D Code
for i = 1:length(fvect)
    tstart=tic;
    [dT_3D(i),npoints(i)]=rombint2D(@(XI,ETA,FF) TDTR_3DAni_getTintegral(XI,ETA,FF,lambda_tensor,C,h,wp,Qp,ws),a,b,c,d,fvect(i));
    tcycle(i)=toc(tstart);
end
ttotal=toc(total)

%% Plot the results
loglog(fvect(:),real(dT_radial(:)),'b-',fvect(:),imag(dT_radial(:)),'g-',fvect(:),real(dT_3D(:)),'bo',fvect(:),-imag(dT_3D(:)),'go');
%axis([1e5 3e9 -2 2])



