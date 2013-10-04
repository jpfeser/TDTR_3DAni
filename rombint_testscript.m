%trap2d_testscript
clear all
%% Define materials properties
k=zeros(3,6);
k(1,:)=[140 140 140 0 0 0];
k(2,:)=[0.2 0.2 0.2 0 0 0];

theta=0*pi/180;
R=[cos(theta) -sin(theta) 0;
    sin(theta) cos(theta) 0;
    0 0 1];
Rtrans=R';
kmat=[1.3 0 0;
    0 1.3 0;
    0 0 1.3];
kmat_rot=R*kmat*Rtrans;
k(3,1)=kmat_rot(1,1);
k(3,2)=kmat_rot(2,2);
k(3,3)=kmat_rot(3,3);
k(3,4)=kmat_rot(1,2);
k(3,5)=kmat_rot(1,3);
k(3,6)=kmat_rot(2,3);

C=ones(3,1)*1e6;
 h=[70 1 1e6]'*1e-9;
% k = [2000 2000 2000 0 0 0];
% C = 1e6;
% h = 1e-3;

ws=14e-6;
wp=14e-6;
Qp=10e-3;

%% Make ETA's, NU's
Nz=2500;
fvect=linspace(-99.99e3,100.01e3,Nz);

oldintegral=zeros(Nz,1);
boundary = 2/sqrt(ws^2+wp^2)
a=-boundary;
b=boundary;
c=-boundary;
d=boundary;
integral=zeros(Nz,1);

%% Vectorized Concept
% chunksize=50;
% for chunk=1:length(fvect)/chunksize
%     integral(((chunk-1)*chunksize+1):(chunk*chunksize))=rombint2D(@(ETA,NU,FF) TDTR_3DAni_getTintegral(ETA,NU,FF,k,C,h,wp,Qp,ws),a,b,c,d,fvect(((chunk-1)*chunksize+1):(chunk*chunksize)));
% end
%% Parallel Concept
% parfor i=1:Nz;
%     integral(i)=rombint2D(@(ETA,NU,FF) TDTR_3DAni_getTintegral(ETA,NU,FF,k,C,h,wp,Qp,ws),a,b,c,d,fvect(i));
% end


%% Serial Concept
for i=1:Nz;
     integral(i)=rombint2D(@(ETA,NU,FF) TDTR_3DAni_getTintegral(ETA,NU,FF,k,C,h,wp,Qp,ws),a,b,c,d,fvect(i));
     
%     integral=trap2D(@(X,Y,Z) 1/pi*exp(-(X.^2+Y.^2)),X,Y,FF)
 
%     eta = linspace(-10e5,10e5,Nx); 
%     nu = linspace(-10e5,10e5,Ny);
%     [ETA,NU,FF]=ndgrid(eta,nu,fvect);
%     integral=trap2D(@(ETA,NU,FF) TDTR_3DAni_getTintegral(ETA,NU,FF,k,C,h,wp,Qp,ws),ETA,NU,FF);
%     integral(1)

end