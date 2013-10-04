k=zeros(3,6);
k(1,:)=[140 140 140 0 0 0];
k(2,:)=[0.2 0.2 0.2 0 0 0];

theta=45*pi/180;
R=[cos(theta) -sin(theta) 0;
    sin(theta) cos(theta) 0;
    0 0 1];
Rtrans=R';
kmat=[5 0 0;
    0 2000 0;
    0 0 2000];
kmat_rot=R*kmat*Rtrans;
k(3,1)=kmat_rot(1,1);
k(3,2)=kmat_rot(2,2);
k(3,3)=kmat_rot(3,3);
k(3,4)=kmat_rot(1,2);
k(3,5)=kmat_rot(1,3);
k(3,6)=kmat_rot(2,3);

C=ones(3,1)*1e6;
h=[70 1 1e6]'*1e-9;

fvect=1e6;%(1:2500)'*1e6;
 eta = linspace(-10e5,10e5,60); 
 nu = linspace(-10e5,10e5,60);
[ETA,NU,FF]=ndgrid(eta,nu,fvect);
 
%eta=rand(2,1);
%nu=rand(2,1);
%[ETA,FF]=meshgrid(eta,fvect);
%[NU,FF]=meshgrid(nu,fvect);

G= TDTR_3DAni_getG(ETA,NU,FF,k,C,h);

% contour(ETA,NU,G)
% figure(gcf)
% 
% Var1_Grid=ETA;
% Var2_Grid=NU;
% Var3_Grid=G;
% 
% filename_out='5-2000-2000_90rot.txt';
% D=size(Var1_Grid)
% %initialize
% Var1_Vect=zeros(D(1)*D(2),1);
% Var2_Vect=Var1_Vect;
% Var3_Vect=Var1_Vect;
% 
% 
% sum=0;
% for i=1:D(1)
%     for j=1:D(2);
%         sum=sum+1;
%         Var1_Vect(sum)=Var1_Grid(i,j);
%         Var2_Vect(sum)=Var2_Grid(i,j);
%         Var3_Vect(sum)=Var3_Grid(i,j);
%     end
% end
% 
% dlmwrite(filename_out,[Var1_Vect,Var2_Vect,real(Var3_Vect),imag(Var3_Vect)],'delimiter','\t');

