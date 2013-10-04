function [G] = TDTR_3DAni_getG(eta,nu,f,k,C,h)
%TDTR_3DAni_getG evaluates the fourier transform of the Green's function G=Ts/qpointsource 
%   (eta,nu): are vectors containing fourier conjugates of (x,y); they must
%   be the same size.
%   f: is the frequency of the oscillating unit heat source (in this
%   program f is a scaler).
%   lambda: is an array containing the thermal conductivity tensor elements
%   (kx,ky,kz,kxy,kxz,kyz) in the column and different layers in the rows.
%   The last row is the substrate, the first row is the surface.
%   C: is the specific heat vector (each layer in a different row)
%   h: the thickness of each layer; the program sets the last layer to
%   infinity, but a value must still be entered.

%% Extract fourier coordinates
eta2=eta.^2;
nu2=nu.^2;
etanu=eta.*nu;

%% Extract Anisotropic Ratios
phi_xx=k(:,1)./k(:,3);
phi_yy=k(:,2)./k(:,3);
phi_xy=k(:,4)./k(:,3);
phi_xz=k(:,5)./k(:,3);
phi_yz=k(:,6)./k(:,3);

%% rewrite in rad/s instead of Hz
ii=sqrt(-1);
omega=2*pi*f;

%% Construct properties of the bottom layer
Nlayers=length(C);
q2=(ii*omega.*C(Nlayers)/k(Nlayers,3));
lambda1= q2+(phi_xx(Nlayers)*eta2+2*phi_xy(Nlayers)*etanu+phi_yy(Nlayers)*nu2);
lambda2= 2*ii*(phi_xz(Nlayers)*eta+phi_yz(Nlayers)*nu);
un_plus = -lambda2/2 + sqrt((lambda2/2).^2 +lambda1);
un_minus = -lambda2/2 - sqrt((lambda2/2).^2 +lambda1);
gamman_plus = k(Nlayers,3)*un_plus + k(Nlayers,5)*eta + k(Nlayers,6)*nu;
gamman_minus = k(Nlayers,3)*un_minus + k(Nlayers,5)*eta + k(Nlayers,6)*nu;

%%Initialize large arrays
Bplus=zeros(size(eta));
Bminus=ones(size(eta));

%% Use transfer matrix math to get to the top layer
if Nlayers~=1
    for n=Nlayers:-1:2
        q2=(ii*omega.*C(n-1)/k(n-1,3));
        
        lambda1= q2+(phi_xx(n-1)*eta2+2*phi_xy(n-1)*etanu+phi_yy(n-1)*nu2);
        lambda2= 2*ii*(phi_xz(n-1)*eta+phi_yz(n-1)*nu);
        unnext_plus = -lambda2/2 + sqrt((lambda2/2).^2 +lambda1);
        unnext_minus = -lambda2/2 - sqrt((lambda2/2).^2 +lambda1);
        gammannext_plus = k(n-1,3)*unnext_plus + k(n-1,5)*eta + k(n-1,6)*nu;
        gammannext_minus = k(n-1,3)*unnext_minus + k(n-1,5)*eta + k(n-1,6)*nu;
        
        % find next B+
        AA = gammannext_minus-gamman_plus;
        BB = gammannext_minus-gamman_minus;
        expterm1=exp(unnext_plus*h(n-1)).*(AA.*Bplus+BB.*Bminus);
        
        % find next B- (recycle AA,BB,expterm to save memory)
        AA = -gammannext_plus+gamman_plus;
        BB = -gammannext_plus+gamman_minus;
        expterm2=exp(unnext_minus*h(n-1)).*(AA.*Bplus+BB.*Bminus);
        
        Bplus = expterm1./(-gammannext_plus+gammannext_minus);
        Bminus = expterm2./(-gammannext_plus+gammannext_minus);
        
        %% These next 3 lines fix a numerical stability issue if one of the
        %% layers is very thick or resistive;
        %penetration_logic=logical(h(n-1)*abs(unminus)>100);  %if pentration is smaller than layer...set to semi-inf
        %Bplus(penetration_logic)=0;
        %Bminus(penetration_logic)=1;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %now the new becomes old, and the process starts over.
        gamman_plus=gammannext_plus;
        gamman_minus=gammannext_minus;
    end
end

G = -(Bplus + Bminus)./(gamman_plus.*Bplus+gamman_minus.*Bminus);

end

