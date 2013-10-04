function [f] = TDTR_3DAni_getG_savespace(xi,eta,f,k,C,h)
%TDTR_3DAni_getG evaluates the fourier transform of the Green's function G=Ts/qpointsource 
%   (xi,eta): are vectors containing fourier conjugates of (x,y); they must
%   be the same size.
%   f: is the frequency of the oscillating unit heat source (in this
%   program f is a scaler).
%   lambda: is an array containing the thermal conductivity tensor elements
%   (kx,ky,kz,kxy,kxz,kyz) in the column and different layers in the rows.
%   The last row is the substrate, the first row is the surface.
%   C: is the specific heat vector (each layer in a different row)
%   h: the thickness of each layer; the program sets the last layer to
%   infinity, but a value must still be entered.

%% Extract Anisotropic Ratios
phi_xx=k(:,1)./k(:,3);
phi_yy=k(:,2)./k(:,3);
phi_xy=k(:,4)./k(:,3);
phi_xz=k(:,5)./k(:,3);
phi_yz=k(:,6)./k(:,3);

%% rewrite in rad/s instead of Hz
ii=sqrt(-1);

%% Construct properties of the bottom layer
Nlayers=length(C);
q2=(ii*2*pi*f.*C(Nlayers)/k(Nlayers,3));
gamman_plus= q2+(phi_xx(Nlayers)*xi.^2+2*phi_xy(Nlayers)*xi.*eta+phi_yy(Nlayers)*eta.^2);
%^ should be lambda1: using gamman_plus to save space
gamman_minus= 2*ii*(phi_xz(Nlayers)*xi+phi_yz(Nlayers)*eta);
%^ should be lambda2: using gamman_plus to save space
un_plus = -gamman_minus/2 + sqrt((gamman_minus/2).^2 +gamman_plus); %using surrogates for lambda1, lambda2
un_minus = -gamman_minus/2 - sqrt((gamman_minus/2).^2 +gamman_plus);%using surrogates for lambda1, lambda2
gamman_plus = k(Nlayers,3)*un_plus + ii*(k(Nlayers,5)*xi + k(Nlayers,6)*eta);
gamman_minus = k(Nlayers,3)*un_minus + ii*(k(Nlayers,5)*xi + k(Nlayers,6)*eta);

%%Initialize large arrays
Bplus=zeros(size(xi));
Bminus=ones(size(xi));

%% Use transfer matrix math to get to the top layer
if Nlayers~=1
    for n=Nlayers:-1:2
        q2=(ii*2*pi*f.*C(n-1)/k(n-1,3));
        
        gammannext_plus= q2+(phi_xx(n-1)*xi.^2+2*phi_xy(n-1)*xi.*eta+phi_yy(n-1)*eta.^2);
        %^ should be lambda1: using gamman_plus to save memory
        gammannext_minus= 2*ii*(phi_xz(n-1)*xi+phi_yz(n-1)*eta);
        %^ should be lambda2: using gamman_plus to save memory
        unnext_plus = -gammannext_minus/2 + sqrt((gammannext_minus/2).^2 +gammannext_plus); %using surrogates for lambda1, lambda2
        unnext_minus = -gammannext_minus/2 - sqrt((gammannext_minus/2).^2 +gammannext_plus);%using surrogates for lambda1, lambda2
        gammannext_plus = k(n-1,3)*unnext_plus + ii*(k(n-1,5)*xi + k(n-1,6)*eta);
        gammannext_minus = k(n-1,3)*unnext_minus + ii*(k(n-1,5)*xi + k(n-1,6)*eta);
        
        q2 =Bplus; %prevents me accidentally overwriting and subsequently misusing Bplus in the following lines; re-using q2 to save memory
        
        % find next B+
        % find next B-
        Bplus = exp(-unnext_plus*h(n-1)).*(((gammannext_minus-gamman_plus).*q2+(gammannext_minus-gamman_minus).*Bminus))./(-gammannext_plus+gammannext_minus);
        Bminus = exp(-unnext_minus*h(n-1)).*(((-gammannext_plus+gamman_plus).*q2+(-gammannext_plus+gamman_minus).*Bminus))./(-gammannext_plus+gammannext_minus);
        
        %% These next 3 lines fix a numerical stability issue if one of the
        % layers is very thick or resistive;
        penetration_logic=logical(h(n-1)*abs(unnext_plus)>100);  %if pentration is smaller than layer...set to semi-inf
        Bplus(penetration_logic)=0;
        Bminus(penetration_logic)=1;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %% now the new becomes old, and the process starts over.
        gamman_plus=gammannext_plus;
        gamman_minus=gammannext_minus;
    end
end

f = -(Bplus + Bminus)./(gamman_plus.*Bplus+gamman_minus.*Bminus); %reusing f, instead of G to save memory

end

