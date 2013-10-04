function [ktensor] = MakeRadial_ktensor(k,eta)
%MakeIsotropic_ktensor(k)
%   is a shortcut to make an radially symmetric thermal conductivity tensor for
%   multiples layers.  The output is an N x 6 matrix where the six values
%   in the row are ktensor = [kxx,kyy,kzz,kxy,kxz,kyz] = [kr,kr,kz,0,0,0].  k
%   is a vector of length N, where N is also the number of layers in the
%   material stack. eta is the ratio kr/kz which is set equal to one
%   (isotropic) if it is not entered.

k=k(:); %forces k to be a column vector
if nargin==1
    eta = ones(size(k));
else
    eta = eta(:);
end
ktensor = zeros(length(k),6);
ktensor(:,1:3)=[k.*eta,k.*eta,k];
end

