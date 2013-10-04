function [ktensor] = MakeIsotropic_ktensor(k)
%MakeIsotropic_ktensor(k)
%   is a shortcut to make an isotropic thermal conductivity tensor for
%   multiples layers.  The output is an N x 6 matrix where the six values
%   in the row are ktensor = [kxx,kyy,kzz,kxy,kxz,kyz] = [k,k,k,0,0,0].  k
%   is a vector of length N, where N is also the number of layers in the
%   material stack.

k=k(:); %forces k to be a column vector
ktensor = zeros(length(k),6);
ktensor(:,1:3)=[k,k,k];
end

