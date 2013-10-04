function [Kernal] = TDTR_3DAni_getTintegrand_v1(xi,eta,f,k,C,h,wp,Qp,ws,xoffset,yoffset)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%[G] = TDTR_3DAni_getG(xi,eta,f,k,C,h);
[G]=TDTR_3DAni_getG_v1(xi,eta,f,k,C,h);
P = 1/(2*pi)*Qp*exp(-1/8*wp^2*(xi.^2+eta.^2)); 
S = 1/(2*pi)*exp(-1/8*ws^2*((-xi).^2+(-eta).^2));
Kernal = G.*P.*S.*exp(-sqrt(-1)*(xoffset*xi + yoffset*eta));
%[G] = ones(size(P));
end

