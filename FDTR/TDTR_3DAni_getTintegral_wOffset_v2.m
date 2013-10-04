function [Kernal] = TDTR_3DAni_getTintegral_wOffset_v2(xi,eta,f,k,C,h,wp_x,wp_y,Qp,ws_x,ws_y,xoffset,yoffset)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%[G] = TDTR_3DAni_getG(xi,eta,f,k,C,h);
[G] = TDTR_3DAni_getG_savespace(xi,eta,f,k,C,h);
%P = 1/(2*pi)*Qp*exp(-1/8*wp^2*(xi.^2+eta.^2)); 
%S = 1/(2*pi)*exp(-1/8*ws^2*((-xi).^2+(-eta).^2));

P = 1/(2*pi)*Qp*exp(-1/8*(wp_x^2*xi.^2+wp_y^2*eta.^2)); 
S = 1/(2*pi)*exp(-1/8*(ws_x^2*(-xi).^2+wp_y^2*(-eta).^2));
Kernal = G.*P.*S.*exp(-sqrt(-1)*(xoffset*xi + yoffset*eta));
end

