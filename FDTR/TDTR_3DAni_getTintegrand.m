function [Kernal] = TDTR_3DAni_getTintegrand(xi,eta,f,k,C,h,wp_x,wp_y,Qp,ws_x,ws_y,xoffset,yoffset)
%TDTR_3DAni_getTintegrand calculates the integrand of the FDTR integral
%   Detailed explanation goes here


G = TDTR_3DAni_getG_savespace(xi,eta,f,k,C,h); %FDTR point response
P = 1/(2*pi)*Qp*exp(-1/8*(wp_x^2*xi.^2+wp_y^2*eta.^2)); %Elliptical Pump
S = 1/(2*pi)*exp(-1/8*(ws_x^2*(-xi).^2+wp_y^2*(-eta).^2)); %Elliptical Probe
S = S.*exp(-sqrt(-1)*(xoffset*xi + yoffset*eta)); %offsets the beam

Kernal = G.*P.*S; %integrand
end

